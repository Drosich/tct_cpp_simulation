#include "charge_injection.hh"
#include "charge_carrier.hh"
#include "detector.hh"
#include "utility.hh"

#include <iostream>
#include <random>
#include <math.h>
#include <algorithm>

#define H_BAR 1.0546e-34

Charge_injection::Charge_injection(float focus,
                                   float wavelength, 
                                   float numerical_aperture,
                                   float refractive_index,
                                   Detector* det,
                                   int type)
{
    _focus = focus;
    _wavelength = wavelength;
    _numerical_aperture = numerical_aperture;
    _refractive_index = refractive_index;
    _type = type;
    _det = det;

    _charges_per_point_init = _compute_xy_beam(10000, -32*3.66e-6, 32*3.66e-6, std::random_device{}(), 1000);
    _create_injection();
}

Charge_injection::~Charge_injection()
{
    for(size_t i = 0; i < _charges.size(); ++i)
    {
        delete _charges.at(i);
    }
}

float Charge_injection::_compute_beam_width(float y)
{
    return std::sqrt(std::pow(_wavelength/(M_PI*_numerical_aperture), 2.) + 
                     std::pow((y - _focus)*_numerical_aperture/(_refractive_index), 2.));
}

std::vector<std::pair<float, float>> Charge_injection::_compute_xy_beam(int N,
                                                                        float y_min,
                                                                        float y_max,
                                                                        unsigned seed = std::random_device{}(),
                                                                        int grid_for_max_search = 100)
{
    if (N <= 0) throw std::invalid_argument("N must be > 0");
    if (!(y_min < y_max)) throw std::invalid_argument("y_min < y_max required");

    // RNGs
    std::mt19937_64 gen(seed);
    std::uniform_real_distribution<float> unif_y(y_min, y_max);
    std::uniform_real_distribution<float> unif01(0.0, 1.0);

    // 1) Find min_w on [y_min,y_max] to compute upper bound of p(y) = 1/w(y).
    //    p(y) ∝ 1/w(y) so max_p = 1/min_w.
    float min_w = std::numeric_limits<float>::infinity();
    for (int i = 0; i < grid_for_max_search; ++i) {
        float frac = (float)i / (grid_for_max_search - 1);
        float y = y_min + frac * (y_max - y_min);
        float w = _compute_beam_width(y);
        if (!(w > 0.0)) throw std::runtime_error("w_of_y must return positive values");
        if (w < min_w) min_w = w;
    }
    if (!std::isfinite(min_w) || min_w <= 0.0) throw std::runtime_error("could not determine min w");

    float max_py = 1.0 / min_w; // p(y) ∝ 1/w(y), so max over interval ≤ 1/min_w

    std::vector<std::pair<float,float>> samples;
    samples.reserve(N);

    // Rejection sampling for y, then sample x | y
    while ((int)samples.size() < N) {
        float y = unif_y(gen);
        float w = _compute_beam_width(y);
        float py = 1.0 / w;               // unnormalized p(y)
        float u = unif01(gen) * max_py;   // sample in [0, max_py)
        if (u <= py) {
            // accept y
            float sigma = w / std::sqrt(2.0); // stddev for x
            // Normal distribution for x
            std::normal_distribution<float> gauss_x(0.0, sigma);
            float x = gauss_x(gen);
            samples.emplace_back(x, y);
        }
        // else: reject and try again
    }

    return samples;
}

void Charge_injection::_create_injection()
{
    for(size_t i = 0; i < _charges_per_point_init.size(); ++i)
    {
        _charges.push_back(new Charge_carrier(_charges_per_point_init.at(i).first, _charges_per_point_init.at(i).second, _type));
    }
}

void Charge_injection::update_speeds()
{
    float E = 0.;
    float x_lim = _det->get_depleted_width();
    if(_det->get_depleted_width() > _det->get_physical_width())
    {
        x_lim = _det->get_physical_width();
    }
    for(size_t i = 0; i < _charges.size(); ++i)
    {
        if(_charges.at(i)->get_position().second > x_lim || _charges.at(i)->get_position().second < 0.)
        {
            _charges.at(i)->set_velocity(0., 0.);
        }
        else
        {
            E = linear_field(_charges.at(i)->get_position().first, _charges.at(i)->get_position().second, _det);
            if(_type == 0)
            {
                _charges.at(i)->set_velocity(0., 450e-4*1e7);
            }
            else
            {
                _charges.at(i)->set_velocity(0., 150e-4*1e7);
            }
            
        }
        
    }
}

std::vector<Charge_carrier*> Charge_injection::get_charges()
{
    return _charges;
}