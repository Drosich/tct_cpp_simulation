#include "charge_injection.hh"
#include "charge_carrier.hh"
#include "detector.hh"
#include "utility.hh"

#include <iostream>
#include <random>
#include <cmath>
#include <algorithm>
#include <filesystem>

#define H_BAR 1.0546e-34 // Planck constant over 2pi

/**
 * @brief class constructor
 * 
 * Initializes the initial position of the charges and fills
 * the charge carrier array.
 * 
 * @param focus depth at which the laser is focused (m)
 * @param wavelength laser wavelength (m)
 * @param numerical_aperture laser numerical aperture
 * @param refractive_index detector material refractive index
 * @param det detector geometry
 * @param type type of the carriers. 0->electrons, 1->holes
 * @param N number of charges
 */
Charge_injection::Charge_injection(float focus,
                                   float wavelength, 
                                   float numerical_aperture,
                                   float refractive_index,
                                   Detector* det,
                                   int type,
                                   int N)
{
    _focus = focus;
    _wavelength = wavelength;
    _numerical_aperture = numerical_aperture;
    _refractive_index = refractive_index;
    _type = type;
    _det = det;
    _n_of_charges = N;

    _charges_per_point_init = _compute_xy_beam(_n_of_charges, -100.e-6, 100.e-6, std::random_device{}(), 200000);
    _create_injection();
    std::cout << "Simulating " << _charges.size() << " charges" << std::endl;

    std::filesystem::path cwd = std::filesystem::current_path().parent_path();
    if(_type == 0)
        readCSV(cwd.string() + "/exp_data/electron_drift_velocity.csv", _E_field_experimental_range, _velocity_exp);
    else
        readCSV(cwd.string() + "/exp_data/hole_drift_velocity.csv", _E_field_experimental_range, _velocity_exp);   
}

/**
 * @brief calculate the beam width of a gaussian beam
 * 
 * calculates the width of a laser with a gaussian emission profile
 * 
 * @param y position in the y axis (longitudinal) (m)
 * 
 * @returns width of the beam at point y (m)
 */
float Charge_injection::_compute_beam_width(float y)
{
    return std::sqrt(std::pow(_wavelength/(M_PI*_numerical_aperture), 2.) + 
                     std::pow((y - _focus)*_numerical_aperture/(_refractive_index), 2.));
}

/**
 * @brief calculate the position of the charges
 * 
 * distributes N charges in space according to the TPA-TCT charge carrier
 * density profile (gaussian^2)
 * 
 * @param N number of charges
 * @param y_min lower limit of the distribution on the y axis (m)
 * @param y_max upper limit of the distribution on the x axis (m)
 * @param seed seed for the random distribution
 * @param grid_for_max_search required for the rejection sampling algorithm
 * 
 * @returns vector of pairs with the x and y coordinates of the N charges
 */
std::vector<std::pair<float, float>> Charge_injection::_compute_xy_beam(int N,
                                                                        float y_min,
                                                                        float y_max,
                                                                        unsigned seed,
                                                                        int grid_for_max_search)
{
    if (N <= 0) throw std::invalid_argument("N must be > 0");
    if (!(y_min < y_max)) throw std::invalid_argument("y_min < y_max required");

    std::mt19937_64 gen(seed);
    std::uniform_real_distribution<float> unif_y(y_min, y_max);
    std::uniform_real_distribution<float> unif01(0.0, 1.0);

    float min_w = std::numeric_limits<float>::infinity();
    for (int i = 0; i < grid_for_max_search; ++i) {
        float frac = (float)i / (grid_for_max_search - 1);
        float y = y_min + frac * (y_max - y_min);
        float w = _compute_beam_width(y);
        if (!(w > 0.0)) throw std::runtime_error("w_of_y must return positive values");
        if (w < min_w) min_w = w;
    }
    if (!std::isfinite(min_w) || min_w <= 0.0) throw std::runtime_error("could not determine positive min w");

    float max_py = 1.0 / (min_w * min_w * min_w);

    std::vector<std::pair<float,float>> samples;
    samples.reserve(N);

    while ((int)samples.size() < N) {
        float y = unif_y(gen);
        float w = _compute_beam_width(y);
        float py = 1.0 / (w*w*w);
        float u = unif01(gen) * max_py;
        if (u <= py) {
            float sigma = w / std::sqrt(8.0);
            std::normal_distribution<float> gauss_x(0.0, sigma);
            float x = gauss_x(gen);
            samples.emplace_back(x, y);
        }
    }

    return samples;
}

/**
 * @brief change the carrier type
 * 
 * changes the type of the carrier and updates the v(E) data
 * 
 * @param type new type. 0->electrons, 1->holes
 */
void Charge_injection::set_type(int type)
{
    _type = type;
    _E_field_experimental_range.clear();
    _velocity_exp.clear();
    std::filesystem::path cwd = std::filesystem::current_path().parent_path();
    if(_type == 0)
    {
        readCSV(cwd.string() + "/exp_data/electron_drift_velocity.csv", _E_field_experimental_range, _velocity_exp);
        std::cout << "Initializing electron injection" << std::endl;
    }
    else
    {
        readCSV(cwd.string() + "/exp_data/hole_drift_velocity.csv", _E_field_experimental_range, _velocity_exp);   
        std::cout << "Initializing hole injection" << std::endl;
    }
}

/**
 * @brief creates the charges
 * 
 * initializes the charges and fills the charge injection array
 */
void Charge_injection::_create_injection()
{
    for(const auto& p : _charges_per_point_init)
    {
        _charges.emplace_back(p.first, p.second, _type);
    }
}

/**
 * @brief Updates the speeds of the carriers
 * 
 * Updates the drift velocities of the charge carriers according to the local
 * electric field at their respective positions. If the charge exits the edges
 * of the detector, the velocity is set to 0
 */
void Charge_injection::update_speeds()
{
    float E = 0.;
    float v = 0.;
    float x_lim = _det->get_depleted_width();
    if (_det->get_depleted_width() > _det->get_physical_width())
        x_lim = _det->get_physical_width();

    for (auto& charge : _charges)
    {
        auto pos = charge.get_position();
        if (pos.second > x_lim || pos.second < 0.)
            charge.set_velocity(0., 0.);
        else
        {
            E = linear_field(pos.first, pos.second, _det);
            v = linear_interpolation(E/1e8, _E_field_experimental_range, _velocity_exp)*1e-2;
            // if(_type == 0)
            // {
            //     v = 450e-4*E;
            // }
            // else
            // {
            //     v = 300e-4*E;
            // }
            
            charge.set_velocity(0., v);
        }
    }
}

/**
 * @brief get the charge injection array
 * 
 * acesses the charge injection array
 * 
 * @returns the charge carrier vector
 */
std::vector<Charge_carrier>& Charge_injection::get_charges()
{
    return _charges;
}
