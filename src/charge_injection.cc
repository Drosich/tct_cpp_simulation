#include "charge_injection.hh"
#include "charge_carrier.hh"
#include "detector.hh"
#include "utility.hh"

#include <iostream>
#include <math.h>
#include <algorithm>

#define H_BAR 1.0546e-34

Charge_injection::Charge_injection(float focus,
                                   float power,
                                   float TPA, 
                                   float pulse_duration, 
                                   float wavelength, 
                                   float numerical_aperture,
                                   float refractive_index,
                                   Detector* det,
                                   int type)
{
    _focus = focus;
    _power = power;
    _TPA = TPA;
    _pulse_duration = pulse_duration;
    _wavelength = wavelength;
    _numerical_aperture = numerical_aperture;
    _refractive_index = refractive_index;
    _type = type;
    _det = det;

    _compute_initial_positions();
    _compute_charges_per_point();
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

void Charge_injection::_compute_initial_positions()
{
    size_t N = 500;
    float dx = _det->get_physical_length() / N;
    float dy = _det->get_physical_width() / N; 
    for(size_t i = 0; i < N; ++i)
    {
        _x_init.push_back(-_det->get_physical_length()/2 + i*dx);
    }
    for(size_t i = 0; i < N; ++i)
    {
        _y_init.push_back(i*dy);
    }
}

void Charge_injection::_compute_charges_per_point()
{
    // TODO YOU ARE HERE
    size_t N = 500;
    float beam_width = 0.;
    float coef = 0.;
    float charge = 0.;
    for(size_t i = 0; i < N; ++i)
    {
        for(size_t j = 0; j < N; ++j)
        {
            beam_width = _compute_beam_width(_y_init.at(i));
            coef = std::pow(_power, 2.)*_TPA*4.*std::log(2.)*_wavelength / 
            (_pulse_duration*H_BAR*2*M_PI*std::pow(M_PI, 5./2.)*std::pow(beam_width, 4.)*std::sqrt(std::log(4)));
            charge = -4*std::pow(_x_init.at(j), 2.)/std::pow(beam_width, 2.);
            if(charge < -10)
            {
                _charges_per_point_init.push_back(0.1);
            }
            else
            {
                _charges_per_point_init.push_back(coef*std::exp(charge));
            }
        }
    }
    float max_val = *std::max_element(_charges_per_point_init.begin(), _charges_per_point_init.end());
    if (max_val == 0.0f)
        throw std::runtime_error("Cannot normalize when maximum is zero");

    for (auto& val : _charges_per_point_init)
        val /= (max_val/10.);
}

void Charge_injection::_create_injection()
{
    int counter = 0;
    for(size_t i = 0; i < _x_init.size(); ++i)
    {
        for(size_t j = 0; j < _y_init.size(); ++j)
        {
            for(int k = 0; k < static_cast<int>(_charges_per_point_init.at(counter)); ++k)
            {
                _charges.push_back(new Charge_carrier(_x_init.at(j), _y_init.at(i), 1));
            }
            counter++;
        }
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