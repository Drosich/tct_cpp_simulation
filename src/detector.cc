#include "detector.hh"

#include <iostream>
#include <math.h>

#define QE 1.602e-19
#define EPS_0 8.854e-12

Detector::Detector(float Nd, float W, float L, float V_bi, float V_bias, float R, std::string material)
{
    _material = material;
    _initialize_material();

    _doping_concentration = Nd;
    _physical_width = W;
    _physical_length = L;
    _built_in_voltage = V_bi;
    _bias_voltage = V_bias;
    _resistance = R;
    _capacitance = 1.6111e-12;

    _depleted_width = _calculate_depleted_width();
    _depletion_voltage = _calculate_depletion_voltage();
}

float Detector::_calculate_depleted_width()
{
    return std::sqrt(2 * _eps * EPS_0 * _bias_voltage / (QE * _doping_concentration));
}

float Detector::_calculate_depletion_voltage()
{
    return QE * _doping_concentration * std::pow(_physical_width, 2.) / (2 * EPS_0 * _eps) - _built_in_voltage;
}

void Detector::_initialize_material()
{
    if(_material == "SiC")
    {
        _eps = 9.72;
        _e_diffusion_constant = 22.e-4;
        _h_diffusion_constant = 3.e-4;
        _e_lifetime = 1.e-9;
        _h_lifetime = 6.e-7;
    }
    else
    {
        std::cout << "Detector::_find_eps: error: unrecognised material. Defaulting to SiC" << std::endl;
        _eps = 9.72;
        _e_diffusion_constant = 22.e-4;
        _h_diffusion_constant = 3.e-4;
        _e_lifetime = 1.e-9;
        _h_lifetime = 6.e-7;
    }
}

void Detector::set_doping_concentration(float Nd)
{
    _doping_concentration = Nd;
    _detector_has_been_modified();
}

void Detector::set_material(const std::string& mat)
{
    _material = mat;
    _detector_has_been_modified();
}

void Detector::set_physical_width(float W)
{
    _physical_width = W;
    _detector_has_been_modified();
}

void Detector::set_bias_voltage(float V_bias)
{
    _bias_voltage = V_bias;
    _detector_has_been_modified();
}

void Detector::set_built_in_voltage(float V_bi)
{
    _built_in_voltage = V_bi;
    _detector_has_been_modified();
}

void Detector::set_resistance(float R)
{
    _resistance = R;
}

void Detector::set_physical_length(float L)
{
    _physical_length = L;
}

void Detector::_detector_has_been_modified()
{
    _initialize_material();
    _depleted_width = _calculate_depleted_width();
    _depletion_voltage = _calculate_depletion_voltage();
}