#include "charge_carrier.hh"

Charge_carrier::Charge_carrier(float pos_x, float pos_y, int type)
{
    _pos_x = pos_x;
    _pos_y = pos_y;
    _type = type;
    _vel_x = 0.;
    _vel_y = 0.;
}

void Charge_carrier::set_position(float x, float y)
{
    _pos_x += x;
    _pos_y += y;
}

void Charge_carrier::set_velocity(float vx, float vy)
{
    _vel_x = vx;
    _vel_y = vy;
}

std::pair<float, float> Charge_carrier::get_position()
{
    return std::pair<float, float> {_pos_x, _pos_y}; 
}

std::pair<float, float> Charge_carrier::get_velocity()
{
    return std::pair<float, float> {_vel_x, _vel_y}; 
}

int Charge_carrier::get_type()
{
    return _type;
}