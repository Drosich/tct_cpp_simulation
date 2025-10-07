#include "charge_carrier.hh"

/**
 * @brief class constructor
 * 
 * class constructor. Initializes positions and drift velocities.
 * Drift velocities in both axes are initialized at 0.
 * 
 * @param pos_x position of the carrier in the x axis (horizontal)
 * @param pos_y position of the carrier in the y axis (vertical)
 * @param type carrier type. 0->electron, 1->hole
 */
Charge_carrier::Charge_carrier(float pos_x, float pos_y, int type)
{
    _pos_x = pos_x;
    _pos_y = pos_y;
    _type = type;
    _vel_x = 0.;
    _vel_y = 0.;
}

/**
 * @brief position setter
 * 
 * change the value of the position relative to the previous one. It is a
 * relative movement, it is NOT absolute
 * 
 * @param x relative movement in the x axis
 * @param y relative movement in the y axis
 */
void Charge_carrier::set_position(float x, float y)
{
    _pos_x += x;
    _pos_y += y;
}

/**
 * @brief drift velocity setter
 * 
 * change the drift speed of the carrier
 * 
 * @param vx drift speed in the x axis
 * @param vy drift speed in the y axis
 */
void Charge_carrier::set_velocity(float vx, float vy)
{
    _vel_x = vx;
    _vel_y = vy;
}

/**
 * @brief position getter
 * 
 * access the carrier position
 * 
 * @returns pair containing the x and y positions, in that order
 */
std::pair<float, float> Charge_carrier::get_position()
{
    return std::pair<float, float> {_pos_x, _pos_y}; 
}

/**
 * @brief drift velocity getter
 * 
 * access the carrier drift velocity
 * 
 * @returns pair containing the x and y drift velocities, in that order
 */
std::pair<float, float> Charge_carrier::get_velocity()
{
    return std::pair<float, float> {_vel_x, _vel_y}; 
}

/**
 * @brief carrier type getter
 * 
 * access the type of carrier
 * 
 * @returns carrier type as an integer (see constructor)
 */
int Charge_carrier::get_type()
{
    return _type;
}