#ifndef _CHARGECARRIER_HH_
#define _CHARGECARRIER_HH_

/**
 * @class Charge_carrier
 * @author D. Rosich
 * 
 * Class representing a charge carrier inside a semiconductor. Holds position,
 * drift velocity and type (electron or hole)
*/

#include <utility>

class Charge_carrier
{
    public:
        Charge_carrier(float, float, int);
        ~Charge_carrier() = default;

        void set_position(float, float);
        void set_velocity(float, float);

        std::pair<float, float> get_position();
        std::pair<float, float> get_velocity();
        int get_type();
        
    private:
        float _pos_x;
        float _pos_y;
        float _vel_x;
        float _vel_y;
        int _type;
};

#endif