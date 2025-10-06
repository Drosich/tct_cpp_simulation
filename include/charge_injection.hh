#ifndef _CHARGEINJECTION_HH_
#define _CHARGEINJECTION_HH_

class Charge_carrier;
class Detector;

#include <vector>

class Charge_injection
{
    public:
        Charge_injection(float, float, float, float, float, float, float, Detector*, int);
        ~Charge_injection();

        std::vector<float> get_initial_pos_x(){return _x_init;}
        std::vector<float> get_initial_pos_y(){return _y_init;}
        std::vector<int> get_initial_charges(){return _charges_per_point_init;}

    private:
        float _power;
        float _TPA;
        float _pulse_duration;
        float _wavelength;
        float _numerical_aperture;
        float _refractive_index;
        float _focus;
        int _type;
        Detector* _det = nullptr;
        std::vector<Charge_carrier*> charges;

        std::vector<float> _x_init;
        std::vector<float> _y_init;
        std::vector<int> _charges_per_point_init;

        void _compute_initial_positions();
        void _compute_charges_per_point();
        void _create_injection();

        std::vector<float> _compute_speeds();
        float _compute_beam_width(float);

        void update_positions();
        void update_speeds();
};

#endif