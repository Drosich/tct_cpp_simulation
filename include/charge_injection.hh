#ifndef _CHARGEINJECTION_HH_
#define _CHARGEINJECTION_HH_

class Charge_carrier;
class Detector;

#include <vector>
#include <utility>

class Charge_injection
{
    public:
        Charge_injection(float, float, float, float, Detector*, int);
        ~Charge_injection();

        std::vector<std::pair<float, float>> get_initial_charges(){return _charges_per_point_init;}
        std::vector<Charge_carrier*> get_charges();

        void update_speeds();

    private:
        float _wavelength;
        float _numerical_aperture;
        float _refractive_index;
        float _focus;
        int _type;
        Detector* _det = nullptr;
        std::vector<Charge_carrier*> _charges;

        std::vector<float> _x_init;
        std::vector<float> _y_init;
        std::vector<std::pair<float,float>> _charges_per_point_init;
        
        float _compute_beam_width(float);
        std::vector<std::pair<float, float>> _compute_xy_beam(int, float, float, unsigned, int);
        void _create_injection();
};

#endif