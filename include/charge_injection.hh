#ifndef _CHARGEINJECTION_HH_
#define _CHARGEINJECTION_HH_

/**
 * @class Charge_injection
 * @author D. Rosich
 * 
 * Handles the injection of charge carriers into a detector. Stores and manages 
 * electrons or holes as Charge_carrier objects.
 */

#include "charge_carrier.hh"
#include "detector.hh"
#include <vector>
#include <random>
#include <utility>
#include <filesystem>

class Charge_injection
{
    public:
        Charge_injection(float, float, float, float, Detector*, int, int);
        ~Charge_injection() = default;

        void set_type(int);
        void update_speeds();

        std::vector<Charge_carrier>& get_charges();

    private:
        int _type;
        int _n_of_charges;
        float _focus;
        float _wavelength;
        float _numerical_aperture;
        float _refractive_index;
        Detector* _det;

        std::vector<Charge_carrier> _charges;
        std::vector<std::pair<float, float>> _charges_per_point_init;

        std::vector<float> _E_field_experimental_range;
        std::vector<float> _velocity_exp;

        float _compute_beam_width(float);
        std::vector<std::pair<float, float>> _compute_xy_beam(int, float, float, unsigned seed = std::random_device{}(),
                                                            int grid_for_max_search = 2000);
        void _create_injection();
};

#endif