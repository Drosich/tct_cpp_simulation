#ifndef _DETECTOR_HH_
#define _DETECTOR_HH_

#include <string>

class Detector
{
    public:
        Detector(float, float, float, float, float, float, std::string);
        ~Detector() = default;

        inline float get_doping_concentration(){return _doping_concentration;}
        inline float get_physical_width(){return _physical_width;}
        inline float get_physical_length(){return _physical_length;}
        inline float get_built_in_voltage(){return _built_in_voltage;}
        inline float get_depleted_width(){return _depleted_width;}
        inline float get_bias_voltage(){return _bias_voltage;}
        inline std::string get_material(){return _material;}
        inline float get_depletion_voltage(){return _depletion_voltage;}
        inline float get_e_diffusion_constant(){return _e_diffusion_constant;}
        inline float get_h_diffusion_constant(){return _h_diffusion_constant;}
        inline float get_e_lifetime(){return _e_lifetime;}
        inline float get_h_lifetime(){return _h_lifetime;}
        inline float get_eps(){return _eps;}
        inline float get_resistance(){return _resistance;}
        inline float get_capacitance(){return _capacitance;}

        void set_doping_concentration(float);
        void set_physical_width(float);
        void set_physical_length(float);
        void set_bias_voltage(float);
        void set_built_in_voltage(float);
        void set_material(const std::string&);
        void set_resistance(float);

    private:
        float _doping_concentration;
        float _physical_width;
        float _physical_length;
        float _bias_voltage;
        float _built_in_voltage;
        std::string _material;
        float _resistance;

        float _depleted_width;
        float _depletion_voltage;
        float _e_diffusion_constant;
        float _h_diffusion_constant;
        float _e_lifetime;
        float _h_lifetime;
        float _eps;

        float _capacitance;

        float _calculate_depleted_width();
        float _calculate_depletion_voltage();
        void _initialize_material();
        void _detector_has_been_modified();
};

#endif