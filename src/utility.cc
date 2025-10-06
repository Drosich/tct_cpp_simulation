#include "utility.hh"
#include "detector.hh"

float linear_field(float x, float y, Detector* det)
{
    float E0 = 0., E1 = 0.;
    float V_bias = det->get_bias_voltage();
    float V_d = det->get_depletion_voltage();
    float V_bi = det->get_built_in_voltage();
    float y_lim = det->get_depleted_width();
    float diode_w = det->get_physical_width();
    float diode_l = det->get_physical_length();

    if(V_bias >= V_d)
    {
        if (y > diode_w || y < 0 || x < -diode_l/2. || x > diode_l/2.)
        {
            return 0.;
        }
        E0 = 2*(V_d + V_bi)/diode_w;
        E1 = (V_bias - V_d - V_bi)/diode_w; 
        return E0*(1 - y/diode_w) + E1;
    }
    else
    {
        if (y > y_lim || y < 0 || x < -diode_l/2. || x > diode_l/2.)
        {
            return 0.;
        }
        return 2*(V_bias+V_bi)*(1 - x/y_lim)/y_lim;
    }
}

