#ifndef _CONFIG_HH_
#define _CONFIG_HH_

/**
 * @class Config
 * @author D. Rosich
 * 
 * json reader. Loads the configuration file passed by the user
 */

#include <string>
#include <nlohmann/json.hpp>

class Config {
public:
    explicit Config(const std::string& filepath);
    ~Config() = default;

    // Detector parameters
    float get_Nd() const;
    float get_width() const;
    float get_length() const;
    float get_V_bi() const;
    float get_V_bias() const;
    float get_R() const;
    std::string get_material() const;

    // Injection parameters
    float get_focus() const;
    float get_wavelength() const;
    float get_NA() const;
    float get_refractive_index() const;
    int get_type() const;
    int get_N() const;

    // Simulation
    int get_steps() const;
    float get_dt() const;
    float get_t_pc() const;
    std::string get_sim_type() const;

private:
    nlohmann::json _data;
    void _load_json(const std::string& filepath);
};

#endif
