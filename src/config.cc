#include "config.hh"
#include <fstream>
#include <iostream>

using json = nlohmann::json;

Config::Config(const std::string& filepath) {
    _load_json(filepath);
}

void Config::_load_json(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open configuration file: " + filepath);
    }
    file >> _data;
}

// --- Detector ---
float Config::get_Nd() const { return _data["detector"]["Nd"]; }
float Config::get_width() const { return _data["detector"]["width"]; }
float Config::get_length() const { return _data["detector"]["length"]; }
float Config::get_V_bi() const { return _data["detector"]["V_bi"]; }
float Config::get_V_bias() const { return _data["detector"]["V_bias"]; }
float Config::get_R() const { return _data["detector"]["R"]; }
std::string Config::get_material() const { return _data["detector"]["material"]; }

// --- Injection ---
float Config::get_focus() const { return _data["injection"]["focus"]; }
float Config::get_wavelength() const { return _data["injection"]["wavelength"]; }
float Config::get_NA() const { return _data["injection"]["NA"]; }
float Config::get_refractive_index() const { return _data["injection"]["refractive_index"]; }
int Config::get_type() const { return _data["injection"]["type"]; }

// --- Simulation ---
int Config::get_steps() const { return _data["simulation"]["steps"]; }
float Config::get_dt() const { return _data["simulation"]["dt"]; }
