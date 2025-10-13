#include "utility.hh"
#include "detector.hh"

#include <fstream>
#include <iostream>
#include <sstream>

/**
 * @brief calculate linear electric field
 * 
 * given the position (x,y) return the corresponding linear electric field
 * confined inside the detector. For now only a rectangular detector is
 * considered
 * 
 * @param x x coordinate (m)
 * @param y y coordinate (m)
 * @param det detector geometry
 * 
 * @returns value of the electric field a point (x,y) (V/m)
 */
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
        return 2*(V_bias+V_bi)*(1 - y/y_lim)/y_lim;
    }
}

/**
 * @brief plot electric field in the x-y plane
 * 
 * this function can be invoked to visualize the electric field inside the 
 * detector as a 2D map
 * 
 * @param x_min lower x limit (m)
 * @param x_max upper x limit (m)
 * @param nx number of points along x axis
 * @param y_min lower y limit (m)
 * @param y_max upper y limit (m)
 * @param ny number of points along y axis
 * @param det detector geometry
 * 
 * @returns ROOT TH2F histogram with the plot
 */
TH2F* plot_E_field(float x_min, float x_max, int nx, float y_min, float y_max, int ny, Detector* det)
{
    std::vector<float> x_vals;
    std::vector<float> y_vals;
    for(int i = 0; i < nx; ++i)
    {
        x_vals.push_back(x_min + i*x_max/nx);
    }
    for(int i = 0; i < ny; ++i)
    {
        y_vals.push_back(y_min + i*y_max/ny);
    }

    TH2F* hE = new TH2F("hE", "Electric Field Magnitude; x [m]; y [m]; E[V/m]",
                        nx, x_min, x_max,
                        ny, y_min, y_max);

    // Fill the histogram
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            float x = x_vals[i];
            float y = y_vals[j];
            float E = linear_field(x, y, det);
            hE->SetBinContent(i + 1, j + 1, E);
        }
    }

    return hE;
}


/**
 * @brief Read a CSV
 * 
 * Loads a csv file into memory. The csv file is assumed to lack a header and 
 * it should contain two colums x and y. For example, a file containing
 * the electron mobility as a function of the electric field.
 * 
 * @param filename: absolute path of the file
 * @param x: x variable
 * @param y: function. y = f(x)
 * 
 * @return true if the file could be opened, false otherwise
 */
bool readCSV(const std::string& filename, std::vector<float>& x, std::vector<float>& y) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return false;
    }
    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string val_x, val_y;
        if (!std::getline(ss, val_x, ',')) continue;
        if (!std::getline(ss, val_y, ',')) continue;
        try {
            float dx = std::stod(val_x);
            float dy = std::stod(val_y);
            x.push_back(dx);
            y.push_back(dy);
        } catch (...) {
            std::cerr << "Warning: Skipping invalid line: " << line << std::endl;
            continue;
        }
    }
    return true;
}

/**
 * @brief interpolates the drift velocity curve
 * 
 * Performs a linear interpolation of the drift velocity curve and returns the 
 * value of the function at the queried point
 * 
 * @param E: Query point. Electric field at which to obtain the drift velocity (MV/cm)
 * @param x: x axis of the E-v curve. Electric field (MV/cm)
 * @param y: y axis of the E-v curve. Drift speed (cm/s)
 * 
 * @return drift speed at the queried point (cm/s)
 */
float linear_interpolation(float E, std::vector<float>& x, std::vector<float>& y) {
    // If E is out of range, clamp to ends
    if (E <= x.front()) {
        return y.front();
    }
    if (E >= x.back()) {
        return y.back();
    }
    // Find interval [E_value_exp.at(i), E_value_exp[i+1]] containing E
    // std::lower_bound finds first element >= E
    auto it = std::lower_bound(x.begin(), x.end(), E);
    int idx = std::distance(x.begin(), it);

    // idx points to the element >= E, so interval is idx-1 to idx
    int i0 = idx - 1;
    int i1 = idx;

    float x0 = x[i0];
    float x1 = x[i1];
    float y0 = y[i0];
    float y1 = y[i1];

    // Linear interpolation formula
    float t = (E - x0) / (x1 - x0);
    return y0 + t * (y1 - y0);
}