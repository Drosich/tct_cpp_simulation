#ifndef _UTILITY_HH_
#define _UTILITY_HH_

/**
 * @brief Helper general purpose functions
 * 
 * @author D. Rosich
 */

#include <string>
#include <vector>

#include <TStyle.h>
#include <TAxis.h>
#include <TSystem.h>
#include <TH2F.h>

class Detector;

float linear_field(float, float, Detector*);
TH2F* plot_E_field(float, float, int, float, float, int, Detector*);
bool readCSV(const std::string&, std::vector<float>&, std::vector<float>&);
float linear_interpolation(float E, std::vector<float>& x, std::vector<float>& y);

#endif