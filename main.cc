#include <iostream>
#include <vector>

#include "detector.hh"
#include "charge_injection.hh"
#include "charge_carrier.hh"

#include <TApplication.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TStyle.h>

int main()//(int argc, char** argv)
{
    // TApplication app("app", &argc, argv);

    float Nd = 1.7e20;
    float W = 50e-6;
    float L = 50e-6;
    float V_bi = 3.;
    float V_bias = 450.;
    float R = 50.;
    std::string material = "SiC";
    Detector* det = new Detector(Nd, W, L, V_bi, V_bias, R, material);
    
    float focus = 25e-6;
    float power = 44e-12;
    float TPA = 1.5e-11;
    float pulse_duration = 430e-15;
    float wavelength = 400e-9;
    float NA = 0.15;
    float refractive_index = 2.55;
    int type = 1;
    Charge_injection* injection = new Charge_injection(focus,
                                                       power,
                                                       TPA,
                                                       pulse_duration,
                                                       wavelength,
                                                       NA,
                                                       refractive_index,
                                                       det,
                                                       type);
    
    std::cout << injection->get_charges().at(0)->get_position().second << std::endl;
    // std::vector<float> x = injection->get_initial_pos_x();
    // std::vector<float> y = injection->get_initial_pos_y();
    // std::vector<float> charges = injection->get_initial_charges();

    // float x_min = *std::min_element(x.begin(), x.end());
    // float x_max = *std::max_element(x.begin(), x.end());
    // float y_min = *std::min_element(y.begin(), y.end());
    // float y_max = *std::max_element(y.begin(), y.end());
    
    // int nx = 500;
    // int ny = 500;

    // // Create 2D histogram
    // TH2F* h2 = new TH2F("h2", "Charge distribution;X;Y;Charge", nx, x_min, x_max, ny, y_min, y_max);

    // // Fill histogram: x, y, weight = charge
    // int counter = 0;
    // for (size_t i = 0; i < x.size(); ++i) {
    //     for(size_t j = 0; j < y.size(); ++j)
    //     {
    //         h2->Fill(x.at(j), y.at(i), charges.at(counter)); // z-axis = charge
    //         counter++;
    //     }
    // }

    // // Create canvas
    // TCanvas* c1 = new TCanvas("c1", "Charge distribution", 800, 600);
    // gStyle->SetOptStat(0);

    // // Draw histogram with color map
    // h2->Draw("COLZ"); // COLZ = color-coded z-axis

    // c1->Update();

    // // Keep the window open for interaction
    // app.Run();

    return 0;
}