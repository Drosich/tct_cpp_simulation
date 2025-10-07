#include <iostream>
#include <vector>
#include <random>
#include <filesystem>

#include "detector.hh"
#include "charge_injection.hh"
#include "charge_carrier.hh"
#include "config.hh"

#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TSystem.h>

#define QE 1.602e-19

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <path_to_config.json>" << std::endl;
        return 1;
    }
    std::string config_path = argv[1];
    std::filesystem::path cwd = std::filesystem::current_path().parent_path();
    Config cfg(cwd.string() + "/"  + config_path);

    
    TApplication app("", nullptr, nullptr);

    // Detector initialization
    Detector* det = new Detector(cfg.get_Nd(), cfg.get_width(), cfg.get_length(), 
                                 cfg.get_V_bi(), cfg.get_V_bias(), cfg.get_R(), 
                                 cfg.get_material());
    
    // For simple diffusion
    std::default_random_engine generator;
    std::normal_distribution<float> gaussian(0.0, 1.0);
    
    // simulation time array
    int steps = cfg.get_steps();
    float dt = cfg.get_dt();
    std::vector<float> t(steps);
    std::vector<float> signal_total(steps, 0.0f);
    for(int i = 0; i < steps; ++i) 
    {
        t.at(i) = i * dt;
    }

    // check limits
    float x_lim = det->get_depleted_width();
    if(det->get_depleted_width() > det->get_physical_width())
    {
        x_lim = det->get_physical_width();
    }
    
    std::vector<float> z_array;
    for(int i = 0; i < 50; ++i)
    {
        z_array.push_back(-20e-6 + i*90e-6/50);
    }
    std::vector<float> int_charge;
    for(auto z : z_array)
    {
        std::cout << "=== SIMULATING z = " << z/1.e-6 << std::endl;
        Charge_injection* injection_e = new Charge_injection(z,
                                                            cfg.get_wavelength(),
                                                            cfg.get_NA(),
                                                            cfg.get_refractive_index(),
                                                            det,
                                                            0);
        Charge_injection* injection_h = new Charge_injection(z,
                                                            cfg.get_wavelength(),
                                                            cfg.get_NA(),
                                                            cfg.get_refractive_index(),
                                                            det,
                                                            1);
        float sum = 0.;
        for(int i = 0; i < steps; ++i) //SIMULATION LOOP
        {
            std::cout << "Processing: " << i << " th step" << std::endl;
            float sigma_e = std::sqrt(2.0 * det->get_e_diffusion_constant() * dt);
            float sigma_h = std::sqrt(2.0 * det->get_h_diffusion_constant() * dt);
            sum = 0.;

            // update drift velocities
            injection_e->update_speeds();
            injection_h->update_speeds();
            
            // update positions
            for(size_t j = 0; j < injection_e->get_charges().size(); ++j)
            {
                injection_e->get_charges().at(j)->set_position(dt*injection_e->get_charges().at(j)->get_velocity().first,
                                                            dt*injection_e->get_charges().at(j)->get_velocity().second);
                injection_h->get_charges().at(j)->set_position(-dt*injection_h->get_charges().at(j)->get_velocity().first,
                                                            -dt*injection_h->get_charges().at(j)->get_velocity().second);
                if(injection_h->get_charges().at(j)->get_position().first > -det->get_physical_length()/2. &&
                injection_h->get_charges().at(j)->get_position().first < det->get_physical_length()/2. &&
                injection_h->get_charges().at(j)->get_position().second > 0.)
                {
                    injection_h->get_charges().at(j)->set_position(sigma_h*gaussian(generator),
                                                                sigma_h*gaussian(generator));
                }
                if(injection_e->get_charges().at(j)->get_position().first > -det->get_physical_length()/2. &&
                injection_e->get_charges().at(j)->get_position().first < det->get_physical_length()/2. &&
                injection_e->get_charges().at(j)->get_position().second < det->get_physical_width())
                {
                    injection_e->get_charges().at(j)->set_position(sigma_e*gaussian(generator),
                                                                sigma_e*gaussian(generator));
                }
            }

            // calculate signal. Ramo
            for(size_t j = 0; j < injection_e->get_charges().size(); j++)
            {
                sum += (injection_e->get_charges().at(j)->get_velocity().second +
                    injection_h->get_charges().at(j)->get_velocity().second);
            }
            signal_total.at(i) = sum * QE / x_lim;
        }
        TGraph graph(t.size(), t.data(), signal_total.data());
        float Q = graph.Integral();
        int_charge.push_back(Q);
        delete injection_e;
        delete injection_h;
    }

    TCanvas* c = new TCanvas("c", "Z-Scan", 800, 600);
    c->cd();
    TGraph* z_scan = new TGraph(z_array.size(), z_array.data(), int_charge.data());
    z_scan->SetTitle("z-scan;z [um];Charge [a.u.]");
    z_scan->Draw("APL");

    app.Run();

    return 0;
}