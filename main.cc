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

    
    TApplication app("ParticleAnim", nullptr, nullptr);
    TCanvas* c = new TCanvas("c", "Particle Motion", 800, 600);
    gStyle->SetOptStat(0);

    // Create scatter plot (TGraph)
    TGraph* graph = new TGraph();
    TGraph* graph_h = new TGraph();
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.5);
    graph->SetMarkerColor(kBlue);
    graph_h->SetMarkerStyle(20);
    graph_h->SetMarkerSize(0.5);
    graph_h->SetMarkerColor(kRed);

    // Detector initialization
    Detector* det = new Detector(cfg.get_Nd(), cfg.get_width(), cfg.get_length(), 
                                 cfg.get_V_bi(), cfg.get_V_bias(), cfg.get_R(), 
                                 cfg.get_material());
    
    // For simple diffusion
    std::default_random_engine generator;
    std::normal_distribution<float> gaussian(0.0, 1.0);

    // Electron and hole injection
    Charge_injection* injection_e = new Charge_injection(cfg.get_focus(),
                                                         cfg.get_wavelength(),
                                                         cfg.get_NA(),
                                                         cfg.get_refractive_index(),
                                                         det,
                                                         0);
    Charge_injection* injection_h = new Charge_injection(cfg.get_focus(),
                                                         cfg.get_wavelength(),
                                                         cfg.get_NA(),
                                                         cfg.get_refractive_index(),
                                                         det,
                                                         1);
    
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
        
        graph->Set(0); // clear previous points
        graph_h->Set(0); // clear previous points
        for (size_t i = 0; i < injection_e->get_charges().size(); ++i) {
            auto pos = injection_e->get_charges()[i]->get_position();
            auto pos_h = injection_h->get_charges()[i]->get_position();
            graph->SetPoint(i, pos.first, pos.second);
            graph_h->SetPoint(i, pos_h.first, pos_h.second);
        }

        c->cd();
        graph->Draw("AP");
        graph->GetXaxis()->SetLimits(-25e-6, 25e-6);
        graph->GetYaxis()->SetRangeUser(-20e-6, 70e-6);
        graph_h->Draw("P SAME");
        c->Modified();
        c->Update();

        // Small delay so the animation is visible (~20–30 FPS)
        gSystem->ProcessEvents();
        gSystem->Sleep(30);
    }

    
    TGraph* graph_pulse = new TGraph(t.size());
    for (size_t i = 0; i < t.size(); ++i) {
        graph_pulse->SetPoint(i, t[i], signal_total[i]);
    }

    TCanvas* c_pulse = new TCanvas("c_pulse", "TPA pulse", 800, 600);
    c_pulse->cd();
    graph_pulse->Draw("APL");

    c_pulse->Update();

    app.Run();

    return 0;
}