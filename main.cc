#include <iostream>
#include <vector>
#include <random>

#include "detector.hh"
#include "charge_injection.hh"
#include "charge_carrier.hh"

#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TSystem.h>

#define QE 1.602e-19

int main()
{
    std::default_random_engine generator;
    float D = 3e-3;
    std::normal_distribution<float> gaussian(0.0, 1.0);
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
    Charge_injection* injection_e = new Charge_injection(focus,
                                                       power,
                                                       TPA,
                                                       pulse_duration,
                                                       wavelength,
                                                       NA,
                                                       refractive_index,
                                                       det,
                                                       0);
    Charge_injection* injection_h = new Charge_injection(focus,
                                                        power,
                                                        TPA,
                                                        pulse_duration,
                                                        wavelength,
                                                        NA,
                                                        refractive_index,
                                                        det,
                                                        1);
    
    int steps = 2000;
    float dt = 0.005e-9;
    std::vector<float> t(steps);
    std::vector<float> signal_total(steps, 0.0f);
    for(int i = 0; i < steps; ++i) 
    {
        t.at(i) = i * dt;
    }
    float x_lim = det->get_depleted_width();
    if(det->get_depleted_width() > det->get_physical_width())
    {
        x_lim = det->get_physical_width();
    }
    
    float sum = 0.;
    float prev_x_e = 0.;
    float prev_y_e = 0.;
    float prev_x_h = 0.;
    float prev_y_h = 0.;

    for(int i = 0; i < steps; ++i) //SIMULATION LOOP
    {
        std::cout << "Processing: " << i << " th step" << std::endl;
        float sigma = std::sqrt(2.0 * D * dt);
        sum = 0.;
        injection_e->update_speeds();
        injection_h->update_speeds();
        for(size_t j = 0; j < injection_e->get_charges().size(); ++j)
        {
            
            prev_x_e = injection_e->get_charges().at(j)->get_position().first;
            prev_y_e = injection_e->get_charges().at(j)->get_position().second;
            prev_x_h = injection_h->get_charges().at(j)->get_position().first;
            prev_y_h = injection_h->get_charges().at(j)->get_position().second;
            injection_e->get_charges().at(j)->set_position(prev_x_e + dt*injection_e->get_charges().at(j)->get_velocity().first + sigma*gaussian(generator),
                                                         prev_y_e + dt*injection_e->get_charges().at(j)->get_velocity().second + sigma*gaussian(generator));
            injection_h->get_charges().at(j)->set_position(prev_x_h - dt*injection_h->get_charges().at(j)->get_velocity().first + sigma*gaussian(generator),
                                                         prev_y_h - dt*injection_h->get_charges().at(j)->get_velocity().second + sigma*gaussian(generator));
        }
        for(size_t j = 0; j < injection_e->get_charges().size(); j++)
        {
            sum += (injection_e->get_charges().at(j)->get_velocity().first + injection_e->get_charges().at(j)->get_velocity().second);
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
        graph->GetYaxis()->SetRangeUser(0, 50e-6);
        graph_h->Draw("P SAME");
        c->Modified();
        c->Update();

        // Small delay so the animation is visible (~20–30 FPS)
        gSystem->ProcessEvents();
        gSystem->Sleep(30);
    }

    app.Run();
    return 0;
}