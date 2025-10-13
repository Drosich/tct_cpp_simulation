#include <iostream>
#include <vector>
#include <random>
#include <filesystem>

#include "detector.hh"
#include "charge_injection.hh"
#include "charge_carrier.hh"
#include "config.hh"
#include "utility.hh"

#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TSystem.h>
#include <TH2F.h>

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

    std::default_random_engine generator;
    std::normal_distribution<float> gaussian(0.0, 1.0);

    TApplication app("", nullptr, nullptr);

    Detector det(cfg.get_Nd(), cfg.get_width(), cfg.get_length(), 
                 cfg.get_V_bi(), cfg.get_V_bias(), cfg.get_R(), 
                 cfg.get_material());

    int steps = cfg.get_steps();
    float dt = cfg.get_dt();
    std::vector<float> t(steps);
    std::vector<float> signal_e(steps, 0.0f);
    std::vector<float> signal_h(steps, 0.0f);
    std::vector<float> signal_total(steps, 0.0f);
    for(int i = 0; i < steps; ++i) t[i] = i * dt;
    float x_lim = (det.get_depleted_width() > det.get_physical_width()) ? det.get_physical_width() : det.get_depleted_width();

    if(cfg.get_sim_type() == "visualization")
    {
        TGraph* graph_e = new TGraph();
        TGraph* graph_h = new TGraph();
        graph_e->SetMarkerStyle(20);
        graph_e->SetMarkerSize(0.5);
        graph_e->SetMarkerColor(kBlue);
        graph_h->SetMarkerStyle(20);
        graph_h->SetMarkerSize(0.5);
        graph_h->SetMarkerColor(kRed);

        TCanvas* c = new TCanvas("c", "Particle Motion", 800, 600);
        gStyle->SetOptStat(0);

        Charge_injection injection_e(cfg.get_focus(),
                                     cfg.get_wavelength(),
                                     cfg.get_NA(),
                                     cfg.get_refractive_index(),
                                     &det,
                                     0,
                                     cfg.get_N());
        Charge_injection injection_h = injection_e;
        injection_h.set_type(1);

        for(int step = 0; step < cfg.get_steps(); ++step)
        {
            std::cout << "Processing: " << step << " th step" << std::endl;

            injection_e.update_speeds();
            injection_h.update_speeds();

            auto& charges_e = injection_e.get_charges();
            auto& charges_h = injection_h.get_charges();

            for(size_t j = 0; j < charges_e.size(); ++j)
            {
                auto vel_e = charges_e[j].get_velocity();
                charges_e[j].set_position(cfg.get_dt()*vel_e.first, cfg.get_dt()*vel_e.second);

                auto vel_h = charges_h[j].get_velocity();
                charges_h[j].set_position(-cfg.get_dt()*vel_h.first, -cfg.get_dt()*vel_h.second);
            }

            float sum_e = 0.;
            float sum_h = 0.;
            for(size_t j = 0; j < charges_e.size(); j++)
            {
                sum_e += charges_e[j].get_velocity().second;
                sum_h += charges_h[j].get_velocity().second;
            }
            signal_e[step] = sum_e * QE / x_lim;
            signal_h[step] = sum_h * QE / x_lim;
            signal_total[step] = (signal_e[step] + signal_h[step]);

            graph_e->Set(0);
            graph_h->Set(0);
            for (size_t j = 0; j < charges_e.size(); ++j) {
                auto pos_e = charges_e[j].get_position();
                auto pos_h = charges_h[j].get_position();
                graph_e->SetPoint(j, pos_e.first, pos_e.second);
                graph_h->SetPoint(j, pos_h.first, pos_h.second);
            }

            c->cd();
            graph_e->Draw("AP");
            graph_e->GetXaxis()->SetLimits(-25e-6, 25e-6);
            graph_e->GetYaxis()->SetRangeUser(0, 50e-6);
            graph_h->Draw("P SAME");
            c->Modified();
            c->Update();
            gSystem->ProcessEvents();
            gSystem->Sleep(30);
        }
        TCanvas* c_pulse = new TCanvas("c_pulse", "pulse", 800, 600);
        c_pulse->cd();
        TGraph* gr_pulse_e = new TGraph(t.size(), t.data(), signal_e.data());
        gr_pulse_e->SetLineColor(kBlue);
        TGraph* gr_pulse_h = new TGraph(t.size(), t.data(), signal_h.data());
        gr_pulse_h->SetLineColor(kRed);
        TGraph* gr_pulse_total = new TGraph(t.size(), t.data(), signal_total.data());
        gr_pulse_total->SetLineColor(kBlack);
        gr_pulse_total->Draw("APL");
        gr_pulse_e->Draw("PL SAME");
        gr_pulse_h->Draw("PL SAME");
    }
    else if(cfg.get_sim_type() == "z_scan")
    {
        std::vector<float> z_array(50);
        for(int i = 0; i < 50; ++i)
            z_array[i] = -20.e-6 + i*90.e-6/50;

        std::vector<float> int_charge_e;
        std::vector<float> int_charge_h;
        std::vector<float> int_charge_t;
        std::vector<float> WPC;

        for(auto z : z_array)
        {
            std::cout << "=== SIMULATING z = " << z/1.e-6 << std::endl;

            Charge_injection injection_e(z,
                                         cfg.get_wavelength(),
                                         cfg.get_NA(),
                                         cfg.get_refractive_index(),
                                         &det,
                                         0,
                                         cfg.get_N());
            Charge_injection injection_h = injection_e;
            injection_h.set_type(1);

            auto &ce = injection_e.get_charges();
            auto &ch = injection_h.get_charges();
            size_t inside_e = 0, inside_h = 0;
            float y_min = 0.0f;
            float y_max = std::min(det.get_depleted_width(), det.get_physical_width());
            for (auto &p : ce) if (p.get_position().second >= y_min && p.get_position().second <= y_max) ++inside_e;
            for (auto &p : ch) if (p.get_position().second >= y_min && p.get_position().second <= y_max) ++inside_h;
            std::cerr << "z="<<z<<" injected_e inside="<<inside_e<<"/"<<ce.size()
                    <<" injected_h inside="<<inside_h<<"/"<<ch.size()
                    <<" y_range=["<<y_min<<","<<y_max<<"]\n";
            
            for(int step = 0; step < cfg.get_steps(); ++step)
            {
                injection_e.update_speeds();
                injection_h.update_speeds();

                auto& charges_e = injection_e.get_charges();
                auto& charges_h = injection_h.get_charges();

                for(size_t j = 0; j < charges_e.size(); ++j)
                {
                    auto vel_e = charges_e[j].get_velocity();
                    charges_e[j].set_position(cfg.get_dt()*vel_e.first, cfg.get_dt()*vel_e.second);

                    auto vel_h = charges_h[j].get_velocity();
                    charges_h[j].set_position(-cfg.get_dt()*vel_h.first, -cfg.get_dt()*vel_h.second);
                }

                float sum_e = 0.;
                float sum_h = 0.;
                for(size_t j = 0; j < charges_e.size(); j++)
                {
                    sum_e += charges_e[j].get_velocity().second;
                    sum_h += charges_h[j].get_velocity().second;
                }

                signal_e[step] = sum_e * QE / x_lim;
                signal_h[step] = sum_h * QE / x_lim;
                signal_total[step] = (signal_e[step] + signal_h[step]);
            }

            TGraph graph_e(t.size(), t.data(), signal_e.data());
            TGraph graph_h(t.size(), t.data(), signal_h.data());
            TGraph graph_t(t.size(), t.data(), signal_total.data());
            float Q_e = graph_e.Integral();
            float Q_h = graph_h.Integral();
            float Q_t = graph_t.Integral();
            float WPC_val = linear_interpolation(cfg.get_t_pc(), t, signal_total);
            int_charge_e.push_back(Q_e);
            int_charge_h.push_back(Q_h);
            int_charge_t.push_back(Q_t);
            WPC.push_back(WPC_val);
        }

        std::cout << "END" << std::endl;

        TCanvas* c = new TCanvas("c", "Z-Scan", 800, 600);
        c->cd();
        TGraph* z_scan_e = new TGraph(int_charge_e.size(), z_array.data(), int_charge_e.data());
        z_scan_e->SetLineColor(kBlue);
        TGraph* z_scan_h = new TGraph(int_charge_h.size(), z_array.data(), int_charge_h.data());
        z_scan_h->SetLineColor(kRed);
        TGraph* z_scan_t = new TGraph(int_charge_t.size(), z_array.data(), int_charge_t.data());
        z_scan_t->SetLineColor(kBlack);
        z_scan_e->SetTitle("z-scan;z [um];Charge [a.u.]");
        z_scan_t->Draw("APL");
        z_scan_e->Draw("PL SAME");
        z_scan_h->Draw("PL SAME");
        c->Update();

        TCanvas* c2 = new TCanvas("c2", "WPC", 800, 600);
        c2->cd();
        TGraph* z_scan_WPC = new TGraph(WPC.size(), z_array.data(), WPC.data());
        z_scan_WPC->SetTitle("WPC;z [um];WPC [a.u.]");
        z_scan_WPC->Draw("APL");
        c2->Update();
    }
    else
    {
        std::cout << "Unrecognised sim mode. Exiting" << std::endl;
    }

    app.Run();
    return 0;
}
