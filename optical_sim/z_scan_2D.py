from diffractio import degrees, mm, plt, sp, um, np
from diffractio.scalar_masks_X import Scalar_mask_X
from diffractio.scalar_masks_XZ import Scalar_mask_XZ
from diffractio.scalar_sources_X import Scalar_source_X

import time
from scipy.integrate import trapezoid
import ROOT

N = 256
x0 = np.linspace(-35 * um, 35 * um, N)
z0 = np.linspace(-150 * um, 150 * um, N)
#[0]/(1+((x-[2])*(x-[2]))/([1]*[1]))
wavelength = .4 * um

charges = []
z_s = np.linspace(60*um, -100.*um, 100, dtype=np.float64)

graphs = []
NAs = [0.1,0.2,0.3,0.4,0.5]

for NA in NAs:
    w0 = wavelength/(np.pi*NA)
    u0 = Scalar_source_X(x=x0, wavelength=wavelength)
    u0.gauss_beam(A=1, x0=0 * um, z0 = 150 * um, w0=w0, theta=0 * degrees)
    for i in z_s:
        u1 = Scalar_mask_XZ(x=x0, z=z0, wavelength=wavelength, n_background=1)
        u1.incident_field(u0)

        u1.semi_plane(r0=(0, i), refractive_index=1, angle=0 * degrees)#2.7592
        u1.WPM(verbose=True)

        tpa_charge = np.rot90((u1.intensity())**2, k=1)
        z_min = np.argmin(np.abs((z0 - i)))
        z_max = np.argmin(np.abs(z0 - (i+50*um)))
        
        print(f"Integrating from ({i}){z_min} to ({i+50*um}){z_max}")
        print(f"Charge = {np.mean(tpa_charge[(1024-200):(1024+200),z_min:z_max])}")

        charges.append(np.mean(tpa_charge[:,z_min:z_max]))
        # # --- Visualization ---
        # fig, ax = plt.subplots(figsize=(6, 5))
        # im = ax.imshow(
        #     tpa_charge,
        #     extent=[z0[0] / um, z0[-1] / um, x0[0] / um, x0[-1] / um],
        #     origin='lower',
        #     aspect='auto',
        #     cmap='inferno'
        # )
        # plt.colorbar(im, ax=ax, label='(Intensity)^2')

        # # Draw integration window
        # ax.axvline(z0[z_min] / um, color='cyan', linestyle='--', label='Integration Start')
        # ax.axvline(z0[z_max] / um, color='lime', linestyle='--', label='Integration End')

        # ax.set_xlabel('z [µm]')
        # ax.set_ylabel('x [µm]')
        # ax.set_title(f"Integration window for i = {i/um:.1f} µm")
        # ax.legend(loc='upper right')
        # plt.tight_layout()
        # plt.show()
        u1.clear_field()
    graphs.append(ROOT.TGraph(len(z_s), z_s, np.array(charges, dtype=np.float64)/np.max(charges)))
    charges.clear()

c1 = ROOT.TCanvas()
for id, gr in enumerate(graphs):
    gr.SetLineWidth(2)
    gr.SetTitle(f"Simulation (NA = {NAs[id]})")
    gr.GetXaxis().SetTitle("z [#mum]")
    gr.GetYaxis().SetTitle("Normalized Charge")
    if id == 0:
        gr.Draw("APL")
    else:
        gr.Draw("PL")

# # u1.draw(kind="phase", scale="scaled", draw_borders=True, percentage_intensity=0.01)
#     u1.draw(kind="intensity", logarithm=1, scale="scaled", draw_borders=True)
#     plt.show()