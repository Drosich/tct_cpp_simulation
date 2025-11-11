from diffractio import degrees, eps, mm, no_date, np, um #type: ignore
from diffractio.scalar_fields_XYZ import Scalar_field_XYZ #type: ignore
from diffractio.scalar_masks_XY import Scalar_mask_XY #type: ignore 
from diffractio.scalar_masks_XYZ import Scalar_mask_XYZ #type: ignore
from diffractio.scalar_sources_XY import Scalar_source_XY #type: ignore

import ROOT

N = 512
x0 = np.linspace(-25 * um, 25 * um, N)
y0 = np.linspace(-25 * um, 25 * um, N)
z0 = np.linspace(-100 * um, 100 * um, N)

#FIT: [0]*atan(([1] - (x-[2]))/[3]) + [0]*atan(([1] + (x-[2]))/[3]) + [4]

wavelength = .4 * um

z_s = np.linspace(60*um, -100.*um, 100, dtype=np.float64)
size = 1000*um
u1 = Scalar_mask_XYZ(x=x0, y=y0, z=z0, wavelength=wavelength, n_background=1.)

NA = 0.328
w0 = wavelength/(np.pi*NA)
charges = []

u0 = Scalar_source_XY(x=x0, y=y0, wavelength=wavelength)
u0.gauss_beam(A=1, r0=(0, 0), z0 = 100 * um, w0=(w0,w0), theta=0 * degrees)
size=1000*um

for id, i in enumerate(z_s):
    print(f"Iteration {id}/{len(z_s)}")
    u1.incident_field(u0)
    u1.cube(r0=(0*um, 0*um, i + size/2), size=(size, size, size), refractive_index=2.759)
    u1.WPM(verbose=True)

    # Get the intensity array in float32 to cut memory in half
    tpa = u1.intensity().astype(np.float32)
    tpa **= 2  # in-place squaring

    # Compute charges
    z_min = np.argmin(np.abs(z0 - i))
    z_max = np.argmin(np.abs(z0 - (i + 50.)))
    charges.append(np.sum(tpa[:, :, z_min:z_max], dtype=np.float64))

    # Free memory explicitly
    del tpa
    u1.clear_field()