from diffractio import degrees, eps, mm, no_date, np, um
from diffractio.scalar_fields_XYZ import Scalar_field_XYZ
from diffractio.scalar_masks_XY import Scalar_mask_XY
from diffractio.scalar_masks_XYZ import Scalar_mask_XYZ
from diffractio.scalar_sources_XY import Scalar_source_XY
from diffractio.vector_fields_XYZ import Vector_field_XYZ

from py_pol.jones_vector import Jones_vector  # type: ignore

import ROOT

# Spatial grid
N = 360
x0 = np.linspace(-25 * um, 25 * um, N)
y0 = np.linspace(-25 * um, 25 * um, N)
z0 = np.linspace(-100 * um, 100 * um, N)

# Vector source
wavelength = .4 * um
Exyz = Vector_field_XYZ(x=x0, y=y0, z=z0, wavelength=wavelength)
#Polarization
j0 = Jones_vector()
j0.general_azimuth_ellipticity(azimuth=0 * degrees, ellipticity=0 * degrees)

charges = [];focus_x = [];focus_y = [];focus_z = [];maxes = []
z_s = np.linspace(60*um, -100.*um, 100, dtype=np.float64)

#Initial light source
NA = 0.1
w0 = wavelength/(np.pi*NA)
u0 = Scalar_source_XY(x=x0, y=y0, wavelength=wavelength)
u0.gauss_beam(A=1, r0=(0, 0), z0 = 100 * um, w0=(w0,w0), theta=0 * degrees)

# Incident field and polarization
Exyz.incident_field(u0=u0, j0=j0)

size=1000*um
u1 = Scalar_mask_XYZ(x=x0, y=y0, z=z0, wavelength=wavelength, n_background=1.)

for i in z_s:
    Exyz.incident_field(u0=u0, j0=j0)
    # Surface
    u1.cube(r0=(0*um, 0*um, i+size/2), size=(size,size,size), refractive_index=2.759)
    Exyz.refractive_index_from_scalarXYZ(u1)

    # Propagation
    Exyz.FP_WPM(has_edges=True, verbose=True)

    # Data
    tpa = np.array(Exyz.get("intensity", is_matrix=True))**2
    np.save(f"gaussianbeam_z={i}.npy", tpa)
    ix, iy, iz = np.unravel_index(tpa.argmax(), tpa.shape)
    focus_x.append(x0[ix])
    focus_y.append(y0[iy])
    focus_z.append(z0[iz])
    maxes.append(np.max(tpa))
    z_min = np.argmin(np.abs(z0 - i))
    z_max = np.argmin(np.abs(z0 - (i + 50.)))
    charges.append(np.sum(tpa[:,:,z_min:z_max]))
    
    Exyz.clear_field()
