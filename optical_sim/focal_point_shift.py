from diffractio import np, um, degrees, plt
from diffractio.scalar_sources_X import Scalar_source_X
from diffractio.scalar_fields_XZ import Scalar_field_XZ
from diffractio.scalar_masks_XZ import Scalar_mask_XZ

x0 = np.linspace(-5 * um, 5 * um, 2048)
z0 = np.linspace(-0*um, 300 * um, 2048)
max_pos = []
max_pos_noMask = []
wavelength = .8 * um
zR = 3.66 * um
w0 = np.sqrt(wavelength*zR/np.pi)


for z in np.linspace(0*um, 50*um, 50):
    u0 = Scalar_source_X(x=x0, wavelength=wavelength)
    u0.gauss_beam(A=1, x0=0*um, z0=z, w0=w0, theta=0*degrees)

    u1 = Scalar_mask_XZ(x=x0, z=z0, wavelength=wavelength)
    u1.incident_field(u0, z0=150*um)
    u1.layer(r0=(0*um, 150*um), depth=300*um, refractive_index=2.6)

    u1.WPM(verbose=True)
    # u1.draw(draw_borders = True)

    u0_noMask = Scalar_source_X(x=x0, wavelength=wavelength)
    u0_noMask.gauss_beam(A=1, x0=0*um, z0=z, w0=w0, theta=0*degrees)

    u1_noMask = Scalar_mask_XZ(x=x0, z=z0, wavelength=wavelength)
    u1_noMask.incident_field(u0_noMask, z0=150*um)
    

    u1_noMask.WPM(verbose=True)

    intensity = np.abs(u1.u)**2
    profile = np.mean([intensity[:, 1023], intensity[:, 1024],intensity[:, 1025]], axis=0)
    profile /= np.max(profile)
    max_index = np.argmax(profile)


    intensity = np.abs(u1_noMask.u)**2
    profile = np.mean([intensity[:, 1023], intensity[:, 1024],intensity[:, 1025]], axis=0)
    profile /= np.max(profile)

    max_index_noMask = np.argmax(profile)
    max_pos.append(z0[max_index])
    max_pos_noMask.append(z0[max_index_noMask])

plt.plot(np.array(max_pos_noMask)-150*um, abs(np.array(max_pos)-np.array(max_pos_noMask)))
plt.show()