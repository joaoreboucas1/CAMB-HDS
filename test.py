import camb
h = 0.6756
omegam = 0.27
omegabh2 = 0.02238280
omegab = omegabh2/h**2
omegac = omegam - omegab
omegach2 = omegac*h**2

cosmo = camb.set_params(# Background
	H0=100*h, ombh2=omegabh2, omch2=omegach2, TCMB=2.7255,
	# Dark Energy
	dark_energy_model = 'HybridQuintessence', phi_i = 10,
	# Neutrinos
	omnuh2=0, num_nu_massless=3.044, num_nu_massive = 0,
	nu_mass_degeneracies=[0], nu_mass_numbers = [0],
	# Initial Power Spectrum
	As = 2.100549e-09, ns = 0.9660499, 
	YHe = 0.246, WantTransfer=True
)
results = camb.get_results(cosmo)