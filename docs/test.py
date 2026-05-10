import faulthandler
import sys
import camb

faulthandler.enable(file=sys.stderr, all_threads=True, c_stack=True)

# Standard cosmological parameters
h = 0.6756
omegabh2 = 0.022
omegach2 = 0.12
As = 2.1e-9
ns = 0.962
phi_i = 23

cosmo = camb.set_params(
    # Background
    H0=100*h, ombh2=omegabh2, omch2=omegach2, TCMB=2.7255,
    # Dark Energy
    dark_energy_model='HybridQuintessence', phi_i=phi_i, alpha=1.0, potential_type=3, beta=1.68,
    # Neutrinos
    omnuh2=0.06/93.14,# nnu=3.044,
    nu_mass_degeneracies=[0], nu_mass_numbers=[0], share_delta_neff=True,
    # Initial Power Spectrum
    As=As, ns=ns,
    tau=0.0544,
    YHe=0.246, WantTransfer=True
)
cosmo.NonLinear = camb.model.NonLinear_none
print(cosmo)
result = camb.get_results(cosmo)