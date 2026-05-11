# Evaluating at phi_i = 12.189121325035616, alpha = 4.036521593124252, beta = 0.1152552290104486

import faulthandler
import sys
from random import uniform
import matplotlib.pyplot as plt
import camb

faulthandler.enable(file=sys.stderr, all_threads=True, c_stack=True)

# Standard cosmological parameters
h = 0.6756
omegabh2 = 0.022
omegach2 = 0.12
As = 2.1e-9
ns = 0.962
phi_i = 12.189121325035616
alpha = 4.036521593124252
beta = 0.1152552290104486


cosmo = camb.set_params(
    # Background
    H0=100*h, ombh2=omegabh2, omch2=omegach2, TCMB=2.7255,
    # Dark Energy
    dark_energy_model='HybridQuintessence', potential_type=3, phi_i=phi_i, alpha=alpha, beta=beta,
    # Neutrinos
    omnuh2=0.06/93.14, nnu=3.044,
    nu_mass_degeneracies=[0], nu_mass_numbers=[0], share_delta_neff=True,
    # Initial Power Spectrum
    As=As, ns=ns,
    tau=0.0544,
    YHe=0.246, WantTransfer=True
)
cosmo.NonLinear = camb.model.NonLinear_none
try: result = camb.get_results(cosmo)
except Exception as e: print(f"Test raised Python exception: {e}")
else: print("Test ran successfully")