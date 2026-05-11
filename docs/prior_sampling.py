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

priors = {
    'phi_i': (7, 30),
    'alpha': (0, 5),
    'beta' : (0, 5),
}

num_repeats = 100
points = []
statuses = []
for _ in range(num_repeats):
    phi_i = uniform(*priors['phi_i'])
    alpha = uniform(*priors['alpha'])
    beta  = uniform(*priors['beta'])
    points.append((phi_i, alpha, beta))

    print(f"Evaluating at {phi_i = }, {alpha = }, {beta = }")
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
    except Exception as e: statuses.append(1)
    else: statuses.append(0)

phi_is = [point[0] for point in points]
alphas = [point[1] for point in points]
betas  = [point[2] for point in points]
plt.scatter(phi_is, alphas, c=statuses)
plt.xlabel("$\\phi_i$")
plt.ylabel("$\\alpha$")
plt.savefig("prior.pdf")