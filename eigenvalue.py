import dedalus.public as de2
from eigentools import Eigenproblem, CriticalFinder
import numpy as np
import argparse
import matplotlib.pyplot as plt
from mpi4py import MPI
import time

import logging

logger = logging.getLogger(__name__.split(".")[-1])

comm = MPI.COMM_WORLD

Nz = 128
z_basis = de2.Chebyshev("z", Nz, interval=(0, 1))
d = de2.Domain([z_basis], np.complex128, comm=MPI.COMM_SELF)
z = z_basis.grid()

pert = de2.EVP(d, ["u", "v", "w", "P", "T", "uz", "vz", "wz", "Tz"], eigenvalue="omega")
pert.parameters["Rf"] = 1e5
pert.parameters["Pr"] = 1
pert.parameters["Ta"] = 1e4**0.5
theta = np.radians(5)
pert.parameters["cos_phi"] = np.cos(theta)
pert.parameters["sin_phi"] = np.sin(theta)
pert.parameters["k"] = 2  # horizontal wavenumber
pert.substitutions["dt(A)"] = "omega*A"
pert.substitutions["dy(A)"] = "1j*k*A"

# Equilibrium Temperature
l = 0.1
a = 1 / (l * (1 - np.exp(-1 / l)))
beta = 1.0
C = a * l * l + 1  # Arbitrary constant, C=1+al^2 makes T_eq(z=0)=1


T_eq = d.new_field(name="T_eq")
T_eq["g"] = -a * l * l * np.exp(-z / l) + beta * z * z * 0.5 - a * l * z + C
pert.parameters["T_eq"] = T_eq

T_eq_z = d.new_field(
    name="T_eq_z",
)
T_eq_z["g"] = T_eq.differentiate("z")["g"]
pert.parameters["T_eq_z"] = T_eq_z

# Mass Conservation
pert.add_equation("dy(v) + wz = 0")

# x-component of Navier Stokes equation
pert.add_equation("dt(u) - (dy(dy(u)) + dz(uz)) + Ta*(w*sin_phi - v*cos_phi) = 0")

# y-component of Navier Stokes equation
pert.add_equation("dt(v) - (dy(dy(v)) + dz(vz)) + dy(P) + Ta*(u*cos_phi) = 0")

# z-component of Navier Stokes equation
pert.add_equation(
    "dt(w) - (dy(dy(w)) + dz(wz)) + dz(P) - Ta*(u*sin_phi) - (Rf/Pr) * T = 0"
)

# Temperature equation
pert.add_equation("dt(T) + (v*dy(T_eq) - w*T_eq_z) - (1/Pr) * (dy(dy(T)) + dz(Tz)) = 0")

# dz substitutions
pert.add_equation("uz - dz(u) = 0")
pert.add_equation("vz - dz(v) = 0")
pert.add_equation("wz - dz(w) = 0")
pert.add_equation("Tz - dz(T) = 0")

# Boundary conditions
# Stress Free
pert.add_bc("left(uz) = 0")
pert.add_bc("right(uz) = 0")
pert.add_bc("left(vz) = 0")
pert.add_bc("right(vz) = 0")
# Impermeable
pert.add_bc("left(w) = 0")
pert.add_bc("right(w) = 0")
# Insulating
pert.add_bc("left(Tz) = 0")
pert.add_bc("right(Tz) = 0")

EP = Eigenproblem(pert)

cf = CriticalFinder(EP, ("k", "Rf"), comm, find_freq=True)

if comm.rank == 0:
    print("\n### Generating Grid ###\n")
start = time.time()

direc = "grid"

nRa, nk = 50, 50
Rapoints = np.linspace(500, 3000, nRa)
kpoints = np.linspace(2, 10, nk)

try:
    cf.load_grid(f"{direc}.h5")
except:
    cf.grid_generator((kpoints, Rapoints), sparse=True)
    if comm.rank == 0:
        cf.save_grid(direc)


end = time.time()
if comm.rank == 0:
    print("### Grid Generation Complete ###")
    print(f"Time taken: {(end-start)/60:10.5f} m ({end-start:10.5f} s) \n")

logger.info("Beginning critical finding with root polishing...")
begin = time.time()
crit = cf.crit_finder(polish_roots=True, tol=1e-5)
end = time.time()
logger.info("critical finding/root polishing time: {:10.5f} sec".format(end - start))
