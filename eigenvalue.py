from eigentools import Eigenproblem, CriticalFinder
import dedalus.public as de2
import numpy as np
import argparse
import matplotlib.pyplot as plt
from mpi4py import MPI
import time
from scipy import interpolate

import logging

logger = logging.getLogger(__name__.split(".")[-1])

comm = MPI.COMM_WORLD


parser = argparse.ArgumentParser()
parser.add_argument("--Ra", type=float, help="Rayleigh number")

args = parser.parse_args()
Ra = 1e3 if args.Ra is None else args.Ra
# ? Why do we set Ra if we want to sweep through Ra for the root?
# ? Matt's example sets Ra=1000 but sweeps through Ra=200-700
Ly = 4
n = 1
ky_restrict = []
Lz = 1
Pr = 1
# ? Relatively low Ta to begin with, Ra_crit should be same OoM as non-rotating
Ta = 1e4
phi = np.radians(5)
Ny = 128
Nz = 256

# Create parameter space
nRa, nk = 50, 50
Rapoints = np.linspace(500, 3000, nRa)
kpoints = np.linspace(2, 10, nk)
# print(kpoints)

zbasis = de2.Chebyshev("z", Nz, interval=(0, Lz), dealias=3 / 2)
domain = de2.Domain([zbasis], comm=MPI.COMM_SELF)
z = zbasis.grid()

string = de2.EVP(
    domain,
    ["u", "v", "w", "P", "T", "uz", "vz", "wz", "Tz"],
    eigenvalue="omega",
    tolerance=1e-10,
)
string.parameters["Ra"] = Ra
string.parameters["Pr"] = Pr
string.parameters["Ta"] = np.sqrt(Ta)
string.parameters["cos_phi"] = np.cos(phi)
string.parameters["sin_phi"] = np.sin(phi)
string.parameters["k"] = 5  # horizontal wavenumber #? Also why set k to 1 if sweeping?
string.substitutions["dt(A)"] = "omega*A"
string.substitutions["dy(A)"] = "1j*k*A"

# non-constant co-efficients
# internal heating function
heat = domain.new_field(name="heat")
#! Heat Function from Kazemi et al. 2022
l = 0.1
beta = 1.0
a = 1 / (l * (1 - np.exp(-1 / l)))
heat_func = a * (np.exp(-z / l)) - beta
heat["g"] = a * (np.exp(-z / l)) - beta
string.parameters["heat"] = heat
# print("max heat value = ", np.max(heat["g"]))

string.add_equation("uz - dz(u) = 0")
string.add_equation("vz - dz(v) = 0")
string.add_equation("wz - dz(w) = 0")
string.add_equation("Tz - dz(T) = 0")

#! Mass continuity equation
string.add_equation("dy(v) + wz = 0")

#! 2.5D so all dx terms have been set = 0, but u and uz still exist.
#! x-component of Navier Stokes equation
string.add_equation(
    "dt(u) - (dy(dy(u)) + dz(uz))           + Ta*(w*cos_phi - v*sin_phi)       = -(v*dy(u) + w*uz)"
)

#! y-component of Navier Stokes equation
string.add_equation(
    "dt(v) - (dy(dy(v)) + dz(vz))  + dy(P)  + Ta*u*sin_phi                     = -(v*dy(v) + w*vz)"
)

#! z-component of Navier Stokes equation
string.add_equation(
    "dt(w) - (dy(dy(w)) + dz(wz))  + dz(P)  - Ta*u*cos_phi    - (Ra/Pr)*T      = -(v*dy(w) + w*wz)"
)
#! Temperature evolution equation
string.add_equation(
    "dt(T) - ( (dy(dy(T)) + dz(Tz) )) * (1/Pr)                               = -(v*dy(T) + w*Tz)"
)
# ? The heating term should be included in here - ideally in - ( (dy(dy(T)) + dz(Tz) ) + heat ) * (1/Pr) = ...
# ? - if I add it there, I get error that cannot add dependent and independent terms. From
# ? https://groups.google.com/g/dedalus-users/c/PzXIfJNj8b8/m/EdSWv1vmFAAJ, it seems that
# ? maybe constant terms should go on the RHS of the equation (although this isn't true
# ? of the NCC terms in Matt's ROT_EVP_newsyntax.py code), but adding it on the RHS
# ? gives me an error that RHS is no longer homogeneous, as the heat term is not
# ? equal to 0. I am unsure how to add the heating term into the equation.

#! Stress-Free horizontal boundaries
#! d(u)/dz = 0 at top and bottom
string.add_bc("left(uz) = 0")
string.add_bc("right(uz) = 0")

#! d(v)/dz = 0 at top and bottom
string.add_bc("left(vz) = 0")
string.add_bc("right(vz) = 0")

#! Impermeable top and bottom boundaries
#! w = 0 at top and bottom
string.add_bc("left(w) = 0")
string.add_bc("right(w) = 0")

#! Insulating horizontal boundaries
string.add_bc("left(Tz) = 0")
string.add_bc("right(Tz) = 0")


EP = Eigenproblem(string, factor=10)

cf = CriticalFinder(EP, ("k", "Ra"), comm, find_freq=True)

if comm.rank == 0:
    print("\n### Generating Grid ###\n")
start = time.time()

direc = "grid"

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
exit(0)
#! Will look at this bit once the grid is generated correctly. At the moment, am
#! plotting the grid with plot_grid.ipynb


logger.info("Beginning critical finding with root polishing...")
begin = time.time()
crit = cf.crit_finder(polish_roots=True, tol=1e-5)
end = time.time()
logger.info("critical finding/root polishing time: {:10.5f} sec".format(end - start))

if comm.rank == 0:
    print("crit = {}".format(crit))
    print("critical wavenumber k = {:10.5f}".format(crit[0]))
    print("critical Ra = {:10.5f}".format(crit[1]))
    print("critical freq = {:10.5f}".format(crit[2]))

    mask = np.isfinite(cf.roots)
    yy_root = cf.parameter_grids[0][0, mask]
    rroot = cf.roots[mask]
    root_fn = interpolate.interp1d(
        yy_root, rroot, kind="cubic"
    )  # Interpolating over masked roots
    y_fn = np.linspace(
        yy_root[0], yy_root[-1], 100
    )  # Constructing a kx array for plotting
    #
    print(yy_root, rroot, n * 2 * np.pi / Ly, yy_root[-1])
    # Restrict values of kx to those which fit in our box shape Lx
    while n * 2 * np.pi / Ly < yy_root[-1]:
        if n * 2 * np.pi / Ly < yy_root[0]:
            pass
        else:
            ky_restrict.append(n * 2 * np.pi / Ly)
        n += 1
    print(f"Restricted ky values for Ly = {Ly}:")
    print(ky_restrict)
    Ra_restrict = root_fn(ky_restrict)

    index = np.where(Ra_restrict == np.min(Ra_restrict))[0][0]
    #
    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(111)
    xx = cf.parameter_grids[0]
    yy = cf.parameter_grids[1]
    grid = cf.evalue_grid.real
    biggest_val = 2 * np.abs(grid).std()
    plt.pcolormesh(xx, yy, grid, cmap="RdBu_r", vmin=-biggest_val, vmax=biggest_val)
    plt.colorbar()
    x = cf.parameter_grids[0][0, :]
    y = cf.roots[:]
    plt.scatter(x, y)
    plt.plot(y_fn, root_fn(y_fn), ":")
    plt.ylim(yy.min(), yy.max())
    plt.xlim(xx.min(), xx.max())
    plt.xlabel("ky")
    plt.ylabel("Ra")

    plt.title(
        "All ky: Ra_c = {:.2f}, ky_c = {:.3f} \n Restricted ky: Lx = {}, Ra_c = {:.2f}, ky_c = {:.3f} \n ".format(
            crit[1], crit[0], Ly, Ra_restrict[index], ky_restrict[index]
        )
    )
    plt.tight_layout()
