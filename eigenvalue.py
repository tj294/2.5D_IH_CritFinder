import dedalus.public as de2
from eigentools import Eigenproblem, CriticalFinder
import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib import transforms
from mpi4py import MPI
import time
from scipy import interpolate, optimize
import argparse
import os
import logging

logger = logging.getLogger(__name__.split(".")[-1])


def plot_grid(cf, Ta):
    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(111)
    yy = cf.parameter_grids[0]
    zz = cf.parameter_grids[1]
    grid = cf.evalue_grid.real
    biggest_val = 2 * np.abs(grid).std()
    plt.pcolormesh(yy, zz, grid, cmap="RdBu_r", vmin=-biggest_val, vmax=biggest_val)
    plt.colorbar()
    plt.xlabel("ky")
    plt.ylabel("Ra")
    plt.savefig(f"Ta{args.Ta}_grid.png")


parser = argparse.ArgumentParser(description="Critical Finder")
parser.add_argument(
    "--Ra", type=str, nargs=2, default=["1e2", "1e4"], help="Limits of Ra Sweep"
)
parser.add_argument(
    "--k", type=str, nargs=2, default=["2", "10"], help="Limits of k Sweep"
)
parser.add_argument("--n", type=str, default="50", help="Resolution of Sweep")
parser.add_argument("--Ta", type=str, default="1e4", help="Taylor number to simulate")
parser.add_argument(
    "--overwrite", "-o", action="store_true", help="Overwrite existing grid"
)

args = parser.parse_args()
# print(args)
comm = MPI.COMM_WORLD

Nz = 256
z_basis = de2.Chebyshev("z", Nz, interval=(0, 1))
d = de2.Domain([z_basis], np.complex128, comm=MPI.COMM_SELF)
z = z_basis.grid()

pert = de2.EVP(d, ["u", "v", "w", "P", "T", "uz", "vz", "wz", "Tz"], eigenvalue="omega")
pert.parameters["Rf"] = 1e7
pert.parameters["Pr"] = 1
pert.parameters["Ta"] = float(args.Ta) ** 0.5
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

direc = f"Ta{args.Ta}_grid"

nRa = nk = int(args.n)
Rapoints = np.linspace(float(args.Ra[0]), float(args.Ra[1]), nRa)
kpoints = np.linspace(float(args.k[0]), float(args.k[1]), nk)

if args.overwrite:
    if os.path.isfile(f"{direc}.h5"):
        os.remove(f"{direc}.h5")

try:
    cf.load_grid(f"{direc}.h5")
except:
    cf.grid_generator((kpoints, Rapoints), sparse=True)
    if comm.rank == 0:
        cf.save_grid(direc)

plot_grid(cf, args.Ta)
# cf.plot_crit(xlabel="k", ylabel="Ra", cmap="RdBu_r")

end = time.time()
if comm.rank == 0:
    print("### Grid Generation Complete ###")
    print(f"Time taken: {(end-start)/60:10.5f} m ({end-start:10.5f} s) \n")

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
    n = 128
    Ly = 4
    ky_restrict = []

    # print(yy_root, rroot, n * 2 * np.pi / Ly, yy_root[-1])
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

    # index = np.where(Ra_restrict == np.min(Ra_restrict))[0][0]
    index = 0
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

    # plt.title(
    #     "All ky: Ra_c = {:.2f}, ky_c = {:.3f} \n Restricted ky: Lx = {}, Ra_c = {:.2f}, ky_c = {:.3f} \n ".format(
    #         crit[1], crit[0], Ly, Ra_restrict[index], ky_restrict[index]
    #     )
    # )
    plt.title(
        "All ky: Ra_c = {:.2f}, ky_c = {:.3f} \n ".format(
            crit[1],
            crit[0],
        )
    )
    plt.tight_layout()
    plt.savefig("crit_finder.png")
