# Equation Derivation for viscous time

The non-dimensional equations for convection using $\tau_\nu$ are:

$\bold{\nabla} \cdot \bold{u} = 0$,

$\partial_t \bold{u} + (\bold{u}\cdot\bold{\nabla})\bold{u} = -\nabla P + \frac{R_f}{Pr}T\hat{\bold{z}} + \bold{\nabla}^2\bold{u} - Ta^\frac{1}{2} \bold{\Omega}\times\bold{u}$,

$\partial_t T + (\bold{u}\cdot\bold{\nabla})T = \frac{1}{Pr}\left[\nabla^2 T + \mathcal{H}\right]$,

## Applying perturbation theory

From here we can decompose to an equilibrium state and a small order perturbation:

$\bold{u} = \bold{u}_\text{eq} + \epsilon \bold{u}'$, $\quad T = T_\text{eq} + \epsilon T'$, $\quad P = P_\text{eq} + \epsilon P'$,

where $\epsilon$ is a small parameter. We can substitute these into the non-dimensionalised equations, gather terms of the same order in epsilon, ignore terms of order $\epsilon^2$ and higher, and remember that all $\epsilon^0$ terms are a steady-state equilibrium, so $\partial_t$ terms = 0 and $\bold{u}_\text{eq} = 0$. With this, we get:

### $\epsilon^0$ (equilibrium equations):

$\bold{\nabla} \cdot \bold{u}_\text{eq} = 0$,

$\nabla P_\text{eq} = \frac{R_f}{Pr}\ T_\text{eq}\hat{\bold{z}}$,

$-\nabla^2 T_\text{eq} = \mathcal{H}$

### $\epsilon^1$ (perturbation equations):

$\bold{\nabla}\cdot\bold{u}' = 0$,

$\partial_t \bold{u}' + \nabla P' - \frac{R_f}{Pr} T' \hat{\bold{z}} - \bold{\nabla}^2 \bold{u}' + Ta^\frac{1}{2} \bold{\Omega}\times\bold{u}' = 0$,

$\partial_t T' + (\bold{u'} \cdot \bold{\nabla})T_\text{eq} - \frac{1}{Pr} \nabla^2 T' = 0$

In the case of the Kazemi+22 heating function, we can solve the third $\epsilon^0$ equation to calculate $T_\text{eq}$. Since

$\mathcal{H} = a e^\frac{-z}{\ell} - \beta$,

we can show that

$\partial_z T_\text{eq} = a\ell e^\frac{-z}{\ell} + \beta z + C$,

and we can use the boundary condition that $\partial_z T_\text{eq}|_{z=0} = 0$ to show that $C = -a\ell$. We can then find

$T_\text{eq} = -a\ell^2e^\frac{-z}{\ell} + \frac{\beta z^2}{2} - a\ell z + C$,

where $C$ is an arbitrary constant. It can be shown that choosing $C = 1 + a\ell^2$ will set $T_\text{eq}|_{z=0} = 1$, so we can deduce our equilibrium temperature profile as

$T_\text{eq} = -a\ell^2e^\frac{-z}{\ell} + \frac{\beta z^2}{2} - a\ell z + 1 + a \ell^2$.

## Boundary Conditions

The boundary conditions for the perturbation equations are:

Impermeable top and bottom: $\\w'(z=0) = w'(z=L) = 0$,

Free-Slip top and bottom: $\\\partial_z u'(z=0) = \partial_z u'(z=L) = 0$; $\\ \partial_z v'(z=0) = \partial_z v'(z=L) = 0$,

Insulating top and bottom: $\\\partial_z T'(z=0) = \partial_z T'(z=L) = 0$.

## Solving the equations

The perturbation equations and boundary conditions are inputted into eigenvalue.py and solved eigentools, which should recover the critical Rayleigh number and critical wave number for a given Taylor number.
