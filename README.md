# Equation Derivation for Freefall Time

Free-Fall time is the crossing time if buoyancy acts unimpeded, so the advection term will balance with the buoyancy term.

$(\bold{u}\cdot\bold{\nabla})\bold{u} \sim \alpha g_0 T$

$\therefore U^2/L \sim \alpha g T$ where $U$, $L$, $T$ are the velocity, length and temperature scales respectively.

From this, we can see that the free-fall velocity,

$U_{ff} \sim \sqrt{\alpha g T L}$,

and so the free-fall time,

$\tau_{ff} \sim \frac{L}{U_{ff}} \sim \frac{L}{\sqrt{\alpha g T L}} \sim \sqrt{\frac{L}{\alpha g T}}$.

**BUT:** In my case the temperature scale is more complicated than the traditional $T\sim \Delta T$, due to the heating function $\mathcal{H}$. Following Kazemi+22, as $\mathcal{H} \sim Q$ with units temperature / time, the temperature scale becomes $T \sim \frac{L^2\mathcal{H}}{\kappa}$, which makes $u_{ff} \sim \sqrt{\frac{\alpha g L^3 Q}{\kappa}}$ and $\tau_{ff} \sim \sqrt{\frac{\kappa}{\alpha g L Q}}$.

From the dimensional equations:

$\bold{\nabla} \cdot \bold{u} = 0$,

$\partial_t\bold{u} + (\bold{u} \cdot \bold{\nabla})\bold{u} = -\frac{1}{\rho_0}\nabla P + \alpha g T \hat{\bold{z}} + \nu \bold{\nabla}^2\bold{u} - 2\bold{\Omega}\times\bold{u}$

$\partial_t T + (\bold{u}\cdot \bold{\nabla})T = \kappa \nabla^2 T + \mathcal{H}$

Using the scalings (where $\hat\cdot$ represents a non-dimensional variable):

$\bold{\nabla} \rightarrow \frac{1}{L} \hat{\bold{\nabla}}$; $\quad\bold{u} \rightarrow \sqrt{\frac{\alpha g L^3 Q}{\kappa}}\hat{\bold{u}}$; $\quad\partial_t \rightarrow \sqrt{\frac{\alpha g L Q}{\kappa}} \partial_{\hat{t}}$; $\quad P \rightarrow \frac{\rho_o \kappa^2}{L^2} \hat{P}$; $\quad T \rightarrow \frac{L^2 Q}{\kappa} \hat{T}$; $\quad \bold{\Omega} \rightarrow \Omega_0 \hat{\bold{\Omega}}$; $\quad \mathcal{H} \rightarrow Q\hat{\mathcal{H}}$,

and re-arranging, the equations become:

$\hat{\bold{\nabla}} \cdot \hat{\bold{u}} = 0$,

$\partial_{\hat{t}}\hat{\bold{u}} + (\hat{\bold{u}} \cdot \hat{\bold{\nabla}})\hat{\bold{u}} = \frac{\kappa^3}{\alpha g L^5 Q}\hat{{\nabla}}\hat{P} + \hat{T}\bold{\hat{z}} + \sqrt{\frac{\nu^2\kappa}{\alpha g L^5 Q}} \hat{\bold{\nabla}}^2\hat{\bold{u}} - \sqrt{\frac{4 \Omega_0^2 \kappa}{\alpha g L Q}} \hat{\bold{\Omega}} \times \hat{\bold{u}}$

$\partial_{\hat{t}}\hat{T} + (\hat{\bold{u}}\cdot \hat{\bold{\nabla}})\hat{T}  = \sqrt{\frac{\kappa^3}{\alpha g L^5 Q}} \left[ \hat\nabla^2\hat T + \hat{\mathcal{H}} \right]$,

and since $R_F = \frac{\alpha g L^5 Q}{\nu \kappa^2}$, $Pr = \frac{\nu}{\kappa}$, $Ta = \frac{4 \Omega_0^2 L^4}{\nu^2}$, the equations become (dropping the $\hat\cdot$ notation):

$\bold{\nabla} \cdot \bold{u} = 0$,

$\partial_t \bold{u} + (\bold{u}\cdot\bold{\nabla})\bold{u} = \frac{-1}{R_f Pr}\nabla P + T\hat{\bold{z}} + \sqrt{\frac{Pr}{R_F}} \bold{\nabla}^2\bold{u} - \sqrt{\frac{Ta Pr}{R_f}} \bold{\Omega}\times\bold{u}$,

$\partial_t T + (\bold{u}\cdot\bold{\nabla})T = \frac{1}{\sqrt{R_fPr}}\left[\nabla^2 T + \mathcal{H}\right]$,

## Applying perturbation theory

From here we can decompose to an equilibrium state and a small order perturbation:

$\bold{u} = \bold{u}_\text{eq} + \epsilon \bold{u}'$, $\quad T = T_\text{eq} + \epsilon T'$, $\quad P = P_\text{eq} + \epsilon P'$,

where $\epsilon$ is a small parameter. We can substitute these into the non-dimensionalised equations, gather terms of the same order in epsilon, ignore terms of order $\epsilon^2$ and higher, and remember that all $\epsilon^0$ terms are a steady-state equilibrium, so $\partial_t$ terms = 0 and $\bold{u}_\text{eq} = 0$. With this, we get:

### $\epsilon^0$ (equilibrium equations):

$\bold{\nabla} \cdot \bold{u}_\text{eq} = 0$,

$\frac{1}{R_f Pr}\nabla P_\text{eq} = T_\text{eq}\hat{\bold{z}}$,

$-\nabla^2 T_\text{eq} = \mathcal{H}$

### $\epsilon^1$ (perturbation equations):

$\bold{\nabla}\cdot\bold{u}' = 0$,

$\partial_t \bold{u}' + \frac{1}{R_f Pr}\nabla P' - T' \hat{\bold{z}} - \sqrt{\frac{Pr}{R_f}} \bold{\nabla}^2 \bold{u}' + \sqrt{\frac{Ta Pr}{R_f}} \bold{\Omega}\times\bold{u}' = 0$,

$\partial_t T' + (\bold{u'}\bold{\nabla})T_\text{eq} - \frac{1}{\sqrt{R_f Pr}} (\nabla^2 T') = 0$

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
