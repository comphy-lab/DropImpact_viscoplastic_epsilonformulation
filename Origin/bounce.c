/**
 * Version 2.0
 * Author: Vatsal Sanjay
 * Last updated: Oct 11, 2024

# Introduction:
We investigate the classical problem of VP drop impacting a solid surface.
# Numerical code
Id 1 is for the Viscoplastic liquid drop, and Id 2 is Newtonian gas.
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED // Smear density and viscosity jumps
/**
To model Viscoplastic liquids, we use a modified version of [two-phase.h](http://basilisk.fr/src/two-phase.h). [two-phaseVP.h](two-phaseVP.h) contains these modifications.
*/
#include "two-phaseVP.h"
/**
 You can use: conserving.h as well. Even without it, I was still able to conserve the total energy (also momentum?) of the system if I adapt based on curvature and vorticity/deformation tensor norm (see the adapt even). I had to smear the density and viscosity anyhow because of the sharp ratios in liquid (Bingham) and the gas.
*/
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"
#include "distance.h"

#define tsnap (0.01)

// Error tolerancs
#define fErr (1e-3)                                 // error tolerance in f1 VOF
#define VelErr (1e-2)                               // error tolerances in velocity -- Use 1e-2 for low Oh and 1e-3 to 5e-3 for high Oh/moderate to high J
#define D2Err (1e-2)
#define KAPPAErr (1e-2)

// gas properties!
#define RHO21 (1e-3)
#define MU21 (1e-2)

// domain properties!
#define Ldomain 4

// Distance and radius of drop calculations
#define Xdist (1.02)
#define R2Drop(x,y) (sq(x - Xdist) + sq(y))

// boundary conditions
u.t[left] = dirichlet(0.0);
f[left] = dirichlet(0.0);
// all symmetry planes

int MAXlevel;
double We, Oh, J, tmax;
char nameOut[80], dumpFile[80];

int  main(int argc, char const *argv[]) {

  // Ensure that all the variables were transferred properly from the terminal or job script.
  // if (argc < 6){
  //   fprintf(ferr, "Lack of command line arguments. Check! Need %d more arguments, Level, tauy, We, Oh, tmax\n",6-argc);
  //   return 1;
  // }

  L0 = Ldomain;
  origin (0., 0.);
  init_grid (1 << 8);
  // Values taken from the terminal
  MAXlevel = 9; // atoi(argv[1]);
  J = 1e-1; // atof(argv[2]); // plasto-capillary number
  We = 1e1; // atof(argv[3]); // Weber number
  Oh = 1e-2; // atof(argv[4]); // Ohnesorge number
  tmax = 1e0; // atof(argv[5]);

  fprintf(ferr, "Level %d, We %2.1e, Oh %2.1e, J %4.3f\n", MAXlevel, We, Oh, J);

  // Create a folder named intermediate where all the simulation snapshots are stored.
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
  // Name of the restart file. See writingFiles event.
  sprintf (dumpFile, "dump");

  /**
  We consider the drop impacting a solid surface. To nondimensionalise the governing equations, we use the initial drop radius $R_0$, and inertia-capillary velocity $V_\gamma = \sqrt{\gamma/(\rho_lR_0)}$, respectively. Pressure and stresses are scaled with the characteristic capillary pressure, $\gamma/R_0$. The dimensionless equations for mass and momentum conservation, for the liquid phase, then read

  $$
  \nabla \cdot \boldsymbol{u} = 0
  $$
  $$
  \frac{\partial\boldsymbol{u}}{\partial t} + \nabla\boldsymbol{\cdot}\left(\boldsymbol{uu}\right) = -\nabla p + \nabla\boldsymbol{\cdot}\boldsymbol{\tau} - \mathcal{B}o\,\hat{\boldsymbol{e}}_{\boldsymbol{\mathcal{Z}}},
  $$
  where $\boldsymbol{u}$ is the velocity vector, $t$ is time, $p$ is the pressure and $\boldsymbol{\tau}$ represents the deviatoric stress tensor. We use the regularized Bingham model with
  $$
  \boldsymbol{\tau} = 2\,\text{min}\left(\frac{\mathcal{J}}{2\|\boldsymbol{\mathcal{D}}\|+\epsilon} + \mathcal{O}h\right)\boldsymbol{\mathcal{D}}
  $$
  where $\|\boldsymbol{\mathcal{D}}\|$ is the second invariant of the deformation rate tensor, $\boldsymbol{\mathcal{D}}$, and $\mathcal{O}h_{max}$ is the viscous regularisation parameter. The three dimensionless numbers controlling the equations above are the capillary-Bingham number $\left(\mathcal{J}\right)$, which accounts for the competition between the capillary and yield stresses, the Ohnesorge number $\left(\mathcal{O}h\right)$ that compares the inertial-capillary to inertial-viscous time scales, and the Bond number $\left(\mathcal{B}o\right)$, which compares gravity and surface tension forces:
  $$
  \mathcal{J} = \frac{\tau_yR_0}{\gamma},\,\,\mathcal{O}h = \frac{\mu_l}{\sqrt{\rho_l\gamma R_0}},\,\,\mathcal{B}o = \frac{\rho_l gR_o^2}{\gamma}.
  $$

  Here, $\gamma$ is the liquid-gas surface tension coefficient, and $\tau_y$ and $\rho_l$ are the liquid's yield stress and density, respectively. Next, $\mu_l$ is the constant viscosity in the Bingham model. Note that in our simulations, we also solve the fluid's motion in the gas phase, using a similar set of equations (Newtonian). Hence, the further relevant non-dimensional groups in addition to those above are the ratios of density $\left(\rho_r = \rho_g/\rho_l\right)$ and viscosity $\left(\mu_r = \mu_g/\mu_l\right)$. In the present study, these ratios are kept fixed at $10^{-3}$ and $2 \times 10^{-2}$, respectively (see above). 
  */

  // epsilon = t < tsnap ? 1e-1 : 1e-3;  // epsilon regularisation value of effective viscosity
  epsilon = 1e-2;  // epsilon regularisation value of effective viscosity
  rho1 = 1., rho2 = RHO21;
  mu1 = Oh/sqrt(We), mu2 = MU21*Oh/sqrt(We);
  f.sigma = 1.0/We;
  tauy = J/We;
  CFL = 1e-1;
  run();
}

event init (t = 0) {
  if (!restore (file = dumpFile)){
    refine((R2Drop(x, y) < 1.05) && (level < MAXlevel));
    fraction(f, 1. - R2Drop(x, y));
    foreach() {
      u.x[] = -1.0 * f[];
      u.y[] = 0.0;
    }
  }
}

/**
## Adaptive Mesh Refinement
*/
event adapt(i++){
  if (t < 1e-2){
    adapt_wavelet ((scalar *){f, u.x, u.y},
    (double[]){fErr, VelErr, VelErr},
    MAXlevel);
  } else {
    scalar KAPPA[], D2c[];
    curvature(f, KAPPA);
    foreach() {
      double D11 = (u.y[0,1] - u.y[0,-1])/(2*Delta);
      double D22 = (u.y[]/y);
      double D33 = (u.x[1,0] - u.x[-1,0])/(2*Delta);
      double D13 = 0.5*( (u.y[1,0] - u.y[-1,0] + u.x[0,1] - u.x[0,-1])/(2*Delta) );
      double D2 = (sq(D11)+sq(D22)+sq(D33)+2.0*sq(D13));
      D2c[] = f[]*(D2);
    }
    adapt_wavelet ((scalar *){f, u.x, u.y, KAPPA, D2c},
    (double[]){fErr, VelErr, VelErr, KAPPAErr, D2Err},
    MAXlevel);
  }
}

/**
## Dumping snapshots
*/
event writingFiles (t = 0; t += tsnap; t <= tmax+tsnap) {
  dump (file = dumpFile);
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file=nameOut);
}

/**
## Ending Simulation
*/
event end (t = end) {
  fprintf(ferr, "Done: Level %d, We %2.1e, Oh %2.1e, J %4.3f\n", MAXlevel, We, Oh, J);
}

/**
## Log writing
*/
event logWriting (i++) {
  double ke = 0.;
  foreach (reduction(+:ke)){
    ke += (2*pi*y)*(0.5*rho(f[])*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
  }
  if (pid() == 0){
    static FILE * fp;
    if (i == 0) {
      fprintf (ferr, "i dt t ke\n");
      fp = fopen ("log", "w");
      fprintf (fp, "i dt t ke\n");
      fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
      fclose(fp);
    } else {
      fp = fopen ("log", "a");
      fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
      fclose(fp);
    }
    fprintf (ferr, "%d %g %g %g\n", i, dt, t, ke);
  }
  // Ensure that the cut-off Kinetic energy is smaller than or equal to 1e-6 times the maximum kinetic energy of the system.
  assert(ke > -1e-10);
  assert(ke < 1e3);

  if (ke < 1e-6){
    if (i > 1e2){
      fprintf(ferr, "Kinetic energy is too small. Exiting...\n");
      return 1;
    }
  }

}

/**
## Running the code
~~~bash
#!/bin/bash
qcc -fopenmp -Wall -O2 bounce_VP.c -o bounce_VP -lm -disable-dimensions
export OMP_NUM_THREADS=8
./bounce_VP 10 0.25 1e-2 5.0
~~~
**/