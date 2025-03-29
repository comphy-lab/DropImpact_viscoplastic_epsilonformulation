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
#include "two-phaseVP-HB.h"
/**
 You can use: conserving.h as well. Even without it, I was still able to conserve the total energy (also momentum?) of the system if I adapt based on curvature and vorticity/deformation tensor norm (see the adapt even). I had to smear the density and viscosity anyhow because of the sharp ratios in liquid (Bingham) and the gas.
*/
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"
#include "distance.h"
#include "tag.h"
#include "adapt_wavelet_limited.h"

#define MINlevel 3

// Error tolerancs
#define fErr (1e-3)    // error tolerance in VOF
#define KAPPAErr (1e-3)    // error tolerance in KAPPA
#define VelErr (1e-2)  // error tolerances in velocity
#define D2Err (1e-3) // error tolerances in dissipation

// gas properties!
#define RHO21 (1e-3)
#define MU21 (1e-2)

// Distance and radius of drop calculations
#define SPdist (0.02)
// a0 is the aspect ratio of the ellipse: radial to vertical. 
#define R2Drop(x,y,a0) (sq((x - pow(a0,-0.6666666667) - SPdist)/pow(a0,-0.6666666667)) + sq(y/pow(a0,0.3333333333)))

// boundary conditions
u.t[left] = dirichlet(0.0);
f[left] = dirichlet(0.0);

p[right] = dirichlet(0.);
u.n[right] = neumann(0.);

p[top] = dirichlet(0.);
u.n[top] = neumann(0.);

int MAXlevel;
double We, Oh, J, Bo, a0;

double tmax, Ldomain;
double tsnap = 0.01;
char nameOut[80], resultsName[80], dumpFile[80];
int main(int argc, char const *argv[])
{
  origin(0., 0.);
  init_grid(1 << 5);
  MAXlevel = atoi(argv[1]); // 10
  J = atof(argv[2]);        // 0 for Newtonian.
  We = atof(argv[3]);       // 10
  Oh = atof(argv[4]);       // 0.01
  Bo = atof(argv[5]);       // 0 without considering density
  tmax = atof(argv[6]);     // 10
  Ldomain = atof(argv[7]);  // 8
  DT = atof(argv[8]);       // 1e-4
  sprintf(resultsName, "%s", argv[9]);
  a0 = atof(argv[10]); // non-sphericity of drop

  L0 = Ldomain;
  NITERMAX = 500;
  TOLERANCE = 1e-4;

  char comm[80];
  sprintf(comm, "mkdir -p intermediate");
  system(comm);

  // convert non-dimension number normalized by D to that by R
  // We = We / 2.0;
  // Bo = Bo / 4.0;
  // Oh = Oh * 1.414;
  // J = J / 2;

  CFL = 0.05;
  epsilon = t < tsnap ? 1e-1 : 1e-3;
  rho1 = 1., rho2 = RHO21;
  mu1 = Oh / sqrt(We), mu2 = MU21 * Oh / sqrt(We);
  f.sigma = 1.0 / We;
  tauy = J / We;
  G.x = -Bo / We; // uncomment only if Gravity is needed!
  run();
}

event init(t = 0)
{
  if (!restore(file = "dump", list = all))
  {
    fprintf(ferr, "We: %g, Oh: %g, Bo:%g, J:%g, MAXlevel:%d,Ldomain:%g,tmax:%g\n", We*2.0, Oh*1.414, Bo*4.0, J*2, MAXlevel, Ldomain, tmax);
    fprintf(ferr, "t,i,time\n");
    refine((R2Drop(x, y, a0) < 1.05) && (level < MAXlevel));
    fraction(f, 1. - R2Drop(x, y,a0));
    foreach ()
    {
      u.x[] = -1.0 * f[];
      u.y[] = 0.0;
    }
    boundary((scalar *){f, u.x, u.y});
  }
}

int refRegion(double x, double y, double z)
{
  return (((y < 3.0 && x < 0.002) || (y < 0.002)) ? MAXlevel : ((y < 3.0 && x < 0.005) || (y < 0.005)) ? MAXlevel-1 : MAXlevel-2);
}

/**
## Adaptive Mesh Refinement
*/
event adapt(i++)
{
  if (t < 1e-2)
  {
    adapt_wavelet((scalar *){f, u.x, u.y},
                  (double[]){fErr, VelErr, VelErr},
                  MAXlevel);
  }
  else
  {
    scalar KAPPA[];
    curvature(f, KAPPA);
    scalar D2c[];
    foreach ()
    {
      double D11 = (u.y[0, 1] - u.y[0, -1]) / (2 * Delta);
      double D22 = (u.y[] / max(y, 1e-20));
      double D33 = (u.x[1, 0] - u.x[-1, 0]) / (2 * Delta);
      double D13 = 0.5 * ((u.y[1, 0] - u.y[-1, 0] + u.x[0, 1] - u.x[0, -1]) / (2 * Delta));
      double D2 = (sq(D11) + sq(D22) + sq(D33) + 2.0 * sq(D13));
      D2c[] = f[] * D2;
    }
    // adapt_wavelet_limited((scalar *){f, u.x, u.y, KAPPA, D2c},
    //                       (double[]){fErr, VelErr, VelErr, KAPPAErr, D2Err},
    //                       refRegion, MINlevel);
    adapt_wavelet((scalar *){f, u.x, u.y, KAPPA, D2c},
                          (double[]){fErr, VelErr, VelErr, KAPPAErr, D2Err},
                          MAXlevel);
  }
  double time_now = perf.t / 60.0;
  fprintf(ferr, "%g,%d,%g\n", t, i, time_now);
}

double t_last = 0.0;
double DeltaT = 0.0;
int count_run = 0;
double ke_last = 100;
double Rmax_last=0.;
event postProcess(t = 0, t += 0.001; t <= 20)
{
  double ke = 0., xMin = Ldomain, mvx = 0., mvy = 0., Rmax=0.;
  foreach (reduction(+ : ke) reduction(+ : mvx) reduction(+ : mvy) reduction(min : xMin) reduction(max : Rmax))
  {
    mvx += (2 * pi * y) * (clamp(f[], 0., 1.) * u.x[]) * sq(Delta);
    mvy += (2 * pi * y) * (clamp(f[], 0., 1.) * u.y[]) * sq(Delta);
    ke += (2 * pi * y) * (sq(u.x[]) + sq(u.y[])) * rho(f[]) / 2. * sq(Delta);
    if (f[] > 1e-6)
    {
      if ((x < xMin))
      {
        xMin = x;
      }
      if ((y > Rmax))
      {
        Rmax = y;
      }
    }
  }

  double pdatum = 0, wt = 0;
  foreach_boundary(top, reduction(+ : pdatum), reduction(+ : wt))
  {
    pdatum += 2 * pi * y * p[] * (Delta);
    wt += 2 * pi * y * (Delta);
  }
  // only run 5 step for the fine step
  if (DT < 5e-5 && count_run > 2)
  {
    return 1;
  }
  count_run = count_run + 1;

  if (ke < 0 || ke > 5)
  {
    fprintf(ferr, "Ke_Error, Exit...\n");
    return 1;
  }
  else
  {
    p.nodump = true;
    dump(file = "dump");
  }
  // log
  DeltaT = perf.t / 60.0 - t_last;
  t_last = perf.t / 60.0;
  static FILE *fp1;
  if (pid() == 0)
  {
    if (i == 0)
    {
      fp1 = fopen("log_run", "w");
      fprintf(fp1, "t,i,Cell,Wallclocktime(min),CPUtime(min),ke,mvx,mvy\n");
      fflush(fp1);
    }
    fp1 = fopen("log_run", "a");
    fprintf(fp1, "%g,%d,%ld,%g,%g,%g,%g,%g\n", t, i, grid->tn, perf.t / 60.0, DeltaT, ke, mvx, mvy);
    fflush(fp1);
  }
  // stop condition
  // if ((t > tmax - tsnap) || (t > 0.5 && (ke < 1e-4 || (xMin > 0.04))))
  if ((t > tmax - tsnap) || (t > 0.5 && (Rmax < Rmax_last)))
  {
    char comm[256];
    fprintf(ferr, "Reach Max time. Or Kinetic energy is too small or droplet bounce off. Exiting...\n");
    return 1;
  }
  Rmax_last=Rmax;
}

event snapshot(t = 0, t += tsnap; t <= 20)
{
  p.nodump = true;
  sprintf(nameOut, "intermediate/restart_snapshot-%5.4f", t);
  dump(file = nameOut);
  p.nodump = false;
  sprintf(nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file = nameOut);
}
