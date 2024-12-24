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
#include "tag.h"

#define MINlevel 3
#define tsnap (0.01)

// Error tolerancs
#define fErr (1e-3)   // error tolerance in f1 VOF
#define VelErr (1e-2) // error tolerances in velocity -- Use 1e-2 for low Oh and 1e-3 to 5e-3 for high Oh/moderate to high J
#define D2Err (1e-2)
#define KAPPAErr (1e-6)

// gas properties!
#define RHO21 (1e-3)
#define MU21 (1e-2)

// Distance and radius of drop calculations
#define Xdist (1.02)
#define R2Drop(x, y) (sq(x - Xdist) + sq(y))

// boundary conditions
u.t[left] = dirichlet(0.0);
f[left] = dirichlet(0.0);

p[right] = dirichlet(0.);
u.n[right] = neumann(0.);

p[top] = dirichlet(0.);
u.n[top] = neumann(0.);

int MAXlevel;
double We, Oh, J, Bo;
double tmax, Ldomain;
char nameOut[80], resultsName[80], dumpFile[80];

int main(int argc, char const *argv[])
{

  origin(0., 0.);
  init_grid(1 << 6);
  MAXlevel = atoi(argv[1]);
  J = atof(argv[2]);  // plasto-capillary number
  We = atof(argv[3]); // Weber number
  Oh = atof(argv[4]); // Ohnesorge number
  Bo = atof(argv[5]);
  epsilon = atof(argv[6]); // 1e-2
  tmax = atof(argv[7]);    // 10
  Ldomain = atof(argv[8]); // 6
  sprintf(resultsName, "%s", argv[9]);
  fprintf(ferr, "Ldomain %4.3e, tmax %4.3e, Level %d, We %2.1e, Oh %2.1e, Bo %2.1e, J %4.3f, epsilon %4.3e\n", Ldomain, tmax, MAXlevel, We, Oh, Bo, J, epsilon);

  L0 = Ldomain;
  DT = 1e-5;
  NITERMAX = 1000;

  G.x = -Bo / We; // uncomment only if Gravity is needed!
  rho1 = 1., rho2 = RHO21;
  mu1 = Oh / sqrt(We), mu2 = MU21 * Oh / sqrt(We);
  f.sigma = 1.0 / We;
  tauy = J / We;
  CFL = 1e-1;
  run();
}

event init(t = 0)
{
  // Create a folder named intermediate where all the simulation snapshots are stored.
  char comm[80];
  sprintf(comm, "mkdir -p intermediate");
  system(comm);
  // Name of the restart file. See writingFiles event.
  sprintf(dumpFile, "dump");

  if (!restore(file = dumpFile))
  {
    refine(R2Drop(x, y) < 1.1 && (level < MAXlevel));
    fraction(f, 1. - R2Drop(x, y));
    foreach ()
    {
      u.x[] = -1.0 * f[];
      u.y[] = 0.0;
    }
  }
}

/**
## Adaptive Mesh Refinement
*/
event adapt(i++)
{
  if (t < 1e-2)
  {
    adapt_wavelet((scalar *){f, u.x, u.y}, (double[]){fErr, VelErr, VelErr}, MAXlevel);
  }
  else
  {
    scalar KAPPA[], D2c[];
    curvature(f, KAPPA);
    foreach ()
    {
      double D11 = (u.y[0, 1] - u.y[0, -1]) / (2 * Delta);
      double D22 = (u.y[] / y);
      double D33 = (u.x[1, 0] - u.x[-1, 0]) / (2 * Delta);
      double D13 = 0.5 * ((u.y[1, 0] - u.y[-1, 0] + u.x[0, 1] - u.x[0, -1]) / (2 * Delta));
      double D2 = (sq(D11) + sq(D22) + sq(D33) + 2.0 * sq(D13));
      D2c[] = f[] * (D2);
    }
    adapt_wavelet((scalar *){f, u.x, u.y, KAPPA, D2c}, (double[]){fErr, VelErr, VelErr, KAPPAErr, D2Err}, MAXlevel);
  }
  // if (t>1){
  //   DT=1e-4;
  // }
}

/**
## Dumping snapshots
*/
event writingFiles(t = 0; t += tsnap; t <= tmax)
{
  dump(file = dumpFile);
  sprintf(nameOut, "intermediate/snapshot-%5.4f", t);
  p.nodump = false;
  dump(file = nameOut);
}

double t_last = 0.0;
double DeltaT = 0.0;
event logstats(t += tsnap)
{
  DeltaT = perf.t / 60.0 - t_last;
  t_last = perf.t / 60.0;
  // Output i, timestep, no of cells, real time elapsed, cpu time
  static FILE *fp;
  static FILE *fp_state;
  if (pid() == 0)
  {
    if (i == 0)
    {
      fp = fopen("log_run", "w");
      fprintf(fp, "t i Cell Wallclocktime(min) CPUtime(min) \n");
      fflush(fp);
    }
    fprintf(fp, "%g %i %ld %g %g\n", t, i, grid->tn, perf.t / 60.0, DeltaT);
    fflush(fp);
  }
}

event postProcess(t = 0; t += tsnap; t <= tmax)
{
  // calculate the kinetic energy
  double ke = 0.;
  foreach (reduction(+ : ke))
  {
    ke += sq(Delta) * (2 * pi * y) * (sq(u.x[]) + sq(u.y[])) * rho(f[]) / 2.;
  }

  scalar d[];
  double threshold = 1e-4;
  // d array stores if the liquid is higher than threshold
  foreach ()
  {
    d[] = (f[] > threshold);
  }

  // Any connected region for f > threshold is given a unique tag from 0 to n-1
  int n = tag(d), size[n];
  for (int i = 0; i < n; i++)
  {
    size[i] = 0;
  }

  // size array stores the size of each connected region
  foreach_leaf()
  {
    if (d[] > 0)
    {
      size[((int)d[]) - 1]++;
    }
  }
  // MainPhase is the tag of the largest connected region
  int MaxSize = 0;
  int MainPhase = 0;
  for (int i = 0; i < n; i++)
  {
    // fprintf(ferr, "%d %d\n",i, size[i]);
    if (size[i] > MaxSize)
    {
      MaxSize = size[i];
      MainPhase = i + 1;
    }
  }
  // X coord of the bottom point of the biggest droplet
  double x_min = Ldomain;
  face vector s[];
  s.x.i = -1;
  foreach ()
  {
    if (f[] > 1e-6 && f[] < 1. - 1e-6 && d[] == MainPhase)
    {
      coord n1 = facet_normal(point, f, s);
      double alpha1 = plane_alpha(f[], n1);
      coord segment1[2];
      if (facets(n1, alpha1, segment1) == 2)
      {
        double x1 = x + (segment1[0].x + segment1[1].x) * Delta / 2.;
        if (x1 < x_min)
        {
          x_min = x1;
        }
      }
    }
  }

  static FILE *fp1;
  if (pid() == 0)
  {
    if (i == 0)
    {
      fp1 = fopen("result", "w");
      fprintf(fp1, "t,n_d,ke,x_min\n");
      fclose(fp1);
    }
    fp1 = fopen("result", "a");
    fprintf(fp1, "%g,%d,%g,%g\n", t, n, ke, x_min);
    fclose(fp1);
  }

  // if ((t > 1 && (ke < 1e-6 || x_min > 0.1)))
  if (t==tmax)
  {
    if (pid() == 0) {
      char comm[256];
      sprintf(comm, "cp result ../%s", resultsName);
      system(comm);
    }
    fprintf(ferr, "Kinetic energy is too small or droplet bounce off. Exiting...\n");
    return 1;
  }
}
