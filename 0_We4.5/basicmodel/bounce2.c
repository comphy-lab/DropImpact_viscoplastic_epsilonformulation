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
#define fErr (1e-3)     // error tolerance in VOF
#define KErr (1e-6)     // error tolerance in KAPPA
#define VelErr (1e-2)   // error tolerances in velocity
#define DissErr (1e-5)  // error tolerances in dissipation
#define OmegaErr (1e-2) // error tolerances in vorticity

// gas properties!
#define RHO21 (1e-3)
#define MU21 (1e-2)

// Distance and radius of drop calculations
#define Xdist (1.02)
#define R2Drop(x, y) (sq(x - Xdist) + sq(y))

// boundary conditions
u.t[left] = dirichlet(0.0);
u.n[left] = dirichlet(0.0);
uf.t[left] = dirichlet(0.0);
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
  init_grid(1 << 5);
  MAXlevel = atoi(argv[1]); //10
  J = atof(argv[2]);  // 0 for Newtonian.
  We = atof(argv[3]); // 10
  Oh = atof(argv[4]); // 0.01
  Bo = atof(argv[5]); // 0 without considering density
  epsilon = atof(argv[6]); // 1e-2
  tmax = atof(argv[7]);    // 10
  Ldomain = atof(argv[8]); // 8
  sprintf(resultsName, "%s", argv[9]); 
  fprintf(ferr, "Ldomain %4.3e, tmax %4.3e, Level %d, We %2.1e, Oh %2.1e, Bo %2.1e, J %4.3f, epsilon %4.3e\n", Ldomain, tmax, MAXlevel, We, Oh, Bo, J, epsilon);

  L0 = Ldomain;
  DT = 1e-6;
  NITERMAX = 1000;

  char comm[80];
  sprintf(comm, "mkdir -p intermediate");
  system(comm);

  rho1 = 1., rho2 = RHO21;
  mu1 = Oh / sqrt(We), mu2 = MU21 * Oh / sqrt(We);
  f.sigma = 1.0 / We;
  tauy = J / We;
  G.x = -Bo / We; // uncomment only if Gravity is needed!
  run();
}

event init(t = 0)
{
  // sprintf(dumpFile, "dump");
  // if (!restore(file = dumpFile))
  // {
  refine(R2Drop(x, y) < 1.05 && (level < MAXlevel));
  fraction(f, 1. - R2Drop(x, y));
  foreach ()
  {
    u.x[] = -1.0 * f[];
    u.y[] = 0.0;
  }
  // }
}

/**
## Adaptive Mesh Refinement
*/
event adapt(i++)
{
  if (t < 1e-2){
    adapt_wavelet ((scalar *){f, u.x, u.y},
    (double[]){fErr, VelErr, VelErr},
    MAXlevel);
  }else{
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
    adapt_wavelet((scalar *){f, KAPPA, u.x, u.y, D2c},
                  (double[]){fErr, KErr, VelErr, VelErr, DissErr},
                  MAXlevel, MINlevel);
    unrefine(x > 0.95 * Ldomain); // ensure there is no backflow from the outflow walls!
  }
  if (t>1.0){
    DT=1e-4;
  }
  if (t>5.0){
    DT=1e-3;
  }
}

/**
## Dumping snapshots
*/
event writingFiles(t = 0; t += tsnap; t <= tmax+tsnap)
{
  p.nodump = false;
  dump(file = dumpFile);
  sprintf(nameOut, "intermediate/snapshot-%5.4f", t);
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

double x_min_min = HUGE;
event postProcess(t += tsnap)
{
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

  double ke = 0., vcm = 0., vcm1 = 0.;
  face vector s[];
  s.x.i = -1;
  double yMax = -HUGE;
  double xMax = -HUGE;
  double xMin = HUGE;
  foreach (reduction(+ : ke) reduction(+ : vcm) reduction(+ : vcm1) reduction(max : yMax) reduction(max : xMax) reduction(min : xMin))
  {
    //vcm and vcm1 are used for the calculation of the normal force.
    vcm += (2 * pi * y) * (clamp(f[], 0., 1.) * u.x[]) * sq(Delta);
    ke += sq(Delta) * (2 * pi * y) * (sq(u.x[]) + sq(u.y[])) * rho(f[]) / 2.;
    if (d[] == MainPhase)
    {
      vcm1 += (2 * pi * y) * (clamp(f[], 0., 1.) * u.x[]) * sq(Delta);      
    }
    //yMax,xMax and xMin anre the boundary of the main droplet
    if (f[] > 1e-6 && f[] < 1. - 1e-6 && d[] == MainPhase)
    {
      coord n1 = facet_normal(point, f, s);
      double alpha1 = plane_alpha(f[], n1);
      coord segment1[2];
      if (facets(n1, alpha1, segment1) == 2)
      {
        double x1 = x + (segment1[0].x + segment1[1].x) * Delta / 2.;
        double y1 = y + (segment1[0].y + segment1[1].y) * Delta / 2.;
        if (y1 > yMax)
        {
          yMax = y1;
        }
        if (y1 < pow(0.5, MAXlevel) * Ldomain*1.5)
        {
          if (x1 > xMax)
          {
            xMax = x1;
          }
        }
        if (x1 < xMin)
        {
          xMin = x1;
        }
      }
    }    
  } 
  if (xMin < x_min_min)
  {
    x_min_min = xMin;
  }

  double pdatum = 0, wt = 0;
  foreach_boundary(top, reduction(+ : pdatum), reduction(+ : wt))
  {
    pdatum += 2 * pi * y * p[] * (Delta);
    wt += 2 * pi * y * (Delta);
  }
  if (wt > 0)
  {
    pdatum /= wt;
  }

  double pforce = 0.;
  foreach_boundary(left, reduction(+ : pforce))
  {
    pforce += 2 * pi * y * (p[] - pdatum) * (Delta);
  }

  static FILE *fp1;
  if (pid() == 0)
  {
    if (i == 0)
    {
      fp1 = fopen("result", "w");
      fprintf(fp1, "t,n,ke,vcm,vcm1,Rmax,Zmax,Zmin,pforce\n");
      fclose(fp1);
    }
    fp1 = fopen("result", "a");
    fprintf(fp1, "%g,%d,%g,%g,%g,%g,%g,%g,%g\n", t,n,ke,vcm,vcm1,yMax,xMax,xMin,pforce);
    fclose(fp1);
  }

  // if ((t == tmax) || (t > 1 && (ke < 1e-6 || (x_min - x_min_min > pow(0.5, MAXlevel) * Ldomain))))
  if ((t == tmax) || (t > 1 && (ke < 1e-6 || (xMin - x_min_min > 0.02))))
  {
    char comm[256];
    sprintf(comm, "cp result ../%s", resultsName);
    system(comm);
    fprintf(ferr, "Kinetic energy is too small or droplet bounce off. Exiting...\n");
    return 1;
  }
}
