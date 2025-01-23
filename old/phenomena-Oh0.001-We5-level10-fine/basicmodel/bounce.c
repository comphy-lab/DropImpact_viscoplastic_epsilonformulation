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
#define tsnap (0.01)

// Error tolerancs
#define fErr (1e-3)     // error tolerance in VOF
#define KErr (1e-2)     // error tolerance in KAPPA
#define VelErr (1e-2)   // error tolerances in velocity
#define DissErr (1e-2)  // error tolerances in dissipation
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
  DT = atof(argv[9]); // 1e-3
  CFL = atof(argv[10]); // 1e-3
  sprintf(resultsName, "%s", argv[11]); 

  //log
  fprintf(ferr, "We,Oh,J,MAXlevel,epsilon,MU21,DT,Ldomain,tmax,Bo\n");
  fprintf(ferr, "%g,%g,%g,%d,%4.3e,%g,%g,%g,%g,%g\n", We, Oh, J, 13, epsilon, MU21, DT, Ldomain, tmax, Bo);
  
  L0 = Ldomain;
  NITERMAX = 1000;

  rho1 = 1., rho2 = RHO21;
  mu1 = Oh / sqrt(We), mu2 = MU21 * Oh / sqrt(We);
  f.sigma = 1.0 / We;
  tauy = J / We;
  G.x = -Bo / We; // uncomment only if Gravity is needed!
  run();
}

event init(t = 0)
{
  refine(R2Drop(x, y) < 1.05 && (level < MAXlevel));
  fraction(f, 1. - R2Drop(x, y));
  foreach ()
  {
    u.x[] = -1.0 * f[];
    u.y[] = 0.0;
  }
}

int refRegion(double x, double y, double z){
  return (((y < 2.0 && x < 0.002) || (y < 0.002)) ? (MAXlevel < 13 ? 13 : MAXlevel):
          ((y < 2.0 && x < 0.005) || (y < 0.005)) ? (MAXlevel < 12 ? 12 : MAXlevel):
          ((y < 2.0 && x < 0.01) || (y < 0.01)) ? (MAXlevel < 11 ? 11 : MAXlevel):
          ((y < 2.0 && x < 0.02) || (y < 0.02)) ? (MAXlevel < 10 ? 10 : MAXlevel):
          MAXlevel);
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
  }
  else{
    scalar KAPPA[];
    curvature(f, KAPPA);
    scalar D2c[];
    foreach (){
      double D11 = (u.y[0, 1] - u.y[0, -1]) / (2 * Delta);
      double D22 = (u.y[] / max(y, 1e-20));
      double D33 = (u.x[1, 0] - u.x[-1, 0]) / (2 * Delta);
      double D13 = 0.5 * ((u.y[1, 0] - u.y[-1, 0] + u.x[0, 1] - u.x[0, -1]) / (2 * Delta));
      double D2 = (sq(D11) + sq(D22) + sq(D33) + 2.0 * sq(D13));
      D2c[] = f[] * D2;
    }
    adapt_wavelet_limited((scalar *){f, KAPPA, u.x, u.y, D2c},
              (double[]){fErr, KErr, VelErr, VelErr, DissErr},
              refRegion, MINlevel);
    unrefine(x > 0.95 * Ldomain);
  }
  if (t>5){
    DT=1e-3;
  }
}

/**
## Dumping snapshots
*/
event writingFiles(t += tsnap)
{
  if (i==0){
    char comm[80];
    sprintf(comm, "mkdir -p intermediate");
    system(comm);
  }
  p.nodump = false;
  dump(file = "dump");
  sprintf(nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file = nameOut);
}

double t_last = 0.0;
double DeltaT = 0.0;
double x_min_min = HUGE;
event postProcess(t += tsnap)
{
  scalar d[];
  double threshold = 1e-4;
  // d array stores if the liquid is higher than threshold
  foreach (){
    d[] = (f[] > threshold);
  }

  // Any connected region for f > threshold is given a unique tag from 0 to n-1
  int n_d = tag(d), size[n_d];
  for (int i = 0; i < n_d; i++){
    size[i] = 0;
  }

  foreach (serial)
    if (d[] > 0)
      size[((int)d[]) - 1]++;
  #if _MPI
    MPI_Allreduce(MPI_IN_PLACE, size, n_d, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  #endif

  // MainPhase is the tag of the largest connected region
  int MaxSize = 0;
  int MainPhase = 0;
  for (int i = 0; i < n_d; i++){
    if (size[i] > MaxSize){
      MaxSize = size[i];
      MainPhase = i + 1;
    }
  }

  double ke = 0., xMin = HUGE;
  foreach (reduction(+ : ke) reduction(min : xMin)){
    ke += sq(Delta) * (2 * pi * y) * (sq(u.x[]) + sq(u.y[])) * rho(f[]) / 2.;
    if (d[] == MainPhase){
      if ((x < xMin)){
        xMin = x;
      }
    }               
  }
  if (xMin < x_min_min){
    x_min_min = xMin;
  }

  DeltaT = perf.t / 60.0 - t_last;
  t_last = perf.t / 60.0;
  static FILE *fp1;
  if (pid() == 0){
    if (i == 0){
      fp1 = fopen("log_run", "w");   
      fprintf(fp1, "t,i,Cell,Wallclocktime(min),CPUtime(min),ke,Zmin,Zmin1\n");fflush(fp1);
    }
    fp1 = fopen("log_run", "a");  
    fprintf(fp1, "%g,%d,%d,%g,%g,%g,%g,%g\n", t,i,grid->tn,perf.t / 60.0, DeltaT,ke,xMin,xMin-x_min_min);fflush(fp1);
  }

  if ((t > tmax-tsnap) || (ke < 1e-6) || (xMin - x_min_min > 0.02))
  {
    char comm[256];
    sprintf(comm, "cp log_run ../Results_Running/log_%s.csv", resultsName);
    system(comm);
    fprintf(ferr, "Kinetic energy is too small or droplet bounce off. Exiting...\n");
    return 1;
  }
}

event end (t = tmax) 
{

}