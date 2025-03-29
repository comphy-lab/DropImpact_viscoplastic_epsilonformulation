/**
 * @file bounce.c
 * @author Vatsal Sanjay (vatsalsanjay@gmail.com)
 * vatsalsanjay.com
 * Physics of Fluids
 * @date Aug 24, 2024
 * @version 2.0 
 * 
*/

// 1 is drop
#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"

// Error tolerances
#define fErr (1e-3)                                 // error tolerance in VOF
#define KErr (1e-6)                                 // error tolerance in KAPPA
#define VelErr (1e-2)                            // error tolerances in velocity
#define DissErr (1e-5)                            // error tolerances in dissipation

// air-water
#define Rho21 (1e-3)
// Calculations!
#define SPdist (0.02)
// a0 is the aspect ratio of the ellipse: radial to vertical. 
#define R2Drop(x,y,a0) (sq((x - pow(a0,-0.6666666667) - SPdist)/pow(a0,-0.6666666667)) + sq(y/pow(a0,0.3333333333)))

/*
To describe the equation of an ellipse with a given aspect ratio \( a_0 \), where the horizontal axis is along the \( y \)-axis and the vertical axis is along the \( x \)-axis, we start with the standard form of the ellipse equation:

\[
\frac{y^2}{b^2} + \frac{x^2}{a^2} = 1
\]

Given the aspect ratio \( a_0 \), which is the ratio of the horizontal diameter (2b) to the vertical diameter (2a), we have:

\[
a_0 = \frac{2b}{2a} = \frac{b}{a}
\]

Thus, \( b = a a_0 \). Substituting \( b \) in the ellipse equation, we get:

\[
\frac{y^2}{(a a_0)^2} + \frac{x^2}{a^2} = 1
\]

Simplifying, the equation of the ellipse becomes:

\[
\frac{y^2}{a^2 a_0^2} + \frac{x^2}{a^2} = 1
\]

or equivalently:

\[
\frac{y^2}{a^2 a_0^2} + \frac{x^2}{a^2} = 1
\]

When this ellipse is rotated around the \( x \)-axis, it forms an ellipsoid. The volume \( V \) of an ellipsoid with semi-axes \( a \), \( b \), and \( c \) is given by:

\[
V = \frac{4}{3} \pi a b c
\]

For our specific case, the semi-axes are \( a \) (along \( x \)-axis), \( a a_0 \) (along \( y \)-axis), and \( a a_0 \) (along \( z \)-axis). Therefore, the volume \( V \) of the resulting ellipsoid is:

\[
V = \frac{4}{3} \pi a (a a_0) (a a_0) = \frac{4}{3} \pi a^3 a_0^2
\]

Thus, the volume of the ellipsoid is:

\[
V = \frac{4}{3} \pi a^3 a_0^2
\]

Now, if we need to keep the volume of the ellipsoid constant, we can set \( V = \frac{4}{3} \pi a^3 a_0^2 = \text{constant} \). This implies that the semi-axis \( a \) is equal to \( a_0^{-2/3} \). Therefore, we can write:

\[ \left(\frac{y}{a_0^{1/3}}\right)^2 + \left(\frac{x}{a_0^{-2/3}}\right)^2 = 1 \]

*/

// boundary conditions
u.t[left] = dirichlet(0.0);
f[left] = dirichlet(0.0);

u.n[right] = neumann(0.);
p[right] = dirichlet(0.0);

u.n[top] = neumann(0.);
p[top] = dirichlet(0.0);

int MAXlevel;
double tmax, We, Ohd, Ohs, Bo, a0, Ldomain;
#define MINlevel 2                                            // maximum level
#define tsnap (0.01)

int main(int argc, char const *argv[]) {
  if (argc < 8){
    fprintf(ferr, "Lack of command line arguments. Check! Need %d more arguments\n",8-argc);
    return 1;
  }

  MAXlevel = atoi(argv[1]);
  tmax = atof(argv[2]);
  We = atof(argv[3]); // We is 1 for 0.22 m/s <1250*0.22^2*0.001/0.06>
  Ohd = atof(argv[4]); // <\mu/sqrt(1250*0.060*0.001)>
  Ohs = atof(argv[5]); //\mu_r * Ohd
  Bo = atof(argv[6]);
  Ldomain = atof(argv[7]); // size of domain. must keep Ldomain \gg 1
  a0 = atof(argv[8]); // non-sphericity of drop
  
  fprintf(ferr, "Level %d tmax %g. We %g, Ohd %3.2e, Ohs %3.2e, Bo %g, a %g, Lo %g\n", MAXlevel, tmax, We, Ohd, Ohs, Bo, a0, Ldomain);

  L0=Ldomain;
  X0=0.; Y0=0.;
  init_grid (1 << (4));

  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);

  rho1 = 1.0; mu1 = Ohd/sqrt(We);
  rho2 = Rho21; mu2 = Ohs/sqrt(We);
  f.sigma = 1.0/We;
  G.x = -Bo/We; // Gravity
  run();
}

event init(t = 0){
  if(!restore (file = "dump")){
    refine((R2Drop(x,y,a0) < 1.05) && (level < MAXlevel));
    fraction (f, 1. - R2Drop(x,y,a0));
    foreach () {
      u.x[] = -1.0*f[];
      u.y[] = 0.0;
    }
  }

  // this is a workaround to address: https://github.com/Computational-Multiphase-Physics/basilisk-C_v2024_Jul23/issues/2
  //dump (file = "dump");
  //return 1;

}

scalar KAPPA[], D2c[];
event adapt(i++){
  curvature(f, KAPPA);
  foreach(){
    double D11 = (u.y[0,1] - u.y[0,-1])/(2*Delta);
    double D22 = (u.y[]/max(y,1e-20));
    double D33 = (u.x[1,0] - u.x[-1,0])/(2*Delta);
    double D13 = 0.5*( (u.y[1,0] - u.y[-1,0] + u.x[0,1] - u.x[0,-1])/(2*Delta) );
    double D2 = (sq(D11)+sq(D22)+sq(D33)+2.0*sq(D13));
    D2c[] = f[]*D2;
  }
  adapt_wavelet ((scalar *){f, KAPPA, u.x, u.y, D2c},
     (double[]){fErr, KErr, VelErr, VelErr, DissErr},
      MAXlevel, MINlevel);
  unrefine(x>0.95*Ldomain); // ensure there is no backflow from the outflow walls!
}

// Outputs
// static
event writingFiles (t = 0, t += tsnap; t <= tmax) {
  p.nodump = false; // dump pressure to calculate force in post-processing: see getEpsNForce.c
  dump (file = "dump");
  char nameOut[80];
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);
}

event logWriting (i+=10) {
  double ke = 0., vol = 0.;
  foreach (reduction(+:ke), reduction(+:vol)){
    ke += 2*pi*y*(0.5*rho(f[])*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
    vol += 2*pi*y*f[]*sq(Delta);
  }

  static FILE * fp;

  if (pid() == 0){
    if (i == 0) {
      fprintf (ferr, "i dt t ke p\n");
      fp = fopen ("log", "w");
      fprintf(fp, "Level %d tmax %g. We %g, Ohd %3.2e, Ohs %3.2e, Bo %g, a %g, Lo %g\n", MAXlevel, tmax, We, Ohd, Ohs, Bo, a0, Ldomain);
      fprintf (fp, "i dt t ke vol\n");
      fprintf (fp, "%d %g %g %g %g\n", i, dt, t, ke, vol);
      fclose(fp);
    } else {
      fp = fopen ("log", "a");
      fprintf (fp, "%d %g %g %g %g\n", i, dt, t, ke, vol);
      fclose(fp);
    }
    fprintf (ferr, "%d %g %g %g %g\n", i, dt, t, ke, vol);
  }

}
