#include "vof.h"

scalar f[], * interfaces = {f};
scalar D2[];
face vector D2f[];
double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0.;
double epsilon = 1e-6, tauy = 0., n = 1.;

face vector alphav[];
scalar rhov[];

event defaults (i = 0) {
  alpha = alphav;
  rho = rhov;
  mu = new face vector;
}

#ifndef rho
# define rho(f) (clamp(f,0.,1.)*(rho1 - rho2) + rho2)
#endif
#ifndef mu
# define mu(muTemp, mu2, f)  (clamp(f,0.,1.)*(muTemp - mu2) + mu2)
#endif

#ifdef FILTERED
scalar sf[];
#else
# define sf f
#endif

event tracer_advection (i++) {

  /**
  When using smearing of the density jump, we initialise *sf* with the
  vertex-average of *f*. */

#ifndef sf
#if dimension <= 2
  foreach()
    sf[] = (4.*f[] +
	    2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +
	    f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;
#else // dimension == 3
  foreach()
    sf[] = (8.*f[] +
	    4.*(f[-1] + f[1] + f[0,1] + f[0,-1] + f[0,0,1] + f[0,0,-1]) +
	    2.*(f[-1,1] + f[-1,0,1] + f[-1,0,-1] + f[-1,-1] +
		f[0,1,1] + f[0,1,-1] + f[0,-1,1] + f[0,-1,-1] +
		f[1,1] + f[1,0,1] + f[1,-1] + f[1,0,-1]) +
	    f[1,-1,1] + f[-1,1,1] + f[-1,1,-1] + f[1,1,1] +
	    f[1,1,-1] + f[-1,-1,-1] + f[1,-1,-1] + f[-1,-1,1])/64.;
#endif
#endif

#if TREE
  sf.prolongation = refine_bilinear;
  sf.dirty = true; // boundary conditions need to be updated
#endif
}

event properties (i++) {
  foreach_face(x) {
    double ff = (sf[] + sf[-1])/2.;
    alphav.x[] = fm.x[]/rho(ff);
    double muTemp = mu1;
    face vector muv = mu;

    double dudx = (u.x[] - u.x[-1])/Delta;
    double dudy = (u.x[0,1] + u.x[-1,1] - u.x[0,-1] - u.x[-1,-1])/(4.*Delta);
    double dudz = (u.x[0,0,1] + u.x[-1,0,1] - u.x[0,0,-1] - u.x[-1,0,-1])/(4.*Delta);

    double dvdx = (v[] - v[-1])/Delta;
    double dvdy = (v[0,1] - v[0,-1])/(2.*Delta);
    double dvdz = (v[0,0,1] - v[0,0,-1])/(2.*Delta);

    double dwdx = (w[] - w[-1])/Delta;
    double dwdy = (w[0,1] - w[0,-1])/(2.*Delta);
    double dwdz = (w[0,0,1] - w[0,0,-1])/(2.*Delta);

    double D11 = dudx;
    double D22 = dvdy;
    double D33 = dwdz;
    double D12 = 0.5*(dudy + dvdx);
    double D13 = 0.5*(dudz + dwdx);
    double D23 = 0.5*(dvdz + dwdy);

    double D2temp = sqrt((sq(D11) + sq(D22) + sq(D33) + 
                         2*(sq(D12) + sq(D13) + sq(D23))))/sqrt(2.);

    if (tauy > 0.)
      muTemp = tauy/(2.*D2temp + epsilon) + mu1*pow(2.*D2temp + epsilon, n-1);
    
    muv.x[] = fm.x[]*mu(muTemp, mu2, ff);
    D2f.x[] = D2temp;
  }

  foreach_face(y) {
    double ff = (sf[0,0] + sf[0,-1])/2.;
    alphav.y[] = fm.y[]/rho(ff);
    double muTemp = mu1;
    face vector muv = mu;

    double dvdy = (u.y[] - u.y[0,-1])/Delta;
    double dvdx = (u.y[1] + u.y[1,-1] - u.y[-1] - u.y[-1,-1])/(4.*Delta);
    double dvdz = (u.y[0,0,1] + u.y[0,-1,1] - u.y[0,0,-1] - u.y[0,-1,-1])/(4.*Delta);

    double dudy = (u.x[] - u.x[0,-1])/Delta;
    double dudx = (u.x[1] - u.x[-1])/(2.*Delta);
    double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);

    double dwdy = (w[] - w[0,-1])/Delta;
    double dwdx = (w[1] - w[-1])/(2.*Delta);
    double dwdz = (w[0,0,1] - w[0,0,-1])/(2.*Delta);

    double D11 = dudx;
    double D22 = dvdy;
    double D33 = dwdz;
    double D12 = 0.5*(dudy + dvdx);
    double D13 = 0.5*(dudz + dwdx);
    double D23 = 0.5*(dvdz + dwdy);

    double D2temp = sqrt((sq(D11) + sq(D22) + sq(D33) + 
                         2*(sq(D12) + sq(D13) + sq(D23))))/sqrt(2.);

    if (tauy > 0.)
      muTemp = tauy/(2.*D2temp + epsilon) + mu1*pow(2.*D2temp + epsilon, n-1);

    muv.y[] = fm.y[]*mu(muTemp, mu2, ff);
    D2f.y[] = D2temp;
  }

  foreach_face(z) {
    double ff = (sf[0,0,0] + sf[0,0,-1])/2.;
    alphav.z[] = fm.z[]/rho(ff);
    double muTemp = mu1;
    face vector muv = mu;

    double dwdz = (u.z[] - u.z[0,0,-1])/Delta;
    double dwdx = (u.z[1] + u.z[1,0,-1] - u.z[-1] - u.z[-1,0,-1])/(4.*Delta);
    double dwdy = (u.z[0,1] + u.z[0,1,-1] - u.z[0,-1] - u.z[0,-1,-1])/(4.*Delta);

    double dudz = (u.x[] - u.x[0,0,-1])/Delta;
    double dudx = (u.x[1] - u.x[-1])/(2.*Delta);
    double dudy = (u.x[0,1] - u.x[0,-1])/(2.*Delta);

    double dvdz = (v[] - v[0,0,-1])/Delta;
    double dvdx = (v[1] - v[-1])/(2.*Delta);
    double dvdy = (v[0,1] - v[0,-1])/(2.*Delta);

    double D11 = dudx;
    double D22 = dvdy;
    double D33 = dwdz;
    double D12 = 0.5*(dudy + dvdx);
    double D13 = 0.5*(dudz + dwdx);
    double D23 = 0.5*(dvdz + dwdy);

    double D2temp = sqrt((sq(D11) + sq(D22) + sq(D33) + 
                         2*(sq(D12) + sq(D13) + sq(D23))))/sqrt(2.);

    if (tauy > 0.)
      muTemp = tauy/(2.*D2temp + epsilon) + mu1*pow(2.*D2temp + epsilon, n-1);

    muv.z[] = fm.z[]*mu(muTemp, mu2, ff);
    D2f.z[] = D2temp;
  }

  foreach(){
    rhov[] = cm[]*rho(sf[]);
    D2[] = (D2f.x[] + D2f.y[] + D2f.z[] +
            D2f.x[1] + D2f.y[0,1] + D2f.z[0,0,1])/6.;
  }

#if TREE
  sf.prolongation = fraction_refine;
  sf.dirty = true;
#endif
}