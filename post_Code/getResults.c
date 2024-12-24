/* Title: Calculating the jet height, and tip location
# For estimating the theta till the very end, condition on y is eased
# Authors: Vatsal Sanjay & Ayush Dixit
# vatsalsanjay@gmail.com
# Physics of Fluids
*/
#include "axi.h"
#include "navier-stokes/centered.h"
#include "vof.h"
#include "tag.h"
#include "heights.h"
#include "curvature.h"

// #include "droplet_stat.h"
double rho1 = 1., rho2 = 1e-3;
#define rho(f) (clamp(f,0.,1.)*(rho1 - rho2) + rho2)

char filename[80], nameTrack[80];
scalar f[], * interfaces = {f};

// First argument should be filename.
int main(int a, char const *arguments[])
{
  sprintf (filename, "%s", arguments[1]);
  // sprintf(nameTrack, "%s", arguments[2]);

  restore (file = filename);
  boundary((scalar *){f, u.x, u.y});

  // tag all liquid parts starts
  // to tag all the liquid parts
  scalar d[];
  double threshold = 1e-4;
  // d array stores if the liquid is higher than threshold
  foreach(){
    d[] = (f[] > threshold);
  }

  // Any connected region for f > threshold is given a unique tag from 0 to n-1
  int n = tag (d), size[n];
  for (int i = 0; i < n; i++){
    size[i] = 0;
  }

  // size array stores the size of each connected region
  foreach_leaf(){
    if (d[] > 0){
      size[((int) d[]) - 1]++;
    }
  }

  // MainPhase is the tag of the largest connected region
  int MaxSize = 0;
  int MainPhase = 0;
  for (int i = 0; i < n; i++){
     // fprintf(ferr, "%d %d\n",i, size[i]);
    if (size[i] > MaxSize){
      MaxSize = size[i];
      MainPhase = i+1;
    }
  }

  // top droplet 
  double v[n];
  double px[n];
  double ux[n];
  for (int ii=0;ii<n;ii++){
    v[ii]=px[ii]=ux[ii]=0.0;
  }
  foreach (serial){
    if (d[] > 0) {
      int j = d[] - 1;
      v[j] += sq(Delta)*(2*pi*y)*f[];
      px[j] += sq(Delta)*(2*pi*y)*f[]*x;
      ux[j] += sq(Delta)*(2*pi*y)*f[]*u.x[];
    }
  }
  // To get the upper droplet and the corresponding coord, volume and velocity.
  double x_top_droplet,v_top_droplet,u_top_droplet;
  int max_index = -1;
  for (int i_px = 0; i_px < n; ++i_px) {
    if (px[i_px]/v[i_px] > x_top_droplet) {
        max_index = i_px;
    }
  }
  x_top_droplet = px[max_index]/v[max_index];
  v_top_droplet = v[max_index];
  u_top_droplet = ux[max_index]/v[max_index];
  double r_top_droplet=pow(v_top_droplet/(4.*pi/3.), 1./3.);
  // tag all liquid parts ends

  
  face vector s[];
  s.x.i = -1;
  double yMin = 0.2;
  double vTip, xTip = -30;
  double kappaTip=0;
  double vTop, xTop = -30;
  double R_max = 0; 
  double R_max_bottom = 0; 
  scalar kappa[];
  curvature(f, kappa);
  double xcmax = 0.0, ycmax = 0.0;
  double kappamax = -HUGE;

  // Finding the tip position
  foreach(){
    
    // if (kappa[] != nodata && x < 0. && f[] < 1-1e-6 && y > 0.1 && x > -1.5 && y < 1)
    if (kappa[] != nodata && f[] < 1-1e-6)
    {   
        if (kappa[] > kappamax)
        {
            kappamax = kappa[];
            xcmax = x;
            ycmax = y;
        }
      //fprintf(ferr, "%f %f %f %f\n", x, y, f[], kappa[]);
    }

    if (f[] > 1e-6 && f[] < 1. - 1e-6 && d[] == MainPhase) {
      coord n1 = facet_normal (point, f, s);
      double alpha1 = plane_alpha (f[], n1);
      coord segment1[2];
      if (facets (n1, alpha1, segment1) == 2){
        double x1 = x + (segment1[0].x+segment1[1].x)*Delta/2.;
        double y1 = y + (segment1[0].y+segment1[1].y)*Delta/2.;
        if (y1 > R_max){
          R_max = y1;
        }
      }
    }
  }
  foreach_boundary(left){
    if (f[] > 1e-6 && f[] < 1. - 1e-6 && d[] == MainPhase) {
      coord n1 = facet_normal (point, f, s);
      double alpha1 = plane_alpha (f[], n1);
      coord segment1[2];
      if (facets (n1, alpha1, segment1) == 2){
        double x1 = x + (segment1[0].x+segment1[1].x)*Delta/2.;
        double y1 = y + (segment1[0].y+segment1[1].y)*Delta/2.;
        if (y1 > R_max_bottom && x1 < yMin){
          R_max_bottom = y1;
        }
      }
    }
  }

  foreach_boundary(bottom){
    if (f[] > 1e-6 && f[] < 1. - 1e-6) {
      coord n1 = facet_normal (point, f, s);
      double alpha1 = plane_alpha (f[], n1);
      coord segment1[2];
      if (facets (n1, alpha1, segment1) == 2){
        double x1 = x + (segment1[0].x+segment1[1].x)*Delta/2.;
        double y1 = y + (segment1[0].y+segment1[1].y)*Delta/2.;
        if (x1 > xTip && d[] == MainPhase){
          xTip = x1;
          vTip = u.x[];
          if (kappa[] != nodata){
            kappaTip=kappa[];
          }
          
        }
        if (x1 > xTop){
          xTop = x1;
          vTop = u.x[];
        }
      }
    }
  }

  double rr = sqrt(sq(xcmax) + sq(ycmax));
  double thetap=0;
  if (rr > 0){
    thetap = acos(xcmax / rr)/pi;
  }

  // Beigin:velocity of the center mass normal force F
  double ke = 0., vb=0., sb=0., xb=0;
  foreach (reduction(+:ke) reduction(+:sb) reduction(+:vb)){
    double dv = clamp(f[], 0., 1.) * sq(Delta)*(2*pi*y);
    sb += dv;
    vb += u.x[] * dv;
    xb += x * dv;
    ke += sq(Delta)*(2*pi*y)*(sq(u.x[]) + sq(u.y[]))*rho(f[])/2.;
  }
  vb = vb/sb;
  xb = xb/sb;
  double pdatum = 0, wt = 0;
  foreach_boundary(top, reduction(+:pdatum), reduction(+:wt)){
    pdatum += 2*pi*y*p[]*(Delta);
    wt += 2*pi*y*(Delta);
  }
  if (wt >0){
    pdatum /= wt;
  }

  double pforce = 0.;
  foreach_boundary(left, reduction(+:pforce)){
    pforce += 2*pi*y*(p[]-pdatum)*(Delta);
  }
  // End

  FILE * fp = ferr;
  fprintf(ferr, "%f %i %7.6e %7.6e %7.6e %7.6e %7.6e %7.6e %7.6e %7.6e %7.6e %7.6e %7.6e %7.6e %7.6e %7.6e %7.6e %7.6e %7.6e %7.6e\n", 
              t, n, xTip, vTip,kappaTip, xTop, vTop, kappamax,
              xcmax, ycmax, thetap, R_max,R_max_bottom,
              x_top_droplet,r_top_droplet,u_top_droplet,
              vb,xb,ke,pforce);
  //fprintf(ferr, "%7.6e %7.6e\n", xcmax, ycmax);
  fflush (fp);
  fclose (fp);

//  FILE *fp2;
//  fp2 = fopen (nameTrack, "a");
  //fprintf(fp2, "%f %7.6e %7.6e %7.6e\n", t, xTip, yMin, vTip);
//  fclose(fp2);

}
