@if _XOPEN_SOURCE < 700
  @undef _XOPEN_SOURCE
  @define _XOPEN_SOURCE 700
@endif
@if _GNU_SOURCE
@include <stdint.h>
@include <string.h>
@include <fenv.h>
@endif
#define _CATCH
#define dimension 2
#define BGHOSTS 2
#include "common.h"
#include "grid/quadtree.h"
#ifndef BASILISK_HEADER_0
#define BASILISK_HEADER_0
#line 1 "getEpsForce.c"
/* Title: Getting Energy
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
# Last Update: Sep 06 2021
*/
#include "axi.h"
#include "navier-stokes/centered.h"
#include "fractions.h"
#include "tag.h"
#include "heights.h"

scalar f[];
double Oh, mu1, J, tauy, epsOh, epsJ, We, pForce;

char filename[80], nameEnergy[80];

int main(int a, char const *arguments[])
{
  sprintf(filename, "%s", arguments[1]);
  Oh = atof(arguments[2]);
  We = atof(arguments[3]);
  J = atof(arguments[4]);

  sprintf(filename, "%s", arguments[1]);
  restore(file = filename);

  // boundary conditions
  u.t[left] = dirichlet(0.);
  f[left] = dirichlet(0.0);
  u.n[right] = neumann(0.);
  p[right] = dirichlet(0.0);
  u.n[top] = neumann(0.);
  p[top] = dirichlet(0.0);

  mu1 = Oh / sqrt(We);
  tauy = J / We;
  f.prolongation = refine_bilinear;
  boundary((scalar *){f, u.x, u.y, p});

  // fprintf(ferr, "Oh %3.2e, We %g\n", Oh, We);
  // return 1;

  // tag all liquid parts starts
  scalar d[];
  double threshold = 1e-4;
  foreach ()
  {
    d[] = (f[] > threshold);
  }
  int n = tag(d), size[n];
  for (int i = 0; i < n; i++)
  {
    size[i] = 0;
  }
  foreach_leaf()
  {
    if (d[] > 0)
    {
      size[((int)d[]) - 1]++;
    }
  }
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
  // tag all liquid parts ends

  scalar sf[];
  foreach ()
    sf[] = (4. * f[] +
            2. * (f[0, 1] + f[0, -1] + f[1, 0] + f[-1, 0]) +
            f[-1, -1] + f[1, -1] + f[1, 1] + f[-1, 1]) /
           16.;
  sf.prolongation = refine_bilinear;
  boundary({sf});
  /*
  Do calculations start
  */
  epsOh = 0., epsJ = 0., pForce = 0.;

  face vector s[];
  s.x.i = -1;
  double yMax = -HUGE;
  double xMax = -HUGE;
  double vTip = 0., xTP = 0.;
  double uTip = 0., yTP = 0.;

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
        double y1 = y + (segment1[0].y + segment1[1].y) * Delta / 2.;
        if (y1 > yMax)
        {
          yMax = y1;
          xTP = x1;
          vTip = interpolate(u.y, xTP, yMax);
        }
        if (y1 < 0.01)
        {
          if (x1 > xMax)
          {
            xMax = x1;
            yTP = y1;
            uTip = interpolate(u.x, xMax, yTP);
          }
        }
      }
    }
  }
  double vcm = 0., vol = 0.;
  double vcm1 = 0., vol1 = 0.;
  foreach ()
  {
    // zcm += (2*pi*y)*(clamp(sf[], 0., 1.)*x)*sq(Delta);
    vcm += (2 * pi * y) * (clamp(sf[], 0., 1.) * u.x[]) * sq(Delta);
    vol += (2 * pi * y) * clamp(sf[], 0., 1.) * sq(Delta);
    if (d[] == MainPhase)
    {
      vcm1 += (2 * pi * y) * (clamp(sf[], 0., 1.) * u.x[]) * sq(Delta);
      vol1 += (2 * pi * y) * clamp(sf[], 0., 1.) * sq(Delta);
    }

    double D11 = (u.y[0, 1] - u.y[0, -1]) / (2 * Delta);
    double D22 = (u.y[] / max(y, 1e-20));
    double D33 = (u.x[1, 0] - u.x[-1, 0]) / (2 * Delta);
    double D13 = 0.5 * ((u.y[1, 0] - u.y[-1, 0] + u.x[0, 1] - u.x[0, -1]) / (2 * Delta));
    double D2 = (sq(D11) + sq(D22) + sq(D33) + 2.0 * sq(D13));
    epsOh += (2 * pi * y) * (2 * mu1 * clamp(sf[], 0., 1.) * D2) * sq(Delta);
    epsJ += (2 * pi * y) * (tauy * clamp(sf[], 0., 1.) * sqrt(D2)) * sq(Delta);
  }

  // calculate the force on the substrate
  double pdatum = 0, wt = 0;
  foreach_boundary(top)
  {
    pdatum += 2 * pi * y * p[] * (Delta);
    wt += 2 * pi * y * (Delta);
  }
  if (wt > 0)
  {
    pdatum /= wt;
  }
  foreach_boundary(left)
  {
    pForce += 2 * pi * y * (p[] - pdatum) * (Delta);
  }

  /*
  Do calculations end
  */

  fprintf(ferr, "%6.5e %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e %6.5e\n",
          t, vcm, vcm / vol, vcm1, vcm1 / vol1, epsOh, epsJ, pForce, yMax, vTip, xMax, uTip);
}

#endif
