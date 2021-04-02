// keplers_goat_herd.cpp -- Oliver Philcox, 2021.
// Solve Kepler's equation with a variety of numerical techniques for an array of mean anomalies
//
// Timing is performed using the `chrono' package.

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <complex>
#include <chrono>
#include <cmath>
using namespace std::chrono;

// Could swap between single and double precision here.
typedef double Float;

Float mToE(double m, double e);
Float mToE(double m, double e, int N_it);

Float mToE(double m, double e){
  if(e <= 0.5)
    return mToE(m,e,10);

  if(e <= 0.9)
    return mToE(m,e,25);

  if(e <= 0.95)
    return mToE(m,e,50);

  if(e <= 0.99)
    return mToE(m,e,128);

  return mToE(m,e,256);
}


Float mToE(double m, double e, int N_it){
  // Solve Kepler's equation via the contour integration method of Philcox et al. (2021)
  // This uses techniques described in Ullisch (2020) to solve the `geometric goat problem'.
  // N_it specifies the number of grid-points.

  Float ft_gx2, ft_gx1, this_ell, freq, zR, zI, cosC, sinC, esinRadius, ecosRadius, center;
  Float fxR, fxI, ftmp, tmpcosh, tmpsinh, tmpcos, tmpsin;

  // Define sampling points (actually use one more than this)
  int N_points = N_it-2;
  int N_fft = (N_it-1)*2;

  // Define contour radius
  Float radius = e/2;

  // Generate e^{ikx} sampling points and precompute real and imaginary parts
  Float cf, sf, exp2R[N_points], exp2I[N_points], exp4R[N_points], exp4I[N_points], coshI[N_points], sinhI[N_points], ecosR[N_points], esinR[N_points];
  for(int jj=0;jj<N_points;jj++){
    // NB: j = jj+1
    freq = 2.0*M_PI*(jj+1)/N_fft;
    cf = cos(freq);
    sf = sin(freq);
    exp2R[jj] = cf;
    exp2I[jj] = sf;
    exp4R[jj] = cf*cf-sf*sf;
    exp4I[jj] = 2.0*cf*sf;
    coshI[jj] = cosh(radius*exp2I[jj]);
    sinhI[jj] = sinh(radius*exp2I[jj]);
    ecosR[jj] = e*cos(radius*exp2R[jj]);
    esinR[jj] = e*sin(radius*exp2R[jj]);
  }

  // Precompute e sin(e/2) and e cos(e/2)
  esinRadius = e*sin(radius);
  ecosRadius = e*cos(radius);

  this_ell = m;

  // Define contour center for each ell and precompute sin(center), cos(center)
  if(this_ell<M_PI) center = this_ell+e/2;
  else center = this_ell-e/2;
  sinC = sin(center);
  cosC = cos(center);
  double output = center;

  // Accumulate Fourier coefficients
  // NB: we halve the range by symmetry, absorbing factor of 2 into ratio

  ///////////////
  // Separate out j = 0 piece, which is simpler

  // Compute z in real and imaginary parts (zI = 0 here)
  zR = center + radius;

  // Compute e*sin(zR) from precomputed quantities
  tmpsin = sinC*ecosRadius+cosC*esinRadius; // sin(zR)

  // Compute f(z(x)) in real and imaginary parts (fxI = 0)
  fxR = zR - tmpsin - this_ell;

  // Add to array, with factor of 1/2 since an edge
  ft_gx2 = 0.5/fxR;
  ft_gx1 = 0.5/fxR;

  ///////////////
  // Compute for j = 1 to N_points
  // NB: j = jj+1
  for(int jj=0;jj<N_points;jj++){

    // Compute z in real and imaginary parts
    zR = center + radius*exp2R[jj];
    zI = radius*exp2I[jj];

    // Compute f(z(x)) in real and imaginary parts
    // can use precomputed cosh / sinh / cos / sin for this!
    tmpcosh = coshI[jj]; // cosh(zI)
    tmpsinh = sinhI[jj]; // sinh(zI)
    tmpsin = sinC*ecosR[jj]+cosC*esinR[jj]; // e sin(zR)
    tmpcos = cosC*ecosR[jj]-sinC*esinR[jj]; // e cos(zR)

    fxR = zR - tmpsin*tmpcosh-this_ell;
    fxI = zI - tmpcos*tmpsinh;

    // Compute 1/f(z) and append to array
    ftmp = fxR*fxR+fxI*fxI;
    fxR /= ftmp;
    fxI /= ftmp;

    ft_gx2 += (exp4R[jj]*fxR+exp4I[jj]*fxI);
    ft_gx1 += (exp2R[jj]*fxR+exp2I[jj]*fxI);
  }

  ///////////////
  // Separate out j = N_it piece, which is simpler

  // Compute z in real and imaginary parts (zI = 0 here)
  zR = center - radius;

  // Compute sin(zR) from precomputed quantities
  tmpsin = sinC*ecosRadius-cosC*esinRadius; // sin(zR)

  // Compute f(z(x)) in real and imaginary parts (fxI = 0 here)
  fxR = zR - tmpsin-this_ell;

  // Add to sum, with 1/2 factor for edges
  ft_gx2 += 0.5/fxR;
  ft_gx1 += -0.5/fxR;

  ///////////////
  // Compute E(ell)
  output += radius*ft_gx2/ft_gx1;

  return output;
}

// ========================== Timing Comparison ================

int main(int argc, char *argv[]) {

  // PARAMETERS
  int N_ell = 10000; // ell array size
  Float e = 0.5; // Eccentricity
  Float tol = 1e-12; // tolerance for error acceptance

  // Print parameters
  printf("N_ell = %d\n",N_ell);
  printf("e = %.2f\n",e);
  printf("tolerance = %.2e\n",tol);

  // Define ell array from a linearly spaced grid of E
  Float* E_exact = new Float[N_ell];
  Float* ell_arr = new Float[N_ell];
  for(int i=0;i<N_ell;i++){
    E_exact[i] = 2.0*M_PI*(i+0.5)/N_ell;
    ell_arr[i] = E_exact[i]-e*sin(E_exact[i]);
  }

  // Initialize timers
  auto start = high_resolution_clock::now();
  auto stop = high_resolution_clock::now();

  // Output estimates
  Float* E_contour = new Float[N_ell];

  // Compute contour estimate
  int N_contour = 2; // number of integration steps
  Float err_contour;
  long long int duration_contour;

  while (N_contour<256){ // max limit!
    start = high_resolution_clock::now(); // starting time
    for(int i=0; i<N_ell; i++){
      E_contour[i] = mToE(ell_arr[i],e,N_contour);
    }
    stop = high_resolution_clock::now(); // ending time
    duration_contour = duration_cast<microseconds>(stop - start).count(); // duration

    err_contour = 0;
    for(int i=0;i<N_ell;i++) err_contour += fabs(E_exact[i]-E_contour[i])/N_ell;
    if(err_contour<tol) break;
    N_contour+=1;
  }
  printf("Computed contour estimate in %d steps after %.1f ms with mean-error %.2e\n",N_contour,float(duration_contour/1000.),err_contour);

}
