// keplers_goat_herd.cpp -- Oliver Philcox, 2021.
// Solve Kepler's equation with a variety of numerical techniques for an array of mean anomalies
//
// Methods implemented are:
//   - Newton-Raphson: The quadratic Newton-Raphson root finder.
//   - Danby: The quartic root finder described in Danby (1988).
//   - Series: An elliptic series method, as described in Murray & Dermott.
//   - Contour: A method based on contour integration, described in Philcox et al. (2021).
//
// Given an array of mean anomalies, an eccentricity and a desired precision, the code will estimate the eccentric anomaly using each method.
// We generate a grid of mean anomalies from a uniformly spaced grid in eccentric anomaly.
// The hyperparameter of each approach is increased until the desired precision is reached, relative to the true result.
// Timing is performed using the `chrono' package.

// ported to java by aholinch 

public class keplers_goat_herd
{
	public static final double M_PI = Math.PI;

	// ========================== Class to hold basic attributes ================

	public static class Approximations {
		private int N_ell; // number of ell points
		private double e; // eccentricity
		private double[] ell_arr; // array of mean anomalies

		public Approximations(int _N_ell, double _e, double[] _ell_arr){
			N_ell = _N_ell;
			e = _e;

			// Create and copy ell array
			ell_arr = new double[N_ell];
			for(int i=0;i<N_ell;i++) ell_arr[i] = _ell_arr[i];

		}

		// ========================== Approximation Methods ================

		public void compute_newton_raphson(int N_it, double output[]){
			/// Compute Newton-Raphson method with given number of steps
			// This has quadratic convergence.
			// The initial step is defined as E_0 = ell + sgn(sin(ell))*e*k following Danby (1988)

			double k = 0.85;
			double f_E, fP_E, this_ell, old_E;

			for(int i=0;i<N_ell;i++){
				this_ell = ell_arr[i];

				// Define initial estimate
				if((Math.sin(this_ell))<0) old_E = this_ell - k*e;
				else old_E = this_ell + k*e;

				// Perform Newton-Raphson estimate
				for(int j=0;j<N_it;j++) {

					// Compute f(E) and f'(E)
					f_E = old_E - e*Math.sin(old_E)-this_ell;
					fP_E = 1. - e*Math.cos(old_E);

					// Update E
					old_E -= f_E/fP_E;
				}

				// Add to array
				output[i] = old_E;
			}
		}

		public void compute_danby(int N_it, double output[]){
			/// Compute Danby (1988) fourth-order root-finding method with given number of steps
			// This has quartic convergence.
			// The initial step is defined as E_0 = ell + sgn(Math.sin(ell))*e*k following Danby (1988)

			double k = 0.85;
			double f_E, fP_E, fPP_E, fPPP_E, this_ell, old_E, delta_i1, delta_i2, delta_i3, esinE, ecosE;

			for(int i=0;i<N_ell;i++){
				this_ell = ell_arr[i];

				// Define initial estimate
				if((Math.sin(this_ell))<0) old_E = this_ell - k*e;
				else old_E = this_ell + k*e;

				// Perform Newton-Raphson estimate
				for(int j=0;j<N_it;j++) {

					// Compute f(E), f'(E), f''(E) and f'''(E), avoiding recomputation of sine and cosine.
					esinE = e*Math.sin(old_E);
					ecosE = e*Math.cos(old_E);
					f_E = old_E - esinE-this_ell;
					fP_E = 1. - ecosE;
					fPP_E = esinE;
					fPPP_E = ecosE;

					delta_i1 = -f_E/fP_E;
					delta_i2 = -f_E/(fP_E+1./2.*delta_i1*fPP_E);
					delta_i3 = -f_E/(fP_E+1./2.*delta_i2*fPP_E+1./6.*fPPP_E*delta_i2*delta_i2);

					// Update E
					old_E += delta_i3;
				}

				// Add to array
				output[i] = old_E;
			}
		}

		public void compute_series(int N_it, double output[]){
			// Solve Kepler's equation via the series method described in Murray & Dermot.
			// We use a specifies maximum number of iterations and precompute the coupling coefficients.

			double coeff = 0;

			// Take initial guess
			for(int i=0; i<N_ell;i++){
				output[i] = ell_arr[i];
			}

			// Iterate over number of iterations
			for(int s=1; s<=N_it; s++){
				// TODO find a Bessel function implementation compatible with MIT
				//coeff = 2*jn(s,s*e)/s; // define coefficient

				for(int i=0;i<N_ell;i++){
					output[i] += coeff*Math.sin(s*ell_arr[i]);
				}
			}
		}

		public void compute_contour(int N_it, double output[]){
			// Solve Kepler's equation via the contour integration method of Philcox et al. (2021)
			// This uses techniques described in Ullisch (2020) to solve the `geometric goat problem'.
			// N_it specifies the number of grid-points.

			double ft_gx2, ft_gx1, this_ell, freq, zR, zI, cosC, sinC, esinRadius, ecosRadius, center;
			double fxR, fxI, ftmp, tmpcosh, tmpsinh, tmpcos, tmpsin;

			// Define sampling points (actually use one more than this)
			int N_points = N_it-2;
			int N_fft = (N_it-1)*2;

			// Define contour radius
			double radius = e/2;

			// Generate e^{ikx} sampling points and precompute real and imaginary parts
			double cf, sf;
			double exp2R[] = new double[N_points];
			double exp2I[] = new double[N_points];
			double exp4R[] = new double[N_points];
			double exp4I[] = new double[N_points];
			double coshI[] = new double[N_points];
			double sinhI[] = new double[N_points];
			double ecosR[] = new double[N_points];
			double esinR[] = new double[N_points];
			for(int jj=0;jj<N_points;jj++){
				// NB: j = jj+1
				freq = 2.0*M_PI*(jj+1)/N_fft;
				cf = Math.cos(freq);
				sf = Math.sin(freq);
				exp2R[jj] = cf;
				exp2I[jj] = sf;
				exp4R[jj] = cf*cf-sf*sf;
				exp4I[jj] = 2.0*cf*sf;
				coshI[jj] = Math.cosh(radius*exp2I[jj]);
				sinhI[jj] = Math.sinh(radius*exp2I[jj]);
				ecosR[jj] = e*Math.cos(radius*exp2R[jj]);
				esinR[jj] = e*Math.sin(radius*exp2R[jj]);
			}

			// Precompute e sin(e/2) and e cos(e/2)
			esinRadius = e*Math.sin(radius);
			ecosRadius = e*Math.cos(radius);

			// Iterate over array of mean anomalies
			for(int i=0;i<N_ell;i++){
				this_ell = ell_arr[i];

				// Define contour center for each ell and precompute sin(center), cos(center)
				if(this_ell<M_PI) center = this_ell+e/2;
				else center = this_ell-e/2;
				sinC = Math.sin(center);
				cosC = Math.cos(center);
				output[i] = center;

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
				output[i] += radius*ft_gx2/ft_gx1;
			}
		}
	}

	// ========================== Timing Comparison ================

	public static void main(String args[]) {

		// PARAMETERS
		int N_ell = 1000000; // ell array size
		double e = 0.5; // Eccentricity
		double tol = 1e-12; // tolerance for error acceptance

		// Print parameters
		System.out.printf("N_ell = %d\n",N_ell);
		System.out.printf("e = %.2f\n",e);
		System.out.printf("tolerance = %.2e\n",tol);

		// Define ell array from a linearly spaced grid of E
		double[] E_exact = new double[N_ell];
		double[] ell_arr = new double[N_ell];
		for(int i=0;i<N_ell;i++){
			E_exact[i] = 2.0*M_PI*(i+0.5)/N_ell;
			ell_arr[i] = E_exact[i]-e*Math.sin(E_exact[i]);
		}

		// Create output class to hold methods
		Approximations approx = new Approximations(N_ell, e, ell_arr);

		// Initialize timers
		long tstart = 0;
		long tstop = 0;

		// Output estimates
		double[] E_newton_raphson = new double[N_ell];
		double[] E_Danby = new double[N_ell];
		double[] E_series = new double[N_ell];
		double[] E_contour = new double[N_ell];

		// Compute Newton-Raphson quadratic estimate
		int N_NR = 0; // Newton-Raphson iterations
		double err_NR = 0;
		double duration_NR = 0;

		// Increase N_NR until we reach tolerance!
		while (N_NR<100){ // max limit!
			tstart = System.nanoTime(); // starting time
			approx.compute_newton_raphson(N_NR,E_newton_raphson);
			tstop =  System.nanoTime(); // ending time
			duration_NR = (tstop-tstart)/1000.0d; // duration

			err_NR = 0;
			for(int i=0;i<N_ell;i++) err_NR += Math.abs(E_exact[i]-E_newton_raphson[i])/N_ell;
			if(err_NR<tol) break;
			N_NR ++;

		}
		System.out.printf("Computed Newton-Raphson estimate in %d steps after %.1f ms with mean-error %.2e\n",N_NR,duration_NR/1000,err_NR);

		// Compute Danby quartic estimate
		int N_Danby = 0; // Danby iterations
		double err_Danby = 0;
		double duration_Danby = 0;

		while (N_Danby<100){ // max limit!
			tstart = System.nanoTime(); // starting time
			approx.compute_danby(N_Danby,E_Danby);
			tstop =  System.nanoTime(); // ending time
			duration_Danby = (tstop-tstart)/1000.0d; // duration

			err_Danby = 0;
			for(int i=0;i<N_ell;i++) err_Danby += Math.abs(E_exact[i]-E_Danby[i])/N_ell;
			if(err_Danby<tol) break;
			N_Danby++;

		}
		System.out.printf("Computed Danby estimate in %d steps after %.1f ms with mean-error %.2e\n",N_Danby,duration_Danby/1000,err_Danby);

		// Compute series estimate
		int N_series = 0; // Series iterations
		double err_series = 0;
		double duration_series = 0;

		if(1>0)
		{
			System.out.println("Need a Bessel implementation in order to compute series");
		}
		else
		{
			if(e>0.6627434){
				System.out.printf("### Series method is non-convergent; skipping!\n");
			}
			else{
				while (N_series<100){ // max limit!
					tstart = System.nanoTime(); // starting time
					approx.compute_series(N_series,E_series);
					tstop =  System.nanoTime(); // ending time
					duration_series = (tstop-tstart)/1000.0d; // duration

					err_series = 0;
					for(int i=0;i<N_ell;i++) err_series += Math.abs(E_exact[i]-E_series[i])/N_ell;
					if(err_series<tol) break;
					N_series++;
				}
				System.out.printf("Computed series estimate in %d steps after %.1f ms with mean-error %.2e\n",N_series,duration_series/1000,err_series);
			}
		}

		// Compute contour estimate
		int N_contour = 2; // number of integration steps
		double err_contour = 0;
		double duration_contour = 0;

		while (N_contour<256){ // max limit!
			tstart = System.nanoTime(); // starting time
			approx.compute_contour(N_contour,E_contour);
			tstop =  System.nanoTime(); // ending time
			duration_contour = (tstop-tstart)/1000.0d; // duration

			err_contour = 0;
			for(int i=0;i<N_ell;i++) err_contour += Math.abs(E_exact[i]-E_contour[i])/N_ell;
			if(err_contour<tol) break;
			N_contour+=1;
		}
		System.out.printf("Computed contour estimate in %d steps after %.1f ms with mean-error %.2e\n",N_contour,duration_contour/1000,err_contour);

	}
}
