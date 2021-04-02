// contour implementation from Oliver Philcox, 2021.
// https://github.com/oliverphilcox/Keplers-Goat-Herd
// java port by aholinch

public class kgh
{
    
    /**
     * Solve keplers equation for E where m = E - e*sinE
     * m is in radians
     * e is unitless between 0 and 1
     */
    public static double mToE(double m, double e){
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
    
    /**
     * Solve keplers equation for E where m = E - e*sinE
     * m is in radians
     * e is unitless between 0 and 1
     * N_it integration steps
     */
    public static double mToE(double m, double e, int N_it){
        // skip trivial cases
        if(e == 0) return m;
        if(m == 0) return 0;
        if(m == 2.0*Math.PI) return 2.0*Math.PI;
        
        
        // Solve Kepler's equation via the contour integration method of Philcox et al. (2021)
        // This uses techniques described in Ullisch (2020) to solve the `geometric goat problem'.
        // N_it specifies the number of grid-points.

        double ft_gx2, ft_gx1, freq, zR, zI, cosC, sinC, esinRadius, ecosRadius, center;
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
            freq = 2.0*Math.PI*(jj+1)/N_fft;
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

        double output = 0;
        
        // Define contour center for each ell and precompute sin(center), cos(center)
        if(m<Math.PI) center = m+e/2;
        else center = m-e/2;
        sinC = Math.sin(center);
        cosC = Math.cos(center);
        output = center;

        // Accumulate Fourier coefficients
        // NB: we halve the range by symmetry, absorbing factor of 2 into ratio

        ///////////////
        // Separate out j = 0 piece, which is simpler

        // Compute z in real and imaginary parts (zI = 0 here)
        zR = center + radius;

        // Compute e*sin(zR) from precomputed quantities
        tmpsin = sinC*ecosRadius+cosC*esinRadius; // sin(zR)

        // Compute f(z(x)) in real and imaginary parts (fxI = 0)
        fxR = zR - tmpsin - m;

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

            fxR = zR - tmpsin*tmpcosh-m;
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
        fxR = zR - tmpsin-m;

        // Add to sum, with 1/2 factor for edges
        ft_gx2 += 0.5/fxR;
        ft_gx1 += -0.5/fxR;

        ///////////////
        // Compute E(ell)
        output += radius*ft_gx2/ft_gx1;

        return output;
    }

    public static void main(String args[])
    {
        // PARAMETERS
        int N_ell = 10000; // ell array size
        double e = 0.5; // Eccentricity
        double tol = 1e-12; // tolerance for error acceptance
      
        // Print parameters
        System.out.printf("N_ell = %d\n",N_ell);
        System.out.printf("e = %.2f\n",e);
        System.out.printf("tolerance = %.2e\n",tol);

        // Define ell array from a linearly spaced grid of E
        double[] E_exact = new double[N_ell];
        double[] ell_arr = new double[N_ell];
        double[] E_contour = new double[N_ell];

        for(int i=0;i<N_ell;i++){
            E_exact[i] = 2.0*Math.PI*(i+0.5)/N_ell;
            ell_arr[i] = E_exact[i]-e*Math.sin(E_exact[i]);
        }

        // Compute contour estimate
        int N_contour = 2; // number of integration steps
        double err_contour = 0;
        double duration_contour = 0;
        long tstart = 0;
        long tstop = 0;

        while (N_contour<256){ // max limit!
            tstart = System.nanoTime(); // starting time
            for(int i=0; i<N_ell; i++)
            {
                E_contour[i] = mToE(ell_arr[i],e,N_contour);
            }
            tstop =  System.nanoTime(); // ending time
            duration_contour = (tstop-tstart)/1000.0d; // duration

            err_contour = 0;
            for(int i=0;i<N_ell;i++) err_contour += Math.abs(E_exact[i]-E_contour[i])/N_ell;
            if(err_contour<tol) break;
            N_contour+=1;
        }
        System.out.printf("Computed contour estimate in %d steps after %.1f ms with mean-error %.2e\n",N_contour,duration_contour/1000,err_contour);

        System.exit(0);
    }
}
