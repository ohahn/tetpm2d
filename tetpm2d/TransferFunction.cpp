//
//  TransferFunction.cpp
//  tetpm2d
//
//  Created by Oliver Hahn on 8/8/12.
//  Copyright (c) 2012 Oliver Hahn. All rights reserved.
//

#include "TransferFunction.hh"
#include "Global.h"
#include <gsl/gsl_integration.h>

/* Convenience from Numerical Recipes in C, 2nd edition */
double sqrarg;
double cubearg;
double pow4arg;

/* cosmological parameters */
double
 g_H0,       // Hubble constant
g_Omega_b,      // baryon content
g_Omega_m,//0.276,       // matter content
g_Omega_L,//0.724,       // dark energy content
g_sigma_8,//0.811,        // spectrum normalization
g_nspect,//0.961,        // spectral index
g_zstart;

const double
g_deltac  = 1.686,// critical density for spherical collapse
g_G       = 43.0117902, // gravitational constant in Mpc km^2/s^2 1.0e+10Msun
g_rhoc    = 3.0*10000.0/8.0/M_PI/g_G; // critical density

#define REL_PRECISION 1.0e-7

TransferFunction_Eisenstein *ptf = NULL;

double cosmo_norm_fac = 1.0;

/**********/

double integrate( double (* func) (double x, void * params), double a, double b )
{
    gsl_function F;
    F.function = func;
    F.params = NULL;
    
    double result;
    double error;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(10000);
    gsl_integration_qags( &F, a, b, 0, REL_PRECISION, 10000, w, &result, &error );
    gsl_integration_workspace_free(w);
    
    
    //if( error > REL_PRECISION )
    //    std::cerr << " - Warning: no convergence in function 'integrate', rel. error=" << error << std::endl;
    
    return result;
}

double integrate( double (* func) (double x, void * params), std::vector<double>& params,
                 double a, double b )
{
    gsl_function F;
    F.function = func;
    F.params = (void*)(&params[0]);
    
    double result;
    double error;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(10000);
    gsl_integration_qags( &F, a, b, 0, REL_PRECISION, 10000, w, &result, &error );
    gsl_integration_workspace_free(w);
    
    
    //if( error > REL_PRECISION )
    //    std::cerr << " - Warning: no convergence in function 'integrate', rel. error=" << error << std::endl;
    
    return result;
}

inline static double GrowthIntegrand( double a, void *Params=NULL )
{
    double eta = sqrt(g_Omega_m/a+g_Omega_L*a*a
                      +1.0-g_Omega_m-g_Omega_L);
    return 2.5/(eta*eta*eta);
}

double Dplus( double a )
{
    double eta =  sqrt(g_Omega_m/a + g_Omega_L*a*a
                       + 1.0 - g_Omega_m-g_Omega_L );
    double integral = integrate( &GrowthIntegrand, 0.0, a );
    return eta/a*integral;
}

double dsigma( double k, void *Params )
{
    if( k<=0.0 )
        return 0.0f;
    
    double x = k*((double*)Params)[0];

    double w = 3.0*(sin(x)-x*cos(x))/(x*x*x);
    
    double tfk = ptf->compute(k);
    return k*k * w*w * pow(k,g_nspect) * tfk*tfk;
}

float ComputeVFact( float a ){
    float fomega, dlogadt, eta;
    float Omega_k = 1.0 - g_Omega_m - g_Omega_L;
    
    float Dp = Dplus( a )/Dplus(1.0);
    
    eta     = sqrt( (double)(g_Omega_m/a+ g_Omega_L*a*a + Omega_k ));
    fomega  = (2.5/Dp-1.5*g_Omega_m/a-Omega_k)/eta/eta;
    dlogadt = a*eta;
    
    //... /100.0 since we would have to multiply by H0 to convert
    //... the displacement to velocity units. But displacement is
    //... in Mpc/h, and H0 in units of h is 100.
    return fomega * dlogadt/a *100.0;
}

void compute_powerspectrum( const std::vector<double>& k, std::vector<double>& Pk, double z )
{
    double D0 = Dplus(1.0);
    double Dz = Dplus(1.0/(1.0+z))/D0;
    
    //.. compute normalization first ..
    double sigma0;
    {
        std::vector<double> scale(1,0.0);
        scale[0] = 8.0;
        double intt(0.0), dint(1.0), dx(2.0*M_PI/scale[0]), atmp(0.0);
        for(;;){
            dint = integrate( &dsigma, scale, atmp, atmp+dx );
            atmp+=dx;
            intt+=dint;
            if( dint/intt < REL_PRECISION )
                break;
        }
        sigma0 = sqrt(4.0*M_PI*intt);
    }
    
    
    //... compute normalized power spectrum
    Pk.clear();
    
    static double nspect = g_nspect;
#if TRANSFERF_TYPE==0
    static TransferFunction_BBKS TF(false);
#elif TRANSFERF_TYPE==1
    static TransferFunction_BBKS TF(true);
#elif TRANSFERF_TYPE==2
    static TransferFunction_Eisenstein TF;
#else
#error "Undefined Transfer Function Type."
#endif
    
    for(unsigned i=0; i<k.size(); ++i ){  
        double tfk = TF.compute(k[i]);
        Pk.push_back( Dz*Dz*pow(k[i],nspect)*tfk*tfk
                     /sigma0/sigma0*g_sigma_8*g_sigma_8 );
        
    }
    
}


/**********/

double cosmo_get_amp_k( double k )
{
    
    if( ptf == NULL )
    {
        ptf = new TransferFunction_Eisenstein();
        
        double D0 = Dplus(1.0);
        double Dz = Dplus(1.0/(1.0+g_zstart))/D0;
        
        //.. compute normalization first ..
        double sigma0;
        {
            std::vector<double> scale(1,0.0);
            scale[0] = 8.0;
            double intt(0.0), dint(1.0), dx(2.0*M_PI/scale[0]), atmp(0.0);
            for(;;){
                dint = integrate( &dsigma, scale, atmp, atmp+dx );
                atmp+=dx;
                intt+=dint;
                if( dint/intt < REL_PRECISION )
                    break;
            }
            sigma0 = sqrt(4.0*M_PI*intt);
            
            cosmo_norm_fac = Dz * g_sigma_8/sigma0;
        }
        
    }
    
    return ptf->compute(k) * pow( k, 0.5*g_nspect ) * cosmo_norm_fac;
}

void cosmo_set_parameters( double zstart, double Omega_m, double Omega_b, double Omega_L, double n_s, double sigma8, double h )
{
    g_zstart  = zstart;
    g_Omega_m = Omega_m;
    g_Omega_b = Omega_b;
    g_Omega_L = Omega_L;
    g_nspect  = n_s;
    g_sigma_8 = sigma8;
    g_H0      = 100.0 * h;
}

/**************/

//create cosmo particle setup
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "fftw3.h"

void gaussran( float& gaussran1, float& gaussran2 )
{
    
}

void cosmo_init_particles( unsigned seed )
{
    fftwf_real *data = new fftwf_real[ nres * (nres+2) ];
    fftwf_complex *cdata = reinterpret_cast<fftwf_complex*>(data);
    
    fftwf_real *data2 = new fftwf_real[ nres * (nres+2) ];
    fftwf_complex *cdata2 = reinterpret_cast<fftwf_complex*>(data2);
    
    gsl_rng	*RNG = gsl_rng_alloc( gsl_rng_mt19937 );
	gsl_rng_set( RNG, seed );
	
    
    
    fftwf_plan plan, iplan, plan2, iplan2;
    plan  = fftwf_plan_dft_r2c_2d( nres, nres, data, cdata, FFTW_MEASURE ),
    iplan = fftwf_plan_dft_c2r_2d( nres, nres, cdata, data, FFTW_MEASURE );
    
    plan2  = fftwf_plan_dft_r2c_2d( nres, nres, data2, cdata2, FFTW_MEASURE ),
    iplan2 = fftwf_plan_dft_c2r_2d( nres, nres, cdata2, data2, FFTW_MEASURE );
    
    /////////////////////////////////
    
    int nresp = nres/2+1;
    float kfac = 2.0*M_PI/boxlength;
    float gaussran1, gaussran2;
    
    float fftnorm = 1.0f / (float)nres * (2.0f*M_PI/boxlength);
    
    for( int i=0; i<nres; ++i )
        for( int j=0; j<nres; ++j )
        {
            int idx = i*(nres+2)+j;
            data[idx] = gsl_ran_ugaussian_ratio_method( RNG ) / nres;
        }
    fftwf_execute( plan );
    
    /////////////////////////////////
    
    
    for( int i=0; i<nres; ++i )
        for( int j=0; j<nresp; ++j )
        {
            float kx = i>=nresp? (float)(i-nres)*kfac : (float)i*kfac;
            float ky = (float)j*kfac;
            
            float kk = sqrtf(kx*kx+ky*ky);
            
            int idx = i*nresp+j;
            
            float ampk = cosmo_get_amp_k( kk ); //*sqrtf(kk);
            
            if( kk >= nresp*kfac )
                ampk = 0.0;
            
            //cdata[idx][0] = gsl_ran_ugaussian_ratio_method( RNG ) * ampk * fftnorm;
            //cdata[idx][1] = gsl_ran_ugaussian_ratio_method( RNG ) * ampk * fftnorm;
            
            cdata[idx][0] *= ampk * fftnorm;
            cdata[idx][1] *= ampk * fftnorm;
            
            
            cdata2[idx][0] = cdata[idx][0];
            cdata2[idx][1] = cdata[idx][1];
            
            
            /*cdata[idx][0] *= ampk * fftnorm * fftnorm;
            cdata[idx][1] *= ampk * fftnorm * fftnorm;*/
            
            
            
        }
    cdata2[0][0] = 0.0f;
    cdata2[0][1] = 0.0f;
      
    
    // insert code to make random numbers independent of resolution (have rectangle outliens)
    
    
    float dx = boxlength / nres;
    float vfact = ComputeVFact( 1.0f/(1.0f+g_zstart));
    
    /////////////////////////////////
    // generate x-component
    for( int i=0; i<nres; ++i )
        for( int j=0; j<nresp; ++j )
        {
            float kx = i>=nresp? (float)(i-nres)*kfac : (float)i*kfac;
            float ky = (float)j*kfac;
            
            float kk = sqrtf(kx*kx+ky*ky);
            
            int idx = i*nresp+j; // (a+ib) * ik = iak -bk
            
            cdata2[idx][0] = kx/kk/kk * cdata[idx][1];
            cdata2[idx][1] = -kx/kk/kk * cdata[idx][0];
        }
    
    cdata2[0][0] = 0.0f;
    cdata2[0][1] = 0.0f;
    
    fftwf_execute( iplan2 );
    
    for( int i=0; i<nres; ++i )
        for( int j=0; j<nres; ++j )
        {
            int idx = i*(nres+2)+j;
            int ii = i*nres+j;
            P[ii].x = (float)i*dx + data2[idx];
            P[ii].vx = data2[idx] * vfact;
            P[ii].id = ii;
            P[ii].acc[0] = 0.0f;
            
        }
    
    /////////////////////////////////
    // generate y-component
    for( int i=0; i<nres; ++i )
        for( int j=0; j<nresp; ++j )
        {
            float kx = i>=nresp? (float)(i-nres)*kfac : (float)i*kfac;
            float ky = (float)j*kfac;
            
            float kk = sqrtf(kx*kx+ky*ky);
            
            int idx = i*nresp+j;
            
            std::cerr << "ky = " << ky << " ,  kx = " << kx << std::endl;
            
            cdata2[idx][0] = ky/kk/kk * cdata[idx][1];
            cdata2[idx][1] = -ky/kk/kk * cdata[idx][0];
        }
    
    cdata2[0][0] = 0.0f;
    cdata2[0][1] = 0.0f;
    
    fftwf_execute( iplan2 );
    
    for( int i=0; i<nres; ++i )
        for( int j=0; j<nres; ++j )
        {
            int idx = i*(nres+2)+j;
            int ii = i*nres+j;
            P[ii].y = (float)j*dx + data2[idx];
            P[ii].vy = data2[idx] * vfact;
            P[ii].acc[1] = 0.0f;
        }
    
    /////////////////////////////////
    
    delete[] data;
    delete[] data2;
    fftwf_destroy_plan(plan);
    fftwf_destroy_plan(iplan);
    fftwf_destroy_plan(plan2);
    fftwf_destroy_plan(iplan2);
    gsl_rng_free( RNG );
}


    
