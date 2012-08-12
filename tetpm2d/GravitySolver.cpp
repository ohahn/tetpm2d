//
//  GravitySolver.cpp
//  tetpm2d
//
//  Created by Oliver Hahn on 8/8/12.
//  Copyright (c) 2012 Oliver Hahn. All rights reserved.
//

#include "GravitySolver.h"

#include <fftw3.h>
#include <omp.h>

#include "Global.h"

extern int timedir;

gravity_solver::gravity_solver( unsigned n )
: n_( n )
{
    //fftwf_init_threads();
	//fftwf_plan_with_nthreads(omp_get_max_threads());
    
    data = new fftwf_real[ n_ * (n_+2) ];
    force = new fftwf_real[ n_ * (n_+2) ];
    cdata = reinterpret_cast<fftwf_complex*>(data);
    
    box_ = boxlength;
    box05_ = 0.5f * boxlength;
    
    plan  = fftwf_plan_dft_r2c_2d( n_, n_, data, cdata, FFTW_MEASURE ),
    iplan = fftwf_plan_dft_c2r_2d( n_, n_, cdata, data, FFTW_MEASURE );
    
    UnitLength_in_cm = 3.08568025e24f; //      ;  1.0 Mpc
    UnitMass_in_g    = 1.989e43f; //           ;  1.0e10 solar masses
    UnitVelocity_in_cm_per_s = 1e5f; //                ;  1 km/sec
    UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;
    GRAVITY = 6.672e-8f;
    G = GRAVITY / pow(UnitLength_in_cm, 3) * UnitMass_in_g * pow(UnitTime_in_s, 2);
    Omega_m = 1.0; //0.276;
    Omega_L = 0.0; //0.724;
    
    aforce = 0.0;
    stepno=0;
}

gravity_solver::~gravity_solver()
{
    fftwf_destroy_plan(plan);
    fftwf_destroy_plan(iplan);
    
    //fftwf_cleanup_threads();
    
    delete[] data;
    delete[] force;
}

void gravity_solver::cic_deploy( float pweight )
{
    unsigned nn2p = n_*(n_+2), nn2 = n_*n_, np=n_+2;
    
    for( unsigned i=0; i<nn2p; ++i )
        data[i] = -1.0f;
    
#define TETPM
#ifdef TETPM
    unsigned nres2 = nres*nres;
    const int vert[4][2] = {{0,0},{1,0},{0,1},{1,1}};
    const int conn[2][3] = {{0,1,2},{1,3,2}};
    
    float gfac = (float)n_ / box_;
    int             slab_x, slab_y;
	int             slab_xx, slab_yy;
    float dx,dy;
    
    pweight *= (float)n_*(float)n_/(nres*nres);//(float)nres*(float)nres/(n_*n_);
    
    for( unsigned i=0; i<nres2; ++i )
    {
        int iy,ix;
        iy = i%nres;
        ix = (i-iy)/nres;
        
        static int      square_vertices[4];
        for (int k = 0; k < 4; ++k)
            square_vertices[k] = ((ix + vert[k][0]) % nres) * nres + ((iy + vert[k][1]) % nres);
        
        for (unsigned j = 0; j < 2; ++j) {
            static float    xc[2*3];
            static int      vertids[3];

            for (int k = 0; k < 3; ++k)
                vertids[k] = square_vertices[conn[j][k]];
            
            get_triangle_approx( vertids, xc );
            
            for( unsigned k=0; k<3; ++k )
            {
                slab_x = gfac * xc[2*k+0];
                dx = gfac * xc[2*k+0] - slab_x;
                slab_x = (slab_x+n_)%n_;
                slab_xx = (slab_x+1+n_)%n_;
                
                
                slab_y = gfac * xc[2*k+1];
                dy = gfac * xc[2*k+1] - slab_y;
                slab_y = (slab_y+n_)%n_;
                slab_yy = (slab_y+1+n_)%n_;
                
                data[slab_x * np + slab_y ]   += pweight * (1.0f - dx) * (1.0f - dy);
                data[slab_xx * np + slab_y ]  += pweight * dx * (1.0f - dy);
                data[slab_x * np + slab_yy ]  += pweight * (1.0f - dx) * dy;
                data[slab_xx * np + slab_yy ] += pweight * dx * dy;
            }
        }
    }
    
#else
   
    unsigned nres2 = nres*nres;
    
    float gfac = (float)n_ / box_;
    int             slab_x, slab_y;
	int             slab_xx, slab_yy;
    float dx,dy;
    
    pweight *= 6.0f * (float)n_*(float)n_/(nres*nres);
    
    for( unsigned i=0; i<nres2; ++i )
    {
        slab_x = gfac * P[i].x;
        dx = gfac * P[i].x - slab_x;
        slab_x = (slab_x+n_)%n_;
        slab_xx = (slab_x+1+n_)%n_;
        
        
        slab_y = gfac * P[i].y;
        dy = gfac * P[i].y - slab_y;
        slab_y = (slab_y+n_)%n_;
        slab_yy = (slab_y+1+n_)%n_;
        
        data[slab_x * np + slab_y ]   += pweight * (1.0f - dx) * (1.0f - dy);
        data[slab_xx * np + slab_y ]  += pweight * dx * (1.0f - dy);
        data[slab_x * np + slab_yy ]  += pweight * (1.0f - dx) * dy;
        data[slab_xx * np + slab_yy ] += pweight * dx * dy;


    }
    
#endif
    
    
}

void gravity_solver::drift_particles( float a, float da )
{
    float H0 = 100.0f;
    float Omega_k = 1.0f - Omega_m - Omega_L;
    float a2 = a*a, a3 = a2*a;
    float adot = H0 * a * sqrt( Omega_L + Omega_k / a2 + Omega_m / a3 );
    float vfac = 1.0f / adot / a2;
    
    unsigned npart = nres*nres;
    
    for( int i=0; i<npart; ++i )
    {
        P[i].x += vfac * da * P[i].vx;
        P[i].y += vfac * da * P[i].vy;
        
        if( P[i].x >= box_ ) P[i].x -= box_;
        else if( P[i].x < 0.0f ) P[i].x += box_;
        
        if( P[i].y >= box_ ) P[i].y -= box_;
        else if( P[i].y < 0.0f ) P[i].y += box_;
        
    }
}

void gravity_solver::kick_particles( float a, float da )
{
    float H0 = 100.0f;
    float Omega_k = 1.0f - Omega_m - Omega_L;
    float a2 = a*a, a3 = a2*a;
    float adot = H0 * a * sqrt( Omega_L + Omega_k / a2 + Omega_m / a3 );
   
    unsigned npart = nres*nres;
    
    for( int i=0; i<npart; ++i )
    {
        P[i].vx += da/adot * P[i].acc[0];
        P[i].vy += da/adot * P[i].acc[1];
    }
    
}

float gravity_solver::step( float a, float da )
{
    if( timedir )
        da = -da;
    
    if( a > aforce )
    {
        cic_deploy(0.25f);
        compute_force( a );
    }
    kick_particles( a, 0.5f * da );
    
    drift_particles( a, da );
    
    cic_deploy(1.0f/6.0f);
    compute_force( a+da ); aforce = a+da;
    
    kick_particles( a + 0.5f*da, 0.5f * da );

    
    //fprintf(stderr,"step %05d : a=%f,  da=%f complete.\n", ++stepno,a+da,da);
    
    return a+da;
}


void gravity_solver::compute_force( float a )
{
    unsigned npart = nres*nres;
    float gfac = (float)n_ / box_;
    int             slab_x, slab_y;
	int             slab_xx, slab_yy;
    float dx,dy;
    
    float Omega_k = 1.0f - Omega_m - Omega_L;
    float H0 = 100.0f;
    float grav_fac = 3.0f/2.0f * H0 * H0 * Omega_m / a;
    
    
    
    float fac;
    
    // 1/n^3 * (pi/L)^2 L^2/L^3
    
    //fac = grav_fac / (M_PI * box_);	/* to get potential */
	
    fac = grav_fac / (n_*n_) / pow(2.0*M_PI/box_,2);
    fac *= 1.0 / (2 * box_ / n_);	/* for finite differencing */

    
    /*fftwf_plan plan, iplan;
    plan = fftwf_plan_dft_r2c_2d( n_, n_, data, cdata, FFTW_MEASURE);
    iplan = fftwf_plan_dft_c2r_2d( n_, n_, cdata, data, FFTW_MEASURE);
    */
    
    fftwf_execute(plan);
    
    float kfac = M_PI/n_, kkx, kky, kkk, kernel;
    int kx, ky, kz, kk, idx;
    int np = n_+2, n05 = n_/2, npp = n_/2+1;
    
    for( int ix=0; ix<n_; ix++ )
        for( int iy=0; iy<npp; iy++ )
        {
            if( ix > n05 )
                kx = ix - n_;
            else
                kx = ix;
            
            if( iy > n05 )
                ky = iy - n_;
            else
                ky = iy;
            
            kk = kx*kx + ky *ky;
    
            if( kk > 0 )
            {
                //kkx = kx * kfac;
                //kky = ky * kfac;
                
                kernel = -1.0f / ((float)kk);
                
                idx = ix * npp + iy;
                
                cdata[idx][0] *= kernel;
                cdata[idx][1] *= kernel;
                
            }
        }
    
    cdata[0][0] = 0.0f;
    cdata[0][1] = 0.0f;
    
    fftwf_execute(iplan);
    
    int xrr,xr,xl,xll,yrr,yr,yl,yll;
    
    for( int dim=0; dim < 2; ++dim )
    {
        for( int ix=0; ix<n_; ix++ )
            for( int iy=0; iy<n_; iy++ )
            {
                xrr = xll = xr = xl = ix;
                yrr = yll = yr = yl = iy;
                
                switch (dim) {
					case 0:
						xr = (ix + 1 + n_)%n_;
						xrr = (ix + 2 + n_)%n_;
						xl = (ix - 1 + n_)%n_;
						xll = (ix - 2 + n_)%n_;
						break;
					case 1:
						yr = (iy + 1 + n_)%n_;
						yl = (iy - 1 + n_)%n_;
						yrr = (iy + 2 + n_)%n_;
						yll = (iy - 2 + n_)%n_;
						break;
                }
                
                
                force[ ix*np +iy ] = fac * ((4.0/3.0) *
                                            (data[xl*np+yl] - data[xr*np+yr])
                                            -(1.0/6.0) *
                                            (data[xll*np+yll] - data[xrr*np+yrr]));
            }
        
        
        for( int i=0; i<npart; ++i )
        {
            slab_x = gfac * P[i].x;
            dx = gfac * P[i].x - slab_x;
            slab_x = (slab_x+n_)%n_;
            slab_xx = (slab_x+1+n_)%n_;
            
            
            slab_y = gfac * P[i].y;
            dy = gfac * P[i].y - slab_y;
            slab_y = (slab_y+n_)%n_;
            slab_yy = (slab_y+1+n_)%n_;
            
            float acc;
            
            acc  = force[ slab_x*np+ slab_y ] * (1.0 - dx) * (1.0 - dy);
            acc += force[ slab_xx*np+ slab_y ] * dx * (1.0 - dy);
            acc += force[ slab_x*np+ slab_yy ] * (1.0 - dx) * dy;
            acc += force[ slab_xx*np+ slab_yy ] * dx * dy;
            
            P[i].acc[dim] = acc;
        }
                                            
        
        
    }
    
    /*fftwf_destroy_plan( plan );
    fftwf_destroy_plan( iplan );*/
    
}
