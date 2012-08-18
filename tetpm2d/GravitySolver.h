//
//  GravitySolver.h
//  tetpm2d
//
//  Created by Oliver Hahn on 8/8/12.
//  Copyright (c) 2012 Oliver Hahn. All rights reserved.
//

#ifndef __tetpm2d__GravitySolver__
#define __tetpm2d__GravitySolver__

#include <iostream>
#include <cmath>
#include "fftw3.h"
#include "Global.h"


class gravity_solver
{
public:
    float UnitLength_in_cm;
    float UnitMass_in_g; //           ;  1.0e10 solar masses
    float UnitVelocity_in_cm_per_s; //                ;  1 km/sec
    float UnitTime_in_s;
    float GRAVITY;
    float G;
    float Omega_m, Omega_L;
    
    float aforce;
    
    int stepno;
    
protected:
    fftwf_plan plan, iplan;
    fftwf_real *data;
    fftwf_real *force;
    fftwf_complex *cdata;
    
    float box_, box05_;
    unsigned n_;
    
    inline void
    get_triangle_approx( const int *vertids, float *xpp )
    {
        const float     x0[2] = {P[vertids[0]].x, P[vertids[0]].y};
        static float xc[2];
        
        //static const float fac_point = 0.4472135955f; // = sqrt(5.0)/5.0;
        //static const float fac_centr = 0.5527864045f; // = 1.0 - sqrt(5.0)/5.0;
        
        const float fac_point = 0.5f;
        const float fac_centr = 0.5f;
        
        int             i;
        
        static float xp[3*2];
        
        xc[0] = xc[1] = 0.0f;
        xp[0] = xp[1] = 0.0f;
        
        for (i = 1; i < 3; ++i) {
            float           dx[] = {P[vertids[i]].x - x0[0], P[vertids[i]].y - x0[1]};
            
            
            if (dx[0] < -box05_)
                dx[0] += box_;
            else if (dx[0] > box05_)
                dx[0] -= box_;
            if (dx[1] < -box05_)
                dx[1] += box_;
            else if (dx[1] > box05_)
                dx[1] -= box_;
            
            xc[0] += dx[0];
            xc[1] += dx[1];
            
            xp[2*i+0] = dx[0];
            xp[2*i+1] = dx[1];
            
        }
        
        xc[0] *= 0.3333333333f;
        xc[1] *= 0.3333333333f;
        
        // xc and xp is relative to x0
        //xxp0 = 0.5*(xc+xp0)+x0;
        
        for( i=0; i<3; ++i )
        {
            xpp[2*i+0] = fmodf( ( fac_centr * xc[0] + fac_point * xp[2*i+0] ) + x0[0] + box_, box_ );
            xpp[2*i+1] = fmodf( ( fac_centr * xc[1] + fac_point * xp[2*i+1] ) + x0[1] + box_, box_ );
        }
    }
    
    inline void
    get_triangle_centroid( const int *vertids, float *xpp )
    {
        const float     x0[2] = {P[vertids[0]].x, P[vertids[0]].y};
        static float xc[2];
        
        int             i;
        
        static float xp[3*2];
        
        xc[0] = xc[1] = 0.0f;
        xp[0] = xp[1] = 0.0f;
        
        for (i = 1; i < 3; ++i) {
            float           dx[] = {P[vertids[i]].x - x0[0], P[vertids[i]].y - x0[1]};
            
            
            if (dx[0] < -box05_)
                dx[0] += box_;
            else if (dx[0] > box05_)
                dx[0] -= box_;
            if (dx[1] < -box05_)
                dx[1] += box_;
            else if (dx[1] > box05_)
                dx[1] -= box_;
            
            xc[0] += dx[0];
            xc[1] += dx[1];
            
            xp[2*i+0] = dx[0];
            xp[2*i+1] = dx[1];
            
        }
        
        xc[0] *= 0.3333333333f;
        xc[1] *= 0.3333333333f;
        
        // xc and xp is relative to x0
        //xxp0 = 0.5*(xc+xp0)+x0;
        
        xpp[0] = fmodf( xc[0] + x0[0] + box_, box_ );
        xpp[1] = fmodf( xc[1] + x0[1] + box_, box_ );
    }
    
    
protected:
    
    void deploy_tet3pm( void );
    void deploy_stdpm( void );
    void deploy_tcm( void );
    
public:
    
    explicit gravity_solver( unsigned nres );
    ~gravity_solver();

    void cic_deploy( int mass_deploy_method );
    void compute_force( float a );
    void kick_particles( float a, float da );
    void drift_particles( float a, float da );
    
    void reset_force( void )
    {
        aforce = 0.0f;
    }
    
    void change_PM_resolution( int newres )
    {
        delete[] data;
        delete[] force;
        
        n_ = newres;
        
        data = new fftwf_real[ n_ * (n_+2) ];
        force = new fftwf_real[ n_ * (n_+2) ];
        cdata = reinterpret_cast<fftwf_complex*>(data);
        
        fftwf_destroy_plan(plan);
        fftwf_destroy_plan(iplan);
        
        
        plan  = fftwf_plan_dft_r2c_2d( n_, n_, data, cdata, FFTW_MEASURE ),
        iplan = fftwf_plan_dft_c2r_2d( n_, n_, cdata, data, FFTW_MEASURE );
        
        aforce = 0.0;
    }
    
    float step( float a, float da, int mass_deploy_mode );
};


#endif /* defined(__tetpm2d__GravitySolver__) */
