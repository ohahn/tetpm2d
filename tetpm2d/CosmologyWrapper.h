//
//  CosmologyWrapper.h
//  tetpm2d
//
//  Created by Oliver Hahn on 8/8/12.
//  Copyright (c) 2012 Oliver Hahn. All rights reserved.
//

#ifndef tetpm2d_CosmologyWrapper_h
#define tetpm2d_CosmologyWrapper_h

void cosmo_set_parameters( double zstart, double Omega_m, double Omega_b, double Omega_L, double n_s, double sigma8, double h );
double cosmo_get_amp_k( double k );
void cosmo_init_particles( unsigned seed);

#endif
