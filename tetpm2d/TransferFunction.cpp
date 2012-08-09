//
//  TransferFunction.cpp
//  tetpm2d
//
//  Created by Oliver Hahn on 8/8/12.
//  Copyright (c) 2012 Oliver Hahn. All rights reserved.
//

#include "TransferFunction.hh"

/* Convenience from Numerical Recipes in C, 2nd edition */
double sqrarg;
double cubearg;
double pow4arg;


/**********/
TransferFunction_Eisenstein *ptf = NULL;

double power_getk( double k )
{
    
    if( ptf == NULL )
        ptf = new TransferFunction_Eisenstein();
    
    
    return ptf->compute(k);
}
    
