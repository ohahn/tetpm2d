//
//  Global.h
//  tetpm2d
//
//  Created by Oliver Hahn on 8/9/12.
//  Copyright (c) 2012 Oliver Hahn. All rights reserved.
//

#ifndef tetpm2d_Global_h
#define tetpm2d_Global_h

typedef struct
{
    float x,y,vx,vy;
    float acc[2];
    unsigned id;
}particle;

extern particle *P;
extern int nres;
extern float boxlength;


#endif
