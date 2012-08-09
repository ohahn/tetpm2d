//
//  MyOpenGLView.m
//  tetpm2d
//
//  Created by Oliver Hahn on 8/8/12.
//  Copyright (c) 2012 Oliver Hahn. All rights reserved.
//

#include <math.h>
#include <OpenGL/gl.h>
#import "MyOpenGLView.h"
#include "CosmologyWrapper.h"

static int stepcount = 0;

int nres = 128;

const float boxlength = 10.0;

typedef struct 
{
    float x,y,vx,vy;
    unsigned id;
}particle;

particle *P;


void create_particles()
{
    for( int i=0; i<nres; ++i )
        for( int j=0; j<nres; ++j )
        {
            unsigned id = i*nres+j;
            
            P[id].x = (float)j * boxlength/(float)nres + boxlength*0.05 * ((double)rand()/RAND_MAX-0.5);
            P[id].y = (float)i * boxlength/(float)nres + boxlength*0.05 * ((double)rand()/RAND_MAX-0.5);
            P[id].vx = 0.0f;
            P[id].vy = 0.0f;
            P[id].id = id;
            
            power_getk( (float) id );
        }
}


@implementation MyOpenGLView

- (void)awakeFromNib
{
    NSOpenGLPixelFormatAttribute WindowedAttributes[] =
    {
        NSOpenGLPFADoubleBuffer,
        NSOpenGLPFAAccelerated,
        NSOpenGLPFADepthSize, 24,
        NSOpenGLPFAStencilSize, 8,
        NSOpenGLPFASingleRenderer,
        NSOpenGLPFAScreenMask, CGDisplayIDToOpenGLDisplayMask(kCGDirectMainDisplay),
        NSOpenGLPFANoRecovery,
        0
    };
    
    NSOpenGLPixelFormat* windowedPixelFormat = [[[NSOpenGLPixelFormat alloc]
                                                 initWithAttributes:WindowedAttributes] autorelease];
    
    _windowedContext = [[NSOpenGLContext alloc] initWithFormat:windowedPixelFormat
                                                  shareContext:nil];
    if (_windowedContext == nil)
    {
        NSLog(@"Got nil windowed context");
        [self dealloc];
        return;
    }
    
    GLint vsync = 1;
    [_windowedContext setValues:&vsync forParameter:NSOpenGLCPSwapInterval];
    
    _firstFrame = YES;
    
    unsigned i;
    for (i = 0; i < FRAME_DURATION_HISTORY; ++i)
    {
        _frameDurations[i] = 0.0;
    }
    _nextFramePosition = 0;
    
    (void)[NSTimer scheduledTimerWithTimeInterval:0.0
                                           target:self
                                         selector:@selector(drawFrame)
                                         userInfo:nil
                                          repeats:YES];
    
    P = (particle *)malloc( sizeof(particle) *nres*nres );
    create_particles();
}

- (void)dealloc
{
    [_windowedContext release];
    free(P);
    [super dealloc];
}

float get_vertices( int itr, int i, float box, float box05, float *trvert )
{
    //const int numtr = 2;
    
    const int vert[4][2] = {{0,0},{1,0},{0,1},{1,1}};
    const int conn[2][3] = {{0,1,2},{1,3,2}};
    
    float x,y;
    x = P[i].x;
    y = P[i].y;
    
    int ix,iy;
    iy = i % nres;
    ix = (i-iy)/nres;
    
    float dx[3], dy[3];
    
    for( int k=0; k<3; ++k )
    {
        int idxp = ((ix+vert[ conn[itr][k] ][0])%nres)*nres + (iy+vert[ conn[itr][k] ][1])%nres;
        
        dx[k] = P[idxp].x - x;
        dy[k] = P[idxp].y - y;
        
        if( dx[k] < -box05 ) dx[k] += boxlength;
        else if( dx[k] > box05 ) dx[k] -= boxlength;
        if( dy[k] < -box05 ) dy[k] += boxlength;
        else if( dy[k] > box05 ) dy[k] -= boxlength;
        
        trvert[2*k+0] = x+dx[k];
        trvert[2*k+1] = y+dy[k];
    }
    
    float area = 0.5f*(dx[1]*dy[2]-dy[1]*dx[2]);
    return area;
}

void drawAnObject ( NSOpenGLContext *glcontext  )
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   
    unsigned npar = nres*nres;
    float box05 = 0.5f * boxlength;
    
    //glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);     // Turn Blending On
    
    
    const int numtr = 2;
    
    float meanA = 0.5f*boxlength*boxlength/(nres*nres);
    
    float trvert[3*2];

    glBegin(GL_TRIANGLES);
    {
        for( unsigned i=0; i<npar; ++i )
        {
            for( int j=0; j<numtr; ++j )
            {
                float A;
                A = get_vertices( j, i, boxlength, box05, trvert );
                
                float dens  = meanA/fabs(A);
                float col = dens/2.0;///4.0;
                
                if( col > 1.0f ) col = 1.0f;
                
                glColor4f(col, col, col,0.9f);
                glBlendFunc(GL_ONE,GL_ONE);
                
                for( int k=0; k<3; ++k )
                {
                    float dx,dy;
                    dx = 2.0f * (trvert[2*k+0]/boxlength - 0.5f);
                    dy = 2.0f * (trvert[2*k+1]/boxlength - 0.5f);
                    
                    glVertex3f( dx, dy, 0.0f );
                }
            }
        }
        
    }
    glEnd();
    

    float dt  = 1e-3;
    
    P[126].vx = 1.0;
    P[126].vy = 1.0;
    
    
    
    for( unsigned i=0; i<npar; ++i )
    {
        P[i].x += dt * P[i].vx;
        P[i].y += dt * P[i].vy;
    }
    

    glFinish();
}

- (void)drawFrame
{
    /*if (![[[NSRunLoop currentRunLoop] currentMode] isEqualToString:NSDefaultRunLoopMode])
    {
        return;
    }*/
    
    NSOpenGLContext* context;
    
 
    context = _windowedContext;
    if (_firstFrame)
    {
        [_windowedContext setView:self];
    }
    
    if (_firstFrame)
    {
        [context makeCurrentContext];
        NSSize contextSize;
        contextSize = [self bounds].size;
        _firstFrame = NO;
        
        //glClearColor(0, 0, 0, 0);
        //glClear(GL_COLOR_BUFFER_BIT);
        
        glDisable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);     // Turn Blending On
        
        glBlendFunc(GL_SRC_ALPHA,GL_ONE);
        /* glBlendFunc(GL_DST_ALPHA,GL_ONE);*/
        //glBlendFunc(GL_ONE_MINUS_SRC_COLOR,GL_ONE);
        
    }
    
    drawAnObject( [self openGLContext] );
    stepcount++;
    
    
    
    [context flushBuffer];
}

-(void) drawRect: (NSRect) bounds
{
    [self drawFrame];
    
    /*glClearColor(0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT);
    drawAnObject();
    glFlush();*/
}

@end
