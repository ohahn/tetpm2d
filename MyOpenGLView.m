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

static int stepcount = 0;

int nres = 32;

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
            
            P[id].x = i * boxlength;
            P[id].y = j * boxlength;
            P[id].vx = 0.0f;
            P[id].vy = 0.0f;
            P[id].id = id;
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

void drawAnObject ( NSOpenGLContext *glcontext  )
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    //glLoadIdentity();
    
    const int vert[4][3] = {{0,0},{1,0},{0,1},{1,1}};
    const int conn[2][3] = {{0,1,2},{1,3,2}};
    
    unsigned npar = nres*nres;
    
    
    float xc[][3] = {{0.0,0.6,0.0},{-0.2,-0.3,0.0},{0.2,-0.3,0.0}};
    float v[][3] = {{1.0,0.0,0.0},{1.0,0.0,0.0},{1.0,0.1,0.0}};
    float xp[9][3];
    float dt = 1.0;
    
    // update pos
    /*for( int i=0; i<3; ++i ) // all vertices
        for( int j=0; j<3; ++j ) // all coords
            xp[i][j] = xc[i][j] + (float)stepcount * dt * v[i][j];
    */
    
    
    float phi1 = fmod( (float)stepcount *1e-2, 2.0*M_PI );
    float phi2 = fmod( -(float)stepcount *1e-2, 2.0*M_PI );
    float phi3 = fmod( -1.2*(float)stepcount *1e-2, 2.0*M_PI );

    
    float twopithirds = 2.0*M_PI/3.0;
    for( int i=0; i<3; ++i ) // all vertices
    {
        xp[i][0] = cos( phi1 + (float)i*twopithirds);
        xp[i][1] = sin( phi1 + (float)i*twopithirds);
        xp[i][2] = 0.0f;
        
        
        xp[3+i][0] = cos( phi2 + (float)i*twopithirds);
        xp[3+i][1] = sin( phi2 + (float)i*twopithirds);
        xp[3+i][2] = 0.0f;
        
        xp[6+i][0] = cos( phi3 + (float)i*twopithirds);
        xp[6+i][1] = sin( phi3 + (float)i*twopithirds);
        xp[6+i][2] = 0.0f;
    }
    
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);     // Turn Blending On
    
    glBlendFunc(GL_SRC_ALPHA,GL_ONE);
    //glBlendFunc(GL_ONE,GL_ONE);
    
    {
        //glColor4f(1.0f, 1.0f, 1.0f,0.5f);
        glColor4f(0.5f, 0.5f, 0.5f,0.5f);
   
        glBegin(GL_TRIANGLES);
        {
            for( int j=0; j<3; ++j )
            for( int i=0; i<3; ++i )
                glVertex3f( xp[3*j+i][0], xp[3*j+i][1], xp[3*j+i][2] );
        }
        glEnd();
    }
    //glSwapAPPLE();
    //[glcontext flushBuffer];
    //glFlush();
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
