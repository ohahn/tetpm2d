//
//  MyOpenGLView.h
//  tetpm2d
//
//  Created by Oliver Hahn on 8/8/12.
//  Copyright (c) 2012 Oliver Hahn. All rights reserved.
//
#import <Cocoa/Cocoa.h>
#import <Foundation/Foundation.h>
#include "GravitySolver.h"

#define FRAME_DURATION_HISTORY 15

@interface MyOpenGLView : NSOpenGLView
{
    NSTimeInterval   _lastFrameTime;
    NSTimeInterval   _frameDurations[FRAME_DURATION_HISTORY];
    unsigned         _nextFramePosition;

    BOOL             _firstFrame;

    NSOpenGLContext *_windowedContext;
    
    // Used variables and types
    GLuint fbo; // frame buffer object
    GLuint depthbuffer; // depth buffer object
    GLuint img; // FBO "texture"
    GLuint rbo; // render buffer object
    
    GLuint program_postproc, attribute_v_coord_postproc, uniform_fbo_texture;
    GLuint fbo_texture, rbo_depth;
    GLuint vbo_fbo_vertices;
    GLuint vs,fs;
    
    gravity_solver *pgsolve;
}

- (void) drawRect: (NSRect) bounds;
- (void) drawAnObject : (NSOpenGLContext *)glcontext;
- (IBAction)toggle_run:(id)sender;



@end
