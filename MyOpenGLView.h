//
//  MyOpenGLView.h
//  tetpm2d
//
//  Created by Oliver Hahn on 8/8/12.
//  Copyright (c) 2012 Oliver Hahn. All rights reserved.
//
#import <Cocoa/Cocoa.h>
#import <Foundation/Foundation.h>

#define FRAME_DURATION_HISTORY 15

@interface MyOpenGLView : NSOpenGLView
{
    NSTimeInterval   _lastFrameTime;
    NSTimeInterval   _frameDurations[FRAME_DURATION_HISTORY];
    unsigned         _nextFramePosition;

    BOOL             _firstFrame;

    NSOpenGLContext *_windowedContext;
}

- (void) drawRect: (NSRect) bounds;

@end
