//
//  AppDelegate.m
//  tetpm2d
//
//  Created by Oliver Hahn on 8/8/12.
//  Copyright (c) 2012 Oliver Hahn. All rights reserved.
//

#import "AppDelegate.h"

extern int is_running;
extern int timedir;
extern int restart_sim;

@implementation AppDelegate

- (void)dealloc
{
    [super dealloc];
}

- (void)applicationDidFinishLaunching:(NSNotification *)aNotification
{
    // Insert code here to initialize your application
}

- (IBAction)toggle_run:(id)sender
{
    is_running = !is_running;
}

- (IBAction)toggle_timedir:(id)sender
{
    timedir = !timedir;
}

- (IBAction)restart_run:(id)sender
{
    restart_sim = 1;
}

@end
