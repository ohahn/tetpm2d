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
extern int new_nres;
extern int changed_particle_number;
extern int new_PMres;
extern int changed_PMres;
extern int render_mode;
extern int mass_deploy_mode;
extern int sim_type;

@implementation AppDelegate

@synthesize window = _window;

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

- (IBAction)toggle_rendermode:(id)sender
{
    render_mode = !render_mode;
}

- (IBAction)change_numpart:(id)sender
{
    NSSlider *slider = (NSSlider *)sender;
    new_nres = (int)(pow(2,slider.doubleValue)+0.5);
    
    [numpart setStringValue:[NSString stringWithFormat:@"%d",new_nres ]];
    changed_particle_number = 1;
}

- (IBAction)change_numcell:(id)sender
{
    NSSlider *slider = (NSSlider *)sender;
    new_PMres = (int)(pow(2,slider.doubleValue)+0.5);
    
    [numcell setStringValue:[NSString stringWithFormat:@"%d",new_PMres ]];
    changed_PMres = 1;
}

- (IBAction)toggle_massdeploymode:(id)sender
{
    switch ([[sender selectedCell] tag]) {
        case 0:
            mass_deploy_mode = 0;
            break;
            
        case 1:
            mass_deploy_mode = 1;
            break;
            
        case 2:
            mass_deploy_mode = 2;
            break;
            
        default:
            break;
    }
}


- (IBAction) showTheSheet:(id)sender {
    
    [NSApp beginSheet:SimTypeSheet
       modalForWindow:(NSWindow *)_window
        modalDelegate:self
       didEndSelector:nil
          contextInfo:nil];
    
}

-(IBAction)endTheSheet:(id)sender {
    
    NSInteger tag = [[simtype selectedItem] tag];
    
    sim_type = (int)tag;
    restart_sim = 1;
    
    fprintf(stderr,"set simtype to %d\n",tag);
    
    
    [NSApp endSheet:SimTypeSheet];
    [SimTypeSheet orderOut:sender];
    
}

@end
