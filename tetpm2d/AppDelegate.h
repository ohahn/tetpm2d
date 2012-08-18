//
//  AppDelegate.h
//  tetpm2d
//
//  Created by Oliver Hahn on 8/8/12.
//  Copyright (c) 2012 Oliver Hahn. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@interface AppDelegate : NSObject <NSApplicationDelegate>
{
    
    IBOutlet id timedir_button;
    IBOutlet NSSlider *numpart_slider;
    IBOutlet NSSlider *numcell_slider;
    
    
@public
    IBOutlet NSButton *run_button;
    IBOutlet NSTextField *zcurrlabel;
    IBOutlet NSTextField *numpart;
    IBOutlet NSTextField *numcell;
    IBOutlet NSPanel *SimTypeSheet;
    IBOutlet NSPopUpButtonCell *simtype;
    
}

@property (assign) IBOutlet NSWindow *window;

- (IBAction)toggle_run:(id)sender;
- (IBAction)toggle_timedir:(id)sender;
- (IBAction)restart_run:(id)sender;
- (IBAction)change_numpart:(id)sender;
- (IBAction)change_numcell:(id)sender;
- (IBAction)toggle_rendermode:(id)sender;
- (IBAction)toggle_massdeploymode:(id)sender;
- (IBAction)showTheSheet:(id)sender;
- (IBAction)endTheSheet:(id)sender;

@end
