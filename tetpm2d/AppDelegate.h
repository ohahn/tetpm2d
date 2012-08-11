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
    
@public
    IBOutlet NSButton *run_button;
    IBOutlet NSTextField *zcurrlabel;
}

@property (assign) IBOutlet NSWindow *window;

- (IBAction)toggle_run:(id)sender;
- (IBAction)toggle_timedir:(id)sender;
- (IBAction)restart_run:(id)sender;

@end
