//
//  MyOpenGLView.m
//  tetpm2d
//
//  Created by Oliver Hahn on 8/8/12.
//  Copyright (c) 2012 Oliver Hahn. All rights reserved.
//

#include <math.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>

#import "MyOpenGLView.h"
#import "AppDelegate.h"
#include "CosmologyWrapper.h"

static int stepcount = 0;

int is_running = 0;
int timedir = 0;
int restart_sim = 0;
int changed_particle_number = 0;
int render_mode = 0;
int mass_deploy_mode = 0;

int nres = 128;
int new_nres = nres;

int PMres = 128;
int new_PMres = PMres;
int changed_PMres = 0;

float boxlength = 10.0;
float aexp;
float astart = 1.0/101.0;

particle *P;

void ripple_wave( void )
{
    
    double epsilon_a = 0.2;//1.0;//0.2;
	double kp = 1.0;
	double ka = 2.0;
    
    
    double zcross = 7.7-1.0; // expansion factor of 7.7 after shell crossing as in Melott et al. (1997)
    double across = 1.0/(1.0+zcross);
    double waveamp = 1.0/across;// * pow((double)nres,1.5);
    
    double vfact = 100.0 * pow(astart,-1.5) *astart*astart;// * sqrt(astart);
    
    double Aini = waveamp * astart * boxlength;
    double twopi = 2.0*M_PI;
    
    for( int i=0; i<nres; ++i )
        for( int j=0; j<nres; ++j )
        {
            unsigned id = i*nres+j;
            P[id].x = (float)i * boxlength/(float)nres;
            P[id].y = (float)j * boxlength/(float)nres;
            
            P[id].vx = 0.0f;
            P[id].vy = 0.0f;
            
            P[id].acc[0] = 0.0f;
            P[id].acc[1] = 0.0f;
            P[id].id = id;
         
            
            //... perturbation
            float
            qx = ((float)i+0.5)/nres,
            qy = ((float)j+0.5)/nres;
            
            double vy = Aini/(twopi*kp)*sin(twopi*kp*qy + epsilon_a*(kp*kp)/(ka*ka)*cos(twopi*ka*qx));
            double vx = -Aini/(twopi*ka)*epsilon_a*sin(twopi*kp*qy + epsilon_a*(kp*kp)/(ka*ka)*cos(twopi*ka*qx))*sin(twopi*ka*qx);
            P[id].x += vx;
            P[id].y += vy;
            P[id].vx += vfact * vx;
            P[id].vy += vfact * vy;
            
        }
}

void plane_wave( void )
{
    double zcross = 7.7-1.0; // expansion factor of 7.7 after shell crossing as in Melott et al. (1997)
    double across = 1.0/(1.0+zcross);
    double waveamp = 1.0/across;// * pow((double)nres,1.5);
    
    double vfact = 100.0 * pow(astart,-1.5) *astart*astart;// * sqrt(astart);
    
    for( int i=0; i<nres; ++i )
        for( int j=0; j<nres; ++j )
        {
            unsigned id = i*nres+j;
            P[id].x = (float)j * boxlength/(float)nres;
            P[id].y = (float)i * boxlength/(float)nres;
            
            P[id].vx = 0.0f;
            P[id].vy = 0.0f;
            
            P[id].x += astart * waveamp * sin( 2.0*M_PI *((float)j/(float)nres) );
            P[id].vx += astart * waveamp * vfact * sin( 2.0*M_PI *((float)j/(float)nres) );
            
            P[id].acc[0] = 0.0f;
            P[id].acc[1] = 0.0f;
            P[id].id = id;
        }
    
    
}


void init_particles()
{
    //plane_wave();
    ripple_wave();
}

/**
 * Store all the file's contents in memory, useful to pass shaders
 * source code to OpenGL
 */
char* file_read(const char* filename)
{
    FILE* in = fopen(filename, "rb");
    if (in == NULL) return NULL;
    
    int res_size = BUFSIZ;
    char* res = (char*)malloc(res_size);
    int nb_read_total = 0;
    
    while (!feof(in) && !ferror(in)) {
        if (nb_read_total + BUFSIZ > res_size) {
            if (res_size > 10*1024*1024) break;
            res_size = res_size * 2;
            res = (char*)realloc(res, res_size);
        }
        char* p_res = res + nb_read_total;
        nb_read_total += fread(p_res, 1, BUFSIZ, in);
    }
    
    fclose(in);
    res = (char*)realloc(res, nb_read_total + 1);
    res[nb_read_total] = '\0';
    return res;
}

/**
 * Display compilation errors from the OpenGL shader compiler
 */
void print_log(GLuint object)
{
    GLint log_length = 0;
    if (glIsShader(object))
        glGetShaderiv(object, GL_INFO_LOG_LENGTH, &log_length);
    else if (glIsProgram(object))
        glGetProgramiv(object, GL_INFO_LOG_LENGTH, &log_length);
    else {
        fprintf(stderr, "printlog: Not a shader or a program\n");
        return;
    }
    
    char* log = (char*)malloc(log_length);
    
    if (glIsShader(object))
        glGetShaderInfoLog(object, log_length, NULL, log);
    else if (glIsProgram(object))
        glGetProgramInfoLog(object, log_length, NULL, log);
    
    fprintf(stderr, "%s", log);
    free(log);
}

/**
 * Compile the shader from file 'filename', with error handling
 */
GLuint create_shader(const char* filename, GLenum type)
{
    const GLchar* source = file_read(filename);
    if (source == NULL) {
        fprintf(stderr, "Error opening %s: ", filename); perror("");
        return 0;
    }
    
    GLuint res = glCreateShader(type);
    const GLchar* sources[2] = {
#ifdef GL_ES_VERSION_2_0
        "#version 100\n"
        "#define GLES2\n",
#else
        "#version 120\n",
#endif
        source };
    glShaderSource(res, 2, sources, NULL);
    
    free((void*)source);
    
    glCompileShader(res);
    GLint compile_ok = GL_FALSE;
    glGetShaderiv(res, GL_COMPILE_STATUS, &compile_ok);
    if (compile_ok == GL_FALSE) {
        fprintf(stderr, "%s:", filename);
        print_log(res);
        glDeleteShader(res);
        return 0;
    }
    
    
    
    return res;
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
    
    
    aexp = 1.0f/101.0f;
    
    init_particles();
    
    pgsolve = new gravity_solver( PMres );
    
    is_running = 0;
    
}

- (void)dealloc
{
    /* free_resources */
    glDeleteRenderbuffers(1, &rbo_depth);
    glDeleteTextures(1, &fbo_texture);
    glDeleteFramebuffers(1, &fbo);
    
    glDeleteBuffers(1, &vbo_fbo_vertices);
    glDeleteProgram(program_postproc);
    
    [_windowedContext release];
    free(P);
    [super dealloc];
    
    delete pgsolve;
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

#define B2S(x) (2.0f*((x)/boxlength-0.5f))

void drawAnObject( int rendermode )
{
    
    
    //... render
    
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   
    unsigned npar = nres*nres;
    float box05 = 0.5f * boxlength;
    
    //glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);     // Turn Blending On
    glBlendFunc(GL_ONE,GL_ONE);
    
    
    const int numtr = 2;
    
    float meanA = 0.5f*boxlength*boxlength/(nres*nres);
    
    float trvert[3*2];

    if( rendermode == 0 )
    {
    
        glBegin(GL_TRIANGLES);
        {
            for( unsigned i=0; i<npar; ++i )
            {
                for( int j=0; j<numtr; ++j )
                {
                    float A;
                    A = get_vertices( j, i, boxlength, box05, trvert );
                    
                    float dens  = meanA/fabs(A);
                    float col = dens/8.0;///4.0;
                    
                    if( col > 1.0f ) col = 1.0f;
                    
                    glColor4f(col, col, col,0.1f);
                    glBlendFunc(GL_ONE,GL_ONE);
                    
                    if( A< 0)
                    for( int k=0; k<3; ++k )
                    {
                        float dx,dy;
                        dx = B2S(trvert[2*k+0]);
                        dy = B2S(trvert[2*k+1]);
                        
                        glVertex3f( dx, dy, 0.0f );
                    }
                    else
                    {
                        float dx,dy;
                        
                        
                        dx = B2S(trvert[2*2+0]);
                        dy = B2S(trvert[2*2+1]);
                        
                        glVertex3f( dx, dy, 0.0f );
                        
                        dx = B2S(trvert[2*1+0]);
                        dy = B2S(trvert[2*1+1]);
                        
                        glVertex3f( dx, dy, 0.0f );
                        
                        dx = B2S(trvert[2*0+0]);
                        dy = B2S(trvert[2*0+1]);
                        
                        glVertex3f( dx, dy, 0.0f );
                    }
                }
            }

            
        }
        glEnd();
        
    }else{
    
        glBegin(GL_TRIANGLES);
        {
            for( unsigned i=0; i<npar; ++i )
            {
                glColor4f(0.8,0.8,0.8,0.2);
                glBlendFunc(GL_ONE,GL_ONE);
                
                glVertex3f( B2S(P[i].x)-0.002, B2S(P[i].y), 0.0f );
                glVertex3f( B2S(P[i].x)+0.002, B2S(P[i].y), 0.0f );
                glVertex3f( B2S(P[i].x), B2S(P[i].y)+0.004, 0.0f );
                
            }
            
            
        }
        glEnd();
    }
    
    /*float dt  = 1e-3;
    
    P[126].vx = 1.0;
    P[126].vy = 1.0;
    
    
    
    for( unsigned i=0; i<npar; ++i )
    {
        P[i].x += dt * P[i].vx;
        P[i].y += dt * P[i].vy;
    }
    */
    
    

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
        _firstFrame = NO;
        
        [context makeCurrentContext];
        NSSize contextSize;
        contextSize = [self bounds].size;
        
        /* init_resources */
        
        
        //... compile postproc shaders
        /* init_resources */
        /* Post-processing */
        if ((vs = create_shader("/Users/ohahn/Sources/tetpm2d/tetpm2d/postproc_v.glsl", GL_VERTEX_SHADER))   == 0) return;
        if ((fs = create_shader("/Users/ohahn/Sources/tetpm2d/tetpm2d/postproc_f.glsl", GL_FRAGMENT_SHADER)) == 0) return;
        
        GLint link_ok = GL_FALSE, validate_ok = GL_FALSE;
        program_postproc = glCreateProgram();
        glAttachShader(program_postproc, vs);
        glAttachShader(program_postproc, fs);
        glLinkProgram(program_postproc);
        glGetProgramiv(program_postproc, GL_LINK_STATUS, &link_ok);
        if (!link_ok) {
            fprintf(stderr, "glLinkProgram:");
            print_log(program_postproc);
            abort();
        }
        glValidateProgram(program_postproc);
        glGetProgramiv(program_postproc, GL_VALIDATE_STATUS, &validate_ok);
        if (!validate_ok) {
            fprintf(stderr, "glValidateProgram:");
            print_log(program_postproc);
        }
        
        const GLchar attribute_name[] = "v_coord";
        attribute_v_coord_postproc = glGetAttribLocation(program_postproc, attribute_name);
        if (attribute_v_coord_postproc == -1) {
            fprintf(stderr, "Could not bind attribute %s\n", attribute_name);
            abort();
        }
        
        const GLchar uniform_name[] = "fbo_texture";
        uniform_fbo_texture = glGetUniformLocation(program_postproc, uniform_name);
        if (uniform_fbo_texture == -1) {
            fprintf(stderr, "Could not bind uniform %s\n", uniform_name);
            abort();
        }
        
        
        /* Create back-buffer, used for post-processing */
        //NSSize contextSize;
        //contextSize = [self bounds].size;
        
        GLuint screen_width = contextSize.width;
        GLuint screen_height = contextSize.height;

      
#if 0
        glGenTextures(1, &img);
        glBindTexture(GL_TEXTURE_2D, img);
        glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, screen_width, screen_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0); // place multisampling here too!
        glBindTexture(GL_TEXTURE_2D, 0);
        
        
        glGenFramebuffersEXT(1, &fbo);
        glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo);
        glGenRenderbuffersEXT(1, &rbo);
        glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, rbo);
        glRenderbufferStorageEXT(GL_RENDERBUFFER_EXT, GL_DEPTH_COMPONENT, screen_width, screen_height); // one nice place to put multisampling!
        glBindRenderbufferEXT(GL_RENDERBUFFER_EXT, 0);
        
        glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, img, 0);
        
        glFramebufferRenderbufferEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_RENDERBUFFER_EXT, rbo);
        
        GLenum status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
        if(status == GL_FRAMEBUFFER_COMPLETE_EXT)
        {
            printf("FBO Setup complete!\n");
        }
        else
        {
            printf("FBO Error!\n");
        }
        
        glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
        
        glDisable(GL_DEPTH_TEST);
        //glEnable(GL_BLEND);     // Turn Blending On
#endif
        
    }
    int msaa = 1;
    
    NSSize contextSize;
    contextSize = [self bounds].size;
    
    /*GLuint screen_width = contextSize.width;
    GLuint screen_height = contextSize.height;
    

    glViewport(0, 0, screen_width*msaa, screen_height*msaa);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60.0f, (float)(screen_width)/(screen_height), 0.01f, 10000.0f);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    

    //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo);
    glDisable(GL_DEPTH_TEST);*/
    
    
    // render here!!!!
    
    //...
    //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    //glHint (GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
    drawAnObject( render_mode );
    stepcount++;
    
    if( is_running )
    {
        
        if( changed_particle_number )
        {
            free(P);
            nres = new_nres;
            P = (particle *)malloc( sizeof(particle) *nres*nres );
            fprintf(stderr, "changed numpart to %d\n",nres);
            aexp = astart;
            init_particles();
            pgsolve->reset_force();
            changed_particle_number = 0;
        }
        
        if( changed_PMres )
        {
            pgsolve->change_PM_resolution( new_PMres );
            PMres = new_PMres;
            changed_PMres = 0;
        }
        
        if( restart_sim )
        {
            aexp = astart;
            init_particles();
            restart_sim = 0;
        }
        
        aexp = pgsolve->step( aexp, 0.001, mass_deploy_mode );
        
        AppDelegate *app = (AppDelegate *)[[NSApplication sharedApplication] delegate];
        NSTextField *myLabel = app->zcurrlabel;
        [myLabel setStringValue:[NSString stringWithFormat:@"%.3f",1.0f/aexp-1.0f ]];
        
        if( aexp <= astart && timedir == 1 )
        {
            NSButton *run_button = app->run_button;
            [run_button setNextState];
            is_running = 0;
        }
    }
    //...
    
    //... show texture
#if 0
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, screen_width , screen_height , 0, -1, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    glDisable(GL_DEPTH_TEST);
    
    glUseProgram(program_postproc);
    
    glDisable(GL_DEPTH_TEST);
    

    glBindTexture(GL_TEXTURE_2D, img);
    glBegin(GL_QUADS);
    glTexCoord2f(0,1);   glVertex2f(0,0);
    glTexCoord2f(0,0);   glVertex2f(0,screen_height);
    glTexCoord2f(1,0);   glVertex2f(screen_width,screen_height);
    glTexCoord2f(1,1);   glVertex2f(screen_width,0);
    glEnd();

    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0); // unbind
    glEnable(GL_TEXTURE_2D);
    glDisable(GL_DEPTH_TEST);
    
    
    glBindTexture(GL_TEXTURE_2D, img);
    glGenerateMipmapEXT(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, 0);
    
    glDisable(GL_TEXTURE_2D);
    
    glViewport(0, 0, screen_width, screen_height);
#endif
    
    //glutSwapBuffers();
    
    [context flushBuffer];
}

-(void) drawRect: (NSRect) bounds
{
    glClearColor(0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT);
    
    [self drawFrame];
    
    glFlush();
    
    
    
    /*glClearColor(0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT);
    drawAnObject();
    glFlush();*/
}

@end
