#include <GL/gl.h>
#include <GL/glu.h>
#include <stdio.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
static SDL_Window* window;
static SDL_Renderer* renderer;

void f_fall (double *x, double *fx, void *params);

void euler_step (int n, double *x0, double *xh, double h,
		 void (*f)(double *x, double *fx, void *params), void *params);

static int const screen_size[2] = {800,600};

float angle = 0;
float omega = 120.0;
double g = 9.8; // m / s2
double radius = 0.3;

double state[2] = {3.0, 0.0};

Uint32 callback(Uint32 interval, void* name)
{
  double dt = interval / 1000.0;
  double newstate[2];
  
  angle += dt * omega;

  midpoint(2, state, newstate, dt, f_fall, &g);

  memcpy(state, newstate, sizeof(newstate));

  if (state[0] - radius < 0)
  {
    state[1] = -state[1];
  }
  
  return interval;
}


static int get_input(void) {
    SDL_Event event;
    while (SDL_PollEvent(&event)) {
        switch (event.type) {
            case SDL_QUIT: return 0; //The little X in the window got pressed
            case SDL_KEYDOWN:
                if (event.key.keysym.sym==SDLK_ESCAPE) {
                    return 0;
                }
                break;
        }
    }
    return 1;
}
static void draw(void) {
    //Clear the screen's color and depth (default color is black, but can change with glClearColor(...))
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
 
    //Drawing to an area starting at the bottom left, screen_size[0] wide, and screen_size[1] high.
    glViewport(0,0,screen_size[0],screen_size[1]);
    //OpenGL is a state machine.  Tell it that future commands changing the matrix are to change OpenGL's projection matrix
    glMatrixMode(GL_PROJECTION);
    //Reset the projection matrix
    glLoadIdentity();
    //Multiply a perspective projection matrix into OpenGL's projection matrix
    gluPerspective(45.0, (double)(screen_size[0])/(double)(screen_size[1]), 0.1,100.0);
 
    //Tell OpenGL that future commands changing the matrix are to change the modelview matrix
    glMatrixMode(GL_MODELVIEW);
    //Reset the modelview matrix
    glLoadIdentity();
    //Multiply OpenGL's modelview matrix with a transform matrix that simulates a camera at (2,3,4) looking towards the location (0,0,0) with up defined to be (0,1,0)
    gluLookAt(4.0,5.0,1.0, 0.0,0.0,0.0, 0.0,0.0,1.0);
 
    //Begin drawing triangles.  Every subsequent triplet of vertices will be interpreted as a single triangle.
    //  OpenGL's default color is white (1.0,1.0,1.0), so that's what color the triangle will be.
    glColor3f(1.0f,1.0f,1.0f);

    glPushMatrix();
    glTranslatef(0.0, 0.f, state[0]);
    glRotated(angle, 0.0, 1.0, 0.0);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    GLUquadric* quadric = gluNewQuadric();
    gluSphere(quadric, radius, 30, 30);
    gluDeleteQuadric(quadric);
    glPopMatrix();

    
    // glBegin(GL_TRIANGLES);
    // //Three vertices follow, these will form a triangle
    // glVertex3f( 0.0f, 0.1f, 0.0f); //Vertex at ( 0.0, 0.1, 0.0)
    // glVertex3f(-0.1f,-0.1f, 0.7f); //Vertex at (-0.1,-0.1, 0.7)
    // glVertex3f( 1.0f,-0.2f, 0.0f); //Vertex at ( 1.0,-0.2, 0.0)
    // //Done drawing triangles
    // glEnd();
 
    //Now we're going to draw some lines to show the cardinal axes.  Every subsequent pair of vertices
    //  will be a single line.
    glBegin(GL_LINES);
    //All subsequent vertices will be red.
    glColor3f(1.0f,0.0f,0.0f);
    glVertex3f(0.0f,0.0f,0.0f); glVertex3f(1.0f,0.0f,0.0f);
    //All subsequent vertices will be green.
    glColor3f(0.0f,1.0f,0.0f);
    glVertex3f(0.0f,0.0f,0.0f); glVertex3f(0.0f,1.0f,0.0f);
    //All subsequent vertices will be blue.
    glColor3f(0.0f,0.0f,1.0f);
    glVertex3f(0.0f,0.0f,0.0f); glVertex3f(0.0f,0.0f,1.0f);
    //Since OpenGL thinks the color is blue now, all subsequent vertices will be blue.  But, we want the
    //  triangle above to be white the *next* time we call this function!  So, reset the color to white.
    glEnd();
 
    //OpenGL works best double-buffered.  SDL automatically sets that up for us.  This will draw what we have
    //  just drawn to the screen so that we can see it.
    SDL_GL_SwapWindow(window);
}
 
int main(int argc, char* argv[]) {
    //Initialize everything, but don't catch fatal signals; give them to the OS.
    SDL_Init(SDL_INIT_EVERYTHING|SDL_INIT_NOPARACHUTE);
    //Creates the window
    window = SDL_CreateWindow("2381 Dynamics", SDL_WINDOWPOS_UNDEFINED,SDL_WINDOWPOS_UNDEFINED,
			      screen_size[0],screen_size[1], SDL_WINDOW_OPENGL);
    //Create an OpenGL context.  In SDL 1, this was done automatically.
    SDL_GLContext context = SDL_GL_CreateContext(window);
 
    //We now have an OpenGL context, and can call OpenGL functions.
 
    //Objects need to test each other to see which one is in front.  If you don't do this, you'll "see through" things!
    glEnable(GL_DEPTH_TEST);

    SDL_TimerID timerID = SDL_AddTimer(20, callback, "SDL");
    //Main application loop
    while (1) {
        if (!get_input()) break;
        draw();
    }

    SDL_RemoveTimer(timerID);
 
    SDL_GL_DeleteContext(context);
    SDL_DestroyWindow(window);
    SDL_Quit();
 
    //Return success; program exits
    return 0;
}

