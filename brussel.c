#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <GLFW/glfw3.h>
#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glext.h>
// compile: gcc -lm brussel.c -o br -lGL -lGLU -lglfw

double x00, y00, vx, vy, a, b, h;
double *xr, *xv, *yr, *yv;
int N;

double X( double x, double y )
{
    return a - ( b + 1 ) * x + pow( x, 2 ) * y;
}

double Y( double x, double y )
{
    return b * x - pow( x, 2 ) * y;
}

void runge( void )
{
    double k[4], l[4], n[4], m[4];
    int i;
    
    xr[0] = x00;
    xv[0] = vx;
    yr[0] = y00;
    yv[0] = vy;
    for ( i = 0; i < N; i++ ) {
        k[0] = xv[i] * h;
        l[0] = X( xr[i], yr[i] ) * h;
        m[0] = yv[i] * h;
        n[0] = Y( xr[i], yr[i] ) * h;
        k[1] = ( xv[i] + 0.5 * l[0] ) * h;
        l[1] = X( xr[i] + 0.5 * k[0], yr[i] + 0.5 * m[0] ) * h;
        m[1] = ( yr[i] + 0.5 * n[0] ) * h;
        n[1] = Y( xr[i] + 0.5 * k[0], yr[i] + 0.5 * m[0] ) * h;
        k[2] = ( xv[i] + 0.5 * l[1] ) * h;
        l[2] = X( xr[i] + 0.5 * k[1], yr[i] + 0.5 * m[1] ) * h;
        m[2] = ( yv[i] + 0.5 * n[1] ) * h;
        n[2] = Y( xr[i] + 0.5 * k[1], yr[i] + 0.5 * m[1] ) * h;
        k[3] = ( xv[i] + l[2] ) * h;
        l[3] = X( xr[i] + k[2], yr[i] + m[2] ) * h;
        m[3] = ( yv[i] + n[2] ) * h;
        n[3] = Y( xr[i] + k[2], yr[i] + m[2] ) * h;
        xr[i+1] = xr[i] + ( 1.0 / 6.0 ) * ( k[0] + 2*k[1] + 2*k[2] + k[3] );
        xv[i+1] = xv[i] + ( 1.0 / 6.0 ) * ( l[0] + 2*l[1] + 2*l[2] + l[3] );
        yr[i+1] = yr[i] + ( 1.0 / 6.0 ) * ( m[0] + 2*m[1] + 2*m[2] + m[3] );
        yv[i+1] = yv[i] + ( 1.0 / 6.0 ) * ( n[0] + 2*n[1] + 2*n[2] + n[3] );
    }
}

GLFWwindow* initOpenGL() {
    // Initialize GLFW
    if (!glfwInit()) {
        fprintf(stderr, "Failed to initialize GLFW\n");
        return NULL;
    }

    // Create a windowed mode window and its OpenGL context
    GLFWwindow* window = glfwCreateWindow(640, 480, "Brusselator Diffusion", NULL, NULL);
    if (!window) {
        fprintf(stderr, "Failed to open GLFW window\n");
        glfwTerminate();
        return NULL;
    }

    // Make the window's context current
    glfwMakeContextCurrent(window);

    // Initialize GLEW
    glewExperimental = GL_TRUE; // Needed for core profile
    if (glewInit() != GLEW_OK) {
        fprintf(stderr, "Failed to initialize GLEW\n");
        return NULL;
    }

    return window;
}

int main( void )
{
    FILE *f;
    int i;

    printf( ">> Enter N: " );
    scanf( "%d", &N );
    printf( ">> Enter h: " );
    scanf( "%lf", &h );
    printf( ">> Enter a, b: " );
    scanf( "%lf %lf", &a, &b );
    printf( ">> Enter x0, y0: " );
    scanf( "%lf %lf", &x00, &y00 );
    printf( ">> Enter xv, yv: " );
    scanf( "%lf %lf", &vx, &vy );
    xr = (double *) malloc( N * sizeof(double) );
    xv = (double *) malloc( N * sizeof(double) );
    yr = (double *) malloc( N * sizeof(double) );
    yv = (double *) malloc( N * sizeof(double) );
    puts( ">> Calculate..." );
    runge();
    puts( ">> Done" );
    puts( "[DATA]:   step xr xv yr yv" );
    f = fopen( "output.txt", "w" );
    if ( f == NULL ) {
        puts( "[ERROR]:  Can't create output.txt!" );
        puts( "[OUTPUT]: stdout" );
        f = stdout;
    } else {
        puts( "[OUTPUT]: output.txt" );
    }
    for ( i = 0; i < N; i++ ) {
        fprintf( f, "%.16lf %.16lf %.16lf %.16lf %.16lf\n", 
            (double) i * h, xr[i], xv[i], yr[i], yv[i] );
    }
    if ( f != stdout ) {
        fclose( f );
    }
    free( xr );
    free( xv );
    free( yr );
    free( yv );

    //OpenGl
    GLFWwindow* window = initOpenGL();
    if (!window) {
        return -1;
    }
    while (!glfwWindowShouldClose(window)) {
        runge();
        glClear(GL_COLOR_BUFFER_BIT);

        // Render the Brusselator's state here
        // This involves mapping the values of xr, xv, yr, yv to colors and drawing them on the screen
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    glfwTerminate();


    // Gnuplot
    FILE *gnuplotScript = fopen("gnuplotScript.gp", "w");
    if (gnuplotScript == NULL) {
        perror("Error opening file");
        return 1;
    }

    // Запись команд в скрипт
    fprintf(gnuplotScript, "set terminal x11\n");
    fprintf(gnuplotScript, "plot 'output.txt' u 1:2 t 'X(t)' w l, 'output.txt' u 1:4 t 'Y(t)' w l\n");
    fprintf(gnuplotScript, "set terminal x11 1\n");
    fprintf(gnuplotScript, "plot 'output.txt' u 2:3 t 'X(XV)' w l\n");
    fprintf(gnuplotScript, "set terminal x11 2\n");
    fprintf(gnuplotScript, "plot 'output.txt' u 4:5 t 'Y(YV)' w l\n");
    fprintf(gnuplotScript, "set terminal x11 3\n");
    fprintf(gnuplotScript, "plot 'output.txt' u 2:4 t 'X(Y)' w l\n");
    fprintf(gnuplotScript, "pause -1 'Press Enter to close'\n");
    fclose(gnuplotScript);

    // Запуск Gnuplot с временным скриптом
    system("gnuplot gnuplotScript.gp");

    // Удаление временного файла
    remove("gnuplotScript.gp");

    return 0;
}