//
//  juliasets_fractalis.cpp
//  PaPP_Projekt_2
//
//  Created by Richard Krizan on 27/11/2018.
//  Copyright © 2018 Richard Krizan. All rights reserved.
//
#define GL_SILENCE_DEPRECATION

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <omp.h>



#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/freeglut.h>
#include <GL/freeglut_ext.h>
#endif

using namespace std;

typedef struct {
    GLubyte r;
    GLubyte g;
    GLubyte b;
} Pixel;

typedef struct {
    GLfloat real;
    GLfloat imag;
} Complex;

#define width 1280
#define height 720
int wid;
int nproc, rank;

double cRe = -0.79;
double cIm = 0.15;
double xmin   = -2.0;
double ymin   = -2.0;
double nx = 4.0;
double ny = 4.0;
double moveY = 0;
double moveX = 0;
double zoomX = 1;
double zoomY = 1;
//img_min+(real_max-real_min)*height/width; //bottom border
double newZx,newZy,oldZx,oldZy;
Pixel image[height][width];
Pixel mapping[16];
int thread_count = 1,iterations = 1000;
GLuint texture;

void RenderFrame() {
  

    double start_time = omp_get_wtime();
    unsigned x,y;
#pragma omp parallel for num_threads(thread_count) private (x,y, newZx,newZy,oldZx,oldZy)
    for(y=0; y<height; ++y) {
        for(x=0; x<width; ++x) {
          
            newZx = xmin * zoomX + ((x+ moveX) * nx * zoomX ) / width  ;
            newZy = ymin * zoomY + ((y+ moveY) * ny * zoomY ) / height  ;
            int n;
            
            for(n=0; n<iterations; ++n) {
                oldZx = newZx;
                oldZy = newZy;
               
                newZx = oldZx * oldZx - oldZy * oldZy + cRe;
                newZy = 2 * oldZx * oldZy + cIm;
                if((newZx * newZx + newZy * newZy) > 4 )break;
            }
            
            Pixel p;
            if(n < iterations){
                p = mapping[n % 16];
            }
            else{
                p.r = 0;
                p.g = 0;
                p.b = 0;
            }
            image[y][x] = p;
            
        }
    }
    
    double end_time = omp_get_wtime();
    printf("%f\n", end_time-start_time);
}

void init_color_mapping() {
    mapping[0].r = 66; mapping[0].g = 30; mapping[0].b = 15;
    mapping[1].r = 25; mapping[1].g = 7; mapping[1].b = 26;
    mapping[2].r = 9; mapping[2].g = 1; mapping[2].b = 47;
    mapping[3].r = 4; mapping[3].g = 4; mapping[3].b = 73;
    mapping[4].r = 0; mapping[4].g = 7; mapping[4].b = 100;
    mapping[5].r = 12; mapping[5].g = 44; mapping[5].b = 138;
    mapping[6].r = 24; mapping[6].g = 82; mapping[6].b = 177;
    mapping[7].r = 57; mapping[7].g = 125; mapping[7].b = 209;
    mapping[8].r = 134; mapping[8].g = 181; mapping[8].b = 229;
    mapping[9].r = 211; mapping[9].g = 236; mapping[9].b = 248;
    mapping[10].r = 241; mapping[10].g = 233; mapping[10].b = 191;
    mapping[11].r = 248; mapping[11].g = 201; mapping[11].b = 95;
    mapping[12].r = 255; mapping[12].g = 170; mapping[12].b = 0;
    mapping[13].r = 204; mapping[13].g = 128; mapping[13].b = 0;
    mapping[14].r = 153; mapping[14].g = 87; mapping[14].b = 0;
    mapping[15].r = 106; mapping[15].g = 52; mapping[15].b = 3;
}


void display() {
    
    // Call user image generation
    RenderFrame();
    
    // Copy image to texture memory
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, image);
    
    // Clear screen buffer
    glClear(GL_COLOR_BUFFER_BIT);
    
    // Render a quad
    glBegin(GL_QUADS);
    glTexCoord2f(1,0); glVertex2f(1,-1);
    glTexCoord2f(1,1); glVertex2f(1,1);
    glTexCoord2f(0,1); glVertex2f(-1,1);
    glTexCoord2f(0,0); glVertex2f(-1,-1);
    glEnd();
    
    // Display result
    glFlush();
    glutSwapBuffers();
    
}



void init_glut (int argc, char ** argv) {
    
  
    // Texture setup
    glEnable(GL_TEXTURE_2D);
    glGenTextures( 1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    
    // Other
    glClearColor(0,0,0,0);
    gluOrtho2D(-1,1,-1,1);
    glLoadIdentity();
    glColor3f(1,1,1);
    
    init_color_mapping();
}


void keypress(unsigned char key, int x, int y) {
    printf("%c\n",key);
    switch (key) {
        case 'w':
            moveY+=5/zoomY*0.6;
            
            break;
        case 's':
            moveY-=5/zoomY*0.6;
            
            break;
        case 'd':
            moveX+=5/zoomX*0.6;
            
            break;
        case 'a':
            moveX-=5/zoomX*0.6;
            
            break;
        case 'q':
            zoomX*=0.5;
            zoomY*=0.5;
            
            break;
        case 'e':
            zoomX*=1.5;
            zoomY*=1.5;
            break;
        case 27:
            glutDestroyWindow(wid);
            exit(EXIT_SUCCESS);
            break;
        default:
            break;
    }
    
    glutPostRedisplay();
}




int main(int argc, char ** argv) {
    
    glutInit(&argc, argv);
    glutInitWindowSize(width,height);
    glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB|GLUT_DEPTH);
    wid=glutCreateWindow("JuliaSets Fraktalis");
    glutKeyboardFunc(keypress);
    glutDisplayFunc(display);
    
    init_glut(argc,argv);
    //Main exec. loop
    
    glutMainLoop();
    
    return EXIT_SUCCESS;
}

