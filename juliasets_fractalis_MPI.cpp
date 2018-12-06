//
//  juliasets_fractalis.cpp
//  PaPP_Projekt_2
//
//  Created by Richard Krizan on 27/11/2018.
//  Copyright © 2018 Richard Krizan. All rights reserved.
//
#define GL_SILENCE_DEPRECATION

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>
#include <sys/time.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/freeglut.h>
#include <GL/freeglut_ext.h>
#endif



#define width 1280
#define height 720

#define END_TAG    0
#define DATA_TAG    1
#define RESULT_TAG  2
#define VAR_INIT_TAG 3



typedef struct {
    GLubyte r,g,b;
} Pixel;

typedef struct {
    int iterations, color_profile, omp_enabled, thread_count;
} s_variables;


GLuint texture;
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
double newZx,newZy,oldZx,oldZy;
int iterations  = 1000;
int thread_count = 10;
int wid;
Pixel image[height][width];
Pixel mapping[16];
int nproc, rank;
int omp_enabled = 1;



void display();
void calc_mandelbrot(unsigned row, unsigned *row_data);
Pixel calc_color(unsigned n);


void master() {
    
    MPI_Status mpi_status;
    
    timeval t_start, t_end;
    double t_start_omp = 0.0, t_end_omp;
    
    int row = 0;
    int row_data[width+1];
    
    if (omp_enabled) t_start_omp = omp_get_wtime();
    else gettimeofday(&t_start, NULL);
    
    for (int i = 1; i < nproc; ++i) {
        MPI_Send(&row, 1, MPI_INT, i, DATA_TAG, MPI_COMM_WORLD);
        row++;
    }
    
    int complete_rows = 0;
  
    while (complete_rows < height) {
        MPI_Recv(&row_data, width+1, MPI_INT, MPI_ANY_SOURCE, RESULT_TAG, MPI_COMM_WORLD, &mpi_status);
        
        int slave_done = mpi_status.MPI_SOURCE;
        int received_row = row_data[0];
        
        for (int column = 0; column < width; ++column) {
            image[received_row][column] = calc_color(row_data[column+1]);
        }
        
        complete_rows++;
        if (row < height) {
            MPI_Send(&row, 1, MPI_INT, slave_done, DATA_TAG, MPI_COMM_WORLD);
            row++;
        }
    }
    
    for (int i = 1; i < nproc; ++i) {
        MPI_Send(0, 0, MPI_INT, i, END_TAG, MPI_COMM_WORLD);
    }
    
    if (omp_enabled) {
        t_end_omp = omp_get_wtime();
        printf("(DONE) Total computation time %f seconds\n", t_end_omp - t_start_omp);
        //printf("%f\n",t_end_omp-t_start_omp);
        //exit(EXIT_SUCCESS);
    }
    else {
        gettimeofday(&t_end, NULL);
        float total_time = ((t_end.tv_sec - t_start.tv_sec) * 1000000u + t_end.tv_usec - t_start.tv_usec) / 1.e6;
        printf("(DONE) Total computation time %.06lf seconds\n", total_time);
        //printf("%.06lf\n",total_time);
        //exit(EXIT_SUCCESS);
    }
    
    glutDisplayFunc(display);
   
    glutMainLoop();
    
}

void slave() {
    
    MPI_Status mpi_status;
    
    unsigned row = 0;
    unsigned row_data[width+1];
    
    int slave_rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &slave_rank);
    
    MPI_Recv(&row, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);
    
    while (mpi_status.MPI_TAG == DATA_TAG) {
        if (mpi_status.MPI_TAG == END_TAG) exit(EXIT_SUCCESS);
        calc_mandelbrot(row, row_data);
        MPI_Send(row_data, width+1, MPI_INT, 0, RESULT_TAG, MPI_COMM_WORLD);
        MPI_Recv(&row, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);
    }
    
}



void init_mpi(int argc, char ** argv) {

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    if (nproc < 2) {
        printf("ERROR: Aspon dva porcesory su potrebne, mas %d\n",nproc);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    printf("Number of processors: %d\n", nproc);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
}

Pixel calc_color(unsigned n) {
    Pixel p;
    
        if (n < iterations && n > 0) {
            
            p = mapping[n%16];
        }
        else if(n==0 || n==iterations) {
            p.r = 0;
            p.g = 0;
            p.b = 0;
        }
    return p;
}

void calc_mandelbrot(unsigned y, unsigned *row_data) {
    unsigned x;
#pragma omp parallel for if(omp_enabled) private(x,newZx,newZy,oldZx,oldZy) num_threads(thread_count)
    for(x=0; x<width; ++x) {
        newZx = xmin * zoomX + ((x+ moveX) * nx * zoomX ) / width  ;
        newZy = ymin * zoomY + ((y+ moveY) * ny * zoomY ) / height  ;
        int n;
        
        for(n=0; n<iterations; n++) {
            oldZx = newZx;
            oldZy = newZy;
           
            newZx = oldZx * oldZx - oldZy * oldZy + cRe;
            newZy = 2 * oldZx * oldZy + cIm;
            if((newZx * newZx + newZy * newZy) > 4 ) break;
        }
        row_data[x+1] = n;
    }
 row_data[0] = y;
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
    mapping[13].r = 153; mapping[13].g = 87; mapping[13].b = 0;
    mapping[13].r = 106; mapping[13].g = 52; mapping[13].b = 3;
}


void display() {
    
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, image);
    glClear(GL_COLOR_BUFFER_BIT);
    glBegin(GL_QUADS);
    glTexCoord2f(1,0); glVertex2f(1,-1);
    glTexCoord2f(1,1); glVertex2f(1,1);
    glTexCoord2f(0,1); glVertex2f(-1,1);
    glTexCoord2f(0,0); glVertex2f(-1,-1);
    glEnd();
    glFlush();
    glutPostRedisplay();
    glutSwapBuffers();
    
}

void init_glut(int argc, char ** argv) {
    
    glutInit(&argc, argv);
    glutInitWindowSize(width,height);
    glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB|GLUT_DEPTH);
    wid = glutCreateWindow("Julia Sets Fraktalis");
    
    glEnable(GL_TEXTURE_2D);
    glGenTextures( 1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    
    glClearColor(0,0,0,0);
    gluOrtho2D(-1,1,-1,1);
    glLoadIdentity();
    glColor3f(1,1,1);
    init_color_mapping();
    
}



int main(int argc, char ** argv) {
    init_mpi(argc, argv);
    
  
    const int nitems = 4;
    int blocklengths[nitems] = {1,1,1,1};
    MPI_Datatype types[nitems] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT};
    MPI_Datatype mpi_s_variables;
    MPI_Aint offsets[nitems];
    
    offsets[0] = offsetof(s_variables, iterations);
    offsets[1] = offsetof(s_variables, color_profile);
    offsets[2] = offsetof(s_variables, omp_enabled);
    offsets[3] = offsetof(s_variables, thread_count);
    
    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_s_variables);
    MPI_Type_commit(&mpi_s_variables);
    
    if (rank==0) {
        s_variables send;
        send.iterations = iterations;
        send.omp_enabled = omp_enabled;
        send.thread_count = thread_count;
        
        for (int i = 1; i < nproc; ++i) {
            MPI_Send(&send, 1, mpi_s_variables, i, VAR_INIT_TAG, MPI_COMM_WORLD);
        }
        
        init_glut(argc, argv);
        master();
    }
    else {
        
        s_variables recv;
        
        MPI_Status status;
        MPI_Recv(&recv, 1, mpi_s_variables, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        iterations = recv.iterations;
        omp_enabled = recv.omp_enabled;
        thread_count = recv.thread_count;
        
        slave();
    }
    
    MPI_Finalize();
    return EXIT_SUCCESS;
}
