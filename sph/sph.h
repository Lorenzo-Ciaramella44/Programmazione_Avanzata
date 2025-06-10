#ifndef SPH_H
#define SPH_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "time.h"
#define DEBUG 1			/* printf("DEBUG %d\n", DEBUG); */

#define N 3600							/* total particle number */
#define N_L 3200						/* sx particle number */
#define N_R 400							/* dx particle number */
#define PRESSURE_CONVERSION 100000.0  	/* conversion into Pascal */
#define ENERGY_CONVERSION 100000.0    	/* conversion into Joule per kg */
#define GAMMA (5.0/3.0) 				/* coefficient of monoatomic gas */
#define EPSILON 0.01					/* constant to avoid problems */
#define ALPHA 1.0						/* artificial viscosity coefficient */
#define BETA 2.0 						/* artificial viscosity coefficient */
#define K 0.1							/* time step coefficient */
#define MAX_NEIGHBORS 64				/* max number of neighbors*/


/* Struct of particle */
typedef struct {
    double x;        // Position 
    double rho;      // Density 
    double velocity; // Velocity 
    double pressure; // Pression 
    double accel;    // Acceleration 
    double u;        // Internal energy 
    double du_dt;    // Time derivative of internal energy
    double h;        // Smoothing lenght 
} Particle;


/* Struct of neighbor list */
typedef struct {
    int count;                     	// Neighbors number 
    int neighbors[MAX_NEIGHBORS];   // Indices
} NeighborList;


/* Struct to order particles due to distance */
typedef struct {
    int index;         // Particle index
    double distance;   // Distance
} NeighborDistance;


/* Struct for the grid */
typedef struct {
    int particle_indices[N]; 	// Indices of particle in the selected cell
    int count;
} Cell;


double sgn(double x);
void condition_0(double *dim, double *vel1, double *vel2, double *press1, double *press2, double *dens1, double *dens2);
void write_file_particles(int step_counter);
double periodic_boundary_condition(double dx, double dim);
int compare_distances(const void *a, const void *b);
void qfind_neighbors_with_h(double *dim);
void write_file_neighbors(int step_counter);
void pression();
void density(double *dim, double *m);
double sound_speed(int i);
double artificial_viscosity(int a, int b, double dim);
void du_dt(double *dim, double *m);
void acceleration(double *dim, double *m);
double der_W(double x, double h);
double W_func(double x, double h);
double compute_delta_t();
double kick_drift(double *dim);
void kick(double *t_0, double *t);
void find_neighbors_with_h(Cell *grid, double *dim);
Cell* assign_particles_to_cells(double *dim);



extern Particle particles[N];
extern NeighborList neigh[N];

#endif