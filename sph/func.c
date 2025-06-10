#include "sph.h"


/* Define the sign of x */
double sgn(double x) {
    return (x > 0) - (x < 0);  
}


/* Set condition 0 */
void condition_0(double *dim, double *vel1, double *vel2, double *press1, double *press2, double *dens1, double *dens2) {
	double L_velocity, R_velocity;
	double L_pressure, R_pressure;
	double L_density, R_density;
	double tube_lenght;
	
	L_velocity = *vel1;
	R_velocity = *vel2;
	L_pressure = *press1;
	R_pressure = *press2;
	L_density = *dens1;
	R_density = *dens2;
	tube_lenght = *dim;
	
    // Inizialize the list of sx particles
    for (int i = 0; i < N_L; i++) {
        particles[i] = (Particle){
            .x = i * (tube_lenght / (2 * N_L)),
            .rho = L_density,
            .velocity = L_velocity,
            .pressure = L_pressure,
            .u = L_pressure / (L_density * (GAMMA - 1)),
            .h = 32*tube_lenght/N
        };
    }

    // Inizialize the list of dx particles
    for (int i = N_L; i < N; i++) {
        particles[i] = (Particle){
            .x = (i - N_L) * (tube_lenght / (2 * N_R)) + (tube_lenght / 2),
            .rho = R_density,
            .velocity = R_velocity,
            .pressure = R_pressure,
            .u = R_pressure / (R_density * (GAMMA - 1)),
            .h = 32*tube_lenght/N
        };
    }
}


/* Writing particle's data on file */
void write_file_particles(int step_counter) {
    char filename[50];
    sprintf(filename, "particles_step_%d.txt", step_counter);
    FILE *file = fopen(filename, "w");

    if (file == NULL) {
        fprintf(stderr, "Error: impossible opening file %s in writing.\n", filename);
        return;
    }

    // Header
    fprintf(file, "# x [m]   rho [kg/m^3]   velocity [m/s]   pressure [Pa]   accel [m/s^2]   energy [J/kg]   du_dt [J/(kg·s)]   h [m]\n");

    for (int i = 0; i < N; i++) {
        double converted_pressure = particles[i].pressure * PRESSURE_CONVERSION; // Converte la pressione in Pascal
        double converted_energy = particles[i].u * ENERGY_CONVERSION; // Converte l'energia in Joule/kg
        fprintf(file, "%f %f %f %f %f %f %f %f\n", 
                particles[i].x, 
                particles[i].rho, 
                particles[i].velocity, 
                converted_pressure, 
                particles[i].accel, 
                converted_energy, 
                particles[i].du_dt,
                particles[i].h);
    }

    fclose(file);
}


/* periodic boundary condition */
double periodic_boundary_condition(double dx, double dim) {
    if (dx > 0.5 * dim) {  
        dx -= dim;  
    } else if (dx < -0.5 * dim) {  
        dx += dim; 
    }
    return dx;
}



/* Density function */
void density(double *dim, double *m) {
	double tube_lenght, MASS;
	
	tube_lenght = *dim;
	MASS = *m;
	
    for (int i = 0; i < N; i++) {
    	particles[i].rho = 0.0;
        for (int j = 0; j < neigh[i].count; j++) {
        	int b = neigh[i].neighbors[j];
            double dx = particles[i].x - particles[b].x;
            dx = periodic_boundary_condition(dx, tube_lenght);
            
            // Define the smoothing length
            double H = particles[i].h;  
            particles[i].rho += MASS * W_func(dx, H);
            
        } 
    }
}


/* Kernel Cubic Spline for SPH */
double W_func(double x, double h) {
    double r = fabs(x);		// absolute distance between particles is needed
    double q = r / h;
    double norm = 2.0 / (3.0 * h);
    double W;

    if (q <= 1.0) {
        W = norm * (1.0 - 1.5 * q * q + 0.75 * q * q * q);
    } 
    else if (q <= 2.0) {
        double term = 2.0 - q;
        W = norm * 0.25 * term * term * term;
    } 
    else {
    	W = 0.0;
    }
    return W;
}


/* Pressure function */
void pression() {
    for (int i = 0; i < N; i++) {
        particles[i].pressure = (GAMMA - 1) * particles[i].rho * particles[i].u;  
    }
}


/* Assigning particles to cells function */
Cell* assign_particles_to_cells(double *dim) {
	double tube_lenght, cell_size;
	double H = particles[0].h; 
	int grid_size;
	
	tube_lenght = *dim;
	
	/* Maximize h-value finder */
    for (int i = 1; i < N; i++) {
        if (particles[i].h > H) {
            H = particles[i].h;
        }
    }
    
    cell_size = H;
	grid_size = (int)(tube_lenght / cell_size);
	
	Cell *grid = (Cell*)malloc(grid_size * sizeof(Cell));
    if (grid == NULL) {
        fprintf(stderr, "Errore di allocazione della memoria per la griglia!\n");
        exit(1);
    }
	
	for (int i = 0; i < grid_size; i++){
        grid[i].count = 0; 		// Reset
    }

    for (int i = 0; i < N; i++) {
        int cell_index = (int)(particles[i].x / cell_size);
        
        /* re-set over limit accesses */
        if (cell_index >= grid_size) {
            cell_index = grid_size - 1;
        }
        
        int index = grid[cell_index].count;  			// Save count value
        grid[cell_index].particle_indices[index] = i;  	// Assign the correct value
        grid[cell_index].count++;  						// Increase the counter after assign operation
    }
    
    return grid;
}


/* Neighbors function */
void find_neighbors_with_h(Cell *grid, double *dim) {
    NeighborDistance distances[N];
    double tube_lenght, cell_size;
	int grid_size;
	double H = particles[0].h; 
	
    tube_lenght = *dim;
    
    /* Maximize h-value finder */
    for (int i = 1; i < N; i++) {
        if (particles[i].h > H) {
            H = particles[i].h;
        }
    }
    
	cell_size = H;
	grid_size = (int)(tube_lenght/cell_size);
    
	/* Checking validity of grid_size */
    if (grid_size <= 0) {
        fprintf(stderr, "Errore: grid_size non valido (%d)\n", grid_size);
        return;
    }
	
    /*Check Cell count */
	/*
	for (int i = 0; i < grid_size; i++) {
		printf("Cell %d has %d particles\n", i, grid[i].count);
	}
	*/
	
    for (int i = 0; i < N; i++) {
        double x_i = particles[i].x;
        int num_neighbors = 0;

        /* Defining particle's cell */
        int cell_index = (int)(x_i / cell_size);
		
        /* re-set over limit accesses */
		if (cell_index >= grid_size) {
			cell_index = grid_size - 1; 
		}
        /* checking cell_index */
        if (cell_index >= grid_size || cell_index < 0) {
            fprintf(stderr, "Error: Indix cell %d out of boundaries!\n", cell_index);
            continue;
        }

        /* Checking on near cells (with periodicity) */
        for (int offset = -1; offset <= 1; offset++) {
            int neighbor_cell = (cell_index + offset + grid_size) % grid_size;
            
            for (int j = 0; j < grid[neighbor_cell].count; j++) {
                int neighbor_index = grid[neighbor_cell].particle_indices[j];
                if (i == neighbor_index) continue; 			// skip himself

                double dx = particles[neighbor_index].x - x_i;
                dx = periodic_boundary_condition(dx, tube_lenght);
                distances[num_neighbors].index = neighbor_index;
                distances[num_neighbors].distance = fabs(dx);
                num_neighbors++;
            }
        }

        /* Sorting neoghbors by distance */
        qsort(distances, num_neighbors, sizeof(NeighborDistance), compare_distances);

        /* MAX_NEIGHBORS limit */
        for (int j = 0; j < MAX_NEIGHBORS; j++) {
            neigh[i].neighbors[j] = distances[j].index;
        }
        neigh[i].count = MAX_NEIGHBORS;

        /* Getting the new h */
        particles[i].h = 0.5 * distances[MAX_NEIGHBORS - 1].distance;
    }
    
}


/* Sound speed function */
double sound_speed(int i) {
    /* Pressure check */
    if (particles[i].pressure <= 0) {
        //printf("Attention: no valid pressure for particle %d (pressure = %f)\n", i, particles[i].pressure);
        return 0;
    }
    
    /* Density check */
    if (particles[i].rho <= 0) {
        //printf("Attention: no valid density for particle %d (density = %f)\n", i, particles[i].rho);
        return 0;
    }

    double cs = sqrt(GAMMA * particles[i].pressure / particles[i].rho);
    return cs;
}


/* Artificial viscosity */
double artificial_viscosity(int a, int b, double dim) {
    double v_ab = particles[a].velocity - particles[b].velocity;
    
    double dx = particles[a].x - particles[b].x;
    dx = periodic_boundary_condition(dx, dim);

    double cs_ab = (sound_speed(a) + sound_speed(b)) / 2.0;
    
    double H = particles[a].h;  

    if (v_ab * dx < 0) {
        double mu_ab = H * v_ab * dx / (dx * dx + EPSILON * H * H);  
        double rho_ab = (particles[a].rho + particles[b].rho) / 2.0;  // Media della densità

        return (-ALPHA * mu_ab * cs_ab + BETA * mu_ab * mu_ab) / rho_ab;
    }

    return 0.0;  
}


/* Time derivative of internal energy */
void du_dt(double *dim, double *m) {
	double tube_lenght, MASS;
	
	tube_lenght = *dim;
	MASS = *m;
    for (int i = 0; i < N; i++) {
        double sum = 0.0;
        for (int j = 0; j < neigh[i].count; j++) {
            int b = neigh[i].neighbors[j];
           
            double dx = particles[i].x - particles[b].x;
            dx = periodic_boundary_condition(dx, tube_lenght);  

            double H = particles[i].h;
            double gradW = der_W(dx, H);
            double pi_ab = artificial_viscosity(i, b, tube_lenght);
            double v_ab = particles[i].velocity - particles[b].velocity;

            sum += MASS * (particles[i].pressure / (particles[i].rho * particles[i].rho) + pi_ab / 2.0) * v_ab * gradW;
        }
        particles[i].du_dt = sum;
    }
}

/* Acceleration function */
void acceleration(double *dim, double *m) {
	double tube_lenght, MASS;
	
	tube_lenght = *dim;
	MASS = *m;
    for (int i = 0; i < N; i++) {
        double sum = 0.0;  
        for (int j = 0; j < neigh[i].count; j++) {
            int b = neigh[i].neighbors[j];
            double dx = particles[i].x - particles[b].x;
            dx = periodic_boundary_condition(dx, tube_lenght);  

            double H = particles[i].h;
            double gradW = der_W(dx, H);
            double pi_ab = artificial_viscosity(i, b, tube_lenght);
            double term = (particles[i].pressure / (particles[i].rho * particles[i].rho) + 
                           particles[b].pressure / (particles[b].rho * particles[b].rho) + pi_ab);
            
            sum -= MASS * term * gradW;
        }
        particles[i].accel = sum;
    }
}


/* Writing function for neighbor list */
void write_file_neighbors(int step_counter) {
    char filename[50];
    sprintf(filename, "neighbors_step_%d.txt", step_counter);
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        fprintf(stderr, "Error opening file\n");
        return;
    }

    for (int i = 0; i < N; i++) {
        fprintf(file, "Particella id %d: ", i);
        for (int j = 0; j < neigh[i].count; j++) {
            fprintf(file, "%d ", neigh[i].neighbors[j]);
        }
        fprintf(file, "\n");
    }

    fclose(file);
}


/* Derivative of Kernel Cubic Spline */
double der_W(double x, double h) {
    double r = fabs(x);  	// absolute distance between particles is needed
    double q = r / h;    
    double norm = 2.0 / (3.0 * h * h);
    double grad;

    if (q <= 1.0) {
        grad = norm * (-3.0 * q + 2.25 * q * q);
    } else if (q <= 2.0) {
        double term = 2.0 - q;
        grad = -norm * 0.75 * term * term;
    } else {
        grad = 0.0;
    }

    grad *= sgn(x);  		// to get the correct sign in order of distance
    return grad;
}


/* Time step function */
double compute_delta_t() {
    /* Maximize h-value finder */
    double H = particles[0].h;  
    for (int i = 1; i < N; i++) {
        if (particles[i].h > H) {
            H = particles[i].h;
        }
    }

    /* Inizialization of min_delta_t with 1st value of particle */
    double min_delta_t = K * H / sound_speed(0);

    for (int i = 1; i < N; i++) {
        double cs = sound_speed(i);
        double delta_t = K * H / cs;        
        
        if (delta_t < min_delta_t) {
            min_delta_t = delta_t;
        }
    }
    
    return min_delta_t;
}


/* Kick-drift function */
double kick_drift(double *dim) {
    double tube_lenght;
	/* Getting time step */
    double delta_t = compute_delta_t();

    tube_lenght = *dim;
    /* Kick */
    for (int i = 0; i < N; i++) {
        particles[i].velocity += particles[i].accel * delta_t / 2.0;
        particles[i].u += particles[i].du_dt * delta_t / 2.0;
    }

    /* Drift */
    for (int i = 0; i < N; i++) {
        particles[i].x += particles[i].velocity * delta_t;  

        /* Periodic boundary condition */
        if (particles[i].x < 0.0) {
            particles[i].x += tube_lenght;
        } else if (particles[i].x >= tube_lenght) {
            particles[i].x -= tube_lenght;
        }
    }
    return delta_t;
}


/* Kick function */
void kick(double *t_0, double *t) {
	double time_0, delta_t;
	
	time_0 = *t_0;
	delta_t = *t;
	/* Kick */
	for (int i = 0; i < N; i++) {
        particles[i].velocity += particles[i].accel * delta_t / 2.0;
        particles[i].u += particles[i].du_dt * delta_t / 2.0;
    }
}





/* Compare function for qsort */
int compare_distances(const void *a, const void *b) {
    double d1 = ((NeighborDistance *)a)->distance;  // Distance of 1st particle
    double d2 = ((NeighborDistance *)b)->distance;  // Distance of 2nd particle
    return (d1 > d2) - (d1 < d2);  					// Getting -1, 0 or 1
}


/* Neighbor function with h parameter and quick sort */
void qfind_neighbors_with_h(double *dim) {
    NeighborDistance distances[N]; 	// Array to memorize distances and indices
    double tube_lenght;
    
    tube_lenght = *dim;
    
    for (int i = 0; i < N; i++) {
        double x_i = particles[i].x;  
        int num_neighbors = 0;

        for (int j = 0; j < N; j++) {
            if (i == j) continue;  // Skip himself

            double dx = particles[j].x - x_i;  
            dx = periodic_boundary_condition(dx, tube_lenght);

            distances[num_neighbors].index = j;
            distances[num_neighbors].distance = fabs(dx);
            num_neighbors++;
            //printf("asdf %d %lf %d %d \n", distances[num_neighbors].index, distances[num_neighbors].distance, i, j); 
			
        }

        /* Quick sort by distance*/
        qsort(distances, num_neighbors, sizeof(NeighborDistance), compare_distances);

        /* Debug print for uncorrect number of neighbors */
        if (num_neighbors < MAX_NEIGHBORS) {
            printf("Error: Particle %d has just %d neighbors!\n", i, num_neighbors);
            exit(1);
        }

        for (int j = 0; j < MAX_NEIGHBORS; j++) {
            neigh[i].neighbors[j] = distances[j].index;
        }
        neigh[i].count = MAX_NEIGHBORS;

        /* Max distance between 2 particles is 2h */
        particles[i].h = 0.5 * distances[MAX_NEIGHBORS - 1].distance;
    }
}