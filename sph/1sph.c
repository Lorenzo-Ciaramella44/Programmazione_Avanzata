#include "sph.h"


double main() {
	clock_t start, end;								/* variables for seeing computing time required*/
    double cpu_time_used;
    start = clock();
    FILE *fparam;
    
    int index = 0;
    double present_time = 0.0;						/* starting time in s */
    double Max_time = 0.2;							/* limit ending time in s */
    double delta_time;								/* increment time in s */
    double density_L = 1.0, density_R = 0.125;		/* kg/m^3 */
    double dim_tube = 1.0;							/* m */
    double pressure_L = 1.0, pressure_R = 0.1;		/* Pa */
    double velocity_L = 0.0, velocity_R = 0.0;		/* m/s */
    double area = 1.0;								/* m^2 */
    double mass_p = 0.00015625;						/* kg */
    
    condition_0(&dim_tube, &velocity_L, &velocity_R, &pressure_L, &pressure_R, &density_L, &density_R);	// set the initial condition
    write_file_particles(index);		// write the beginning particle file 
    write_file_neighbors(index);		// write the beginning neighbor file
    
    /* starting the evolution process */
    while (present_time < Max_time) {
		index++;
		Cell* grid = assign_particles_to_cells(&dim_tube);  // create the grid
		find_neighbors_with_h(grid, &dim_tube);				// near neighbors method
		//qfind_neighbors_with_h(&dim_tube);				// quick sort method
		free(grid);											// free the grid
		density(&dim_tube, &mass_p);						// compute density
		pression();											// compute pression
		acceleration(&dim_tube, &mass_p);					// compute acceleration
		du_dt(&dim_tube, &mass_p);							// compute energy derivative 
		delta_time = kick_drift(&dim_tube);					// kick drift and computation of time increment
		present_time += delta_time;							// update time
		
		grid = assign_particles_to_cells(&dim_tube);  		// create the grid
		find_neighbors_with_h(grid, &dim_tube);				// near neighbors method
		//qfind_neighbors_with_h(&dim_tube);				// quick sort method
		free(grid);											// free the grid
		density(&dim_tube, &mass_p);						// compute density
		pression();											// compute pression
		acceleration(&dim_tube, &mass_p);					// compute acceleration
		du_dt(&dim_tube, &mass_p);							// compute energy derivative 
		kick(&present_time, &delta_time);					// kick
		write_file_particles(index);						// rewrite the index+1 particle file
		write_file_neighbors(index);						// rewrite the index+1 particle file
		
		printf("t %lf in round %d, progression %lf %%\n", present_time, index, present_time/Max_time*100);	// updated time	
	}
	/* ending of evolution process with printed time*/
    printf("Final time at the end of the simulation %lf \n", present_time);
    
	/* creating of file for muli_plot.gp */
    fparam = fopen("param.txt", "w");
	fprintf(fparam, "a = %d\n", index);
	fclose(fparam);
	
	
    end = clock();								/* time ends */
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Elaboration time: %f seconds\n", cpu_time_used);
}


/* to run the code unify nucleo.c rename.c func_lib.c
gcc 1sph.c -lm
./a.out

!!! before running changing the index in multi_plot.gp to set the correct axes 

gnuplot multi_plot.gp
ffmpeg -framerate 10 -i frame_%03d.png -c:v libx264 -r 30 -pix_fmt yuv420p 1_output.mp4
    
10-06-2025 21:35

*/