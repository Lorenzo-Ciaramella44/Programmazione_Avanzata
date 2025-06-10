#include "nucleo.h"


double main(int argc, char *argv[]) {
	clock_t start, end;		/* variables for seeing time required*/
    double cpu_time_used;
	FILE *f;
    int N_0 = 1000;
	int n_sim;
	double tau = 4.5, T = 25.0;
	const char* sim = "sim.txt";
	const char* sim_joint = "sim_joint.txt";
	
	start = clock();			/* time starts */
	
	srand(time(NULL));
	
	n_sim = atoi(argv[1]);
	
	for (int i = 0; i < n_sim; i++) {
		decay_law(&N_0, &tau, sim, sim_joint);
		sim = rename_file(sim, i);
		sim = "sim.txt";
	}
	
	check_variance(&N_0, &n_sim, &tau, sim_joint, "result.txt", "final.txt");
	
	end = clock();				/* time ends */
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Elaboration time with parallelization: %f seconds\n", cpu_time_used);
}