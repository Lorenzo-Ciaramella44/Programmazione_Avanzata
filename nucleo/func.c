#include "nucleo.h"

void decay_law(int *count, double *t2, const char* nome, const char* nome2) {
	FILE *f, *f2;
	int N_0, index;
	double tau, t;
	int counter_1 = 0, counter_1_5 = 0, counter_2 = 0;
	int counter_4 = 0, counter_8 = 0, counter_20 = 0;
	
	N_0 = *count;
	tau = *t2;
	
	f = fopen(nome, "w");
	f2 = fopen(nome2, "a");
	
	if (f == NULL || f2 == NULL) {
        perror("Errore apertura file");
        exit(EXIT_FAILURE);
    }
    
	printf("Application function 'decay_law'.....\n");
	
	for (int i = 0; i < N_0; i++) {
		t = -tau*log(1.0-rand()/(RAND_MAX+1.0));
		//printf("%lf \n", t);
		if (t < 1.0) {
			counter_1++;
		}
		if (t < 1.5) {
			counter_1_5++;	
		}
		if (t < 2.0) {
			counter_2++;	
		}
		if (t < 4.0) {
			counter_4++;	
		}
		if (t < 8.0) {
			counter_8++;	
		}
		if (t < 20.0) {
			counter_20++;	
		}
	}
	/* decayed particles are registred */
	fprintf(f, "%d %d %d %d %d %d", counter_1, counter_1_5, counter_2, counter_4, counter_8, counter_20);
	fprintf(f2, "%d %d %d %d %d %d\n", counter_1, counter_1_5, counter_2, counter_4, counter_8, counter_20);
	fclose(f);
	fclose(f2);
	printf("Ending function 'decay_law'\n");
}


void check_variance(int *count, int *sim, double *t2, const char* nome, const char* nome2, const char* nome3) {
	FILE *f, *f2, *f3;
	double t[6] = {1.0, 1.5, 2.0, 4.0, 8.0, 20.0};
	double N_th[6];
	double var_th[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	double NUM[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	double NUM2[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	double mu_real[6], var_real[6];
	int N_0, n_sim, N_tot;
	double tau;
	
	N_0 = *count;
	tau = *t2;
	n_sim = *sim;
	N_tot = n_sim*6;
	
	f = fopen(nome, "r");
	f2 = fopen(nome2, "w");
	f3 = fopen(nome3, "w");
	
	int *dati_N = (int *)malloc(N_tot * sizeof(int));
    if (dati_N == NULL) { 
    	printf("Error Memory Allocation\n");
   		fclose(f2);
   		fclose(f3);
   	}
    	
    char line_N[256];	// Buffer per leggere una riga dal file
    size_t k_N = 0;
    //printf("Reading file_N\n");
	while (fgets(line_N, sizeof(line_N), f) != NULL) { 
		char *ptr_N = line_N;
		//printf("Riga letta: %s", line_N); // Stampa di debug
		while (*ptr_N && k_N < N_tot) {
			if (sscanf(ptr_N, "%d", &dati_N[k_N]) == 1) {
				//printf("Dato letto: %f\n", dati[k]); // Stampa di debug
				k_N++; 
			} 
			else {
				int result = sscanf(ptr_N, "%d", &dati_N[k_N]);
				printf("aaaaa %d \n", result);
			}
			// Sposta il puntatore alla posizione successiva dopo il numero letto 
			while (*ptr_N && *ptr_N != ' ' && *ptr_N != '\t' && *ptr_N != '\n') {
				ptr_N++;
			}
			while (*ptr_N && (*ptr_N == ' ' || *ptr_N == '\t' || *ptr_N == '\n')) { 
				ptr_N++;
			}
		} 
	}
	
	/* particle decayed */
	for (int i = 0; i < 6; i++) {
		N_th[i] = N_0*(1.0-exp(-t[i]/tau));
		//printf("%lf %lf\n", t[i], N_th[i]);
		for (int j = 0; j < N_tot; j+=6) { 
			NUM[i] = NUM[i] + dati_N[j+i];
			NUM2[i] = NUM2[i] + (N_th[i]-dati_N[j+i])*(N_th[i]-dati_N[j+i]);
			var_th[i] = var_th[i] + (N_th[i]-dati_N[j+i])*(N_th[i]-dati_N[j+i]);
			 
			//printf("aaa %lf %lf %d %d %d %d\n", NUM[i], NUM2[i], dati_N[j+i], i, j, i+j);
		}
		var_th[i] = var_th[i]*N_th[i]/N_0/n_sim;
		mu_real[i] = NUM[i]/n_sim;
		var_real[i] = NUM2[i]/n_sim;
		printf("\n");
		printf("Expected Value Theoric: %lf, Std Deviation: %lf\n", N_th[i], sqrt(var_th[i]));
		printf("Expected Value Simulated: %lf, Std Deviation: %lf\n",	mu_real[i], sqrt(var_real[i]));
		printf("\n");
		fprintf(f2, "%d %d %lf %lf %lf %lf\n", n_sim/2+50, n_sim/2-50, N_th[i], sqrt(var_th[i]), mu_real[i], sqrt(var_real[i]));
	}
	
	for (int j = 0; j < N_tot; j+=6) { 
		fprintf(f3, "%d %d %d %d %d %d %d\n", j/6, dati_N[j], dati_N[j+1], dati_N[j+2], dati_N[j+3], dati_N[j+4], dati_N[j+5]);
	}
	
	fclose(f3);
	fclose(f2);
	fclose(f);
	free(dati_N);
}

const char* rename_file(const char *filename, int index) {
    static char nuovo_nome[50];
    sprintf(nuovo_nome, "sim_%d.txt", index);
    
    if (rename(filename, nuovo_nome) != 0) {
        perror("Error renaming file");
        return NULL;
    }
    return nuovo_nome;
}