#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define FILENAME_RAW "output/distances_raw_1000.txt"
#define FILENAME_AVG "output/distances_avg_1000.txt"
#define NUM_MONOMERS 1000
#define NUM_STEPS    1000000

double random_angle(){
	return (double) rand() / RAND_MAX * 2 * M_PI;
}

int random_monomer(){
	return rand() % (NUM_MONOMERS - 1);
}

void saveToFile(char* filename, long size, long step, double *distances){

	// Open file
	FILE* fptr = fopen(filename, "w+");
	if (fptr == NULL){
		printf("Unable to open file: %s.\n", filename);
		return;
	}

	// Write to file
	for (long i = 0; i < size; i++)
		fprintf(fptr, "%ld, %lf\n", i * step, distances[i]);

	fclose(fptr);
	return;
}

int main(){
	
	// Seed rand
	srand(time(NULL));
	
	// Allocate position array
	double **positions = (double**)malloc(NUM_MONOMERS * sizeof(double*));
	for (int I = 0; I < NUM_MONOMERS; I++)
		positions[I] = (double*)malloc(2 * sizeof(double));
	
	// Randomly generate polymer I.C.s
	positions[0][0] = 0;
	positions[0][1] = 0;
	
	double a = 1.0;
	double angle;
	
	for (int I = 1; I < NUM_MONOMERS; I++){
		
		angle = random_angle();
		positions[I][0] = positions[I-1][0] + a * cos(angle);
		positions[I][1] = positions[I-1][1] + a * sin(angle);
	}
	
	int I;
	double disp_x, disp_y;
	double cosine, sine;
	double vector_x, vector_y;
	double rot_vec_x, rot_vec_y;
	double mag;
	
	// End to end distance array
	double *distances = (double*)malloc(NUM_STEPS * sizeof(double));
	
	// Metropolis sampling
	for (long step = 0; step < NUM_STEPS; step++){
		
		// Pick random angle and index
		I = random_monomer();
		angle = random_angle();
		cosine = cos(angle);
		sine = sin(angle);
		
		vector_x = positions[I+1][0] - positions[I][0];
		vector_y = positions[I+1][1] - positions[I][1];
		
		// Perform rotation
		rot_vec_x = cosine * vector_x -   sine * vector_y;
		rot_vec_y =   sine * vector_x + cosine * vector_y;
		
		// Compute displacement
		disp_x = rot_vec_x - vector_x;
		disp_y = rot_vec_y - vector_y;
		
		for (int index = I+1; index < NUM_MONOMERS; index++){
			positions[index][0] += disp_x;
			positions[index][1] += disp_y;
		}
		
		// Compute end to end distance
		vector_x = positions[NUM_MONOMERS-1][0] - positions[0][0];
		vector_y = positions[NUM_MONOMERS-1][1] - positions[0][1];
		mag = pow(pow(vector_x, 2) + pow(vector_y, 2), 0.5);
		distances[step] = mag;
	}
	
	// Save distances to file
	saveToFile(FILENAME_RAW, NUM_STEPS, 1, distances);
	
	// Compute the distances every N steps
	long size = NUM_STEPS / NUM_MONOMERS;
	double *distance_every_N = (double*)malloc(size * sizeof(double));
	for (long i = 0; i < size; i++)
		distance_every_N[i] = distances[i * NUM_MONOMERS];
	
	// Cumulate sum and average
	for (long i = size - 1; i > 0; i--){
		for (long j = 0; j < i; j++){
			distance_every_N[i] += distance_every_N[j] * distance_every_N[j];
		}
		distance_every_N[i] =  sqrt(distance_every_N[i] / (i + 1));
	}
	
	// Save distances to file
	saveToFile(FILENAME_AVG, size, NUM_MONOMERS, distance_every_N);
	return 0;
}
