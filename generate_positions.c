#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define FILENAME     "data/monomer_positions.txt" 
#define NUM_MONOMERS 10
#define NUM_STEPS    100

enum coord{X, Y, Z};

double random_theta(){
	return (double) rand() / RAND_MAX * M_PI;
}

double random_phi(){
	return (double) rand() / RAND_MAX * M_PI * 2.0;
}

int random_monomer_index(){
	return rand() % (NUM_MONOMERS - 1);
}

void saveFile(char *filename, double **positions, char *mode){
	
	// Open file
	FILE* fptr = fopen(filename, mode);
	if (fptr == NULL){
		printf("Unable to open file: %s.\n", filename);
		return;
	}

	// Write to file
	for (int i = 0; i < NUM_MONOMERS; i++)
		for (int j = 0; j <= 2; j++)
			fprintf(fptr, "%lf, ", positions[i][j]);
	fprintf(fptr, "0.0 \n");

	// Close file
	fclose(fptr);
	return;
}

void saveFirstLineToFile(char* filename, double **positions){

	saveFile(filename, positions, "w+");
	return;
}

void saveNextLineToFile(char* filename, double **positions){

	saveFile(filename, positions, "a");
	return;
}

int main(){
	
	// Seed rand
	srand(time(NULL));
	
	// Declare my friendly variables
	int index, next_index;
	double theta, phi;
	double cos_t, sin_t;
	double cos_p, sin_p;
	double vec_x, vec_y, vec_z;
	double disp_x, disp_y, disp_z;
	
	// Allocate position array
	double **positions = (double**)malloc(NUM_MONOMERS * sizeof(double*));
	for (index = 0; index < NUM_MONOMERS; index++)
		positions[index] = (double*)malloc(3 * sizeof(double));
	
	// Randomly generate polymer I.C.s
	positions[0][X] = 0.0;
	positions[0][Y] = 0.0;
	positions[0][Z] = 0.0;
	
	for (index = 1; index < NUM_MONOMERS; index++){
		
		// Pick random angles
		theta = random_theta();
		phi   = random_phi();
		
		cos_t = cos(theta);
		sin_t = sin(theta);
		cos_p = cos(phi);
		sin_p = sin(phi);
		
		positions[index][X] = positions[index - 1][X] + sin_t * cos_p;
		positions[index][Y] = positions[index - 1][Y] + sin_t * sin_p;
		positions[index][Z] = positions[index - 1][Z] + cos_t;
	}
	
	// Save positions
	saveFirstLineToFile(FILENAME, positions);
	
	// Random pivot sampling
	for (long step = 0; step < NUM_STEPS; step++){
		
		// Pick index
		index = random_monomer_index();
		
		// Pick random angles
		theta = random_theta();
		phi   = random_phi();
		
		cos_t = cos(theta);
		sin_t = sin(theta);
		cos_p = cos(phi);
		sin_p = sin(phi);
		
		vec_x = positions[index][X] + sin_t * cos_p;
		vec_y = positions[index][Y] + sin_t * sin_p;
		vec_z = positions[index][Z] + cos_t;
		
		// Determine displacement and set new vectors
		disp_x = vec_x - positions[index + 1][X];
		disp_y = vec_y - positions[index + 1][Y];
		disp_z = vec_z - positions[index + 1][Z];
		
		// Propagate displacement
		for (next_index = index + 1; next_index < NUM_MONOMERS; next_index++){
			positions[next_index][X] += disp_x;
			positions[next_index][Y] += disp_y;
			positions[next_index][Z] += disp_z;
			
			// Save positions
			saveNextLineToFile(FILENAME, positions);
		}
	}
}
