#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define FILENAME     "monomer_positions.txt" 
#define NUM_MONOMERS 20
#define NUM_STEPS    1000

double random_angle(){
	return (double) rand() / RAND_MAX * 2 * M_PI;
}

int random_monomer(){
	return rand() % (NUM_MONOMERS - 1);
}

void saveFirstLineToFile(char* filename, double **positions){

	// Open file
	FILE* fptr = fopen(filename, "w+");
	if (fptr == NULL){
		printf("Unable to open file: %s.\n", filename);
		return;
	}

	// Append to file
	for (int i = 0; i < NUM_MONOMERS; i++)
		for (int j = 0; j < 3; j++)
			fprintf(fptr, "%lf, ", positions[i][j]);
	fprintf(fptr, "0.0 \n");

	/*
	NOTE THAT THE LINES END WITH COMMAS, WHICH WE DO NOT WANT!
	FOR THE MOMENT, I THROW IN A 0.0 SO THAT np.loadtxt()
	DOES NOT COMPLAIN!
	*/

	// Close file and exit
	fclose(fptr);
	return;
}

void saveLineToFile(char* filename, double **positions){

	// Open file
	FILE* fptr = fopen(filename, "a");
	if (fptr == NULL){
		printf("Unable to open file: %s.\n", filename);
		return;
	}

	// Append to file
	for (int i = 0; i < NUM_MONOMERS; i++)
		for (int j = 0; j < 3; j++)
			fprintf(fptr, "%lf, ", positions[i][j]);
	fprintf(fptr, "0.0 \n");

	/*
	NOTE THAT THE LINES END WITH COMMAS, WHICH WE DO NOT WANT!
	FOR THE MOMENT, I THROW IN A 0.0 SO THAT np.loadtxt()
	DOES NOT COMPLAIN!
	*/

	// Close file and exit
	fclose(fptr);
	return;
}

int main(){
	
	// Seed rand
	srand(time(NULL));
	
	// Allocate position array
	double **positions = (double**)malloc(NUM_MONOMERS * sizeof(double*));
	for (int I = 0; I < NUM_MONOMERS; I++)
		positions[I] = (double*)malloc(3 * sizeof(double));
	
	// Declare my friendly variables
	int I;
	double a = 1.0;
	double angle;
	double disp_x, disp_y;
	double cosine, sine;
	double vector_x, vector_y;
	double rot_vec_x, rot_vec_y;
	
	// Randomly generate polymer I.C.s
	positions[0][0] = 0.0;
	positions[0][1] = 0.0;
	positions[0][2] = 0.0;
	
	for (I = 1; I < NUM_MONOMERS; I++){
		
		angle = random_angle();
		positions[I][0] = positions[I-1][0] + a * cos(angle);
		positions[I][1] = positions[I-1][1] + a * sin(angle);
		positions[I][2] = 0.0;
	}
	
	// Save positions
	saveFirstLineToFile(FILENAME, positions);
	
	// Random pivot sampling
	for (long step = 0; step < NUM_STEPS; step++){
		
		// Pick random angle and index
		I = random_monomer();
		angle = random_angle();
		cosine = cos(angle);
		sine = sin(angle);
		
		// Retrieve existing monomer-to-monomer vector
		vector_x = positions[I+1][0] - positions[I][0];
		vector_y = positions[I+1][1] - positions[I][1];
		
		// Perform rotation
		rot_vec_x = cosine * vector_x -   sine * vector_y;
		rot_vec_y =   sine * vector_x + cosine * vector_y;
		
		// Compute displacement
		disp_x = rot_vec_x - vector_x;
		disp_y = rot_vec_y - vector_y;
		
		// Propagate displacement
		for (int index = I+1; index < NUM_MONOMERS; index++){
			positions[index][0] += disp_x;
			positions[index][1] += disp_y;
			
			// Save positions
			saveLineToFile(FILENAME, positions);
		}
	}
}
