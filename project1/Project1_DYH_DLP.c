#include <trexio.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>


// Function to dynamically allocate a 4D array of size [size_i][size_j][size_k][size_l]

double ****allocate_two_e_int(int size_i, int size_j, int size_k, int size_l) {

    double ****array = (double ****)malloc(size_i * sizeof(double ***));

    for (int i = 0; i < size_i; i++) {
        array[i] = (double ***)malloc(size_j * sizeof(double **)); // size_j for second dimension
        for (int j = 0; j < size_j; j++) {
            array[i][j] = (double **)malloc(size_k * sizeof(double *)); // size_k for third dimension
            for (int k = 0; k < size_k; k++) {
                array[i][j][k] = (double *)malloc(size_l * sizeof(double)); // size_l for fourth dimension
            }
        }
    }

    return array;
}

// Function to free a 4D array
void free_integral_array(double ****array, int size_i, int size_j, int size_k, int size_l) {
    for (int i = 0; i < size_i; i++) {
        for (int j = 0; j < size_j; j++) {
            for (int k = 0; k < size_k; k++) {
                free(array[i][j][k]);
            }
            free(array[i][j]);
        }
        free(array[i]);
    }
    free(array);
}


int main(int argc, char *argv[]){
	
	if (argc < 2) {
    		printf("Usage: %s molecule file missing\n", argv[0]);
    	return 1;
	}

	if (argc > 2) {
    		printf("Usage: %s more than one molecule file was provided\n", argv[0]);
    	return 1;
	}

	char *molecule_file = argv[1];
	char *molecules_directory = "./data/";
	char filepath[256];

	snprintf(filepath,sizeof(filepath), "%s%s", molecules_directory, argv[1]);

	char output_filename[256];
	snprintf(output_filename,sizeof(output_filename), "%s.out",  argv[1]);
	
	FILE *output = fopen(output_filename, "w");

	if (output == NULL){
	        printf("Error creating output %s.\n",output_filename );
		return 1;
	}


	fprintf(output," ____  _  _  ____  __                  _  _  ____     _     _  _  ____  ____ \n");	
    	fprintf(output,"(    \\( \\/ )(    \\(  )     ___  ___   / )( \\(  __)   ( )   ( \\/ )(  _ \\(___ \\ \n");
    	fprintf(output," ) D ( )  /  ) D (/ (_/\\  (___)(___)  ) __ ( ) _)   (_ _)  / \\/ \\ ) __/ / __/ \n");
    	fprintf(output,"(____/(__/  (____/\\____/              \\_)(_/(__)     (_)   \\_)(_/(__)  (____) \n");	

	fprintf(output,"\n");
	fprintf(output, "----------------------Calculations performed for %s----------------------\n", argv[1]);

	fprintf(output,"\n#All results are presented in atomic units \n");
 
	//OPENING FILE

	trexio_exit_code rc;
	trexio_t* trexio_file = trexio_open(filepath,'r', TREXIO_AUTO, &rc);
	if (rc != TREXIO_SUCCESS) {
		printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
		exit(1);
	}
	printf("Reading file %s succesfully\n", filepath);

	//READING NUCLEUS REPULSIONS

	trexio_exit_code trexio_read_nucleus_repulsion(trexio_t* const trexio_file, 									double* const energy);
	double energy; // Variable where the energy is read
	rc = trexio_read_nucleus_repulsion(trexio_file, &energy);
	// Check the return code to be sure reading was OK
	if (rc != TREXIO_SUCCESS) {
  		printf("TREXIO Error reading nuclear repulsion energy:\n%s\n",
        	trexio_string_of_error(rc));
		exit(1);
	 }
		
	//READING NUMBER OF OCCUPIED MOLECULAR ORBITALS

	trexio_exit_code trexio_read_electron_up_num(trexio_t* const trexio_file,
                                             int32_t* const n_up);
	int n_up; //variable where number of occupied orbitals is read
	rc = trexio_read_electron_up_num(trexio_file, &n_up);
        // Check the return code to be sure reading was OK
        if (rc != TREXIO_SUCCESS) {
        	printf("TREXIO Error reading number of up-spin electrons:\n%s\n",
        	trexio_string_of_error(rc));
        	exit(1);
         }


	//READING NUMBER OF MOLECULAR ORBITALS

	trexio_exit_code trexio_read_mo_num(trexio_t* const trexio_file,
                                             int32_t* const mo_num);
        int mo_num; //variable where number of molecular orbitals is read
        rc = trexio_read_mo_num(trexio_file, &mo_num);
        if (rc != TREXIO_SUCCESS) {
        	printf("TREXIO Error reading numer of molecular orbitals:\n%s\n",
        	trexio_string_of_error(rc));
        	exit(1);
         }



	//READING ONE ELECTRON INTEGRALS
	

	 trexio_exit_code trexio_read_mo_1e_int_core_hamiltonian(trexio_t* const trexio_file, double* const data);

     	 double* one_int = malloc(mo_num * mo_num * sizeof(double));

         rc = trexio_read_mo_1e_int_core_hamiltonian(trexio_file, one_int);
         if (rc != TREXIO_SUCCESS) {
        	printf("TREXIO Error reading one electron integrals:\n%s\n",
        	trexio_string_of_error(rc));
        	exit(1);
         }


	//READING TWO ELECTRON INTEGRALS


	int64_t n_integrals;
	
	trexio_exit_code trexio_read_mo_2e_int_eri_size(trexio_t* const trexio_file, int64_t* const n_integrals);
	rc = trexio_read_mo_2e_int_eri_size(trexio_file, &n_integrals);	
	if (rc != TREXIO_SUCCESS) {
                printf("TREXIO Error reading number of non zero two electron integrals:\n%s\n",
                trexio_string_of_error(rc));
                exit(1);
	}


	trexio_exit_code trexio_read_mo_2e_int_eri(trexio_t* const file,
                                           const int64_t offset_file,
                                           int64_t* const buffer_size,
                                           int32_t* const index,
                                           double* const value);
	
	int64_t buffer_size = n_integrals;
	int32_t* const index = malloc(4 * n_integrals * sizeof(int32_t));
	if (index == NULL) {
  		fprintf(stderr, "Malloc failed for index");
		exit(1); 
	}
	double* const value = malloc(n_integrals * sizeof(double));
	if (value == NULL) {
  		fprintf(stderr, "Malloc failed for value");
		exit(1);
	 }

	
         rc = trexio_read_mo_2e_int_eri(trexio_file, 0, &buffer_size, index, value);
         if (rc != TREXIO_SUCCESS) {
                printf("TREXIO Error reading two electron integrals:\n%s\n",
                trexio_string_of_error(rc));
                exit(1);
         }
		

	//PREPARATION OF A 4D ARRAY TO STORE TWO ELECTRON INTEGRALS
	
	int size_i, size_j, size_k, size_l;
	double ****two_e_int = allocate_two_e_int(mo_num, mo_num, mo_num, mo_num);
	
	for (int n=0; n<n_integrals; n++){
		int i = index[4*n+0];
        	int j = index[4*n+1];
        	int k = index[4*n+2];
        	int l = index[4*n+3];
		
		two_e_int[i][j][k][l] = value[n];
		two_e_int[k][l][i][j] = value[n];
		two_e_int[k][j][i][l] = value[n];
		two_e_int[i][l][k][j] = value[n];
		two_e_int[j][i][l][k] = value[n];
		two_e_int[l][k][j][i] = value[n];
		two_e_int[j][k][l][i] = value[n];
		two_e_int[l][i][j][k] = value[n];

		//prints all two electron integrals
		//printf("\n integral_two_e(%d,%d,%d,%d) = %f \n",i,j,k,l,two_e_int[i][j][k][l]); 		
	}

	
	
	//READING ORBITAL ENERGIES FOR MP2
	
	trexio_exit_code trexio_read_mo_energy(trexio_t* const file,
                                            double* const mo_energy);
	double* const mo_energy = malloc(mo_num * sizeof(double));
		
	rc = trexio_read_mo_energy(trexio_file, mo_energy);
        if (rc != TREXIO_SUCCESS) {
                printf("TREXIO Error reading molecular orbitals:\n%s\n",
                trexio_string_of_error(rc));
                exit(1);
	}

		
	printf("\nrepulsion energy = %lf\n",energy);
	fprintf(output,"\nrepulsion energy = %lf\n",energy);	

	printf("\noccupied orbitals = %d\n",n_up);
        fprintf(output,"\noccupied orbitals = %d\n",n_up);

	printf("\ntotal molecular orbitals = %d\n",mo_num);
        fprintf(output,"\ntotal molecular orbitals = %d\n",mo_num);

	
	//CLOSE FILE	

	rc = trexio_close(trexio_file);
	if (rc != TREXIO_SUCCESS) {
  	printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
	exit(1);
	 }
	trexio_file = NULL;
	
	//END OF READING VARIABLES

	//HARTREE FOCK ENEGRY COMPUTATION


	double etotal;
	double one_e_sum;
	double two_e_sum;

		//ONE ELECTRON TERM

	for (int i=0; i < n_up; i++){
		one_e_sum += one_int[i*mo_num + i];
	}

	printf("\none electron integral term  = %f \n" ,2*one_e_sum);
        fprintf(output,"\none electron integral term  = %f \n" ,2*one_e_sum);

		//TWO ELECTRON TERM
	
	for (int i=0; i < n_up; i++){
		for (int j=0; j<n_up; j++){
			two_e_sum += 2 * two_e_int[i][j][i][j]  - two_e_int[i][j][j][i];		

		}
	}	

	printf("\ntwo electron integral term = %f \n" , two_e_sum);
        fprintf(output,"\ntwo electron integral term = %f \n" , two_e_sum);

/*
	//PRINTS A SPECIFIC TWO ELECTRON INTEGRAL INITIALIZING i, j, k, l

	int i = 1;
	int j = 2;
	int k = 3;
	int l = 4;

	printf("two electron integral (%d,%d,%d,%d) = %f\n", i,j,k,l, two_e_int[i][j][k][l]);
		
*/
	//HF FINAL ENERGY
	
	double e;
	
	e = energy + 2.0*one_e_sum + two_e_sum;
	printf("\nHF ENERGY = %f \n", e);
        fprintf(output,"\nHF ENERGY = %f \n", e);

		
	free(value);
	free(index);
	free(one_int);
	
	//MP2 ENERGIES COMPUTATION
	
	double e_mp2;

	for (int i=0; i < n_up; i++){
        	        for (int j=0; j<n_up; j++){
                	        for (int a=n_up; a < mo_num; a++){
                        	        for (int b=n_up; b < mo_num; b++){
	
        	                                e_mp2 += (two_e_int[i][j][a][b] * (2*two_e_int[i][j][a][b] - two_e_int[i][j][b][a])) / (mo_energy[i] + mo_energy[j] - mo_energy[a] - mo_energy[b]);
 
                	                }
                        	}
			}	
	}
	printf("\nMP2 ENERGY CORRECTION = %.8lf \n", e_mp2);		 	 
        fprintf(output,"\nMP2 ENERGY CORRECTION = %.8lf \n", e_mp2);
	


	free(two_e_int);			
	
	fprintf(output,"\nHARTREEFOCKATION PERFORMED SUCCESFULLY");
	
}


