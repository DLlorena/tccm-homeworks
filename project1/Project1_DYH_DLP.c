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

////////////////////////////////////////////////////////////////////////////////////////////////
/*
This programs intends to calculate Hartree Fock energy and MP2 correction for specific molecules. The data for each molecule is stored in HDF5 format files and the data is extracted using trexio libraries
*/

int main(int argc, char *argv[]){	//program is executed with arguments needed, being the first argument argv[0] the program itself
	
	if (argc < 2) { 	//if less than two terms (executable + file) are passed the programs returns error
    		printf("Usage: %s molecule file missing\n", argv[0]);
    	return 1;
	}

	if (argc > 2) {		//if more than two terms (executable + file) are passed the programs returns error
    		printf("Usage: %s more than one molecule file was provided\n", argv[0]);
    	return 1;
	}
	
	//Defining the path where molecule files are stored

	char *molecule_file = argv[1];			//Declaring pointer to a string variable to store the second argument (name of molecule file)
	char *molecules_directory = "./data/";		//Directory path where molecule files are stores
	char filepath[256]; 				//Declares a variable large enough to store the whole path of molecule files

	snprintf(filepath,sizeof(filepath), "%s%s", molecules_directory, argv[1]);		//Merges the name of the directory and the molecule file creating the complete path

	char output_filename[256];			//Declares a file to write the outputs 
	snprintf(output_filename,sizeof(output_filename), "%s.out",  argv[1]); 			//Defines the output file name by merging input molecule file and ".out"
	
	FILE *output = fopen(output_filename, "w");	//Opens file in writting mode

	if (output == NULL){				//Handles error opening the file
	        printf("Error creating output %s.\n",output_filename );
		return 1;
	}

	
	//Writes program logo and general info in the output file

	fprintf(output," ____  _  _  ____  __                  _  _  ____     _     _  _  ____  ____ \n");	
    	fprintf(output,"(    \\( \\/ )(    \\(  )     ___  ___   / )( \\(  __)   ( )   ( \\/ )(  _ \\(___ \\ \n");
    	fprintf(output," ) D ( )  /  ) D (/ (_/\\  (___)(___)  ) __ ( ) _)   (_ _)  / \\/ \\ ) __/ / __/ \n");
    	fprintf(output,"(____/(__/  (____/\\____/              \\_)(_/(__)     (_)   \\_)(_/(__)  (____) \n");	

	fprintf(output,"\n");
	fprintf(output, "----------------------Calculations performed for %s----------------------\n", argv[1]);

	fprintf(output,"\n#All results are presented in atomic units \n");
 
	////////////////////////////////////////EXTRACTING DATA SECTION//////////////////////////////////////////////////


	//OPENING HDF5 FILE FOR THE MOLECULE OF INTEREST
	
	trexio_exit_code rc; 	//Variable where return code of trexio functions is store

	trexio_t* trexio_file = trexio_open(filepath,'r', TREXIO_AUTO, &rc);		//Calls trexio function to open and read molecule file using the already defined path

	//If returning code does not return succes the program displays an error refering to this function

	if (rc != TREXIO_SUCCESS) {		
		printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
		exit(1);
	}

	printf("Reading file %s succesfully\n", filepath);	//Specifies which file is the program reading

	//READING NUCLEUS REPULSIONS

	//Function declaration for reading nuclear repulsion. The function takes a pointer to the TREXIO file and a pointer to the energy variable where the result will be stored.
	trexio_exit_code trexio_read_nucleus_repulsion(trexio_t* const trexio_file, double* const energy);

	//variabel where repulsion energy is stored
	double energy;

	//Calls function using trexio file and energy adress as arguments 											
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
	int n_up;
	rc = trexio_read_electron_up_num(trexio_file, &n_up);
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

	//Allocates memory for a mo_num x mo_num size 1D array of doubles to store one electron integrals
     	double* one_int = malloc(mo_num * mo_num * sizeof(double));

        rc = trexio_read_mo_1e_int_core_hamiltonian(trexio_file, one_int);
        if (rc != TREXIO_SUCCESS) {
        	printf("TREXIO Error reading one electron integrals:\n%s\n",
        	trexio_string_of_error(rc));
        	exit(1);
        }


	//READING TWO ELECTRON INTEGRALS

	//Reads the number of non-zero two electron integrals
	int64_t n_integrals;
	
	trexio_exit_code trexio_read_mo_2e_int_eri_size(trexio_t* const trexio_file, int64_t* const n_integrals);
	rc = trexio_read_mo_2e_int_eri_size(trexio_file, &n_integrals);	
	if (rc != TREXIO_SUCCESS) {
                printf("TREXIO Error reading number of non zero two electron integrals:\n%s\n",
                trexio_string_of_error(rc));
                exit(1);
	}

	//Function declaration to read two electron integrals. Takes trexio file pointer, offset, amount of integrals, indexes anda values of the integrals.
	trexio_exit_code trexio_read_mo_2e_int_eri(trexio_t* const file,
                                           const int64_t offset_file,
                                           int64_t* const buffer_size,
                                           int32_t* const index,
                                           double* const value);

	//Buffer size is initialized as the number of two electron integrals.	
	int64_t buffer_size = n_integrals;

	//Allocates memory for a 1D array to store the 4 idexes pero integral value. Returns error if allocation fails.
	int32_t* const index = malloc(4 * n_integrals * sizeof(int32_t));
	if (index == NULL) {
  		fprintf(stderr, "Malloc failed for index");
		exit(1); 
	}
	
	//Allcoates memory for a 1D array to store integral values. Return error if allocation fails.
	double* const value = malloc(n_integrals * sizeof(double));
	if (value == NULL) {
  		fprintf(stderr, "Malloc failed for value");
		exit(1);
	 }

	//Calls function to read two electron integrals
        rc = trexio_read_mo_2e_int_eri(trexio_file, 0, &buffer_size, index, value);
        if (rc != TREXIO_SUCCESS) {
                printf("TREXIO Error reading two electron integrals:\n%s\n",
                trexio_string_of_error(rc));
                exit(1);
         }
		

	//PREPARATION OF A 4D ARRAY TO STORE TWO ELECTRON INTEGRALS
	
	//Declares size of each of the 4 dimensions 
	int size_i, size_j, size_k, size_l;

	//Calls function to allocate 4D array with a size of mo_num for each dimension.
	double ****two_e_int = allocate_two_e_int(mo_num, mo_num, mo_num, mo_num);
	
	//Stores integrals in 4D array. Loops through integral values, gets the fours indexes of it and stores the value in the corret place of the 4D matrix following 8-fold simmetry of the two electon integrals
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

	//Displays in terminal and writes in the output the read data from trexio file
	
	printf("\nrepulsion energy = %lf\n",energy);
	fprintf(output,"\nrepulsion energy = %lf\n",energy);	

	printf("\noccupied orbitals = %d\n",n_up);
        fprintf(output,"\noccupied orbitals = %d\n",n_up);

	printf("\ntotal molecular orbitals = %d\n",mo_num);
        fprintf(output,"\ntotal molecular orbitals = %d\n",mo_num);

	
	//CLOSES FILE	

	rc = trexio_close(trexio_file);
	if (rc != TREXIO_SUCCESS) {
  	printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
	exit(1);
	 }
	trexio_file = NULL;
	
	//END OF READING DATA

	//HARTREE FOCK ENEGRY COMPUTATION


	double etotal;
	double one_e_sum;
	double two_e_sum;

		//ONE ELECTRON TERM
	
	//Sums all the diagonal terms for the occupied orbitals. "i*mo_num + i" is the equivalent in a 1D array of size mo_num.
	for (int i=0; i < n_up; i++){
		one_e_sum += one_int[i*mo_num + i];
	}

	//Displays in terminal and writes in output file
	printf("\none electron integral term  = %f \n" ,2*one_e_sum);
        fprintf(output,"\none electron integral term  = %f \n" ,2*one_e_sum);

		//TWO ELECTRON TERM

	//Loop trough occupied orbitals and computes two electron term.	
	for (int i=0; i < n_up; i++){
		for (int j=0; j<n_up; j++){
			two_e_sum += 2 * two_e_int[i][j][i][j]  - two_e_int[i][j][j][i];	//4D array allows indexing by its 4 dimensions.	

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
	
	//Once obtained one and two electron terms, the final energy is computed.
	e = energy + 2.0*one_e_sum + two_e_sum;
	printf("\nHF ENERGY = %f \n", e);
        fprintf(output,"\nHF ENERGY = %f \n", e);

	//Frees memory from 1D arrays
	free(value);
	free(index);
	free(one_int);
	
	//MP2 ENERGIES COMPUTATION
	
	double e_mp2;

	//Loops through occupied orbitals i and j and virtual orbitals a and b. 
	for (int i=0; i < n_up; i++){
        	        for (int j=0; j<n_up; j++){
                	        for (int a=n_up; a < mo_num; a++){
                        	        for (int b=n_up; b < mo_num; b++){
						
						//Computes MP2 energy correction by indexing two electron integrals and orbital energies.
        	                                e_mp2 += (two_e_int[i][j][a][b] * (2*two_e_int[i][j][a][b] - two_e_int[i][j][b][a])) / (mo_energy[i] + mo_energy[j] - mo_energy[a] - mo_energy[b]);
 
                	                }
                        	}
			}	
	}
	printf("\nMP2 ENERGY CORRECTION = %.8lf \n", e_mp2);		 	 
        fprintf(output,"\nMP2 ENERGY CORRECTION = %.8lf \n", e_mp2);
	

	//Calls function to free memory of a 4D array.
	free(two_e_int);			
	
	fprintf(output,"\nHARTREEFOCKATION PERFORMED SUCCESFULLY");
	
}


