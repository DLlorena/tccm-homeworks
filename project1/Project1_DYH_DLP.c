#include <trexio.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>


double** malloc_2d(size_t m, size_t n) { // Allocating function provide by the teachers
        double** a = malloc(m*sizeof(double*));
        if (a == NULL) {
                perror("Error allocating the memory");
        return NULL;
	}
        a[0] = malloc(n*m*sizeof(double));
        if (a[0] == NULL) {
                free(a);
        return NULL;
	}
        for (size_t i=1 ; i<m ; i++) {
                a[i] = a[i-1]+n;
        }
        return a;
}

void free_2d(double** a) { // Freeing function provide by the teachers
        free(a[0]);
        a[0] = NULL;
        free(a);
}




int main(){
	
	//OPENING FILE

	trexio_exit_code rc;
	trexio_t* trexio_file = trexio_open("./data/h2o.h5",'r', TREXIO_AUTO, &rc);
	if (rc != TREXIO_SUCCESS) {
		printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
		exit(1);
	}


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

	//double** one_int = malloc_2d(mo_num, mo_num); //works passing the values of one_int
	//malloc_2d functions reutrns a pointer to a pointer, then one_int is a doueble*** (pointer to pointer to pointer toa  dobule)

     	 double* one_int = malloc(mo_num * mo_num * sizeof(double));
	 //works flattening a introducing directly a pointer to a pointer. malloc function returns a pointer, then one_int is a pointei to a pointer to a double.


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

	//printf("Number of integrals non-zero integrals = %lld", n_integrals);


	trexio_exit_code trexio_read_mo_2e_int_eri(trexio_t* const file,
                                           const int64_t offset_file,
                                           int64_t* const buffer_size,
                                           int32_t* const index,
                                           double* const value);
	
	int64_t buffer_side = n_integrals;
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
         rc = trexio_read_mo_2e_int_eri(trexio_file, 0, &buffer_side, index, value);
         if (rc != TREXIO_SUCCESS) {
                printf("TREXIO Error reading two electron integrals:\n%s\n",
                trexio_string_of_error(rc));
                exit(1);
         }

	// PRINTS INDEXES OF THE INTEGRALS (EACH 4 IS AN INTEGRAL) AND INTEGRALS
/*	
	printf("buffer_side %lld",buffer_side);
	for (int i=0; i<n_integrals*4+1; i++){
		printf(" %d ",index[i]);
	}
*/
	


/*
        for (int i=0; i < n_integrals+1; i++){
                printf("%f   ", value[i]);
        }					
*/


	printf("repulsion energy = %lf\n",energy);
	printf("occupied orbitals = %d\n",n_up);
	printf("total molecular orbitals = %d\n",mo_num);

/*
	// PRINTS ONE ELECTRON INTEGRALS

	printf("\n one electron integrals \n");
		
	for (int i = 0; i < mo_num*mo_num; i++) {
            	printf("%f  ", one_int[i]);
		
			
	}

*/	

	//PRINTS INTEGRAL NUMBER n AND ITS FOUR BASIS FUNCTIONS
	int n;	
	int i = index[4*n+0];
	int j = index[4*n+1];
	int k = index[4*n+2];
	int l = index[4*n+3];
	double integral = value[n];	

	printf("integral_%d (%d, %d, %d ,%d) = %f\n]",n,i+1,j+1,k+1,l+1,value[n] );
	
	rc = trexio_close(trexio_file);
	if (rc != TREXIO_SUCCESS) {
  	printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
	exit(1);
	 }
	trexio_file = NULL;
	
	//END OF READING VARIABLES

	//HARTREE FOCK ALGORITHM


	double etotal;
	double one_e_sum = 0.0;
	double two_e_sum = 0.0;

		//ONE ELECTRON TERM CALCULATION

	for (int i=0; i < n_up; i++){
		
		printf("integral = %f \n", one_int[i*mo_num+i]);
		one_e_sum += one_int[i*mo_num + i];
	}

	printf("\n one electron integral sum = %f \n" , one_e_sum);


		//TWO ELECTRON TERM COMPUTATION
	
	for (int n=0; n < n_integrals; n++){
			
		int i = index[4*n+0];
       	 	int j = index[4*n+1];
        	int k = index[4*n+2];
        	int l = index[4*n+3];
		
					
			if (i < n_up && j < n_up && k < n_up && l < n_up){
			 	
			//	printf("Integral %d: (i,j,k,l) = (%d,%d,%d,%d), Value                                                     = %f\n",n, i, j, k, l, value[n]);				
				if (i == k && j == l){
  					printf("Integral_nosym %d: (i,j,k,l) = (%d,%d,%d,%d), Value 							= %f\n",n, i, j, k, l, value[n]);
				
					two_e_sum += 2.0*value[n];
				}	
				if (i == l && j == k){
					printf("Integral_sym %d: (i,j,k,l) = (%d,%d,%d,%d), Value                                                   = %f\n",n, i, j, k, l, value[n]);
                                	two_e_sum -= value[n];
                        	} 	  		
			}
	}	


	printf("\n two electron integral sum = %f \n" , two_e_sum);

	//FINAL ENERGY
	
	double e;
	
	e = energy + 2.0*one_e_sum + two_e_sum;
	printf("\n FINAL ENERGY = %f \n", e);
		
	



}


