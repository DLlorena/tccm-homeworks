#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Definition of the functions
size_t read_Natoms(FILE* file){   // Reads the number of atoms from the input
	size_t Natoms;
	fscanf(file, "%zu", &Natoms); // %zu = Read the int without the sign
	return Natoms;
}

void read_molecule(FILE* file, // This functions reads the coordinates (x, y, z) and the masses
	size_t Natoms,
	double** coord,
	double* mass){

	for (size_t i = 0; i < Natoms; i++){
		fscanf(file, "%lf %lf %lf %lf",
				&coord[i][0], &coord[i][1], &coord[i][2], &mass[i]);
	}
}

void compute_distances(size_t Natoms, // This function calculates the distance between atoms, calculates the modulus of a vector
	double** coord,
	double** distance){

	for (size_t i = 0; i < Natoms; i++){
		for (size_t j = 0; j < Natoms; j++){
			double dx = coord[i][0] - coord[j][0];
			double dy = coord[i][1] - coord[j][1];
			double dz = coord[i][2] - coord[j][2];
			distance[i][j] = sqrt(dx*dx + dy*dy + dz*dz);
		}
	}
}

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

double V(double epsilon, // Lennard-Jones potential energy function
	double sigma,
	size_t Natoms,
	double** distance){
	// we have splited the equation y 4 terms
	double sigr; // calculates the value of sigma / r of every atom
	double t1; // t1 = (sigma/r)**12
	double t2; // t2 = (sigma/r)**6
	double etemp; // saves the operation 4*epsilon*(t1-t2) 
	double etot = 0.0; // sums up all the potential energies of each iteration

	for (size_t i = 0; i < Natoms; i++){
		for (size_t j = 0; j < Natoms ; j++){
			if (distance[i][j] < 0.00001){;
				continue;
		}
			else{
				sigr = sigma / distance[i][j];
				t1 = pow(sigr, 12);
				t2 = pow(sigr, 6);
				etemp=4*epsilon*(t1-t2);
			//	printf("Etemp = %lf\n", etemp); Used to track if the loop is correct
				etot += etemp;
			//	printf("Etot = %lf\n", etot); Used to track if the loop is correct
			}
		}
	}
return etot;
}

double T(size_t Natoms,
	double* velocity,
	double* mass){

	double ektot= 0.0;

	for (size_t i = 0; i < Natoms; i++){
	double	veltemp = velocity[i];
	double	t3 = pow(veltemp, 2);
		double	t4 = mass[i] * t3;
		ektot += 0.5 * t4;
	//	printf("EKtot = %lf\n", ektot);
		
	}
	return ektot;
}

double** compute_acc(double epsilon,
	double sigma,
	size_t Natoms,
	double** coord,
	double* mass,
	double** distance){
	double** acceleration = malloc_2d(Natoms, 3);
	double U;
	double t5, t6, epsr, sigr;
	double ax, ay, az;


	for (size_t i = 0; i < Natoms ; i++){
		for (size_t j = 0; j < Natoms; j++){
			if ( distance[i][j] < 0.001){
				continue;}
			sigr = sigma / distance[j][i];
			epsr = epsilon / distance [j][i];
			t5 = pow(sigr,6);
			t6 = pow(sigr,12);
			U = 24 * epsr * (t5 - 2 * t6);
	//		printf ("i= %zu\n", i);
	//		printf ("j = %zu\n", j);
	//		printf ("U = %lf\n", U);
			ax = (-1.0/mass[i]) * U * ((coord[i][0] - coord[j][0]) / distance[i][j]) ;
	//		printf ("ax = %fl\n", ax);
			ay = (-1.0/mass[i]) * U * ((coord[i][1] - coord[j][1]) / distance[i][j]) ;
	//		printf ("ay = %fl\n", ay);
			az = (-1.0/mass[i]) * U * ((coord[i][2] - coord[j][2]) / distance[i][j]) ;
	//		printf ("az = %fl\n", az);
	  		acceleration[i][0] += ax;
	  		acceleration[i][1] += ay;
	  		acceleration[i][2] += az;
	//		printf ("axi = %lf\n", acceleration[i][0]);
	//		printf ("ayi = %lf\n", acceleration[i][1]);
	//		printf ("azi = %lf\n", acceleration[i][2]);
		        
		}
	}
	return acceleration;
}


 double** compute_pos(size_t Natoms,
		double** coord,
		double* velocity,
		double** acceleration,
		double tstep){
	double** coordnext = malloc_2d(Natoms, 3);
	for (size_t i = 0; i < Natoms ; i++){
		for (size_t j = 0; j < Natoms ; j++){			
			coordnext[i][j] = coord[i][j] + (velocity[j] * tstep) + acceleration[i][j] * (tstep / 2.0) ;
		}
	}
	return coordnext;
 }


 double* compute_vel(size_t Natoms,
                double* velocity,
                double** acceleration,
		double** coordnext,
                double tstep,
		double epsilon,
		double sigma,
		double* mass,
		double** distance){
	double* velnext = (double*)malloc(Natoms * sizeof(double));
 	double** accnext = malloc_2d(Natoms, 3);
	accnext = compute_acc(epsilon, sigma, Natoms, coordnext, mass, distance);
        for (size_t i = 0; i < Natoms ; i++){
                for (size_t j = 0; j < Natoms ; j++){
                        velnext[i] = velocity[j] + 0.5 * (acceleration[i][j] + accnext[i][j]) * tstep ;
                }
        }
          return velnext;
 }

void write_xyz(FILE* file, size_t it, size_t Natoms, double** coord, const char* atom_symbol, 
               double ektot, double etot) {
    // Calcula la energía total
    double e = ektot + etot;

    fprintf(file, "Iteration number: %zu\n", it);

    // Escribe el número de átomos
    fprintf(file, "Natoms = %zu\n", Natoms);

    // Escribe el comentario con las energías
    fprintf(file, "Kinetic Energy: %.3f, Potential Energy: %.3f, Total Energy: %.3f\n", 
            ektot, etot, e);

    // Escribe las coordenadas de cada átomo
    for (size_t i = 0; i < Natoms; i++) {
        fprintf(file, "%s %.5f %.5f %.5f\n", 
                atom_symbol, coord[i][0], coord[i][1], coord[i][2]);
    }
}



int main(){
	FILE* file = fopen("inp.txt", "r"); //Opens the file
	if (file == NULL){
		perror("Error opening the file");
		return 1;
	}

	size_t Natoms = read_Natoms(file); // Read the number of atoms with the previous function
	
	double* mass = (double*)malloc(Natoms * sizeof(double)); // Here we allocate memory for mass, coordinates and distance variables
	double** coord = malloc_2d(Natoms, 3);
	double** distance = malloc_2d(Natoms,3);

	read_molecule(file, Natoms, coord, mass); // Obtaining the variables
	
	compute_distances(Natoms, coord, distance);

	printf("Nat = %zu\n", Natoms); // Printing Natoms and distance to check if they are correct
	printf("Distances Matrix: \n");
	for (size_t i = 0; i < Natoms; i++) {
	       for (size_t j = 0; j < Natoms; j++) {
        	   printf("%lf ", distance[i][j]);
		}
	printf("\n");}


	double epsilon = 0.0661; //j/mol 
	double sigma = 0.3345; //nm
	double etot = V(epsilon, sigma, Natoms, distance);
	printf("V = %lf\n", etot);

	double* velocity = (double*)malloc(Natoms * sizeof(double)); // Here we allocate memory for mass, coordinates and distance variables
	double ektot = T(Natoms, mass, velocity);
	printf("T = %lf\n", ektot);

	double** acceleration = malloc_2d(Natoms, 3);
	acceleration = compute_acc(epsilon, sigma, Natoms, coord, mass, distance);
	printf("Acceleration Matrix: \n");
	for (size_t i = 0; i < Natoms; i++) {
               for (size_t j = 0; j < Natoms; j++) {
                   printf("%lf ", acceleration[i][j]);
                }
        printf("\n");
	}
	
	double tstep = 0.2; //fs
	
 	double** coordnext = malloc_2d(Natoms, 3);
 	coordnext = compute_pos(Natoms, coord, velocity, acceleration, tstep);
        printf("Postion n+1 Matrix: \n");
        for (size_t i = 0; i < Natoms; i++) {
               for (size_t j = 0; j < Natoms; j++) {
                   printf("%lf ", coordnext[i][j]);
                }
        printf("\n");}

	double** accnext = malloc_2d(Natoms, 3);
	double* velnext = (double*)malloc(Natoms * sizeof(double));
	velnext = compute_vel(Natoms, velocity, acceleration, coordnext, tstep, epsilon, sigma, mass, distance);
  	printf("Vel n+1 Array: \n");
        for (size_t i = 0; i < Natoms; i++) {
	     	printf("%lf", velnext[i]);
        printf("\n");}


	// VERLET ALGORITHM
	int steps = 1000;
	int i, j, k ;
	FILE* output = fopen("trajectory.xyz", "w");
	    if (output == NULL) {
	        perror("Error opening file");
	    return 1;
	    }
	for (i = 0 ; i < steps ; i++){
		double etot = V(epsilon, sigma, Natoms, distance);
//	        printf ("V = %lf\n", etot);
		double ektot = T(Natoms, mass, velocity);
//		printf ("T = %lf\n", ektot);
//	        printf ("E = %lf\n", ektot + etot);
		coordnext = compute_pos(Natoms, coord, velocity, acceleration, tstep);
		acceleration = compute_acc(epsilon, sigma, Natoms, coord, mass, distance);
		accnext = compute_acc(epsilon, sigma, Natoms, coordnext, mass, distance);
		velnext = compute_vel(Natoms, velocity, acceleration, coordnext, tstep, epsilon, sigma, mass, distance);

		if ((i + 1) % 50 == 0){ 
//	        printf ("E = %lf\n", ektot + etot);
//		printf ("Iteration number = %zu\n", i+1);
//		printf ("Position Matrix: \n");
//		for (j = 0; j < Natoms ; j++){
//                        for (k = 0 ; k < Natoms ; k++){
//				 printf("%lf ", coordnext[j][k]);
//	                }
//		       	printf("\n");
//		}
		write_xyz(output, (i+1) ,  Natoms, coordnext, "Ar", ektot, etot);
		}
		coord = coordnext;
		velocity = velnext;

	}














free_2d(coord); // Freeing the memory where were the variables
free_2d(distance);
free(mass);
free(velocity);

}
