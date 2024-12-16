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

	for (size_t i = 0; i < Natoms; i++){ // Distance between the same atonm (diagonal in the matrix) will be 0!!
		for (size_t j = 0; j < 3; j++){ // Functions who need divide by the interatomic distance will have an if condition!!
			double dx = coord[i][0] - coord[j][0];
			double dy = coord[i][1] - coord[j][1];
			double dz = coord[i][2] - coord[j][2];
			distance[i][j] = sqrt(dx*dx + dy*dy + dz*dz);
		}
	}
}

double** malloc_2d(size_t m, size_t n) { //Dynamical allocation function provide by the teachers
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

void free_2d(double** a) { // Freeing dynamical allocation function provide by the teachers
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
		for (size_t j = 0; j < 3 ; j++){ //In all concatenated for loops, the i goes over all the atoms in the input and j over the 3 cartesian axis
			if (distance[i][j] < 0.00001){; //Avoiding divide by 0
				continue;
		}
				sigr = sigma / distance[i][j];
				t1 = pow(sigr, 12);
				t2 = pow(sigr, 6);
				etemp= 4 * epsilon * (t1-t2);
				etot += etemp; //Sums up the potential energy between i and j atoms
		}
	}
return etot;
}

double T(size_t Natoms, // Kinetic energy function
	double** velocity,
	double* mass){
	// we follow the same structure as Lennard-Jones potential
	double ektot= 0.0; //total kinetic energy
	double t4 = 0.0; // t4 is the summatory of t3 multiplied by the mass of each atom
	for (size_t i = 0; i < Natoms; i++){
		double t3 = pow(velocity[i][0], 2) + pow(velocity[i][1],2) + pow(velocity[i][2], 2);
		// t3 is the modulus of velocity vector of each atom. Due the modulus is the sqrt of the sum of the coordinates squared
		// and the formula is velocity squared; sqrt and square cancells each other and it is only the sum of coordinates squared
		t4 += mass[i] * t3; 
		ektot = 0.5 * t4;
		
	}
	return ektot;
}

// Now we will create the functions for the calculation of the acceleration, next position and next velocity
// We have used double** instead of void due to void gave us some problems when we wanted to calculate
// acceleration and the acceleration at the next step. At the end, we rather use double even though the code
// is not perfect from a programming point of view
//
// Also, when we talk about next positon or next velocity is the position or velocity at the next step!


double** compute_acc(double epsilon, //we have introduced the analytical expression of the acceleration
	double sigma,
	size_t Natoms,
	double** coord,
	double* mass,
	double** distance){
	double** acceleration = malloc_2d(Natoms, 3); //dynamical allocation due to it will be need along the verleth algorithm
	double U;
	double t5, t6, epsr, sigr;
	double ax, ay, az; // we have explicit the equation for each axis.

	for (size_t i = 0; i < Natoms ; i++){
		for (size_t j = 0; j < 3; j++){
			if ( distance[i][j] < 0.001){ //We introduce this if to avoid divisions by 0 (self-interaction)
				continue;}
			sigr = sigma / distance[i][j]; 
			epsr = epsilon / distance [i][j];
			t5 = pow(sigr,6); //t5 = (sigma / r) ** 6
			t6 = pow(sigr,12); //t6 = (sigma / r)** 12
			U = 24 * epsr * (t5 - 2 * t6); // Expression for the U
			ax = (-1.0 / mass[i]) * U * ((coord[i][0] - coord[j][0]) / distance[i][j]) ; //acceleration terms for each axis
			ay = (-1.0 / mass[i]) * U * ((coord[i][1] - coord[j][1]) / distance[i][j]) ;
			az = (-1.0 / mass[i]) * U * ((coord[i][2] - coord[j][2]) / distance[i][j]) ;
	  		acceleration[i][0] += ax; //save the acceleration for each axis
	  		acceleration[i][1] += ay;
	  		acceleration[i][2] += az;
		}
	}
	return acceleration;
}


 double** compute_pos(size_t Natoms,
		double** coord,
		double** velocity,
		double** acceleration,
		double tstep){
	double** coordnext = malloc_2d(Natoms, 3); //dynamical allocation due to it will be used in the calculation of the next velocity
	for (size_t i = 0; i < Natoms ; i++){
		for (size_t j = 0; j < 3 ; j++){			
			coordnext[i][j] = coord[i][j] + (velocity[i][j] * tstep) + acceleration[i][j] * (pow(tstep, 2) / 2.0); //function for position at the next step
		}
	}
	return coordnext;
 }


 double** compute_vel(size_t Natoms,
                double** velocity,
                double** acceleration,
		double** coordnext,
                double tstep,
		double epsilon,
		double sigma,
		double* mass,
		double** distance){
		double** velnext = malloc_2d(Natoms, 3); //Allocating the memory for next velocities matrix
	 	double** accnext = malloc_2d(Natoms, 3); //Allocating the memory for next acceleration matrix
		// We have include the calculation of the next acceleration here because it is the step before
		// of the calculation of the velocity at the next step

		accnext = compute_acc(epsilon, sigma, Natoms, coordnext, mass, distance); 
		
        for (size_t i = 0; i < Natoms ; i++){
                for (size_t j = 0; j < 3 ; j++){
                        velnext[i][j] = velocity[i][j] + 0.5 * (acceleration[i][j] + accnext[i][j]) * tstep ;
                }
        }
          return velnext;
 }

void write_xyz(FILE* file, size_t Natoms, double** coord, const char* atom_symbol, 
               double ektot, double etot) {
// we have created a function which prints an output with the correct format to load in programs like Molden    
   
       	fprintf(file,"%zu\n", Natoms);
	// Molden reads the number of atoms (Only the integer)
	double e = ektot + etot;
	// Molden does not read the second line, ergo we have decided to print the kinetic, potential and total energies.
	// We also have chosen that we will use 3 decimal numbers to describe the energy and the coordinates
    fprintf(file, "Kinetic Energy: %.3f, Potential Energy: %.3f, Total Energy: %.3f\n",
            ektot, etot, e);

    for (size_t i = 0; i < Natoms; i++) { //a loop to print the atomic symbol, and the 3 cartesian coordinates
        fprintf(file, "%s %.3f %.3f %.3f\n", 
                atom_symbol, coord[i][0], coord[i][1], coord[i][2]);
    }
}



int main(){
	// VERLET ALGORITHMi
	FILE* file = fopen("inp.txt", "r"); //Opens the file
        if (file == NULL){
                perror("Error opening the file");
                return 1;
        }

        size_t Natoms = read_Natoms(file); // Read the number of atoms with the previous function

        double* mass = (double*)malloc(Natoms * sizeof(double)); // Here we allocate memory for mass, coordinates and distance variables
        double** coord = malloc_2d(Natoms, 3);
        double** distance = malloc_2d(Natoms, 3);

        read_molecule(file, Natoms, coord, mass); // Obtaining the variables

	compute_distances(Natoms, coord, distance); // Calculate the distances of the 1st iteration 

	double epsilon = 0.0661; // Values of epsilon in j/mol 
	double sigma = 0.3345; // Values of sigma in nm
	double etot = V(epsilon, sigma, Natoms, distance); //First calculation of potential energy
	double** velocity = malloc_2d(Natoms, 2); // Creation of the velocity matrix initalise in zero
	double ektot = T(Natoms, velocity, mass); //First calculation of the kinetic energy
	double** acceleration = malloc_2d(Natoms, 3);

	double tstep = 0.2; // SELECTION OF THE TIME STEP
	int steps = 1000; // SELECTION OF THE TOTAL STEPS

 	double** coordnext = compute_pos(Natoms, coord, velocity, acceleration, tstep); // Calculation of positions in the first iteration
	double** accnext = malloc_2d(Natoms, 3);
	double** velnext = compute_vel(Natoms, velocity, acceleration, coordnext, tstep, epsilon, sigma, mass, distance); // Calculation of velocity in the first iteration
	int i, j, k ;

	FILE* output = fopen("outputDYH_DLP.xyz", "w"); //Here we write the name of the output and we open a new file with this name for printing the output
	    if (output == NULL) {
	        perror("Error opening file");
	    return 1;
	    }
	for (i = 0 ; i < steps ; i++){
		compute_distances(Natoms, coord, distance); //First step of the algorithm is calculate the distance due to the coordinates will change in every step,
	       // thus the distance will do the same
		double etot = V(epsilon, sigma, Natoms, distance); //Calculation of the potential energy
		double ektot = T(Natoms, velocity, mass); // Calculation of the potential energy with the actual velocities
		acceleration = compute_acc(epsilon, sigma, Natoms, coord, mass, distance); // Calculation of the actual acceleration
		coordnext = compute_pos(Natoms, coord, velocity, acceleration, tstep); // Calculation of the position at the next step
//		accnext = compute_acc(epsilon, sigma, Natoms, coordnext, mass, distance); // Calculatio
		velnext = compute_vel(Natoms, velocity, acceleration, coordnext, tstep, epsilon, sigma, mass, distance); //Calculation of the velocity at the next step
		coord = coordnext; // Saving the position at the next step to the actual step
		velocity = velnext; // Saving the velocity at the next step to the actual step

		if ((i + 1) % 50 == 0){ //HERE WE SELECT AT M STEPS WE PRINT THE INFORMATION. M = 50
		write_xyz(output,  Natoms, coordnext, "Ar", ektot, etot);
		}
	}

free_2d(coord); // Freeing the memory where were the variables
free_2d(distance);
free_2d(velocity);
free_2d(acceleration);
free(mass);

}
