#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

	for (size_t i = 0; i < Natoms; i++){ // Distance between the same atom (diagonal in the matrix) will be 0!!
		for (size_t j = 0; j < Natoms; j++){ // Functions who need divide by the interatomic distance will have an if condition!!
			double dx = coord[i][0] - coord[j][0];
			double dy = coord[i][1] - coord[j][1];
			double dz = coord[i][2] - coord[j][2];
			distance[i][j] = sqrt(dx * dx + dy * dy + dz * dz);
		}
	}
}

double V(double epsilon, // Lennard-Jones potential energy function
	double sigma,
	size_t Natoms,
	double** distance){
	// we have splited the equation in 4 terms
	double sigr; // calculates the value of sigma / r of every atom
	double t1; // t1 = (sigma/r)**12
	double t2; // t2 = (sigma/r)**6
	double etemp; // saves the operation 4*epsilon*(t1-t2) 
	double potential_energy = 0.0; // sums up all the potential energies of each iteration

	for (size_t i = 0; i < Natoms; i++){
		for (size_t j = i + 1; j < Natoms ; j++){ //In the second loop, j must be bigger than i, thats why it has a +1
			if (distance[i][j] > 0){; //Avoiding divide by 0
				sigr = sigma / distance[i][j];
				t1 = pow(sigr, 12);
				t2 = pow(sigr, 6);
				etemp= 4.0 * epsilon * (t1 - t2);
				potential_energy += etemp; //Sums up the potential energy between i and j atoms
			}
		}
	}
return potential_energy;
}

double T(size_t Natoms, // Kinetic energy function
	double** velocity,
	double* mass){
	// we follow the same structure as Lennard-Jones potential, dividing the equation in 2 terms
	double kinetic_energy= 0.0; //total kinetic energy
	for (size_t i = 0; i < Natoms; i++){
		double t3 = 0.0;
		for (size_t j = 0; j < 3 ; j++){
			t3 += pow(velocity[i][j], 2);
		}
		// t3 is the modulus of velocity vector of each atom. Due the modulus is the sqrt of the sum of the coordinates squared
		// and the formula is velocity squared; sqrt and square cancells each other and it is only the sum of coordinates squared
		double t4 = mass[i] * t3; 
		kinetic_energy += t4;
	}

	kinetic_energy *= 0.5;	//Multiply once the 1/2 factor
	return kinetic_energy;
}

double E(double V, double T){ //Calculation of the total energy, that it is the sum of potential and kinetic energy
	double E_tot = V + T;
	return E_tot;
	}

void compute_acc(double epsilon, //We introduce the analytical expression of the acceleration
	double sigma,
	size_t Natoms,
	double** coord,
	double* mass,
	double** distance, 
	double** acceleration){

	for (size_t i = 0; i < Natoms ; i++){
		double ax = 0.0; // We have explicit the equation for each axis.
		double ay = 0.0;
		double az = 0.0; 
		for (size_t j = 0; j < Natoms; j++){
			if (i != j){
				if ( distance[i][j] > 0){ //This condition avoids divisions by 0 (self-interaction)
			double	sigr = sigma / distance[i][j]; //Also, we follow the same structure as potential energy;
			double	epsr = epsilon / distance [i][j]; //Splitting the expression into some terms
			double	t5 = pow(sigr, 6); //t5 = (sigma / r) ** 6
			double	t6 = pow(sigr, 12); //t6 = (sigma / r)** 12
			double	U = 24.0 * epsr * (t5 - 2.0 * t6); // Expression for the U
				ax += (-1.0 / mass[i]) * U * ((coord[i][0] - coord[j][0]) / distance[i][j]) ; //Acceleration terms for each axis
				ay += (-1.0 / mass[i]) * U * ((coord[i][1] - coord[j][1]) / distance[i][j]) ;
				az += (-1.0 / mass[i]) * U * ((coord[i][2] - coord[j][2]) / distance[i][j]) ;
				}
			}
		}
	  		acceleration[i][0] = ax; //Save the acceleration for each axis
	  		acceleration[i][1] = ay;
	  		acceleration[i][2] = az;
	}
}

void write_xyz(FILE* file, size_t Natoms, double** coord, const char* atom_symbol, 
               double kinetic_energy, double potential_energy, double E_tot) {
// We have created a function which prints an output with the correct format to load in programs like Molden    
   
       	fprintf(file,"%zu\n", Natoms);
	// Molden reads the number of atoms (Only the integer)
	// Molden does not read the second line, ergo we have decided to print the kinetic, potential and total energies.
	// We also have chosen that we will use 6 decimal numbers to describe the energy and the coordinates
    fprintf(file, "Kinetic Energy: %.6f J/mol, Potential Energy: %.6f J/mol, Total Energy: %.6f J/mol \n",
            kinetic_energy, potential_energy, E_tot);

    for (size_t i = 0; i < Natoms; i++) { //a loop to print the atomic symbol, and the 3 cartesian coordinates
        fprintf(file, "%s %.6f %.6f %.6f\n", 
                atom_symbol, coord[i][0], coord[i][1], coord[i][2]);
    }
}


//--------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------

int main(){
	FILE* file = fopen("inp.txt", "r"); //Opens the file
        if (file == NULL){
                perror("Error opening the file");
                return 1;
        }

        size_t Natoms = read_Natoms(file); // Read the number of atoms with the previous function

        double** coord = malloc_2d(Natoms, 3);
        double* mass = (double*)malloc(Natoms * sizeof(double)); // Here we allocate memory for coordinates and mass
        
	read_molecule(file, Natoms, coord, mass); // Obtaining the variables
	fclose(file); // Closing the input file

	double** distance = malloc_2d(Natoms, Natoms); // Allocate memory for distance matrix
	compute_distances(Natoms, coord, distance); // Calculate the distances of the 1st iteration 

	double epsilon = 0.0661; // Values of epsilon in j/mol 
	double sigma = 0.3345; // Values of sigma in nm
	double potential_energy = V(epsilon, sigma, Natoms, distance); //First calculation of potential energy

	double** velocity = malloc_2d(Natoms, 3); // Creation of the velocity matrix initalise in zero
  	for (size_t i = 0; i < Natoms; i++) {
	    for (size_t j = 0; j < 3; j++){
		    velocity[i][j] = 0.0;
	    }
	}

        double kinetic_energy = T(Natoms, velocity, mass); //First calculation of the kinetic energy
	double E_tot = E(potential_energy, kinetic_energy); // First calculation of the total energy

	double** acceleration = malloc_2d(Natoms, 3); // Allocate memory for acceleration matrix
	compute_acc(epsilon, sigma, Natoms, coord, mass, distance, acceleration); // Calculation of the acceleration

	FILE* output = fopen("outputDYH_DLP.xyz", "w"); //Here we write the name of the output and we open a new file with this name for printing the output

	// VERLET ALGORITHM

	double tstep = 0.2; // SELECTION OF THE TIME STEP
	size_t steps = 1000; // SELECTION OF THE NUMBER OF STEPS
	
	for (size_t n = 0 ; n < steps ; n++){ 
		for (size_t i = 0 ; i < Natoms ; i++){
			for (size_t j = 0 ; j < 3 ; j++){ //First we calculate r^(n+1), then the part of v^(n+1) which depends on the actual acceleration
				coord[i][j] += velocity[i][j] * tstep + acceleration[i][j] * pow(tstep, 2) / 2.0;
				velocity[i][j] += 0.5 * acceleration[i][j] * tstep;
			}
		}
		compute_distances(Natoms, coord, distance); // Calculate the acceleration in the next step with the new coordingates
		compute_acc(epsilon, sigma, Natoms, coord, mass, distance, acceleration);
		for (size_t i = 0 ; i < Natoms ; i++){
			for (size_t j = 0 ; j < 3 ; j++){ // Adding the part of a^(n+1) into the velocity of the next step
				velocity[i][j] += 0.5 * acceleration[i][j] * tstep;
			}
		}

		double potential_energy = V(epsilon, sigma, Natoms, distance); //Calculation of potential energy	
		double kinetic_energy = T(Natoms, velocity, mass); //Calculation of the kinetic energy
		double E_tot = E(potential_energy, kinetic_energy); // Calculation of the total energy

		if ((n + 1) % 10 == 0){ //HERE WE SELECT AT M STEPS WE PRINT THE INFORMATION. M = 10
		write_xyz(output,  Natoms, coord, "Ar", kinetic_energy, potential_energy, E_tot);
		}
	}

	fclose(output); // When all the steps were done, we close the output file


	free_2d(coord); // Freeing the memory where the variables were
	free_2d(distance);
	free_2d(velocity);
	free_2d(acceleration);
	free(mass);

	return 0;
}
