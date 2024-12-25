
README FILE

"Project1_DYH_DLP.c" should be compiled linking trexio libraries using a C compiler, an example is shown next with gcc:

>>   gcc -I/usr/local/include -L/usr/local/lib -ltrexio Project1_DYH_DLP.c -o HF_MP2\

HF_MP2 executable is used in presence of the molecule file name as argument. Molecule file should be a ".h5" trexio file and stored in a directory named "data" within the same directory as the executable. An example of the execution of the program is shown next:

>> ./HF_MP2 h2o.h5

If the path where molecule files are stored needs to be modified, it could be done by changing variable "molecules_directory" in line 54 in "Project1_DYH_DLP.c". Defined as default as "./data/".

Once the program is executed results will be shown in the terminal and also will generate an output file in the form of:

>> h2o.h5.out

This output file contains HF energy and MP2 correction, also selected data involved in Hartree Fock and MP2 calculations; repulsion energy, total and occupied molecular orbitals, one and two electron integral terms, and orbital energies.

Program also gives the possibility of displaying in terminal specific two electron integrals. The lines to do it remain commented as default between 275-285 lines. Indexes i, j, k, l have to be initialized (i=1, j=2 , k=3, l=4 as default).

