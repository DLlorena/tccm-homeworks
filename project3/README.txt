In project3 folder, you will found Project3_DYH_DLP.c which is the code, Project3.exe which is the executable and outputDYH_DLP.xyz which is the results of the code.

Inside the code you will find comments along the file explaining who we did the projet.

To compile the project, you will have to use the next command:
>> gcc -lm Project3_DYH_DLP.c -o Project3.exe (or the name you want for the file)
                                               
To execute the executable you will use this command:
>> ./Project3.exe

It will generate an output fille named: "outputDYH_DLP.xyz". 
The output shows the number of atoms, kinetic, potential and total energy and the 3 coordinates of each atom in every step.
                                               
To visualize the results we recommend use Molden. You can load the output file by this way:
>> molden outputDYH_DLP.xyz                                              
                                               
When you are in Molden, we recommend to use the second option of the label menu (atoms + number) and disable the shade option.
In addition, to have a better view of the dynamics, we also recommend to rotate a bit the system to see the three argon atoms.
In order to see a smooth dynamics, you should decrease M value, we believe it is affordable to use M = 1 and print all the steps because the dynamics are smoother and can be seen in better detail.

                                              
