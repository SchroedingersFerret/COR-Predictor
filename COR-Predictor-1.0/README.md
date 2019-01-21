COR Predictor version 1.0
=========================
By: J. Ball (SchroedingersFerret)

README

COR Predictor is a command line tool for developing a function  that approximates the coefficient of restitution (COR) of two colliding bodies from their physical properties, _without_ knowing the final velocities. 

COR Predictor does this using a simple supervised learning process. The program reads training data from the files `cor_x.csv` and `cor_y.csv` and numerically optimizes a set of 16 parameters to fit the function to the training data. Once the optimization has reached the required accuracy, the parameters are stored in the file `cor_parameters.csv`, which can be called at startup during future optimizations to improve execution times.

---

Files

In addition to the executable, COR Predictor makes use of the following files:
 
`cor_x.csv`:

The file `cor_x.csv` is an n x 7 table containing the independent variables for n training datapoints. The variables represent two sets of four different physical quantities, one for each of two colliding objects. 


* Yield Strength:

The first quantity, stored in `x[i][0]` and `x[i][1]`, is the yield strength in MPa of objects 1 and 2 respectively. _Yield strength_ here is defined as the limit of elastic behavior, beyond which the stress causes the object to deform plastically. 

* Young's Modulus

The second quantities are stored in `x[i][2]` and `x[i][3]` and are the Young's modulus in MPa of objects 1 and 2 respectively. 

* Density:

`x[i][4]` and `x[i][5]` are the densities in Mg/m^3 of objects 1 and 2 respectively. 

* Velocity: 

`x[i][6]` is the magnitude of the objects' relative velocity.


`cor_y.csv`:

The file `cor_y.csv` is an n x 1 table containing the dependent variables for n training datapoints. Each variable is the COR for the collision represented by each datapoint.


`settings.txt`:

File `settings.txt` contains the settings used in the genetic searching algorithm used by the learning program. The various entries and their functions are covered in the section titled "Using COR Predictor".

`cor_parameters.csv`:

This file is generated at the end of optimization and contains the best-fitting parameter configuration. It can be used at startup to reduce the time needed for further optimization.

---

Building COR Predictor

COR Predictor is built by using CMake to generate a makefile which can be used to build the executables. The executables need to be moved to the parent folder with the other files before they can be run.

* Gnu+Linux

In the terminal, navigate to the parent directory `~/COR-Predictor-1.0`. Create the build folder with the command

`~/COR-Predictor-1.0 $ mkdir build`

Next change to the build directory with the command

`~/COR-Predictor-1.0 $ cd build`

Run cmake with the command 

`~/COR-Predictor-1.0/build $ cmake ..`

and make the executable with 

`~/COR-Predictor-1.0/build $ make`

Two executable will be created in the build directory. Move them to the parent folder using a file manager or with

`~/COR-Predictor-1.0/build $ mv COR-Predictor make-CMakeLists ~/COR-Predictor-1.0`

Run COR Predictor by navigating back to the parent directory and entering

`~/COR-Predictor-1.0 $ ./COR-Predictor`

* Windows

In the command line, enter the command

`Run cmake-gui.exe`

or run the program from the start menu, Program Files, or a desktop shortcut.

In the top entry of the GUI, enter the location of the parent directory.
No configuration is necessary. Simply press 'Generate' to create the makefile.

---

Using COR Predictor

COR Predictor uses a genetic algorithm to fit a set of parameters to the training data. The genetic algorithm starts by creating an initial set of parameter configurations and selecting a percentage of them that best fit the training data as a gene pool. These configurations are encoded as 2-d bit arrays. 

An iterative process then starts where the best-fitting half of the configurations is selected to be kept for reproduction. Pairs of these configurations are selected as parents, and children configurations are created by carrying over each bit from one of the two parents selected at random. The children replace the configurations not selected for reproduction.

A percentage of all the bits in the entire population is then selected to be mutated by being changed from 1 to 0 or vise versa. This step is important, as it introduces variation to prevent inevitable bottlenecking. Bottlenecking happens because of the finite size of the population, which will eventually reach a point where all the configurations are very similar and stall the process unless variation is continually introduced. In terms of the optimization problem, mutation allows poorer-fitting configurations to be accepted into the population. Thus, mutation is necessary because the objective function space contains numerous local minima that impede searching.

After mutation occurs, the iterative process starts over and repeats until an error tolerance specified by the user is reached. The process of selecting the best-fitting configurations combined with scrambling the bits using reproduction and mutation efficiently searches the objective function space for an optimal configuration.

Settings

When COR Predictor is run, the program first reads the file `settings.txt` to obtain the settings used to initiate optimization. These settings strongly determine the successfulness of optimization. They must be listed in the file with the syntax

`[Setting_name]=[value]`

If a setting is to be changed, only make changes to the value and leave no space between the value and the `=` sign. The setting name should not be altered, as changing it will have no effect. The setting corresponding to each value is only determined by the order of the listing. 

* `initial_population_size`:

This setting determines the size of the initial population created at initialization. It should be sufficiently large enough to create a wide selection of configurations. If this value is too small, the error of the configurations chosen at startup may be inflated.

* `gene_pool_size`:

This is the number of configurations in the initial population selected to be kept for the iterations. It should be much smaller than the `initial_population_size` setting so that mostly well-fitting configurations will be selected at startup. This setting must also be large enough to slow the speed at which bottlenecking occurs. `gene_pool_size` is limited by the computer's available processing power, as raising its value scales the number of operations performed on each iteration exponentially. 

* `elite_population_size`:

This setting increases the number of configurations kept in the "elite population", which serves several functions. The elite population is kept safe from mutations that decrease their fitness. This provides protection against divergence, especially at high mutation rates. However, large elite population sizes decrease the speed at which the objective function space is searched. Around 25 percent of the gene pool is an appropriate size.

The elite population also plays a role in reducing the time needed for repeated convergences. At startup, the option is given to use a parameter configuration from a previous convergence, which will be used as each configuration in the elite population. This "jump-starts" the convergence process.

When troubleshooting, it is often useful to set the elite population size to 0, since the program only displays the least squares configuration. Removing the effect of the elite population reveals the behavior of the configurations that are not protected from mutations.

* `mutation_rate`:

The `mutation_rate` setting determines the percentage of bits accross the entire population that are changed during mutation. Increasing its value increases the amount of variation introduced on each iteration and aids in the algorithm's ability to search the objective function space. However, a too-large mutation rate will cause the program to diverge or to scan around wildly without converging. 

* `least_squares_error`:

This is the least squares error value at which iteration will conclude. It is better to be conservative (not too small) with this value, as there is no other way to stop the iterative process prematurely at this point. It is better to perform a convergence in multiple steps to avoid the algorithm stalling without reaching the required tolerance. If that happens, the only recourse as of this release is to close out the terminal and restart the program.

Adding Datapoints

COR Predictor reads training datapoints from the files `cor_x.csv` and `cor_y.csv`. To add a datapoint, simply open the files in a spread sheet and enter the data in the next row below the lowest entry in each file. 

Running COR Predictor

COR Predictor is run by entering the command 

`~/COR-Predictor-1.0 $ ./COR-Predictor`

while in the parent directory. 

Upon startup, the program checks if an existing parameter file exists and asks the user whether to initiate using the file if one is found or to start with random configurations. If the user responds "y", COR predictor will continue optimization from where it left off. If the user responds "n", the program will start from scratch. 

Once the program enters the main loop, the program will display the least squares error for each iteration so that the user can track its progress. The main loop will continue until the specified tolerance is reached or the user stops the program. If the required tolerance is reached the program will ask the user whether to write the new parameter configuration to file. If the user responds "y", the new configuration will be saved, overwriting the previous one. If the user responds "n", the new configuration will be discarded. 

Overall, the entire process should take no longer than several minutes for set of good-quality data with little noise. If optimization takes much longer than that, the settings might need to be changed in order for the algorithm to work faster.

Troubleshooting

* `settings.txt`/`cor_x.csv`/`cor_y.csv`/`cor_parameters.csv` was not found.

Each of these files must be in the same folder as the executable.

* `cor_x.csv`/`cor_y.csv` is/are the wrong dimension(s).

Each row of `cor_x.csv` must have 8 columns. Each row of `cor_y.csv` must have 1 column.

* `cor_x.csv` and `cor_y.csv` do not have the same number of entries.

Each set of independent variables in `cor_x.csv` must have a corresponding dependent variable in `cor_y.csv` and vice versa. Both files must have the same number of rows.

* The least squares error is very high at startup. 

This is to some extent normal and usually not a problem since the algorithm works fastest at the beginning of iterations. However, the fitness of the starting population can be improved by increasing the size of the initial population.

* The iterations proceed very slowly.

The gene pool size increases the number of operations performed per iteration. For slower computers, the gene pool size should not be much larger than 500.

* The error fluctuates wildly or diverges. 

If the mutation rate is too high, the error may fluctuate wildly rather than decreasing, or it may increase. A reasonable value for the mutation rate is 0.005, though a higher number may be more appropriate for larger datasets with more noise.

* The error remains the same from the beginning of iterations or decreases slowly and/or inconsistantly.

Set the elite population size to 0 and retry. If the error fluctuates wildly or diverges, decrease the mutation rate. If that is not the case, then the global minimum may already be reached. If that is the case but the error is too high, the quality of the training data may be low. If you are confident in your training data and cannot solve the problem, contact me. The objective function may need to be redesigned.

* The error decreases steadily then remains the same.

Set the elite population size to 0 and retry. If the behavior remains the same, try increasing the mutation rate. If the error decreases then fluctuates, try decreasing the mutation rate or increasing the size of the elite population.

* The error decreases steadily then decreases progressivly slower.

Try decreasing the size of the elite population or increasing the mutation rate.
If this does not help, change the least squares error tolerance to a less-ambitious estimate and use a series of optimizations to reach the required tolerance.

Generating a Function 

Once the parameters have been fitted using the learning program, a library containing a function can be configured using CMake. 

First, run the executable `make-CMakeLists` by entering this command in the parent directory:

`~/COR-Predictor-1.0 $ ./make-CMakeLists`

This will create a CMakeList file in the folder `function`. Navigate to this folder using 

`~/COR-Predictor-1.0 $ cd function`

Create a build folder using 

`~/COR-Predictor-1.0/function $ mkdir build`

Navigate to this folder:

`~/COR-Predictor-1.0/function $ cd build`

Run CMake using 

`~/COR-Predictor-1.0/function/build $ cmake ..`

and run the makefile:

`~/COR-Predictor-1.0/function/build $ make`

A library called `Restitution.hpp` will be configured using the new parameters and placed in the build folder.

---

Credits

The implementation of the genetic algorithm used in COR Predictor is built off of information available in Tao Pang's "An Introduction to Computational Physics," second edition.

---

Contact 

Contact me at <https://github.com/SchroedingersFerret>.

