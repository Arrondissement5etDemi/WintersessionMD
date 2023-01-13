# WintersessionMD
Introduction to Molecular Dynamics with LAMMPS for Princeton Wintersession 2023

## Introduction

Atoms and molecules are like balls and springs, interacting under rules of quantum mechanics. The interactions can be partially approximated by classical effective forces. [Molecular dynamics](https://en.wikipedia.org/wiki/Molecular_dynamics) (MD) solves the equations of motion for you! It is a computer simulation method for analyzing the physical movements of atoms and molecules. The atoms and molecules are allowed to interact for a fixed period of time, giving a view of the dynamic "evolution" of the system. 

Large-scale Atomic/Molecular Massively Parallel Simulator ([LAMMPS](https://www.lammps.org/#gsc.tab=0)) is an efficient and popular molecular dynamics program from Sandia National Laboratories. It's written for serial and parallel computing alike. For this workshop, we will mainly learn its serial version, but will touch upon the parallel version if we have time.

 ![LAMMPS logo](https://docs.lammps.org/_static/lammps-logo.png)

## Installation

### Getting into the cluster

To ensure that everyone has the same environment, let's use Princeton Research Computing's [Adroit](https://researchcomputing.princeton.edu/systems/adroit) cluster, which is a high-performance computer to do our calculations. SSH into your cluster by typing the following command into your terminal. DO NOT type the dollar sign; it just means that we're working with a terminal.

    $ ssh <YourNetID>@adroit-vis.princeton.edu
    
and enter your password. (The password column remains blank when you type, which is expected.)  

> ⚠️ For Windows users, you can use Command Prompt for SSH. But we also recommend Linux-terminal simulator softwares like [CygWin](https://www.cygwin.com/), as well as SSH manager softwares like [WinSCP](https://winscp.net/eng/index.php).

Now the terminal should look like this:

    [<YourNetID>@adroit-vis ~]$
    
Let's create a working directory called `MD` and go into it, by entering these two commands

    $ mkdir MD
    $ cd MD

### Installing Conda and LAMMPS

There are multiple ways to install LAMMPS, see [here](https://docs.lammps.org/Install.html). Unfortunately, Princeton's clusters won't let you install softwares in random places using commands like `sudo apt-get`. So let's use Conda to install LAMMPS. Conda is a package management system and environment management system that runs on a variety of operation systems, and saves us the troubles of figuring out environment dependencies, etc. It's much faster than building and installing using cmake, too.

See what versions of Conda the Adroit cluster has:

    $ module avail
    
and load the latest version:

    $ module load anaconda3/2022.5
    
Now your terminal should look like this:

    (base) [<YourNetID>@adroit-vis MD]$

We are now ready to install LAMMPS by typing the following commands one by one

    $ conda config --add channels conda-forge
    $ conda create -n my-lammps-env
    $ conda activate my-lammps-env
    $ conda install lammps
    
LAMMPS should be installed! The termial should now look like this:

    (my_lammps_env) [<YourNetID>@adroit-vis WintersessionMD]$

Try typing either `$ lmp` or `$ lmp_serial` in your command line to check if it's installed correctly. (Some distributions require you to type `$ lmp_mpi` instead of `$ lmp` .)     
Expected output for `$ lmp` :   

    LAMMPS (23 Jun 2022 - Update 2)
    OMP_NUM_THREADS environment is not set. Defaulting to 1 thread. (src/comm.cpp:98)
    using 1 OpenMP thread(s) per MPI task
    
Expected output for `$ lmp_serial` :  

    LAMMPS (23 Jun 2022 - Update 2)
    

The program is waiting for LAMMPS commands. `ctrl/command + c` to stop the execution for now, as we will provide the commands in the form of a file.

### Installing OVITO

We will use OVITO to visualize the molecules. No command-line is needed here. It can be downloaded and installed to your local machine via [this link](https://www.ovito.org/). You can also build it via GitHub source code following the instructions [here](https://www.ovito.org/manual/development/build_linux.html).
    
## Your first LAMMPS script! A Lennard-Jones fluid

> ⚠️ Actually, the syntax of LAMMPS is notoriously difficult, as it often mixes function names, arugment names and argument values, without an intuitive way to tell which is which. It can be intimidating when you first encounter a LAMMPS script. We feel that the best way to learn LAMMPS is to start with someone else's script and modify it according to your need. 

The Lennard-Jones (LJ) potential is a potential commonly used to describe the intermolecular forces between noble gas molecules, given by 

$$ \frac{v_{LJ}(r)}{\varepsilon} = 4\left(\left(\frac{r}{\sigma}\right)^{-12} - \left({r\over \sigma}\right)^{-6}\right), $$

where $r$ is the distance between two molecules, $\varepsilon$ is an energy parameter, and $\sigma$ is a distance parameter. The following figure shows a plot of the LJ potential. Note the steep repulsion (positive energy) at small $r$ (i.e., when the molecules are close), and the attraction at intermediate $r$. 

![The LJ potential.](https://upload.wikimedia.org/wikipedia/en/thumb/e/e7/Graph_of_Lennard-Jones_potential.png/1920px-Graph_of_Lennard-Jones_potential.png)
 
We want to simulate a fluid in which each two particles are interacting with $v_{LJ}(r)$. Let's write a LAMMPS script for that! Open a file called `lj.lmp` using the command `vi lj.lmp`, press `i` for insert mode, and copy and paste the following script. (To paste, just right click.) Anything after `#` is a comment explaining the purpose of the line. For more information of each command, see the LAMMPS documentation [here](https://docs.lammps.org).

    #This document is named "lj.lmp"
    
    units lj                      #this means we are using arbitrary length and energy units, as in v_LJ: no physical units are assigned to particles.
    atom_style atomic             #means we are dealing with unbonded atoms
    boundary p p p                #periodic boundary conditions
    
    lattice fcc 0.5               #defines a face-centered cubic lattice of number density 0.5
    region box block 0 5 0 5 0 5  #defines our simulation region named "box", which is a cube with side length = 5 fcc unit-cell side length
    create_atoms 1 box            #creates atoms of type "1" (there are only 1 type of particles) on the lattice in the region "box". 
      
    timestep 0.00025              #defines the time step (dt) of the simulation

    mass 1 1.0                    #sets the mass of atoms of type "1" that we just created
    velocity all create 1 2023    #gives all atom of a random velocity, such that the initial temperature is 1. 2023 is a random seed here.

    pair_style lj/cut 2.5         #defines the LJ potential. It is cut off (v(r) set to 0) after r > 2.5 \sigma
    pair_coeff 1 1 1.0 1.0        #sets the LJ coefficients \varepsilon and \sigma between atom types "1" and "1". We are 

    fix 1 all nvt temp 1.0 1.0 $(100.0*dt)  #sets the equation of state: we are fixing number of particles, volume and temperature during each time step. 
                                            #the 3 numbers after "temp" are thermostat parameters: intial and final temperatures, and temperature damping parameter

    thermo 1000                                                                     #prints thermodynamic information every 1000 time steps 
    thermo_modify format line "%d %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e"     #...and specifically, print thermodynamic information in this suggested format

    dump my_dump all custom 100 lj_cut.lammpstrj id type x y z vx vy vz             #dumps information of all atoms every 100 timesteps in the following custom style:
                                                                                    #id of atom, type of atom, x y z coordinates, and x y z components of its velocity

    run 100000                                                                      #runs 100000 timesteps

Save this document (press 'Esc' then type `:wq`), and now we can run the script.

    $ lmp -in lj.lmp

### Bonus skill! extra information that you can compute

    compute myRDF all rdf 200 1 1
    fix 2 all ave/time 50 6 1000 c_myRDF[*] file lj_cut1.rdf mode vector ave window 6

    compute myvacf all vacf
    fix storeMyVacf all vector 1 c_myvacf[4]

    compute msd_1 all msd
    fix store_msd_1 all vector 10 c_msd_1[4]
    variable fitslope_1 equal slope(f_store_msd_1)/6/(10*dt)
    fix 3 all ave/time 100 1 100 c_msd_1[4] v_fitslope_1 c_myvacf[4] file lj_cut1_diffusion.txt
