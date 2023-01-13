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
    
> If you'd like to copy and paste the commands instead, just copy here and right-click in the terminal.

LAMMPS should be installed! The termial should now look like this:

    (my_lammps_env) [<YourNetID>@adroit-vis MD]$

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

```julia
#This document is named "lj.lmp"

units lj                      #this means we are using arbitrary length and energy units, as in v_LJ: no physical units are assigned to particles.
atom_style atomic             #means we are dealing with unbonded atoms
boundary p p p                #periodic boundary conditions

lattice fcc 0.2                     #defines a face-centered cubic lattice of number density 0.2
region my_region block 0 7 0 7 0 7 #defines a geometric simulation region named "my_region", which is a cube with side length = 7 fcc unit-cell side length
create_box 1 my_region              #creates a simulation box with 1 atom types in the region "my_region" just created
create_atoms 1 box                  #creates atoms of type "1" (there are only 1 type of particles) on the lattice in the box just created
  
timestep 0.0025              #defines the time step (dt) of the simulation

mass 1 1.0                    #sets the mass of atoms of type "1" that we just created
velocity all create 0.9 2023    #gives all atom of a random velocity, such that the initial temperature is 0.9. 2023 is a random seed here.

pair_style lj/cut 2.5         #defines the LJ potential. It is cut off (v(r) set to 0) after r > 2.5 \sigma
pair_coeff 1 1 1.0 1.0        #sets the LJ coefficients \varepsilon and \sigma between atom types "1" and "1". We simply set both \varepsilon and \sigma to be 1.

fix 1 all nvt temp 0.9 0.9 $(100.0*dt)  #sets the equation of state: we are fixing number of particles, volume and temperature during each time step. 
                                        #the 3 numbers after "temp" are thermostat parameters: intial and final temperatures, and temperature damping parameter

thermo 1000                                                                     #prints thermodynamic information every 1000 time steps 
thermo_modify format line "%d %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e"     #...and specifically, print thermodynamic information in this suggested format

dump my_dump all custom 100 lj_cut.lammpstrj id type x y z vx vy vz             #dumps information of all atoms every 100 timesteps in the following custom style:
                                                                                #id of atom, type of atom, x y z coordinates, and x y z components of its velocity

run 50000                                                                      #runs 50,000 timesteps
```

Save this document (press 'Esc' then type `:wq`), and now we can run the script.

    $ lmp -in lj.lmp

The expected output should be like this:

```
LAMMPS (23 Jun 2022 - Update 2)
OMP_NUM_THREADS environment is not set. Defaulting to 1 thread. (src/comm.cpp:98)
  using 1 OpenMP thread(s) per MPI task
Lattice spacing in x,y,z = 2.7144176 2.7144176 2.7144176
Created orthogonal box = (0 0 0) to (19.000923 19.000923 19.000923)
  1 by 1 by 1 MPI processor grid
Created 1372 atoms
  using lattice units in orthogonal box = (0 0 0) to (19.000923 19.000923 19.000923)
  create_atoms CPU = 0.000 seconds
Generated 0 of 0 mixed pair_coeff terms from geometric mixing rule
Neighbor list info ...
  update every 1 steps, delay 10 steps, check yes
  max neighbors/atom: 2000, page size: 100000
  master list distance cutoff = 2.8
  ghost atom cutoff = 2.8
  binsize = 1.4, bins = 14 14 14
  1 neighbor lists, perpetual/occasional/extra = 1 0 0
  (1) pair lj/cut, perpetual
      attributes: half, newton on
      pair build: half/bin/atomonly/newton
      stencil: half/bin/3d
      bin: standard
Setting up Verlet run ...
  Unit style    : lj
  Current step  : 0
  Time step     : 0.0025
Per MPI rank memory allocation (min/avg/max) = 4.501 | 4.501 | 4.501 Mbytes
   Step          Temp          E_pair         E_mol          TotEng         Press
0 9.000000e-01 -4.704000e-01 0.000000e+00 8.786160e-01 -4.451195e-03
1000 8.797737e-01 -1.463792e+00 0.000000e+00 -1.450932e-01 5.004053e-02
2000 8.778546e-01 -1.591754e+00 0.000000e+00 -2.759314e-01 3.936154e-02
3000 8.811496e-01 -1.638126e+00 0.000000e+00 -3.173651e-01 3.638044e-02
4000 9.018802e-01 -1.752085e+00 0.000000e+00 -4.002504e-01 4.542555e-02
5000 8.977694e-01 -1.832009e+00 0.000000e+00 -4.863368e-01 3.909381e-02
6000 9.158538e-01 -1.845641e+00 0.000000e+00 -4.728619e-01 4.086728e-02
7000 8.964092e-01 -1.842455e+00 0.000000e+00 -4.988211e-01 6.344970e-02
8000 9.088412e-01 -1.943820e+00 0.000000e+00 -5.815516e-01 1.599530e-02
9000 9.199665e-01 -1.986218e+00 0.000000e+00 -6.072741e-01 3.355576e-02
10000 8.867009e-01 -2.046282e+00 0.000000e+00 -7.172006e-01 1.817911e-02
```

There will be two extra files in the `MD` folder after the simulation. `log.lammps` is identical to the screen output. `lj_cut.lammpstrj` is the trajectory of the simulation. Copy this file to your local machine and let's visualize it in OVITO.

Open Ovito, click `File` -> `Load File`, then select `lj_cut.lammpstrj`. Press the ▶ button at center bottom. Enjoy!

What do you see as the system evolves? Is the phenomena you see consistent with [this paper](https://pure.uva.nl/ws/files/2199981/29999_3595309057smi922.pdf)?

![A snapshot of what you will see](https://github.com/Arrondissement5etDemi/WintersessionMD/blob/main/snapshot.png)

### OVITO features

We will demonstrate some useful features of OVITO, including  
-change color and size of particles  
-visualizing velocities  
-generate particle trajectories  
-create snapshots and videos  

> Excercise: Modify `lj.lmp` so that it simulates a quench of the LJ fluid at number density $\rho = 0.3$, from initial temperature $T_i=1.5$ to final temperature $T_f=0.01$. Set timestep $dt = 0.001$ and run 20000 steps. Visualize the trajectory in OVITO.



### Bonus skill! extra information that you can compute

    compute myRDF all rdf 200 1 1
    fix 2 all ave/time 50 6 1000 c_myRDF[*] file lj_cut1.rdf mode vector ave window 6

    compute myvacf all vacf
    fix storeMyVacf all vector 1 c_myvacf[4]

    compute msd_1 all msd
    fix store_msd_1 all vector 10 c_msd_1[4]
    variable fitslope_1 equal slope(f_store_msd_1)/6/(10*dt)
    fix 3 all ave/time 100 1 100 c_msd_1[4] v_fitslope_1 c_myvacf[4] file lj_cut1_diffusion.txt
