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
    
## Your first LAMMPS job: A Lennard-Jones fluid

    units lj  
    atom_style atomic  
    boundary p p p  
    
    lattice fcc 0.238732  
    region box block 0 5 0 5 0 5  
    create_box 1 box  
    variable rnseed index 10  
    variable number equal ceil(random(1,1000000000,${rnseed})/100000.0)   
    print ${number}  
    create_atoms 1 random 500 ${number} NULL overlap 0.5 maxtry 500  
      
    timestep 0.00025

    mass 1 1.0
    velocity all create 1 ${number}
    read_dump lj_cut1_31.lammpstrj 3300000 x y z vx vy vz box no

    pair_style lj/cut 1
    pair_coeff 1 1 80.0 1.0 1

    fix 1 all nvt temp 1 1 $(100.0*dt)

    thermo 1000
    thermo_modify format line "%d %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e"

    dump my_dump2 all custom 100 lj_cut1_32.lammpstrj id type x y z vx vy vz

    compute myRDF all rdf 200 1 1
    fix 2 all ave/time 50 6 1000 c_myRDF[*] file lj_cut1.rdf mode vector ave window 6

    compute myvacf all vacf
    fix storeMyVacf all vector 1 c_myvacf[4]

    compute msd_1 all msd
    fix store_msd_1 all vector 10 c_msd_1[4]
    variable fitslope_1 equal slope(f_store_msd_1)/6/(10*dt)
    fix 3 all ave/time 100 1 100 c_msd_1[4] v_fitslope_1 c_myvacf[4] file lj_cut1_diffusion.txt

    run 100000

to be continued...
