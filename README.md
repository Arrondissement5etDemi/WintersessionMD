# WintersessionMD
Introduction to Molecular Dynamics with LAMMPS for Princeton Wintersession 2023

## Introduction

Atoms and molecules are like balls and springs, interacting under rules of quantum mechanics. The interactions can be partially approximated by classical effective forces. [Molecular dynamics](https://en.wikipedia.org/wiki/Molecular_dynamics) (MD) solves the equations of motion for you! It is a computer simulation method for analyzing the physical movements of atoms and molecules. The atoms and molecules are allowed to interact for a fixed period of time, giving a view of the dynamic "evolution" of the system. 

Large-scale Atomic/Molecular Massively Parallel Simulator ([LAMMPS](https://www.lammps.org/#gsc.tab=0)) is an efficient and popular molecular dynamics program from Sandia National Laboratories. It's written for serial and parallel computing alike. For this workshop, we will mainly learn its serial version, but will touch upon the parallel version if we have time.

 ![LAMMPS logo](https://docs.lammps.org/_static/lammps-logo.png)

## Installation

### Get into the cluster

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

To be continued ...
