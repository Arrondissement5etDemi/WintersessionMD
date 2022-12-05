# WintersessionMD
Introduction to Molecular Dynamics with LAMMPS for Princeton Wintersession 2023

## Introduction

Atoms and molecules are like balls and springs, interacting under rules of quantum mechanics. The interactions can be partially approximated by classical effective forces. [Molecular dynamics](https://en.wikipedia.org/wiki/Molecular_dynamics) (MD) solves the equations of motion for you! It is a computer simulation method for analyzing the physical movements of atoms and molecules. The atoms and molecules are allowed to interact for a fixed period of time, giving a view of the dynamic "evolution" of the system. 

Large-scale Atomic/Molecular Massively Parallel Simulator ([LAMMPS](https://www.lammps.org/#gsc.tab=0)) is an efficient and popular molecular dynamics program from Sandia National Laboratories. It's written for serial and parallel computing alike. For this workshop, we will mainly learn its serial version, but will touch upon the parallel version if we have time.

## Installation

To ensure that everyone has the same environment, let's use Princeton Research Computing's [Adroit](https://researchcomputing.princeton.edu/systems/adroit) cluster. SSH into your cluster by typing the following into your terminal:

    $ ssh <YourNetID>@adroit-vis.princeton.edu
    
and enter your password. (The password column remains blank when you type, which is expected.)

Now let's use the conda environment to install LAMMPS.
