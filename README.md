# rex-mode
Computes kinetic energy distributions from molecular frequency calculations

The program takes in the results from a TS frequency calculation carried out in the Gaussian software and performs
the following tasks:

1) Calculates the kinetic energy distribution within the reactive mode.
2) Displaces the TS geometry along the reactive mode.
3) Gauges all interatomic displacements and ranks them as a list of potential bond
formation/breaking events.

With these results, the program suggests the user atom pairs which are most likely to be engaged in bond
formation/cleavage. Not all of the suggested pairs will lead to a real bond formation/cleavage and, hence,
these candidates must be examined by the user to select the chemically relevant ones.

Execution of the program requires Python 3, the Pandas module and Fortran 95.

To execute the program with a Gaussian output (exemplified as output.out), the user only needs to enter
the command:

python3 rmcf.py output.out

An exemplary output is presented below:

#####################################################################
RMCF: Reactive Mode Composition Factor Analysis
Kinetic energy distribution at the reactive mode of transition states
Written by: Mauricio Maldonado-Dominguez and Santiago Alonso-Gil
J. Heyrovsky Institute of Physical Chemistry of the CAS
Originally published in https://doi.org/10.1039/D1SC02826J
Run date: 2021 Jun 25 12:31:13
#####################################################################
Reading Gaussian output: TS-rec4-Gas.log
The Gaussian calculation terminated normally.
The number of atoms is: 19
The last geometry for TS-rec4-Gas.log is:
S6
 Atom X Y Z
1 C 0.093550 1.645520 0.043167
2 C -0.719262 0.783153 -0.689887
3 H 0.535988 2.469753 -0.501510
4 H -0.838256 0.793845 -1.759106
5 C 0.833545 -1.405411 0.482231
6 C 0.938055 -1.075590 -0.855154
7 H 0.292889 -2.250545 0.880086
8 H 0.575828 -1.669862 -1.680475
9 C 1.391141 -0.374958 1.247565
10 C 1.802298 0.667562 0.394837
11 H 1.372357 -0.324240 2.326097
12 H 2.465043 1.456700 0.726153
13 C 1.915339 0.050910 -0.986849
14 H 1.748996 0.723783 -1.824986
15 H 2.924411 -0.370363 -1.092421
16 H -0.207862 1.871674 1.055035
17 O -1.694823 -0.023149 1.210520
18 O -2.530058 -0.598291 -0.715933
19 N -1.715597 0.013939 -0.024369

The reactive mode for TS-rec4-Gas.log
and its kinetic energy distribution (KED) are:

 Atom dx dy dz KED
1 C 0.46 -0.28 0.10 0.377315
2 C 0.18 -0.25 -0.02 0.119860
3 H -0.08 -0.01 0.02 0.000723
4 H -0.08 0.14 0.03 0.002819
5 C 0.04 0.01 0.07 0.008301
6 C -0.18 0.19 -0.01 0.086279
7 H 0.14 -0.07 0.04 0.002736
8 H -0.17 0.14 0.01 0.005094
9 C -0.07 -0.04 -0.01 0.008301
10 C -0.44 0.30 -0.10 0.369266
11 H 0.02 -0.05 -0.01 0.000314
12 H -0.13 0.01 -0.04 0.001949
13 C -0.02 0.03 -0.02 0.002138
14 H 0.11 -0.05 -0.12 0.003039
15 H -0.02 -0.03 0.24 0.006173
16 H 0.02 0.01 -0.06 0.000430
17 O 0.02 0.01 -0.00 0.000838
18 O -0.02 0.02 -0.00 0.001342
19 N 0.04 -0.01 -0.02 0.003081

The most relevant atom pairs for TS-rec4-Gas.log are:

Atom1 Atom2 KED(1) KED(2) SUM(KED) dR(%)
 1 10 0.377 0.369 0.747 -27.2
 2 10 0.120 0.369 0.489 -10.8
 1 6 0.377 0.086 0.464 -9.2
 1 5 0.377 0.008 0.386 -5.9
 1 9 0.377 0.008 0.386 -8.8
 1 13 0.377 0.002 0.379 -8.8
 1 17 0.377 0.001 0.378 1.7
 5 10 0.008 0.369 0.378 2.2
 2 6 0.120 0.086 0.206 -11.4
 2 5 0.120 0.008 0.128 -3.9
 2 9 0.120 0.008 0.128 -3.9
 2 13 0.120 0.002 0.122 -4.8
 
Where dR(%) is the interatomic distance change upon perturbation of the TS along the reactive mode

As can be seen in detail in Chem. Sci., 2021,12, 12682-12694, inspection of TS1 for reaction 4 reveals that the 
1-10 pair is common to both products, and that 6-2 is a chemically relevant pair. 
The 9-17 pair is not suggested, which agrees well with the experimental
observation that the corresponding product is not formed at all. The user can still evaluate fragments not
suggested by the program by summation of their atomic KED components from the “The reactive
mode for TS-rec4-Gas.log and its kinetic energy distribution (KED) are:” section
of the output.
