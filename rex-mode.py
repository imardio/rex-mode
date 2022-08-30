"""
Reactive Mode Composition Factor Analysis

A program to calculate the kinetic energy distribution at the reactive mode of transition states
Written by: Mauricio Maldonado-Dominguez and Santiago Alonso-Gil

The theory behind the program is described in: 
Phys. Chem. Chem. Phys., 2019, 21, 24912-24918.

Applications of the method to the analysis of nonequilibrium reactivity:
J. Am. Chem. Soc. 2020, 142, 8, 3947-3958.
Chem. Sci. 2021, 12, 12682-12694.

Points to implement:
Program capabilities:
    - Identify and respect frozen atoms.
    - Give the user the possibility to analyze real modes too.
Program reshaping:
    - Rewrite the fortran subroutine in Python syntax. ímajor update)
    - Make global variables local. (FIrst step toward a cleaner code).
"""


########
#Modules
########

import pandas as pd
import argparse
import os
import sys
import logging
from datetime import datetime
import itertools
import subprocess as sp
import time
from io import StringIO
import re

########
#Dictionaries 
########

atomic_symbol = {
    1:'H', 2:'He', 3:'Li', 4:'Be', 5:'B', 6:'C', 7:'N', 8:'O', 9:'F', 10:'Ne',
    11:'Na', 12:'Mg', 13:'Al', 14:'Si',15:'P',16:'S',17:'Cl',18:'Ar',19:'K',
    20:'Ca', 21 : "Sc" , 22 : "Ti" , 23 : "V"  , 24 : "Cr" , 25 : "Mn", 26 : "Fe" ,
    27: "Co" , 28 : "Ni" , 29 : "Cu" , 30 : "Zn", 31 : "Ga" , 32 : "Ge" , 33 : "As" ,
    34: "Se" , 35 : "Br", 36 : "Kr" , 37 : "Rb" , 38 : "Sr" , 39 : "Y"  , 40 : "Zr",
    41: "Nb" , 42 : "Mo" , 43 : "Tc" , 44 : "Ru" , 45 : "Rh", 46 : "Pd" , 47 : "Ag" ,
    48: "Cd" , 49 : "In" , 50 : "Sn", 51 : "Sb" , 52 : "Te" , 53 : "I"  , 54 : "Xe" ,
    55: "Cs", 56 : "Ba" , 57 : "La" , 58 : "Ce" , 59 : "Pr" , 60 : "Nd", 61 : "Pm" ,
    62: "Sm", 63 : "Eu" , 64 : "Gd" , 65 : "Tb", 66 : "Dy" , 67 : "Ho" , 68 : "Er" ,
    69: "Tm", 70 : "Yb", 71 : "Lu" , 72 : "Hf" , 73 : "Ta" , 74 : "W"  , 75 : "Re",
    76: "Os", 77 : "Ir" , 78 : "Pt" , 79 : "Au" , 80 : "Hg", 81 : "Tl" , 82 : "Pb" ,
    83: "Bi",
    }

atomic_mass = {
    1:1.0, 2:4.0, 3:7.0, 4:9.0, 5:11.0, 6:12.0, 7:14.0, 8:16.0, 9:19.0, 10:20.0,
    11:23.0, 12:24.0, 13:27.0, 14:28.0, 15:31.0, 16:32.0, 17:35.0, 18:40.0, 19:39.0, 
    20:40.0, 21:45.0, 22:48.0, 23:51.0, 24:52.0, 25:55.0, 26:56.0, 27:59.0, 28:59.0, 
    29:64.0, 30:65.0, 31:70.0, 32:73.0, 33:75.0, 34:79.0, 35:80.0, 36:84.0, 37:85.0, 
    38:88.0, 39:89.0, 40:91.0, 41:93.0, 42:96.0, 43:97.0, 44:101.0, 45:104.0, 46:106.0, 
    47:108.0, 48:112.0, 49:115.0, 50:119.0, 51:122.0, 52:128.0, 53:127.0, 54:131.0,
    55:133.0, 56:137.0,
    }

#######
# Program
#######

def rmcf(filename1):

    logging.basicConfig(filename='RMCF_' + filename1[:-4] + '.LOG', level=logging.DEBUG, format='%(message)s', filemode='w')
    
    datetime_now = datetime.now()
    formatted_datetime = datetime_now.strftime("%Y %b %d %H:%M:%S")

    print('Analyzing ' + filename1)

    logging.info('#####################################################################\n')
    logging.info('RMCF: Reactive Mode Composition Factor Analysis\n')
    logging.info('Kinetic energy distribution at the reactive mode of transition states')
    logging.info('Written by: Mauricio Maldonado-Dominguez and Santiago Alonso-Gil')
    logging.info('J. Heyrovsky Institute of Physical Chemistry of the CAS')
    logging.info('Run date: '+formatted_datetime+'\n')
    logging.info('#####################################################################\n')
    logging.info('Reading Gaussian output: ' + filename1)

    check_success()

    if "NAtoms=" in open(args.GaussFile,'rt').read():
        n_atoms_1()
    elif "Input orientation:" in open(args.GaussFile,'rt').read():
        n_atoms_2()
    elif "Standard orientation:" in open(args.GaussFile,'rt').read():
        n_atoms_3()
    else:
        logging.info('The program aborted with Error 1. Unrecognized format in Gaussian output file.' + '\n')
        print("Error 1. Unrecognized format in Gaussian output file. Aborting")
        exit()

    logging.info('The number of atoms is: ' + str(total_atoms) + '\n')

    if "Input orientation" in open(args.GaussFile,'rt').read():  
        last_input()
        logging.info('The last geometry for ' + filename1 + ' is:' + '\n')
        logging.info(last_xyz.to_string(header=True, index=True) + '\n')
    elif "Standard orientation:" in open(args.GaussFile,'rt').read():
        last_standard()
        logging.info('The last geometry for ' + filename1 + ' is:' + '\n')
        logging.info(last_xyz.to_string(header=True, index=True) + '\n')
    else:
        logging.info('The program aborted with Error 2. Unrecognized format in Gaussian output file.' + '\n')
        print("Error 2. Unrecognized format in Gaussian output file. Aborting")
        exit()

    if "atoms frozen in the vibrational analysis" in open(args.GaussFile, 'rt').read():
        logging.info('Frozen atoms found. This is currently unsupported. Aborting the RMCF analysis.' + '\n')
        print('Frozen atoms found. Aborting the RMCF analysis.')
        exit()
    elif "Harmonic frequencies" in open(args.GaussFile, 'rt').read():
        logging.info('Harmonic frequencies found. The program will now perform the RMCF analysis.' + '\n')
        normal_modes()
        split_modes()
        delete_lines('modes_1', [0,1,2,3,4,5,6])
        ked()
        ambimodal()
        find_pairs()
    else:
        print('Error 3. Normal vibrational modes not found. Aborting' + '\n')
        logging.info('The program aborted with Error 3. Normal vibrational modes not found.' + '\n')
        exit()

    logging.info('The reactive mode for ' + filename1)
    logging.info('and its kinetic energy distribution (KED) are:' + '\n')

    sort_pairs()

    logging.info(mode_xyz.to_string(header=True, index=True) + '\n')
    logging.info('The most relevant atom pairs for ' + filename1 + ' are:' + '\n')
    logging.info(sorted.to_string(header=True, index=False) + '\n')
    logging.info('Where dR(%) is the interatomic distance change')
    logging.info('upon perturbation of the TS along the reactive mode' + '\n')
    logging.info('These atomic pairs represent potential bond formation/cleavage events')
    logging.info('and can be used to partition the TS structure for the calculation of')
    logging.info('product distributions in N-furcating reactions.' + '\n')

    remove_aux()

    print('Success!')

########
#Functions
########

def check_success():
    with open(args.GaussFile) as f:
        if (' Normal termination of Gaussian') in f.read():
            logging.info('The Gaussian calculation terminated normally.'+'\n')
        else: 
            logging.info('The Gaussian calculation did not finish correctly. Aborting.'+'\n')
            print('Error. Please check your Gaussian output.')
            exit()

####Extract the total number of atoms as an integer, without resorting to the "NAtoms" identifier.

def n_atoms_1():
    global total_atoms
    if "NAtoms" in open(args.GaussFile,'rt').read():
        with open(args.GaussFile,'rt') as f:
            for line in f.readlines():
                if "NAtoms" in line:
                    if re.split(r'\s',line)[6].rstrip('\n') is not "":
#                        If true, we found the number of atoms at position 6
                        total_atoms = re.split(r'\s',line)[6].rstrip('\n')
                        with open('NAtoms', 'w') as atoms:
                            atoms.write(str(total_atoms))
                        break
                    elif re.split(r'\s',line)[8].rstrip('\n') is not "":
#                        If true, we found the number of atoms at position 8
                        total_atoms = re.split(r'\s',line)[8].rstrip('\n')
                        with open('NAtoms', 'w') as atoms:
                            atoms.write(str(total_atoms))
                        break
              
def n_atoms_2():
    global total_atoms
    with open(args.GaussFile,'rt') as file: 
        for line in file:
            lines.append(str(line.rstrip('\n'))) 
    index1 = lines.index(str("Standard orientation:"))
    aux_index = []
    for line in lines:
        substring = line[:21]
        if substring == " Rotational constants":
            index3 = lines.index(line)
            aux_index.append(index3)
        else:
            pass
    index2 = aux_index[0]
    input_orientation = lines[index1 + 5:index2 - 1]
    total_atoms = len(input_orientation)
    with open('NAtoms', 'w') as atoms:
        atoms.write(str(total_atoms))

def n_atoms_3():
    global total_atoms
    with open(args.GaussFile,'rt') as file: 
        for line in file:
            lines.append(str(line.rstrip('\n'))) 
    index1 = lines.index(str("Input orientation:"))
    aux_index = []
    for line in lines:
        substring = line[:21]
        if substring == " Rotational constants":
            index3 = lines.index(line)
            aux_index.append(index3)
        else:
            pass
    index2 = aux_index[0]
    input_orientation = lines[index1 + 5:index2 - 1]
    total_atoms = len(input_orientation)
    with open('NAtoms', 'w') as atoms:
        atoms.write(str(total_atoms))

#### Extract normal modes

def normal_modes():
    global lines
    lines = []
    with open(args.GaussFile,'rt') as file: 
        for line in file:
            lines.append(str(line.rstrip('\n'))) 
    index1 = lines.index(str(" and normal coordinates:"))
    index2 = lines.index(str(" - Thermochemistry -"))
    list = lines[index1 + 1:index2 - 2]
    with open('modes', 'w') as file:
        for line in list:
            file.writelines(str(line) + '\n')

#### Cut the 'modes' file into 3-mode chunks and remove all headers

def split_modes():
    total_lines = int(total_atoms) + 7
    with open('modes_1','w') as output:
        with open("modes") as file:
            head = list(itertools.islice(file, total_lines))
            output.writelines(head)

def delete_lines(original_file, line_numbers):
    is_skipped = False
    counter = 0
    # Create name of dummy / temporary file
    dummy_file = original_file + '.bak'
    # Open original file in read only mode and dummy file in write mode
    with open(original_file, 'r') as read_obj, open(dummy_file, 'w') as write_obj:
        for line in read_obj:
            # If current line number exist in list then skip copying that line
            if counter not in line_numbers:
                write_obj.write(line)
            else:
                is_skipped = True
            counter += 1

    if is_skipped:
        os.remove(original_file)
        os.rename(dummy_file, original_file)
    else:
        os.remove(dummy_file)

#### Calculate KED. Currently implemented only for the reactive mode. 
#### Future implementations will address real modes.

def ked():
    global mode_xyz
    mode_xyz = pd.read_table('modes_1', usecols=[2,3,4], delim_whitespace=True, header=None)
    mode_atoms = pd.read_table('modes_1', usecols=[1], delim_whitespace=True, header=None)
    mode_symbols = pd.read_table('modes_1', usecols=[1], delim_whitespace=True, header=None)
    mode_symbols.replace(atomic_symbol, inplace=True)
    mode_atoms.replace(atomic_mass, inplace=True)
    mode_masses = pd.to_numeric(mode_atoms[1], errors='coerce')
    norm2 = (mode_xyz**2).sum(axis=1)
    for k, v in atomic_mass.items():
        atomic_mass[k] = float(v)
    preKED = norm2*mode_masses
    normalize = preKED.sum()
    KED = preKED.div(normalize)
    mode_xyz[5] = KED
    mode_xyz.insert(0, 1, mode_symbols)
    mode_xyz.columns=['Atom', 'dx', 'dy', 'dz', 'KED']
    mode_xyz.index += 1

### Extract the last geometry

def unique_geometry():
    with open('last_geometry', 'w') as unique:
        for element in input_orientation:
            unique.write(element + '\n')
    global only_xyz
    only_xyz = pd.read_table('last_geometry', usecols=[3,4,5], delim_whitespace=True, header=None)
    only_symbols = pd.read_table('last_geometry', usecols=[1], delim_whitespace=True, header=None)
    only_symbols.replace(atomic_symbol, inplace=True)
    only_xyz.insert(0, 1, only_symbols)
    only_xyz.columns=['Atom', 'X', 'Y', 'Z']
    only_xyz.index += 1

def last_input():
    lines = []
    with open(args.GaussFile,'rt') as file: 
        for line in file:
            lines.append(str(line.rstrip('\n'))) 
    input = "Input orientation:"
    locations = [] # Here we will compile where the geometry definitions begin throughout the output file
    with open(args.GaussFile,'r') as f:
        for num, line in enumerate(f, 1):
            if input in line:
                a = (int(num) + 5)
                locations.append(str(a) + '\n')
            else:
                pass
    limit = " ---------------------------------------------------------------------\n"

    delimiters = [] # Here we will compile where the geometry definitions end throughout the output file
    for location in locations:
        with open(args.GaussFile,'r') as f:
                skipped = itertools.islice(f, int(location), None)
                for num, line in enumerate(skipped, int(location) + 1):
                    if limit == line:
                        b = (int(num))
                        delimiters.append(str(b) + '\n')
                        break
                    else:
                        pass
    with open('geometry_locations','w') as output:
        for line in locations:
            output.writelines(str(line))
    with open('delimiter_locations','w') as output:
        for line in delimiters:
            output.writelines(str(line))
    with open('geometry_locations') as file:
        for line in file:
            pass
        index1 = int(line)
    with open('delimiter_locations') as file:
        for line in file:
            pass
        index2 = int(line)
    global last_xyz
    last_xyz = lines[index1 - 1:index2 - 1]
    with open('last_geometry', 'w') as file:
        for line in last_xyz:
            file.writelines(str(line) + '\n')
    last_xyz = pd.read_table('last_geometry', usecols=[3,4,5], delim_whitespace=True, header=None)
    geom_symbols = pd.read_table('last_geometry', usecols=[1], delim_whitespace=True, header=None)
    geom_symbols.replace(atomic_symbol, inplace=True)
    last_xyz.insert(0, 1, geom_symbols)
    last_xyz.columns=['Atom', 'X', 'Y', 'Z']
    last_xyz.index += 1

def last_standard():
    lines = []
    with open(args.GaussFile,'rt') as file: 
        for line in file:
            lines.append(str(line.rstrip('\n'))) 
    input = "Standard orientation:"
    locations = [] # Here we will compile where the geometry definitions begin throughout the output file
    with open(args.GaussFile,'r') as f:
        for num, line in enumerate(f, 1):
            if input in line:
                a = (int(num) + 5)
                locations.append(str(a) + '\n')
            else:
                pass
    limit = " ---------------------------------------------------------------------\n"

    delimiters = [] # Here we will compile where the geometry definitions end throughout the output file
    for location in locations:
        with open(args.GaussFile,'r') as f:
                skipped = itertools.islice(f, int(location), None)
                for num, line in enumerate(skipped, int(location) + 1):
                    if limit == line:
                        b = (int(num))
                        delimiters.append(str(b) + '\n')
                        break
                    else:
                        pass
    with open('geometry_locations','w') as output:
        for line in locations:
            output.writelines(str(line))
    with open('delimiter_locations','w') as output:
        for line in delimiters:
            output.writelines(str(line))
    with open('geometry_locations') as file:
        for line in file:
            pass
        index1 = int(line)
    with open('delimiter_locations') as file:
        for line in file:
            pass
        index2 = int(line)
    global last_xyz
    last_xyz = lines[index1 - 1:index2 - 1]
    with open('last_geometry', 'w') as file:
        for line in last_xyz:
            file.writelines(str(line) + '\n')
    last_xyz = pd.read_table('last_geometry', usecols=[3,4,5], delim_whitespace=True, header=None)
    geom_symbols = pd.read_table('last_geometry', usecols=[1], delim_whitespace=True, header=None)
    geom_symbols.replace(atomic_symbol, inplace=True)
    last_xyz.insert(0, 1, geom_symbols)
    last_xyz.columns=['Atom', 'X', 'Y', 'Z']
    last_xyz.index += 1

def ambimodal():
    subroutine = """      program ambimodal
! ===============================================================
! The subroutiune perturbs the TS and calculates dR(%)
!
! Written by Santiago Alonso-Gil and Mauricio Maldonado-Dominguez
! ===============================================================
       implicit none
       integer :: NAtoms,i,j,N
       real*8 :: a,b,c,SUM,angle
       character(5) :: word
       real*8,allocatable,dimension (:) :: x0,y0,z0,xf,yf,zf,KED,vx,vy,vz
       real*8,dimension(82) :: mAN
       real*8,external :: module2
       integer,allocatable,dimension(:) :: AN
       integer,allocatable,dimension(:,:) :: candidate
!Reading the number of atoms from NAtoms file.    
       open(40,file="NAtoms",status="old",action="read")
       read(40,*) NAtoms
       close(40)
!Allocating the variables
       allocate (x0(NAtoms),y0(NAtoms),z0(NAtoms),xf(NAtoms),yf(NAtoms))
       allocate (zf(NAtoms),KED(NAtoms),AN(NAtoms),vx(NAtoms),vy(NAtoms))
       allocate (vz(NAtoms),candidate(NAtoms,NAtoms))
       x0=0.d0
       y0=0.d0
       z0=0.d0
       xf=0.d0
       yf=0.d0
       zf=0.d0
       vx=0.d0
       vy=0.d0
       vz=0.d0
       KED=0.d0
       SUM=0.d0
       mAN=0.d0
       AN=0
       candidate=0
!Assign a mass to each atom
       mAN(1)=1.00
       mAN(3)=7.00
       mAN(5)=11.00
       mAN(6)=12.00
       mAN(7)=14.00
       mAN(8)=16.00
       mAN(9)=19.00
       mAN(11)=23.00
       mAN(12)=24.00
       mAN(13)=27.00
       mAN(14)=28.00
       mAN(15)=31.00
       mAN(16)=32.00
       mAN(17)=35.00
       mAN(19)=39.00
       mAN(20)=40.00
       mAN(23)=51.00
       mAN(24)=52.00
       mAN(25)=55.00
       mAN(26)=56.00
       mAN(27)=59.00
       mAN(28)=58.00
       mAN(29)=63.00
       mAN(30)=64.00
       mAN(35)=79.00
       mAN(44)=102.00
       mAN(47)=107.00
       mAN(53)=127.00
       mAN(76)=192.00
       mAN(78)=195.00
       mAN(79)=197.00
       mAN(82)=208.00
!Open the coord.txt file and calculate the ri values 
       open(41,file="last_geometry",status="old",action="read")
!       do i=1,3
!          read(41,*)
!       enddo
       do i=1,NAtoms
          read(41,*) N,AN(i),N,x0(i),y0(i),z0(i)
       enddo
       close(41)
!Open the freqs.txt file and calculate the KED values AND HERE WE MUST DECLARE FILE modes_1
       open(42,file="modes_1",status="old",action="read")
!       do i=1,11
!          read(42,*)
!       enddo
       do i=1,NAtoms
          read(42,*) N,N,vx(i),vy(i),vz(i)
          KED(i)=mAN(AN(i))*module2(vx(i),vy(i),vz(i))
          SUM=SUM+KED(i)
       enddo
       KED=KED/SUM
       close(44)
!First criterion: KED(i) > 0.0002
       do i=1,NAtoms
          if (KED(i) > 2.d-4) then
             candidate(i,i+1:NAtoms)=1
             candidate(i+1:NAtoms,i)=1
          else 
             candidate(i,:)=0
             candidate(:,i)=0
          endif
       enddo
!Second criterion: if the change in distance (dR) is less than a series of empirical threshdolds, out
       do i=1,NAtoms
          do j=1,NAtoms
             if (candidate(i,j)==1)then
                a=dsqrt((x0(i)-x0(j))*(x0(i)-x0(j))+(y0(i)-y0(j))*(y0(i)-y0(j))+(z0(i)-z0(j))*(z0(i)-z0(j)))
                xf(i)=x0(i)+0.5*vx(i)
                yf(i)=y0(i)+0.5*vy(i)
                zf(i)=z0(i)+0.5*vz(i)
                xf(j)=x0(j)+0.5*vx(j)
                yf(j)=y0(j)+0.5*vy(j)
                zf(j)=z0(j)+0.5*vz(j)
                b=dsqrt((xf(i)-xf(j))*(xf(i)-xf(j))+(yf(i)-yf(j))*(yf(i)-yf(j))+(zf(i)-zf(j))*(zf(i)-zf(j)))
                c=abs(b-a)*100/a
                if (c<1.5d0) then
                   candidate(i,j)=0
                endif
                if (a>3.3d0) then
                   candidate(i,j)=0
                endif
                if (AN(i)+AN(j)==2) then
                   candidate(i,j)=0
                endif
                if (AN(i)==6) then
                   if (AN(j)==1) then
                      if (a<1.15d0) then
                         candidate(i,j)=0
                      endif
                   endif
                   if (AN(j)==6) then
                      if (a<1.53d0) then
                         candidate(i,j)=0
                      endif
                   endif
                   if (AN(j)==7) then
                      if (a<1.5d0) then
                         candidate(i,j)=0
                      endif
                      if (a>2.4d0) then
                         candidate(i,j)=0
                      endif
                   endif
                   if (AN(j)==8) then
                      if (a<1.5d0) then
                         candidate(i,j)=0
                      endif
                   endif
                   if (AN(j)==17) then
                      if (a<1.8d0) then
                         candidate(i,j)=0
                      endif
                   endif
                endif
                if (AN(j)==6) then
                   if (AN(i)==1) then
                      if (a<1.15d0) then
                         candidate(i,j)=0
                      endif
                   endif
                   if (AN(i)==6) then
                      if (a<1.53d0) then
                         candidate(i,j)=0
                      endif
                   endif
                   if (AN(i)==7) then
                      if (a<1.5d0) then
                         candidate(i,j)=0
                      endif
                      if (a>2.4d0) then
                         candidate(i,j)=0
                      endif
                   endif
                   if (AN(i)==8) then
                      if (a<1.5d0) then
                         candidate(i,j)=0
                      endif
                   endif
                   if (AN(i)==17) then
                      if (a<1.8d0) then
                         candidate(i,j)=0
                      endif
                   endif
                endif
                if (AN(i)==7) then
                   if (AN(j)==1) then
                      if (a<1.1d0) then
                         candidate(i,j)=0
                      endif
                   endif
                   if (AN(j)==6) then
                      if (a>2.4d0) then
                         candidate(i,j)=0
                      endif
                   endif
                   if (AN(j)==7) then
                      if (a<1.5d0) then
                         candidate(i,j)=0
                      endif
                      if (a>2.5d0) then
                         candidate(i,j)=0
                      endif
                   endif
                   if (AN(j)==8) then
                      if (a<1.5d0) then
                         candidate(i,j)=0
                      endif
                   endif
                endif
                if (AN(j)==7) then
                   if (AN(i)==1) then
                      if (a<1.1d0) then
                         candidate(i,j)=0
                      endif
                   endif
                   if (AN(i)==7) then
                      if (a<1.5d0) then
                         candidate(i,j)=0
                      endif
                   endif
                   if (AN(i)==8) then
                      if (a<1.5d0) then
                         candidate(i,j)=0
                      endif
                   endif
                endif
                if (AN(i)==8) then
                   if (AN(j)==1) then
                      if (a<1.1d0) then
                         candidate(i,j)=0
                      endif
                   endif
                   if (AN(j)==7) then
                      if (a<1.5d0) then
                         candidate(i,j)=0
                      endif
                   endif
                   if (AN(j)==8) then
                      if (a<1.5d0) then
                         candidate(i,j)=0
                      endif
                   endif
                endif
                if (AN(j)==8) then
                   if (AN(i)==1) then
                      if (a<1.1d0) then
                         candidate(i,j)=0
                      endif
                   endif
                   if (AN(i)==7) then
                      if (a<1.5d0) then
                         candidate(i,j)=0
                      endif
                   endif
                   if (AN(i)==8) then
                      if (a<1.5d0) then
                         candidate(i,j)=0
                      endif
                   endif
                endif
                if ( AN(j)==1 ) then
                   if ((a>1.8d0)) then
                      candidate(i,j)=0
                   endif
                endif
                if ( AN(i)==1 ) then
                   if ((a>1.8d0)) then
                      candidate(i,j)=0
                   endif
                endif
             endif
          enddo
       enddo
!Printing
       write(6,*) "Atom1 Atom2 KED(1) KED(2) SUM(KED) angle  dR dR(%)"
       do i=1,NAtoms
          do j=i+1,NAtoms
             if(candidate(i,j)==1)then
                angle=vx(i)*vx(j)+vy(i)*vy(j)+vz(i)*vz(j)
                angle=angle/(dsqrt(module2(vx(i),vy(i),vz(i)))*dsqrt(module2(vx(j),vy(j),vz(j))))
                angle=acos(angle)*180.d0/3.14159265d0
                a=dsqrt((x0(i)-x0(j))*(x0(i)-x0(j))+(y0(i)-y0(j))*(y0(i)-y0(j))+(z0(i)-z0(j))*(z0(i)-z0(j)))
                b=dsqrt((xf(i)-xf(j))*(xf(i)-xf(j))+(yf(i)-yf(j))*(yf(i)-yf(j))+(zf(i)-zf(j))*(zf(i)-zf(j)))
                c=(b-a)*100/a
                write(6,"(3X,I3,2X,I3,2X,F5.3,2X,F5.3,2X,F5.3,2X,F5.1,2X,F5.3,2X,F5.1)") i,j,KED(i),KED(j),KED(i)+KED(j),angle,b-a,c
             endif
          enddo
       enddo
end program

real*8 function module2(x,y,z)
implicit none
real*8 :: x,y,z
module2 = x*x+y*y+z*z
end function module2
    """
    with open('ambimodal.f90', 'w') as output:
        output.write(subroutine)
    os.popen('f95 ambimodal.f90 -o ambimodal.x')
    time.sleep(0.15)

def find_pairs():
    file = open('atom_pairs', 'w')
    here = os.path.dirname(os.path.abspath(__file__))
    filename = os.path.join(here, 'ambimodal.x')
    sp.call([filename], stdout=file)

#### Sort the atomic pairs

def sort_pairs():
    global sorted
    sorted = pd.read_table('atom_pairs', delim_whitespace=True).drop(columns=["dR","angle"]).sort_values(by="SUM(KED)",ascending=False)

### Remove auxiliary files

def remove_aux():
    os.remove("delimiter_locations")
    os.remove("geometry_locations")
    os.remove("modes")
    os.remove("atom_pairs")
    os.remove("modes_1")
    os.remove("last_geometry")
    os.remove("ambimodal.f90")
    os.remove("ambimodal.x")
    os.remove("NAtoms")

#######################
# Parsing the input and executing the code
#######################

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This program calculates Kinetic Energy Distributions from a Gaussian TS calculation')
    parser.add_argument('GaussFile', help='Gaussian output file (*.out or *.log)')

    args = parser.parse_args()

    if not os.path.isfile(args.GaussFile):
        print(args.GaussFile + " is not a valid Gaussian output file.")
        quit()

    rmcf(args.GaussFile)
