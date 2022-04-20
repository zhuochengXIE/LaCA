"""
Ovito script to analyze Laves phases (C14, C15, C36) and their lattice defects
Based on Voronoi index to identify atomic clusters
Work with Ovito Pro version 3.7.3
@author: Zhuocheng Xie
@date:   Thu Apr 14 11:02:01 CEST 2022
@usage:  ovitos ./ovito_Laves_crystal_analysis_voronoibased.py -f <CONFIG> [-r <AtomicRadii>] [-t <ParticleType>] [-h <HELP>]
"""
import os, sys, subprocess, argparse
from multiprocessing import Array
import numpy as np
from math import *
from ovito.io import *
from ovito.modifiers import *
from ovito.data import *

#from datetime import datetime
#begin_time = datetime.now()

##########################################################################################
# Function to control option parsing in Python
##########################################################################################
def parseCmdLineArgs(argv):
    from argparse import RawTextHelpFormatter
    parser = argparse.ArgumentParser(description='[>>>ovito_Laves_crystal_analysis.py<<< \nPython script to analyze Laves phases (C14, C15, C36) and their defects. \n\n--> @Author: Zhuocheng Xie <--',formatter_class=RawTextHelpFormatter)
    parser.add_argument('-f','--fileConfig',nargs=1,dest='fConfig',help='Name of the configuration file',required=True)
    parser.add_argument('-r','--AtomicRadius',nargs=2,dest='radii',help='Atomic radii of type 1 and 2 particles in Angstrom for computing the poly-disperse Voronoi tessellation; by default a mono-disperse Voronoi tessellation is computed',required=False,type=float)
    parser.add_argument('-t','--type',nargs=1,dest='type',help='Particle type of binary AB2 Laves phase, e.g., 1 --> (A=AtomType1 and B=AtomType2), 2 --> (B=AtomType1 and A=AtomType2); by default particle type is NOT taken into account',required=False,type=int)
    args=parser.parse_args()
    return args

if __name__=="__main__":
    cargs = parseCmdLineArgs(sys.argv[1:])
    fConfig = cargs.fConfig[0]
    radii = [1,1]
    type = [0]
    if cargs.radii is not None: 
        radii = cargs.radii
        pdvoro = True
    else: 
        pdvoro = False
    if cargs.type is not None:
        type = cargs.type
        
##########################################################################################
# Define Laves crystals and defects 
##########################################################################################

# Define similarity between test and reference vector (A, B1, B2) using euclidean distance
def euclidean_distance(x,y): 
    return sqrt(sum(pow(a-b,2) for a, b in zip(x, y)))

# Define C14 neighbor list
def C14_A(i):
    global atom_voro5, atom_voro6
    voro5 = atom_voro5[i]
    voro6 = atom_voro6[i]
    deviation = euclidean_distance([A, B1, B2],[4, 3, 9])
    if voro5 == 12 and voro6 == 4 and deviation <= sqrt(2):
        atom_laves[i] = Laves_C14_A
        return True
    else: 
        return False

def C14_B1(i):
    global atom_voro5, atom_voro6
    voro5 = atom_voro5[i]
    voro6 = atom_voro6[i]
    deviation = euclidean_distance([A, B1, B2],[6, 0, 6])
    if voro5 == 12 and voro6 == 0 and deviation <= sqrt(2):
        atom_laves[i] = Laves_C14_B1
        return True
    else: 
        return False

def C14_B2(i):
    global atom_voro5, atom_voro6
    voro5 = atom_voro5[i]
    voro6 = atom_voro6[i]
    deviation = euclidean_distance([A, B1, B2],[6, 2, 4])
    if voro5 == 12 and voro6 == 0 and deviation < 1:
        atom_laves[i] = Laves_C14_B2
        return True
    else: 
        return False

# Define C15 neighbor list
def C15_A(i):
    global atom_voro5, atom_voro6
    voro5 = atom_voro5[i]
    voro6 = atom_voro6[i]
    deviation = euclidean_distance([A, B1, B2],[4, 12, 0])
    if voro5 == 12 and voro6 == 4 and deviation <= sqrt(2):
        atom_laves[i] = Laves_C15_A
        return True
    else: 
        return False

def C15_B1(i):
    global atom_voro5, atom_voro6
    voro5 = atom_voro5[i]
    voro6 = atom_voro6[i]
    deviation = euclidean_distance([A, B1, B2],[6, 6, 0])
    if voro5 == 12 and voro6 == 0 and deviation <= sqrt(2):
        atom_laves[i] = Laves_C15_B1
        return True
    else: 
        return False

# Define C14-C15 interface neighbor list
def inf_A1(i):
    global atom_voro5, atom_voro6
    voro5 = atom_voro5[i]
    voro6 = atom_voro6[i]
    deviation = euclidean_distance([A, B1, B2],[4, 6, 6])
    if voro5 == 12 and voro6 == 4 and deviation <= sqrt(2):
        atom_laves[i] = Laves_inf_A1
        return True
    else: 
        return False

def inf_A2(i):
    global atom_voro5, atom_voro6
    voro5 = atom_voro5[i]
    voro6 = atom_voro6[i]
    deviation = euclidean_distance([A, B1, B2],[4, 9, 3])
    if voro5 == 12 and voro6 == 4 and deviation <= sqrt(2):
        atom_laves[i] = Laves_inf_A2
        return True
    else: 
        return False

def inf_B1(i):
    global atom_voro5, atom_voro6
    voro5 = atom_voro5[i]
    voro6 = atom_voro6[i]
    if voro5 == 12 and voro6 == 0 and A == 6 and B1 == 3 and B2 == 3:
        atom_laves[i] = Laves_inf_B1
        return True
    else: 
        return False

##########################################################################################
# OUTPUT Style
##########################################################################################

def write_laves(i):
    global atom_nrs,atom_types,atom_pos,atom_voro5,atom_voro6,atom_csp,atom_laves
    # write line to outfile:
    outstr="%10d %3d %12.6f %12.6f %12.6f %2d %2d %12.6f %d\n" % \
            (atom_nrs[i],atom_types[i],atom_pos[i][0],atom_pos[i][1],atom_pos[i][2],atom_voro5[i],atom_voro6[i],atom_csp[i],atom_laves[i])
    f_laves.write(outstr)

def main():
    global atom_nrs,atom_types,atom_pos,atom_csp,atom_voro,atom_laves,else_csp,atom_voro5,atom_voro6,atom_types1,atom_pos1
    global Laves_C14_A,Laves_C14_B1,Laves_C14_B2,Laves_C15_A,Laves_C15_B1,Laves_inf_A1,Laves_inf_A2,Laves_inf_B1
    
Laves_C14_A  = 1
Laves_C14_B1 = 2
Laves_C14_B2 = 3
Laves_C15_A  = 4
Laves_C15_B1 = 5
Laves_inf_A1 = 6
Laves_inf_A2 = 7
Laves_inf_B1 = 8
OTHER = 0
OTHER_Laves = 9

##########################################################################################
# Calculate Voronoi index and CSP
##########################################################################################

pipeline = import_file(fConfig)
data = pipeline.compute()
filename=fConfig

# Set up the Voronoi analysis modifier.
def assign_particle_radii(frame, data):
    atom_types = data.particles_.particle_types_
    atom_types.type_by_id_(1).radius = radii[0]   # atomic radius assigned to atom type 1
    atom_types.type_by_id_(2).radius = radii[1]   # atomic radius assigned to atom type 2
pipeline.modifiers.append(assign_particle_radii)

pipeline.modifiers.append(VoronoiAnalysisModifier(compute_indices = True, use_radii = pdvoro))
data = pipeline.compute()

pipeline.modifiers.append(ComputePropertyModifier(output_property = 'voro5', expressions = ['VoronoiIndex.5']))
pipeline.modifiers.append(ComputePropertyModifier(output_property = 'voro6', expressions = ['VoronoiIndex.6']))

#print('voronoi time:', datetime.now() - begin_time)

# Set up the csp modifier.
def exportELSE_type0(frame, data):
    data.apply(ExpressionSelectionModifier(expression = '(VoronoiIndex.3 == 0 && VoronoiIndex.4 == 0 && VoronoiIndex.5 == 12 && VoronoiIndex.6 == 4) || (VoronoiIndex.3 != 0 || VoronoiIndex.4 != 0 || VoronoiIndex.5 != 12 || VoronoiIndex.6 != 0)'))
    data.apply(ComputePropertyModifier(output_property = 'Centrosymmetry', expressions = ['-1'], only_selected = True))
    
def exportICO_type0(frame, data):
    data.apply(ExpressionSelectionModifier(expression = '(VoronoiIndex.3 == 0 && VoronoiIndex.4 == 0 && VoronoiIndex.5 == 12 && VoronoiIndex.6 == 4) || (VoronoiIndex.3 != 0 || VoronoiIndex.4 != 0 || VoronoiIndex.5 != 12 || VoronoiIndex.6 != 0)'))
    data.apply(InvertSelectionModifier())
    data.apply(CentroSymmetryModifier(num_neighbors = 6, only_selected = True))
    
def exportELSE_type1(frame, data):
    data.apply(ExpressionSelectionModifier(expression = '(VoronoiIndex.3 == 0 && VoronoiIndex.4 == 0 && VoronoiIndex.5 == 12 && VoronoiIndex.6 == 4) || (ParticleType==1 && (VoronoiIndex.3 != 0 || VoronoiIndex.4 != 0 || VoronoiIndex.5 != 12 || VoronoiIndex.6 != 0))'))
    data.apply(ComputePropertyModifier(output_property = 'Centrosymmetry', expressions = ['-1'], only_selected = True))
    
def exportICO_type1(frame, data):
    data.apply(ExpressionSelectionModifier(expression = '(VoronoiIndex.3 == 0 && VoronoiIndex.4 == 0 && VoronoiIndex.5 == 12 && VoronoiIndex.6 == 4) || (ParticleType==1 && (VoronoiIndex.3 != 0 || VoronoiIndex.4 != 0 || VoronoiIndex.5 != 12 || VoronoiIndex.6 != 0))'))
    data.apply(InvertSelectionModifier())
    data.apply(CentroSymmetryModifier(num_neighbors = 6, only_selected = True))
    
def exportELSE_type2(frame, data):
    data.apply(ExpressionSelectionModifier(expression = '(VoronoiIndex.3 == 0 && VoronoiIndex.4 == 0 && VoronoiIndex.5 == 12 && VoronoiIndex.6 == 4) || (ParticleType==2 && (VoronoiIndex.3 != 0 || VoronoiIndex.4 != 0 || VoronoiIndex.5 != 12 || VoronoiIndex.6 != 0))'))
    data.apply(ComputePropertyModifier(output_property = 'Centrosymmetry', expressions = ['-1'], only_selected = True))
    
def exportICO_type2(frame, data):
    data.apply(ExpressionSelectionModifier(expression = '(VoronoiIndex.3 == 0 && VoronoiIndex.4 == 0 && VoronoiIndex.5 == 12 && VoronoiIndex.6 == 4) || (ParticleType==2 && (VoronoiIndex.3 != 0 || VoronoiIndex.4 != 0 || VoronoiIndex.5 != 12 || VoronoiIndex.6 != 0))'))
    data.apply(InvertSelectionModifier())
    data.apply(CentroSymmetryModifier(num_neighbors = 6, only_selected = True))

if type == [0]:
    pipeline.modifiers.append(exportICO_type0)
    pipeline.modifiers.append(exportELSE_type0)
    
elif type == [1]:
    pipeline.modifiers.append(exportICO_type1)
    pipeline.modifiers.append(exportELSE_type1)

elif type == [2]:
    pipeline.modifiers.append(exportICO_type2)
    pipeline.modifiers.append(exportELSE_type2)

data = pipeline.compute()

#print('csp time:', datetime.now() - begin_time)

##########################################################################################
# Build neighbor list
##########################################################################################

nr_atoms = data.particles.count

atom_nrs=data.particles['Particle Identifier'].array
atom_types = data.particles['Particle Type'].array
atom_pos = data.particles['Position'].array
atom_csp = data.particles['Centrosymmetry'].array
atom_voro5 = data.particles["voro5"].array
atom_voro6 = data.particles["voro6"].array
atom_laves = Array('i',[-1]*nr_atoms)

# Output header
filename_Laves = filename + ".voro_csp.laves"
f_laves = open(filename_Laves, 'a+')
header1 = "head -n 8 %s >> %s"
header2 = "echo \"ITEM: ATOMS id type x y z VoronoiIndex_5 VoronoiIndex_6 Centrosymmetry laves\" >> %s"
subprocess.check_call(header1 % (filename, filename_Laves), shell=True)
subprocess.check_call(header2 % (filename_Laves), shell=True)

# Cutoff to differeniate B1 and B2 with different chemical environment
cutoff_csp = 5

# Build neighbor list array
finderICO = NearestNeighborFinder(12, data)
nblICO = finderICO.find_all()

finderZ16 = NearestNeighborFinder(16, data)
nblZ16 = finderZ16.find_all()

# Iterate over all atoms 
for i in range(nr_atoms):
    count_A = 0
    count_B1 = 0
    count_B2 = 0
	# Visit the 12 nearest neighbors of ICO atom.
    if atom_voro5[i] == 12 and atom_voro6[i] == 0 and atom_laves[i] == -1:
            for neigh in range(12):
                nblICO_i = nblICO[0][i][neigh]
                csp_of_neighbor = atom_csp[nblICO_i]
                if (csp_of_neighbor > -1 and csp_of_neighbor < cutoff_csp):
                    count_B1 = count_B1+1
                if (csp_of_neighbor > cutoff_csp):
                    count_B2 = count_B2+1
                if (csp_of_neighbor < 0):
                    count_A = count_A+1
            A = count_A
            B1 = count_B1
            B2 = count_B2
            if C14_B1(i):
                atom_laves[i]=Laves_C14_B1
                write_laves(i)
                pass
            elif C14_B2(i):
                atom_laves[i]=Laves_C14_B2
                write_laves(i)
                pass
            elif C15_B1(i):
                atom_laves[i]=Laves_C15_B1
                write_laves(i)
                pass
            elif inf_B1(i):
                atom_laves[i]=Laves_inf_B1
                write_laves(i)
                pass
            else:
                atom_laves[i]=OTHER_Laves
                write_laves(i)
                pass
	# Visit the 16 nearest neighbors of Z16 atom.
    elif atom_voro5[i] == 12 and atom_voro6[i] == 4 and atom_laves[i] == -1:
            for neigh in range(16):
                nblZ16_i = nblZ16[0][i][neigh]
                csp_of_neighbor = atom_csp[nblZ16_i]
                if (csp_of_neighbor > -1 and csp_of_neighbor < cutoff_csp):
                    count_B1 = count_B1+1
                if (csp_of_neighbor > cutoff_csp):
                    count_B2 = count_B2+1
                if (csp_of_neighbor < 0):
                    count_A = count_A+1
            A = count_A
            B1 = count_B1
            B2 = count_B2
            if C14_A(i):
                atom_laves[i]=Laves_C14_A
                write_laves(i)
                pass
            elif C15_A(i):
                atom_laves[i]=Laves_C15_A
                write_laves(i)
                pass
            elif inf_A1(i):
                atom_laves[i]=Laves_inf_A1
                write_laves(i)
                pass
            elif inf_A2(i):
                atom_laves[i]=Laves_inf_A2
                write_laves(i)
                pass
            else:
                atom_laves[i]=OTHER_Laves
                write_laves(i)
                pass
    else:
                atom_laves[i]=OTHER
                write_laves(i)
                pass            
f_laves.close()
#print('total time:', datetime.now() - begin_time)