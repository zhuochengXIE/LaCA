# LaCA
Method for atomistic identification of Laves phases and their lattice defects

Laves phase Crystal Analysis (LaCA) is a tool to analyze Laves phases and their lattice defects in atomistic simulations.
This document gives a tutorial on how to use LaCA on analzying Laves phase atomistic structures.

The workflow and case studies of the LaCA algorithm can be found in this paper:
Xie, Z., et al. Laves phase crystal analysis (LaCA): Atomistic identification of lattice defects in C14 and C15 topologically close-packed phases. Journal of Materials Research 36, 2010–2024 (2021). https://doi.org/10.1557/s43578-021-00237-y

If you use LaCA in your work, the citation of the above-mentioned paper will be greatly appreciated!

################ Implementation details ################
Characterization of the atomic cluster is the first step of LaCA, Common Neighbor Analysis (CNA) and Voronoi index can both be used to fulfill the task. Compare to the Voronoi index, the CNA is less computational cost however more sensitive to perturbation of atomic coordinates. Therefore, the Voronoi index shows a higher recognition rate at high temperatures.
In order to use the CNA to identify Z12 and Z16 Frank-Kasper clusters, the source code of CNA modifier in Ovito has to be modified. Due to lack of the access to Python API and source code of Ovito at the same time, the workflow of CNA-based LaCA involves some manual steps at the moment. To simplify the usage of LaCA, Voronoi-based LaCA is recommended.

!!The LaCA method is still in development, it can be extended to analyze other TCP phases!!

################  Software and library requirements ###############
ovito-pro 3.7.3 +

################ File structure ################
ovito_state_files                               # ovito state files to load saved modifiers
Case_studies                                    # Case studies including NEB calculation on synchroshear, nanoindentation, and slip transmission
ovito_Laves_crystal_analysis_Voronoibased_v1.1.py    # Script of the Voronoi-based LaCA

################  Perform Voronoi-based LaCA on the lammps dump file ################
ovitos ovito_Laves_crystal_analysis_Voronoibased.py -f <CONFIG> [-r <AtomicRadii>] [-t <ParticleType>] [-h <HELP>]

optional arguments:
  -h, --help            show this help message and exit
  -f FCONFIG,           Name of the configuration file
  -r RADII RADII,       Atomic radii of type 1 and 2 particles in Angstrom for computing the poly-disperse Voronoi tessellation; by default a mono-disperse Voronoi tessellation is computed; the effect of atomic radii of elements in Laves phases is limited on the Voronoi index.
  -t TYPE,              Particle type of binary AB2 Laves phase, e.g., 1 --> (A=AtomType1 and B=AtomType2), 2 --> (B=AtomType1 and A=AtomType2); by default particle type is NOT taken into account; for binary AB2 system it is recommended to consider the type, which can improve the recognition rate.
  
Output *.voro_csp.laves file

################  Visualize *.voro_csp.laves file ################
Load *.voro_csp.laves file in ovito and load the ovito state file LaCA_voronoibased_color_coding.ovito
Blue -- C14 
Green -- C15
Red -- Interface between C14 and C15
Orange -- Other Laves
White -- Other
Brown -- B(A) anti-site          
Yellow -- A(B) anti-site
!!Due to the limitation of Voronoi, the Voronoi index of surface atom has no meaning!!

################  Perform DXA ################ 
Load *.voro_csp.laves file in ovito and load the ovito state file LaCA_voronoibased_DXA.ovito
