# TEMSIM_modified
Modified version of EJ Kirkland's TEMSIM of computem package multislice simualtion software.
Please cite: Kirkland, Earl J. Advanced computing in electron microscopy. New York: Plenum Press, 1998.
This software produce images suitable for cryo-EM multislice simulation without gaussian noise.
# General compilatlion.
You need sfftw-dev. Just make.
1. You need convert pdb file to a xyz file. Use pdb2xyz_3angle. Complie the program by:
g++ -O3 -o pdb2xyz_3angle pdb2xyz_3angle.cpp slicelib.o floatTIFF.o cfpix.o -lfftw3f_threads -lfftw3f
2. For xyz2slic or xyz2slic_with_random_center_oxygen, use this command below.
g++ -O3 -o xyz2slic xyz2slic.cpp autoslic.cpp slicelib.o floatTIFF.o cfpix.o -lfftw3f_threads -lfftw3f
3. For generating final cryo-EM simulation particle images, use multislice-list-image_modified_v4:
g++ -O3 -o multislice-list-image_modified_v4 multislice-list-image_modified_v4.cpp slicelib.o floatTIFF.o cfpix.o -lfftw3f_threads -lfftw3f
# Usage
1. Convert PDB to a .xyz file. You can either use pdb2xyz-3angle or temsim/pdb2xyz_3angle_python_ver.py. The latter is suggested.

Usage: pdb2xyz_3angle_python_ver.py file.pdb 'pixel_size_you_need' 'boxsize_of_the_particles' 'euler.alt' 'euler.az' 'euler.phi' output.xyz 1 

The system use EMAN's zxz euler system. If your pdb is Cn or Dn symmetry, euler.alt=az=phi=0.0. If your pdb is relion I3 symmetry, alt=0.0, az=180.0, phi=0.0. Others are unknown. Incorrect rotation of xyz file will lead to a reconstrution failure.
