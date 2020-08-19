# TEMSIM_modified
Modified version of EJ Kirkland's TEMSIM of computem package multislice simualtion software.
Please cite: Kirkland, Earl J. Advanced computing in electron microscopy. New York: Plenum Press, 1998.
This software produce images suitable for cryo-EM multislice simulation without gaussian noise.
# PS:
I mode this repository open because the release of https://doi.org/10.1016/j.ultramic.2020.113094 , and I think this package could do better.
I did many simulations, the result was that multiple elastic scattering effect is SIGNIFICANTLY WEAKER THAN THE EWALD SPHERE EFFECT. If we could not overcome the Ewald Sphere effect on small proteins, it will be meaningless to even consider multiple elastic scattering effect.
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

2. Run xyz2slic. This program produce a bunch of TIFF format complex multislice files. "Enter jmin jmax" is the start line of eulerfile.txt. e.g. "50 100". "Wavefunction size in pixels, Nx,Ny" is the same as boxsize. For "Slice thickness (in Angstroms)", Kirkland suggest ~100, but others suggest 10-20. From my point of view 20 for > 1MDa protein, 10 for smaller protein.

xyz2slic_with_random_center_oxygen_v2 provide adding random oxygen atoms to emulate glassy water. The density of LDA is 0.94g/cm^3, so the oxygen density is 0.835g/cm^3.

"xyz2slic_with_random_center_oxygen_v2 -h" for usage. You can also run this program without adding oxygen.

3. Run multislice-list-image_modified_v4 to produce particle image for cryo-EM.

First you need a file index. such as 'ls *.tiff >> index.txt'. multislice-list-image_modified_v4 -a 'pixel_size' -i index.txt -f 'first line' -l 'last line'

To convert the output TIFF to mrc, you can either use imod or proc2d of EMAN. A script below user proc2d to convert.

foreach file (*.tif)
rm -f temp
ls $file > temp
set out = `awk -F. '{print $1".mrc"}' temp`
proc2d $file $out
end
rm -f temp

4. Reconstruction. Although CTF with known defocus0 is applied to the image, because of the depth-of-field and multi-scattering effects, the measured defocus is not the given value. In general defocus = defocus0 + half_particle_size . I suggested to use relion and CTFFIND4 to measure every particle image.

To combine the result, use temsim/read_ctffind_star_write_lst.py

temsim/read_ctffind_star_write_lst.py "relion.star" "metadata_line+4/+5" "boxsize" "output_eman.lst"

To reconstruct to the 3D model, use j3dr of JSPR package (https://jiang.bio.purdue.edu/jspr/), or use the Block-based-recontruction/for_straight_forward_relion/jspr_refine_2_relion_nofilename_for_relion.py to convert a relion star file for relion_reconstruction.

5. If you like to compare the FSC of the pdb model, you can compile temsim/pdb2mrc_lorenzian_gaussian_v3.c by:

g++ temsim/pdb2mrc_lorenzian_gaussian_v3.c -o temsim/pdb2mrc_lorenzian_gaussian_v3 -I/home/user/EMAN/include /home/user/EMAN/lib/libEM.so , which require EMAN installed. It is really slow.

BTW, the lorenzian-gaussian approximation used by Kirkland is quite similar to gauss-ball approximation. At all frequency, the FSC dropped never below 0.5 . So you can simply use pdb2mrc from EMAN.

Remember to remove all the "ENDMDL" symbol on pdb, or pdb2mrc will only produce a single asymetrical unit rather than the whole structure. Also to remember to enable "center" in pdb2mrc, this option center the pdb by mass center, which is implemented at pdb2xyz_3angle_python_ver.py .

To calc FSC, use EMAN's proc3d reconstruct.mrc model.mrc apix=pixel_size fsc=fsc.fsc , or e2proc3d.py reconstruct.mrc fsc.fsc --apix pixel_size --calcfsc=model.mrc in EMAN2.

I provided some xyz files in xyz_file folder, you can test yourself.
