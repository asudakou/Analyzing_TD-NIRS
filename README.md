# Analyzing_TD-NIRS
- The contents are Data and MATLAB scripts used for Publication in 2023. Please see file "Content.png"

> A. Sudakou, H. Wabnitz, A. Liemert, M. Wolf, and A. Liebert.  "Two-layered blood-lipid phantom and method to determine absorption and oxygenation employing changes in moments of DTOFs".  DOI  : https://doi.org/10.1364/BOE.492168 

**The 1st LMA** does multi-layered curve-fitting of the DTOF for determining the absolute optical properties (and/or the thicknesses of layers).

**The 2nd LMA** uses the changes in moments to determine the changes in optical properties (and/or the thicknesses of layers).

**The data** includes 3 experiments involving blood and 2 experiments involving ink, in a new two-layered phantom, measured with a multi-wavelength TD-NIRS system.

To use the codes and data, download all into one folder, add this folder to MATLAB path, and lastly, in files 'Ink_main.m' and 'Blood_main.m' set variable 'folder_calc' to the path of a folder where the calculated data will be stored (some calculations can take a whole day). 

For backup, a copy of these files was saved in Google folder: https://drive.google.com/drive/folders/14ckJFZtuB_6qBHPFUp8uAxJ785Oq4XRf?usp=sharing 

- The moments-based method:

> A. Liebert, H. Wabnitz, J. Steinbrink, H. Obrig, M. Möller, R. Macdonald, A. Villringer, and H. Rinneberg.  "Time-Resolved Multidistance Near-Infrared Spectroscopy of the Adult Head: Intracerebral and Extracerebral Absorption Changes from Moments of Distribution of Times of Flight of Photons," Appl. Opt. 43, 3037-3047 (2004).  DOI:  https://doi.org/10.1364/AO.43.003037

- Generating DTOFs:

> A. Liemert and A. Kienle.  "Light diffusion in N-layered turbid media: frequency and time domains," Journal of biomedical optics 15, 025002 (2010).  DOI:  https://doi.org/10.1117/1.3368682

> A. Liemert and A. Kienle.  "Application of the Laplace transform in time-domain optical spectroscopy and imaging," Journal of biomedical optics 20, 110502 (2015).  DOI:  https://doi.org/10.1117/1.JBO.20.11.110502

- If you use the code or the data, or find it helpful, please cite the above publication(s)
- The code and data can be used and modified, provided the original source is clearly acknowledged, following the GNU General Public License v3.0
