**Helsinki BEM Framework LCISA solver for MEG/EEG**

Copyright Matti Stenroos 2007–2020, licensed under CC BY-NC-ND 4.0 license.
Contact matti.stenroos@aalto.fi

This Helsinki BEM Framework (HBF) package contains an application-oriented version of the linear collocation BEM solver formulated using isolated source approach (LCISA),
used as test solver in articles

*  Stenroos, M., Hunold, A., Haueisen, J., 2014. Comparison of three-shell and simplified volume conductor models in magnetoencephalography. Neuroimage 94, 337–348.
*  Stenroos, M., Nummenmaa, A., 2016. Incorporating and compensating cerebrospinal fluid in surface-based forward models of magneto- and electroencephalography. PLoS ONE 11(7):e0159595, 2016.

Please cite the 2014 paper,
if you use this code for basic MEG modeling (3-shell or 1-shell),
or the 2016 paper,
if you use this code for EEG or advanced MEG modeling (4-compartment).
If (when...) you wish to cite also my method papers that describe,
how the computations are done,
please see the references at the end of this document.
The 2007 paper descibes the basic collocation BEM,
and the 2012 paper the isolated source approach as used here.

The solver is aimed for MEG/EEG use,
but it works just as well or better for MCG/ECG.
There are more accurate BEM approaches (for example linear Galerkin BEM and symmetric BEM),
but taking into account the level of approximation in the head and thorax models,
this solver is adequate for typical experimental use,
when some attention is paid to the distance of source from nearest boundaries (see the papers).

The example scripts present all necessary file formats.
Regarding boundary meshes, it is essential to know that
1) meshes are ordered from the innermost to the outermost mesh, 
2) meshes must be closed,
3) triangle orientation is CCW (please run hbf_CheckTriangleOrientation to check this),
4) the resulting models have been verified with meshes of approx 2500 vertices per surface,
5) the recommended distance for dipole sources from the inner skull is half of triangle side length.
6) all physical measures are in SI units --- spatial variables in meters, dipole moments in ampere meters;
7) all sources must be inside the ISA surface --- in the case of three-shell MEG models, inside surface one.

Before using the solver with your own meshes,
please study the example case to see,
what kind of meshings and parameters to use. 
If you get surprising results (like essentially worse performance than that presented in the paper),
please check the above list and plot your model geometry.
If there is something strange, please drop me a line!
And, don't worry about the computational cost with the 2500-vertex meshes of the example model;
my desktop PC builds the whole model in less than 35 seconds,
so there really is no need to use coarser meshes.

**Example data description**

The example data, 'hbf_samplehead_3shell_wh' is to be downloaded from another repository,

https://github.com/MattiStenroos/hbf_sampledata
  
It contains meshed of the  head geometry of Subject 1 of the data set described in

Wakeman, D.G. & Henson, R.N. (2015). A multi-subject, multi-modal
human neuroimaging dataset. Sci. Data 2:150001 doi: 10.1038/sdata.2015.1,

available through 
ftp://ftp.mrc-cbu.cam.ac.uk/personal/rik.henson/wakemandg_hensonrn/

If you use the head geometry in a publication,
please cite the Wakeman--Henson publication and include the above
information.

The meshes were constructed with SPM, using scripts provided with the
data. Sensor geometries of a 306-channel Neuromag MEG system (as described in
MNE-C software, http://mne.tools/) and 256-channel ABC EEG-layout (as
described by Biosemi, http://www.biosemi.com/download/Cap_coords_all.xls,
were added manually by the author. The sensors and their coregistrations
do not correspond to the sensors used in the of the Wakeman–Henson
dataset.


**References**

* Stenroos, M., Nummenmaa, A., 2016. Incorporating and compensating cerebrospinal fluid in surface-based forward models of magneto- and electroencephalography. PLoS ONE 11(7):e0159595, 2016.
* Stenroos, M., Hunold, A., Haueisen, J., 2014. Comparison of three-shell and simplified volume conductor models in magnetoencephalography. Neuroimage 94, 337–348.
* Stenroos, M., Sarvas, J., 2012. Bioelectromagnetic forward problem: isolated source approach revis(it)ed. Phys Med Biol 57, 3517–3535.
* Stenroos, M., Nenonen, J., 2012. On the accuracy of collocation and Galerkin BEM in the EEG/MEG forward problem. Int J Bioelectromagnetism 14, 29–33.
* Stenroos, M., Mäntynen, V., Nenonen, J., 2007. A Matlab library for solving quasi-static volume conduction problems using the boundary element method. Comput Methods Programs Biomed 88, 256–263.

For BEM solution in a general piece-wise homogeneous model, see
* Stenroos M., 2016. Integral equations and boundary-element solution for static potential in a general piece-wise homogeneous volume conductor. Phys Med Biol 61:N606–N617.
