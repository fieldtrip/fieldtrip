This directory contains a standard Boundary Element Method volume
conduction model of the head that can be used for EEG forward and
inverse computations. The geometry is based on the "colin27" template.

The construction of the volume conduction model is detailled in
  Oostenveld R, Stegeman DF, Praamstra P, van Oosterom A.
  Clin Neurophysiol. 2003 Jul;114(7):1194-202.
  Brain symmetry and topographic analysis of lateralized event-related potentials.
Please cite this reference if you use this standard BEM volume
conduction model in your analyses.

Accompanying electrode positions according to the 10-20, the 10-10
and the 10-5 standards (with different naming schemes) can be found
in the template/electrode directory.

--------------------------------------------------------------------
The following detailled description is from
http://imaging.mrc-cbu.cam.ac.uk/imaging/MniTalairach
--------------------------------------------------------------------
One of the MNI lab members, Colin Holmes, was scanned 27 times, and
the scans were coregistered and averaged to create a very high
detail MRI dataset of one brain. This average was also matched to
the MNI305, to create the image known as "colin27". colin27 is used
in the MNI brainweb simulator. SPM96 used colin27 as its standard
template. [...] SPM96 and later contains a 2mm resolution copy of
the same image, in the canonical directory of the SPM distribution.
In SPM96 this is called T1 in later distributions it is called
single_subj_T1.

