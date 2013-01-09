Anatomical Automatic Labeling (AAL) For SPM8
Version v1
Date 14/06/2010
E-mail : aalgin@cyceron.fr.

Anatomical Automatic Labeling (AAL) is a package for the anatomical labeling of 
functional brain mapping experiments. It is an in-house package made by 
Neurofonctional Imaging Group (GIN, UMR6232, CYCERON, Caen ,France), which is 
available to the scientific community as a copyright freeware under the terms of 
the GNU General Public Licence.

I      Goal
II     Reference
III    Distribution
IV     How to install the software
V      How to use the software
VI     Bibliography
VII    Contact
VIII   History

I Goal

This project has been initiated in the nineties with the construction of a set 
of rules to be used for the anatomical parcellation of the brain according to 
major sulci and gyri. We applied this set of rules to built an anatomical 
parcellation of the spatially normalized single subject high resolution T1 
volume provided by the Montreal Neurological Institute (MNI) (Collins et al., 
1998). The MNI single subject main sulci were first delineated and further used 
as landmarks for the 3D definition of 45 Anatomical Volumes Of Interest (AVOI) 
in each hemisphere. This procedure was performed using a dedicated software 
which allowed a 3D following of the sulci course on the edited brain. Regions of 
interest were then drawn manually using the same software every 2 mm on the 
axial slices of the high resolution MNI single subject. The 90 AVOI were 
reconstructed and assigned a label.
Using this parcellation method, three procedures to perform the automated 
anatomical labelling of functional studies are proposed :
1) labelling of an extremum defined by a set of coordinates,
2) percentage of voxels belonging to each of the AVOI intersected by a sphere 
centered by an extremum,
3) percentage of voxels belonging to each of the AVOI intersected by an 
activated cluster.

An interface with the Statistical Parametric Mapping (SPM8) package (Friston et 
al., 1995) is provided as a freeware to researchers of the neuroimaging 
community. We believe that this tool is an improvement for the macroscopical 
labelling of activated area as compared to labelling assessed using the 
Talairach atlas brain in which deformations are well known. However, this tool 
does not alleviate the need for more sophisticated labelling strategies based on 
anatomical or cytoarchitectonic probabilistic maps.

II Reference

If you need to reference this work, please used the following reference :

Tzourio-Mazoyer N, Landeau B, Papathanassiou D, Crivello F, Etard O, Delcroix N, 
et al. Automated anatomical labelling of activations in spm using a macroscopic 
anatomical parcellation of the MNI MRI single subject brain. Neuroimage 2002; 
15: 273-289.

III Distribution

The distribution includes a readme.txt (this file), 4 parcellation definition 
files (the ROI_MNI_V4* files), and 9 program file (*.m files).

III.1 ROI_MNI_V4.img (flat short integer image) and ROI_MNI_V4.hdr (header file in 
ANALYZE-7 format with slight customizations to the header as described in SPM 
help/spm_format.man and in Neurological Orientation (R is R)). Each anatomical region is 
associated a gray level (see III.2). 

Note that aal programs manage the orientation of your own analyze image (see defaults.analyze.flip variable).


III.2 ROI_MNI_V4_list.mat : Matlab format file giving the correspondence between the 
anatomical region name and the gray level information.

Anatomical-region-name(*) Gray-level
Precentral_L 2001
Precentral_R 2002
Frontal_Sup_L 2101
Frontal_Sup_R 2102
Frontal_Sup_Orb_L 2111
Frontal_Sup_Orb_R 2112
Frontal_Mid_L 2201
Frontal_Mid_R 2202
Frontal_Mid_Orb_L 2211
Frontal_Mid_Orb_R 2212
Frontal_Inf_Oper_L 2301
Frontal_Inf_Oper_R 2302
Frontal_Inf_Tri_L 2311
Frontal_Inf_Tri_R 2312
Frontal_Inf_Orb_L 2321
Frontal_Inf_Orb_R 2322
Rolandic_Oper_L 2331
Rolandic_Oper_R 2332
Supp_Motor_Area_L 2401
Supp_Motor_Area_R 2402
Olfactory_L 2501
Olfactory_R 2502
Frontal_Sup_Medial_L 2601
Frontal_Sup_Medial_R 2602
Frontal_Med_Orb_L 2611
Frontal_Med_Orb_R 2612
Rectus_L 2701
Rectus_R 2702
Insula_L 3001
Insula_R 3002
Cingulum_Ant_L 4001
Cingulum_Ant_R 4002
Cingulum_Mid_L 4011
Cingulum_Mid_R 4012
Cingulum_Post_L 4021
Cingulum_Post_R 4022
Hippocampus_L 4101
Hippocampus_R 4102
ParaHippocampal_L 4111
ParaHippocampal_R 4112
Amygdala_L 4201
Amygdala_R 4202
Calcarine_L 5001
Calcarine_R 5002
Cuneus_L 5011
Cuneus_R 5012
Lingual_L 5021
Lingual_R 5022
Occipital_Sup_L 5101
Occipital_Sup_R 5102
Occipital_Mid_L 5201
Occipital_Mid_R 5202
Occipital_Inf_L 5301
Occipital_Inf_R 5302
Fusiform_L 5401
Fusiform_R 5402
Postcentral_L 6001
Postcentral_R 6002
Parietal_Sup_L 6101
Parietal_Sup_R 6102
Parietal_Inf_L 6201
Parietal_Inf_R 6202
SupraMarginal_L 6211
SupraMarginal_R 6212
Angular_L 6221
Angular_R 6222
Precuneus_L 6301
Precuneus_R 6302
Paracentral_Lobule_L 6401
Paracentral_Lobule_R 6402
Caudate_L 7001
Caudate_R 7002
Putamen_L 7011
Putamen_R 7012
Pallidum_L 7021
Pallidum_R 7022
Thalamus_L 7101
Thalamus_R 7102
Heschl_L 8101
Heschl_R 8102
Temporal_Sup_L 8111
Temporal_Sup_R 8112
Temporal_Pole_Sup_L 8121
Temporal_Pole_Sup_R 8122
Temporal_Mid_L 8201
Temporal_Mid_R 8202
Temporal_Pole_Mid_L 8211
Temporal_Pole_Mid_R 8212
Temporal_Inf_L 8301
Temporal_Inf_R 8302
Cerebelum_Crus1_L 9001
Cerebelum_Crus1_R 9002
Cerebelum_Crus2_L 9011
Cerebelum_Crus2_R 9012
Cerebelum_3_L 9021
Cerebelum_3_R 9022
Cerebelum_4_5_L 9031
Cerebelum_4_5_R 9032
Cerebelum_6_L 9041
Cerebelum_6_R 9042
Cerebelum_7b_L 9051
Cerebelum_7b_R 9052
Cerebelum_8_L 9061
Cerebelum_8_R 9062
Cerebelum_9_L 9071
Cerebelum_9_R 9072
Cerebelum_10_L 9081
Cerebelum_10_R 9082
Vermis_1_2 9100
Vermis_3 9110
Vermis_4_5 9120
Vermis_6 9130
Vermis_7 9140
Vermis_8 9150
Vermis_9 9160
Vermis_10 9170

(*) Note that the cerebral AVOI definitions are fully described in the paper by 
Tzourio-Mazoyer et al. (Tzourio-Mazoyer et al., 2002). The cerebellar AVOI 
definitions are based on the cerebellum parcellation proposed by Schmahmann et 
al. (Schmahmann et al., 1999).

III.3 ROI_MNI_V4_Border.mat : Matlab format file listing the border of each region. 
This file is used in the automatic labeling procedure.

III.4 the matlab program files :

aal.m
gin_clusters.m
gin_clusters_plabels.m
gin_det_dlabels.m
gin_det_plabels.m
gin_dlabels.m
gin_list_dlabels.m
gin_list_plabels.m
gin_rclusters.m

IV How to install the software.

Note that aal procedures have been tested only on unix machines.

IV.1 Copy the archive to the chosen location
In this example we choose to install the software directly at the same location 
than the spm distribution (/usr/local/soft/spm8/toolbox)

unix> cp aal_for_spm8.tar.Z /usr/local/soft/spm8/toolbox
unix> cd /usr/local/soft/spm8/toolbox

IV.2 Expand the archive will create an anat_aal_vb1 directory
unix> uncompress aal_for_spm8.tar.Z
unix> tar xvf aal_for_spm8.tar

IV.3 Add this directory to your  matlab path
unix > setenv MATLABPATH ${MATLABPATH}:/usr/local/soft/spm8/toolbox/aal

V How to use the software.

V.1 Do a regular statistical analysis using spm8

V.2 Launch matlab
unix > matlab

V.3 launch aal
>> aal

V.4 Choose a labeling procedure. The 3 choices are explained and documented in 
the Neuroimage paper (Tzourio-Mazoyer et al., 2002):

Local maxima labeling
Extended local maxima labeling
Cluster labeling

V.5 Select the desired contrast, mask, probability and extent threshold like in 
the regular spm_result

V.6 For "Extended local maxima labeling" input the local maxima radius of the 
sphere in millimeters (default 10 mm).

V.7 Select the anatomical parcellation database :
In /usr/local/soft/spm8/toolbox/aal
The file ROI_MNI_V4.nii

V.8 Results

V.8.1 Local maxima labeling
For each local maxima :
-coordinates in mm x,y,z
-anatomical label (see below)
-distance in millimeter to this region. If the local maxima is inside a region 
this distance is null (0.00). If the local maxima is outside the parcellation 
the nearest region name is displayed in the previous column and the shortest 
distance from the local maxima to this region is listed (exp : 2.30 mm)
-anatomical label of the local maxima to the second nearest region
-shortest distance of the local maxima to the second nearest region
-anatomical label of the local maxima to the third nearest region
-shortest distance of the local maxima to the third nearest region

V.8.2 Extended local maxima labeling
Each local maxima is supposed to be a 10 mm (if the default is used) spherical 
region. The intersection of this volume and the AVOI is computed and the result 
sorted in a descending order according the percentage of overlap (exp : a result 
of Postcentral_L 100 % indicates that the 10mm radius region surrouding the 
local maxima is fully included in the Postcentral_L region)
For each local maxima :
-coordinates in mm x,y,z
-list of anatomical label and percentage of overlap. Percentage less than 1% are 
not listed. If part of the region is outside the parcellation the anatomical 
label will list "OUTSIDE".

V.8.3 Cluster labeling
The intersection of each cluster and the AVOI is computed and the result sorted 
in a descending order according the percentage of overlap.
For each local maxima :
-coordinates in mm x,y,z of the most significative local maxima of the cluster
-list of anatomical label and percentage of overlap. Percentage less than 1% are 
not listed. If part of the region is outside the parcellation the anatomical 
label will list "OUTSIDE".
    
Example : a result of
-44 -22 - 56           	Postcentral_L     55.00
          		Precentral_L     31.00
          		OUTSIDE         10.00
          		Parietal_Sup_L      5.00
indicates that :

55% of the cluster volume is included in the Postcentral_L region
31% of the cluster volume is included in the Precentral_L region
10% of the cluster volume is outside the parcellation
5% of the cluster volume is included in the Parietal_Sup_L region

VI Bibliography

Collins DL, Zijdenbos AP, Kollokian V, Sled JG, Kabani NJ, Holmes CJ, et al.
Design and construction of a realistic digital brain phantom. IEEE Transactions 
on Medical Imaging 1998; 17: 463-68.

Friston KJ, Holmes AP, Worsley KJ, Poline JP, Frith CD, Frackowiak RSJ. 
Statistical parametric maps in functional imaging: A general linear approach. 
Human Brain Mapping 1995; 2: 189-210.

Schmahmann JD, Doyon J, McDonald D, Holmes C, Lavoie K, Hurwitz AS, et al. 
Three-dimensional MRI atlas of the human cerebellum in proportional stereotaxic 
space. Neuroimage 1999; 10: 233-260.

Tzourio-Mazoyer N, Landeau B, Papathanassiou D, Crivello F, Etard O, Delcroix N, 
et al. Automated anatomical labelling of activations in spm using a macroscopic 
anatomical parcellation of the MNI MRI single subject brain. Neuroimage 2002; 
15: 273-289.


VII Contact
Any comments could be send to aalgin@cyceron.fr. 

VIII History

Version v1 Date 14/06/2010
aal version for SPM8










