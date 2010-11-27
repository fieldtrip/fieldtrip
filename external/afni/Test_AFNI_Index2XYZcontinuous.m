%script Test_AFNI_Index2XYZcontinuous
%
%
%
%Purpose:
%   
%   
%   
%Input:
%   
%   
%   
%Output:
%
%
%   
%   
%      
%Key Terms:
%   
%More Info :
%   
%   
%   
%
%     Author : Ziad Saad
%     Date : Tue Sep 5 21:48:26 PDT 2000
%     LBC/NIMH/ National Institutes of Health, Bethesda Maryland


%Debug Flag
DBG = 1;

%launch afni
%unix('afni &');

%Load a Brik Info
[err, Info] = BrikInfo('ARzs_CW_avvr+orig.HEAD');

%XYZ indices
	Indx = [24 36 8; 31 12 0];
	
%get coordinates in RAI
	[err,XYZmm] = AFNI_Index2XYZcontinuous (Indx, Info);
	
%goto these voxels in the AFNI window and check that
fprintf (1,'Voxel indices :');
	Indx
fprintf (1,'Voxel Coords [RAI]:');
	XYZmm
fprintf (1,'Voxel Coords [ILA]:');
	 [err, maplocation, mapsign, Mtrans] = AFNI_CoordChange ('RAI', 'ILA', XYZmm);
Mtrans

%or directly,
	[err,XYZmm] = AFNI_Index2XYZcontinuous (Indx, Info, 'ILA');
	fprintf (1,'Voxel Coords [ILA]:');
	XYZmm
