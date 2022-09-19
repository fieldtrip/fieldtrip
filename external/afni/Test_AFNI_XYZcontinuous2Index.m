%script Test_AFNI_XYZcontinuous2Index
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
%     Date : Thu Sep 7 16:50:58 PDT 2000
%     LBC/NIMH/ National Institutes of Health, Bethesda Maryland


%Debug Flag
DBG = 1;

%Load a Brik Info
[err, Info] = BrikInfo('ARzsspgrax+orig.BRIK');

XYZmmRAI = [ -16.8750  -28.1250  -34.0000; 73.1250   -1.8750  -66.0000];
XYZmmILA = [-34.0000   16.8750  -28.1250;  -66.0000  -73.1250   -1.8750];

%change to Index
	[err,Indx] = AFNI_XYZcontinuous2Index (XYZmmRAI, Info);
	fprintf(1,'From RAI, to 3D:\n');
	Indx
%change to Index 1D
	[err,Indx] = AFNI_XYZcontinuous2Index (XYZmmRAI, Info, '', 1);
	fprintf(1,'From RAI, to 1D:\n');
	Indx

%change from ILA	
	[err,Indx] = AFNI_XYZcontinuous2Index (XYZmmILA, Info, 'ILA',3);
	fprintf(1,'From ILA, to 3D:\n');
	Indx

