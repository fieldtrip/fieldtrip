function [err,XYZmm,map] = THD_dicomm_to_3dmm (Info, XYZdic)
%
%   [err,XYZmm,map] = THD_dicomm_to_3dmm (Info, [XYZdic])
%
%Purpose:
%   
%   Change from image-oriented mm coordinates to AFNI dicom (RAI) coordinates
%   
%Input Parameters:
%   Info: Structure of header file, see BrikInfo
%   XYZdic: Nx3 matrix of dicom coordinates
%   
%Output Parameters:
%   err : 0 No Problem
%       : 1  Problems
%   XYZmm: Nx3 matrix of 3dmm coordinates
%   map: 1x3 vector such that XYZmm = XYZdic(:,map)
%      
%Key Terms:
%   
%More Info :
%   Based on THD_dicomm_to_3dmm in thd_coords.c
%   
%   
%
%     Author : Ziad Saad
%     Date : Wed Feb 18 12:01:08 EST 2004
%     SSCC/NIMH/ National Institutes of Health, Bethesda Maryland


%Define the function name for easy referencing
FuncName = 'THD_dicomm_to_3dmm';

%Debug Flag
DBG = 1;

%initailize return variables
err = 1;

for (i=1:1:3),
   switch (Info.ORIENT_SPECIFIC(i)),
      case 0
         map(i) = 1;
      case 1
         map(i) = 1;
      case 2
         map(i) = 2;
      case 3
         map(i) = 2;
      case 4
         map(i) = 3;
      case 5
         map(i) = 3;
      otherwise
         err = 1;
         fprintf(2,'Error %s:\nBad orientation code (%d).\n', FuncName, i);
         return;
   end
end

XYZmm = [];
if (nargin > 1),
   XYZmm = XYZdic;
   if (~isempty(XYZdic)) XYZmm = XYZdic(:,map); end    
end

err = 0;
return;

