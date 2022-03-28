function [err,XYZdic, map] = THD_3dmm_to_dicomm (Info, XYZmm)
%
%   [err, XYZdic, map] = THD_3dmm_to_dicomm (Info, XYZmm)
%
%Purpose:
%   Transform image coordinate (in mm) to dicom (RAI) in mm
%
%
%Input Parameters:
%   Info: Structure of header file, see BrikInfo
%   XYZmm: Nx3 matrix of 3dmm coordinates
%
%
%
%Output Parameters:
%   err : 0 No Problem
%       : 1  Problems
%   XYZdic: Nx3 matrix of dicom coordinates
%   map: 1x3 vector such that XYZdic = XYZmm(:,map)
%
%
%
%Key Terms:
%
%More Info :
%   Based on THD_3dmm_to_dicomm in thd_coords.c
%
%See Also:
%   AFNI_XYZcontinuous2Index
%   AFNI_Index2XYZcontinuous
%   AFNI_CoordChange
%
%     Author : Ziad Saad
%     Date : Wed Feb 18 13:44:06 EST 2004
%     SSCC/NIMH/ National Institutes of Health, Bethesda Maryland


%Define the function name for easy referencing
FuncName = 'THD_3dmm_to_dicomm';

%Debug Flag
DBG = 1;

%initailize return variables
err = 1;

for (i=1:1:3),
   switch (Info.ORIENT_SPECIFIC(i)),
      case 0
         map(1) = i;
      case 1
         map(1) = i;
      case 2
         map(2) = i;
      case 3
         map(2) = i;
      case 4
         map(3) = i;
      case 5
         map(3) = i;
      otherwise
         err = 1;
         fprintf(2,'Error %s:\nBad orientation code (%d).\n', FuncName, i);
         return;
   end
end

XYZdic = [];
if (nargin > 1),
   XYZdic = XYZmm;
   if (~isempty(XYZmm)),
      XYZdic(:,1) = XYZmm(:,map(1));
      XYZdic(:,2) = XYZmm(:,map(2));
      XYZdic(:,3) = XYZmm(:,map(3));
   end
end



err = 0;
return;

