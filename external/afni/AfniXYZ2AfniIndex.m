function [err, Indx] = AfniXYZ2AfniIndex (XYZ, Nxx, Nyy)
%
%    [err, Indx] = AfniXYZ2AfniIndex (XYZ, Nxx, Nyy)
%  
% Returns  the AFNI index of an AFNI voxel with coordinates XYZ   
%
% XYZ : is the XYZ triplets matrix Nx3
%  If the data in XYZ is not integers, it will
%  be rounded to the nearest integer.
%
% Nxx, Nyy  : Number of pixels in the slice in the X and Y directions
%
% Indx : Nx1 vector containing the AFNI indices
% err = 0 No problem
% err = 1 input Matrix size problems
%
%
%      Ziad Saad   Sun Mar 15 19:24:39 CST 1998 

	[n1,m1] = size(XYZ);	
	
	if (m1 ~= 3),
		fprintf (1,'\a\nError in AfniXYZ2AfniIndex : Bad size for XYZ\n\n');
		err = 1;
		return;
	end 
	
	XYZ = round (XYZ);
		
   Indx = XYZ(:,1) + Nxx.* XYZ(:,2) + Nxx .* Nyy .* XYZ(:,3);
	
err = 0;
return;
