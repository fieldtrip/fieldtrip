function [err, XYZ] = AfniIndex2AfniXYZ (Indx, Nxx, Nyy)
%
%    [err, XYZ] = AfniIndex2AfniXYZ (Indx, Nxx, Nyy);
%
% Returns  the XYZ coordinates of an AFNI voxel with
%  an AFNI index of Indx
%
% Indx : Nx1 vector containing the AFNI indices
%   Indx must be an integer vector. If it is not, it's values
%   are rounded to the nearest integer
%
% Nxx, Nyy  : Number of pixels in the slice in the X and Y directions
%
% err = 0 No problem
% err = 1 input Matrix size problems
%
% XYZ : is the XYZ triplets matrix Nx3
%
%      Ziad Saad   Sun Mar 15 19:24:39 CST 1998

	[n1,m1] = size(Indx);	
	
	if (m1 ~= 1),
		fprintf (1,'\a\nError in AfniIndex2AfniXYZ : Bad size for Indx\n\n');
		err = 1;
		return;
	end
	
	Indx = round (Indx);
	
	NxxNyy = Nxx .* Nyy;

	Z = floor (Indx ./ NxxNyy);
	Y = floor ((Indx - Z .* NxxNyy) ./ Nxx);
	X = Indx - ( Y .* Nxx) - (Z .* NxxNyy) ;

	XYZ = [X Y Z];
	
err = 0;
return;
