function [err,ITout,ColMapout] = ScaleToMap (ITvect,ColMap,Opt)
%
%   [err,Mout,ColMapout] = ScaleToMap (M,ColMap,Opt)
%
%Purpose:
%   This function scales a matrix (or vector) to make it fit a
%   colormap. It was originally deisgned to represent functional
%   data onto 3D and 2D meshes, but it works for a lot of other apps.
%
%
%Input Parameters:
%   M is a matrix (KxL) containing the data that needs scaling.
%     you can choose one value in M to be a "mask" so that
%     data points having the "mask" value will have a special colour.
%     (read on for more info)
%   ColMap is a Nx3 matrix that contains the colours that the data in M
%          should be mapped to. This colour map is applied to values in M
%          different from "mask" value
%   Opt is the options structure with the following fields
%       .DataRange (optional) : 1x2 vector. If ITvect contains data to be pseudocolored,
%               you might not want to have the highest and lowest values
%               to be mapped to the last and first colours in the colormap
%               If you do so, then for different ITvect with different values,
%               The colors will look different. So, to keep the color significance
%               constant across data sets you can specify [min max] in
%               .DataRange so that the first color is always min and the
%               last color is always max.If you do not specify a range,
%               The minimum and maximum values of ITvect will be used.
%               NOTE: The limits in DataRange are overrun by values in the
%               Data set that are lower than min or higher than max.
%       .MaskValue (optional) is the value in M considered as "mask"
%                   default is no MaskValue used
%       .MaskColor (optional) is the color to give to points in M
%            having a value of "mask". default is black. This
%            field can't be specified without .MaskValue
%       .Format (optional) specifies the format ot Mout
%            if set to 1 then Mout has the same size as M,
%            and it contains the scaled values of M
%            if set to 3 then Mout has the size of Kx3 if
%            M has a size of Kx1 and KxLx3 if M has a size of
%            KxL. The values in Mout now represent rgb values
%            interpolated from ColMap.
%            When Mout has a size of Kx3, the Mout(:,1) is the
%            red gun value, Mout(:,2) the green and Mout(:,3) the blue.
%            When Mout has a size of KxLx3, Mout(:,:,1) is the red
%            matrix, Mout(:,:,2) is the green etc...
%
%
%Output Parameters:
%   err : 0 No Problem
%       : 1 Mucho Problems
%   Mout : the scaled version of M (either in indexed colour or true colour)
%   ColMapOut : The colormap made out of Opt.MaskColor for the first color
%                and ColMap for the rest of the colours.
%
%
%More Info :
%   	DispIVSurf, MyDataStructures
%
%Opt.MaskValue = -1;
%Opt.MaskColor = [0 1 0];
%Opt.Format = 1;
%
%ColMap = [0 0 0; 0.25 0.25 0.25;0.5 0.5 0.5 ;0.75 0.75 0.75;1 1 1];
%
%M = [4 8 -1;8 0 2; 6 -1 -1];
%
%[err,Mout,ColMapout] = ScaleToMap (M,ColMap,Opt);
%
%figure
%image(Mout);colormap(ColMapout);colorbar%
%
%If you want the true color, try
%Opt.Format = 3;
%
%ColMap = [0 0 0; 0.25 0.25 0.25;0.5 0.5 0.5 ;0.75 0.75 0.75;1 1 1];
%
%M = [4 8 -1;8 0 2; 6 -1 -1];
%
%[err,Mout,ColMapout] = ScaleToMap (M,ColMap,Opt);
%
% Displaying the result is not obvious here, just look at
% Mout it's pretty obvious
%
%	see also, MakeColorMap
%
%     Author : Ziad Saad
%     Date : Wed Apr 08 15:16:05 CDT 1998


%Define the function name for easy referencing
FuncName = 'ScaleToMap';

%initailize return variables
err = 1;
ITout = [];

%stringent input checking
if (nargin == 2),
	Opt.MaskValue = [];
	Opt.MaskColour = [];
	Opt.Format = 1;
end

if (size(ColMap,2) ~= 3),
	err = ErrEval(FuncName,'Err_Color Map  must have a Nx3 dimention.');	return;
end

tmp = find(ColMap < 0 | ColMap >1);
if (~isempty(tmp)),
	err = ErrEval(FuncName,'Err_Color values are either negative or > 1.');	return
end

if (isfield(Opt,'MaskColor') & ~isfield(Opt,'MaskValue')),
	err = ErrEval(FuncName,'Err_Must specify Opt.MaskValue to use Opt.MaskColor.');	return
end

if (~isfield(Opt,'MaskValue')),
	Opt.MaskValue = [];
	Opt.MaskColor = [0 0 0];
elseif (~isfield(Opt,'MaskColor') | isempty(Opt.MaskColor)),
	Opt.MaskColor = [0 0 0];
end

if (size(Opt.MaskColor,2) ~= 3),
	err = ErrEval(FuncName,'Err_Opt.MaskColor  must have a 1x3 dimention.');	return;
end

if (~isfield(Opt,'Format') | isempty(Opt.Format)),
	Opt.Format = 1;
end

if (~isfield(Opt,'DataRange') | isempty(Opt.DataRange)),
	Opt.DataRange = [];
elseif (length(Opt.DataRange) ~= 2),
	err = ErrEval(FuncName,'Err_Size of Opt.DataRange must be 1x2.');	return;
end

if (isempty(ITvect)),
	err = ErrEval(FuncName,'Err_Empty ITvect !');	return
end

if (~isempty(Opt.DataRange) & ~isempty(Opt.MaskValue)),
	if (~isempty(find(Opt.DataRange == Opt.MaskValue))),
		err = ErrEval(FuncName,'Err_You cannot have a MaskValue equal to values in DataRange.');
		return
	end
end

%reshape input to a vector
[Nr,Nc] = size(ITvect);
ITvect = reshape (ITvect,Nr.*Nc,1);

%Add the min and max values of DataRange if specified for correct scaling
if (~isempty(Opt.DataRange)),
	igood = find (ITvect < Opt.DataRange(1));
	if (~isempty(igood))
		%fprintf(1,'fixing min range, %d elements\n', length(igood));
		ITvect(igood) = Opt.DataRange(1);
	end
	igood = find (ITvect > Opt.DataRange(2));
	if (~isempty(igood))
		%fprintf(1,'fixing max range, %d elements\n', length(igood));
		ITvect(igood) = Opt.DataRange(2);
	end
	
	ITvect = [ITvect ; Opt.DataRange(1) ; Opt.DataRange(2)];
end

%initialize
ITout = ITvect;
if (~isempty(Opt.MaskValue)),
	%find values in ITvect that are not masked and scale those to fit the colour map
	iNotMsk = find (ITvect ~= Opt.MaskValue);
	IToutdata = ITout(iNotMsk);
	ITout = ones(size(ITvect));
	ITout(iNotMsk) = zscale (IToutdata,length(ColMap)+1,2);
	
   %get ridd of DataRange values if present
	if (~isempty(Opt.DataRange)),
		ITout = ITout(1:length(ITvect)-2);
		%Take out the last two values of iNotMsk because they point
		%to the values added because of Opt.DataRange
		iNotMsk = iNotMsk(1:length(iNotMsk)-2);
	end

	
	if (Opt.Format == 1),	%bring output back to input size
		ITout = reshape(ITout,Nr,Nc);
	else	%need return true color values
		IToutdata = ITout(iNotMsk);
		IToutInterp = interp1([2:length(ColMap)+1]',ColMap,IToutdata);
		dumm = ones(Nr.*Nc,1);
		ITout = [Opt.MaskColor(1).*dumm ,...
	         	Opt.MaskColor(2).*dumm ,...
					Opt.MaskColor(3).*dumm ];

		ITout(iNotMsk,:) = IToutInterp;

		if (Nr == 1 | Nc == 1),	%input was a vector
			ITout = reshape(ITout,Nr.*Nc,3);
		else	%input was a matrix
			ITout = reshape(ITout,Nr,Nc,3);
		end
	end

	ColMapout = [Opt.MaskColor;ColMap];

else
	ITout = zscale (ITout,length(ColMap),1);
   %get ridd of DataRange values if present
	if (~isempty(Opt.DataRange)),
		ITout = ITout(1:length(ITvect)-2);
	end
	
	if (Opt.Format == 1),	%bring output back to input size
		ITout = reshape(ITout,Nr,Nc);
	else	%need return true color values
		ITout = interp1([1:length(ColMap)]',ColMap,ITout);
		if (Nr == 1 | Nc == 1),	%input was a vector
			ITout = reshape(ITout,Nr.*Nc,3);
		else	%input was a matrix
			ITout = reshape(ITout,Nr,Nc,3);
		end
	end
	
	ColMapout = ColMap;
end

err = 0;
return;

