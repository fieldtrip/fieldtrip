function len = strlen(str)
%  len = strlen(str)
% compute the # of characters in str (ignoring 0s at the end)


%
% strlen.m
%
% Original Author: Bruce Fischl
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

len = length(str) ;
for i=length(str):-1:1
	if (str(i) ~= 0)
		break ;
	end
	
	len = len-1;
end

		
