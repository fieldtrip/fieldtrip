function [lab, pos] = read_besa_sfp(filename, uniqueonly)

% READ_BESA_SFP reads a besa style electrode location file.
%
% Use as:
%   [lab, pos] = read_besa_sfp(filename, uniqueonly)
%
% Input arguments:
%   filename   = the file name
%   uniqueonly = flag to determine behaviour, to return the positions of the
%                unique labels only (default behaviour: uniqueonly=1), or
%                also return double occurrences, which may be useful when
%                headshape information is represented in the file (as is
%                done in SPM)

if nargin==1
	uniqueonly = 1;
end

fid = fopen(filename);
tmp = textscan(fid, ' %[^ \t]%n%n%n');
fclose(fid);
		
lab = tmp{1};
pos = [tmp{2:4}];

if uniqueonly
	[ulab,ix,iy] = unique(lab);
  n = zeros(max(iy),1);
	for k = 1:max(iy)
		n(k) = sum(iy==k);
	end
	sel = iy(ismember(iy,find(n==1))); % selects the labels that occur once
	tmp = ulab(sel);
	
	[i1,i2] = match_str(lab, tmp);
	lab = lab(i1);
	pos = pos(i1,:);
end		
 
