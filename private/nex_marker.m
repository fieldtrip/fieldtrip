function [n, nm, nl, ts, names, m] = nex_marker(filename, varname)
% nex_marker(filename, varname): Read a marker variable from a .nex file
%
% [n, nm, nl, ts, names, m] = nex_marker(filename, varname)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%   varname - variable name
%
%           continuous (a/d) data come in fragments. Each fragment has a timestamp
%           and a number of a/d data points. The timestamp corresponds to
%           the time of recording of the first a/d value in this fragment.
%           All the data values stored in the vector d. 
% OUTPUT:
%   n - number of markers
%   nm - number of fields in each marker
%   nl - number of characters in each marker field
%   ts - array of marker timestamps (in seconds)
%   names - names of marker fields ([nm 64] character array)
%   m - character array of marker values [n nl nm]

% original from Plexon, download from http://www.plexoninc.com (8/4/02)
% modifications by Robert Oostenveld
%
% $Log: nex_marker.m,v $
% Revision 1.1  2009/01/14 09:24:45  roboos
% moved even more files from fileio to fileio/privtae, see previous log entry
%
% Revision 1.3  2008/09/30 07:47:04  roboos
% replaced all occurences of setstr() with char(), because setstr is deprecated by Matlab
%
% Revision 1.2  2008/07/24 11:57:57  roboos
% converted end of line into unix style
%
% Revision 1.1  2005/02/11 07:46:34  roboos
% downloaded from Plexon website, added support for reading nex files on Solaris and Mac OS X (fopen ieee-le), added log-section to the comments after the help
%

n = 0;
nm = 0;
nl = 0;
ts = 0;
m = 0;
names = 0;

if(nargin ~= 2)
   disp('2 input arguments are required')
   return
end

if(ischar(filename) == 0)
   disp('input arguments should be character arrays')
   return
end

if(ischar(varname) == 0)
   disp('input arguments should be character arrays')
   return
end

if(length(filename) == 0)
   [fname, pathname] = uigetfile('*.nex', 'Select a Nex file');
	filename = strcat(pathname, fname);
end

fid = fopen(filename, 'r', 'ieee-le');
if(fid == -1)
	disp('cannot open file');
   return
end

disp(strcat('file = ', filename));
magic = fread(fid, 1, 'int32');
version = fread(fid, 1, 'int32');
comment = fread(fid, 256, 'char');
freq = fread(fid, 1, 'double');
tbeg = fread(fid, 1, 'int32');
tend = fread(fid, 1, 'int32');
nvar = fread(fid, 1, 'int32');
fseek(fid, 260, 'cof');
name = zeros(1, 64);
found = 0;
for i=1:nvar
	type = fread(fid, 1, 'int32');
	var_version = fread(fid, 1, 'int32');
	name = fread(fid, [1 64], 'char');
	offset = fread(fid, 1, 'int32');
	n = fread(fid, 1, 'int32');
	dummy = fread(fid, 32, 'char');
	adfreq = fread(fid, 1, 'double');
	adtomv = fread(fid, 1, 'double');
	npw = fread(fid, 1, 'int32');
	nm = fread(fid, 1, 'int32');
	nl = fread(fid, 1, 'int32');
	dummy = fread(fid, 68, 'char');
	name = char(name);
	name = deblank(name);
	k = strcmp(name, deblank(varname));
	if(k == 1)
		if type ~= 6
			disp(sprintf('%s is not a marker variable', deblank(varname)));
			return;
		end
		found = 1;
		fseek(fid, offset, 'bof');
		ts = fread(fid, [1 n], 'int32');
		names = zeros(1,64);
		m = zeros(n, nl, nm);
		for j=1:nm
			names(j, :) = fread(fid, [1 64], 'char');
			for p = 1:n
				m(p, :, j) = fread(fid, [1 nl], 'char');
			end
		end
		break
	end
end

fclose(fid);

if found == 0
	disp('did not find variable in the file');
else
	names = char(names);
	m = char(m);
	ts = ts/freq;
	disp(strcat('number of markers = ', num2str(n)));
end
