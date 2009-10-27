function [n, ts_left, ts_right] = nex_int(filename, varname)
% nex_int(filename, varname): Read interval variable from a .nex file
%
% [n, ts_left, ts_right] = nex_int(filename, varname)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%   varname - variable name
% OUTPUT:
%   n - number of intervals
%   ts_left - array of left ends of the intervals (in seconds)
%   ts_right - array of right ends of the intervals (in seconds)

% original from Plexon, download from http://www.plexoninc.com (8/4/02)
% modifications by Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

n = 0;
ts_left = 0;
ts_right = 0;

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
	name = char(name);
	name = deblank(name);
	k = strcmp(name, deblank(varname));
	if(k == 1)
		if type ~= 2
			disp(sprintf('%s is not an interval variable', deblank(varname)));
			return;
		end
		found = 1;
		fseek(fid, offset, 'bof');
		ts_left = fread(fid, [1 n], 'int32');
		ts_right = fread(fid, [1 n], 'int32');
		break
	end
	dummy = fread(fid, 128, 'char');
end

fclose(fid);

if found == 0
	disp('did not find variable in the file');
else
	ts_left = ts_left/freq;
	ts_right = ts_right/freq;
	disp(strcat('number of intervals = ', num2str(n)));
end
