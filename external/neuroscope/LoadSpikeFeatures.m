function features = LoadSpikeFeatures(filename,rate)

%LoadSpikeFeatures - Load spike times, cluster IDs and features from disk.
%
%  USAGE
%
%    features = LoadSpikeFeatures(filename,rate)
%
%    filename            spike file name (either .clu, .res or .fet)
%    rate                sampling rate
%
%  OUTPUT
%
%    The output is a list of (timestamp,group,cluster,features...) t-uples.
%
%  SEE
%
%    See also GetSpikeFeatures.

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

[path,basename,extension] = fileparts(filename);
if isempty(path), path = '.'; end

electrodeGroup = str2num(extension(2:end));
[~,basename] = fileparts(basename);

% Load .clu file
filename = [path '/' basename '.clu.' int2str(electrodeGroup)];
if ~exist(filename),
	error(['File ''' filename ''' not found.']);
end
clu = load(filename);
clu = clu(2:end);

% Load .fet file
filename = [path '/' basename '.fet.' int2str(electrodeGroup)];
if ~exist(filename),
	error(['File ''' filename ''' not found.']);
end
file = fopen(filename,'r');
if file == -1,
	error(['Cannot open file ''' filename '''.']);
end
nFeatures = fscanf(file,'%d',1);
fet = fscanf(file,'%f',[nFeatures,inf])';
fclose(file);

features = [fet(:,end)/rate electrodeGroup*ones(size(clu)) clu fet(:,1:end-1)];

