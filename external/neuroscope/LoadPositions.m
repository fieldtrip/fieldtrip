function positions = LoadPositions(filename,rate)

%LoadPositions - Read positions from file.
%
%  USAGE
%
%    positions = LoadPositions(filename,rate)
%
%    filename            position file name
%    rate                optional sampling rate (if timestamps are not
%                        provided in the file)

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

positions = [];
if ~exist(filename),
	error(['File ''' filename ''' not found.']);
end

[~,~,extension] = fileparts(filename);
if strcmp(extension,'.whl') || strcmp(extension,'.pos'),
	positions = load(filename);
	t = (0:length(positions)-1)'/rate;
	positions = [t positions];
elseif strcmp(extension,'.mqa'),
	positions = load(filename);
end


