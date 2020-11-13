function prov = hcp_provenance

% HCP_PROVENANCE returns a structure with provenance information

% Copyright (C) 2011-2014 by the Human Connectome Project, WU-Minn Consortium (1U54MH091657)
%
% This file is part of megconnectome.
%
% megconnectome is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% megconnectome is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with megconnectome.  If not, see <http://www.gnu.org/licenses/>.

prov.version.matlab = ver('MATLAB');

prov.version.megconnectome = ver('megconnectome');
if isempty(prov.version.megconnectome)
  prov.version.megconnectome = 'unknown'; % otherwise struct2xml fails
end

prov.version.fieldtrip = ver('fieldtrip');
if isempty(prov.version.fieldtrip)
  prov.version.fieldtrip = 'unknown'; % otherwise struct2xml fails
end

if isdeployed
  prov.compiled   = 'true';
else
  prov.compiled   = 'false';
end

% it is not desired to keep the date and time of execution of the analysis pipeline
% see https://github.com/Washington-University/megconnectome/issues/126
% prov.time           = datestr(now);

prov.username       = getusername;
prov.hostname       = gethostname;
prov.architecture   = computer('arch');
if isdeployed
  % the buildtimestamp function is created on the fly by the compile script
  prov.buildtimestamp = buildtimestamp;
else
  prov.buildtimestamp = [];
end
prov.pwd            = pwd;
prov.matlabstack    = printstack(dbstack); % using private helper function
prov.Attributes.xmlns_colon_xsi = 'http://www.w3.org/2001XMLSchema-instance';
prov.Attributes.xsi_colon_noNamespaceSchemaLocation = 'megconnectome.xsd';

function str = printstack(s)
str = sprintf('\n');
% skip the first, it is always hcp_provenance and not interesting
for i=2:length(s)
  str = [str sprintf('In %s at %d\n', s(i).name, s(i).line)];
end
