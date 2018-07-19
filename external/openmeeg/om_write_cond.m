function om_write_cond(condfile,c,names)
%   OM_WRITE_COND
%       [] = OM_WRITE_COND(CONDFILE,C)
%
%   Write conductivity file for OpenMEEG
%
%   Authors: Alexandre Gramfort alexandre.gramfort@inria.fr

% Copyright (C) 2010-2017, OpenMEEG developers

ndomains = length(c);

if nargin < 3
    names = {};
    for k=1:ndomains
        names{k} = ['domain', num2str(k)];
    end
end

if ndomains ~= length(names)
    error('Number of conductivities is not equal to the number of domain names.');
end

cfid = fopen(condfile, 'w');
if cfid == -1
    error(['Failed to open file ''',condfile,'''.']);
end

fprintf(cfid,'# Properties Description 1.0 (Conductivities)\n\n');
fprintf(cfid,'air         0.0\n');

for k=1:ndomains
    fprintf(cfid, '%s       %f\n', names{k}, c(k));
end

fclose(cfid);

end %  function