function om_write_geom(geomfile,bndfile,names)
%   OM_WRITE_COND
%       [] = OM_WRITE_GEOM(GEOMFILE,BNDFILE)
%       [] = OM_WRITE_GEOM(GEOMFILE,BNDFILE,NAMES)
%
%   Write geometry file for OpenMEEG
%
%   Authors: Alexandre Gramfort alexandre.gramfort@inria.fr
%            Paul Czienskowski

% Copyright (C) 2010-2017, OpenMEEG developers

numbnd = length(bndfile);

if nargin < 3
    names = {};
    for k=1:numbnd
        names{k} = ['domain', num2str(k)];
    end
end

if numbnd ~= length(names)
    error('Number of boundary files is not equal to the number of domain names.');
end

gfid = fopen(geomfile, 'w');

if gfid == -1
    error(['Failed to open file ''',geomfile,'''.']);
end

fprintf(gfid,'# Domain Description 1.0\n');
fprintf(gfid,'                        \n');
fprintf(gfid,'Interfaces %d Mesh      \n', numbnd);
fprintf(gfid,'                        \n');

for i=1:numbnd
    fprintf(gfid,'%s                  \n', bndfile{i});
end

fprintf(gfid,'                        \n');
fprintf(gfid,'Domains %d              \n', numbnd+1);
fprintf(gfid,'                        \n');

fprintf(gfid,'Domain air %d           \n', 1);
for i=1:numbnd
    if i < numbnd
        fprintf(gfid,'Domain %s %d %d\n', names{i}, i+1, -i);
    else
        fprintf(gfid,'Domain %s %d   \n', names{i}, -i);
    end
end

fclose(gfid);

end %  function