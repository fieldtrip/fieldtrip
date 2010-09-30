function om_write_geom(geomfile,bndfile,names)
%   OM_WRITE_COND
%       [] = OM_WRITE_GEOM(GEOMFILE,BNDFILE)
%       [] = OM_WRITE_GEOM(GEOMFILE,BNDFILE,NAMES)
%
%   Write geometry file for OpenMEEG
%
%   Authors: Alexandre Gramfort alexandre.gramfort@inria.fr
%            Paul Czienskowski

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
    ind = mod(i,numbnd)+1;
    fprintf(gfid,'%s                  \n', bndfile{ind});
end

fprintf(gfid,'                        \n');
fprintf(gfid,'Domains %d              \n', numbnd+1);
fprintf(gfid,'                        \n');

for i=1:numbnd
    ind = mod(i-2,numbnd) + 1;
    if i < numbnd
        fprintf(gfid,'Domain %s %d -%d\n', names{i}, i, ind);
    else
        fprintf(gfid,'Domain %s -%d   \n', names{i}, ind);
    end
end

fprintf(gfid,'Domain Air %d           \n', numbnd);
fclose(gfid);

end %  function