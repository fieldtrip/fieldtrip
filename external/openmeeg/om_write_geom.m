function om_write_geom(geomfile,bndfile)
%   OM_WRITE_COND
%       [] = OM_WRITE_GEOM(GEOMFILE,BNDFILE)
%
%   Write geometry file for OpenMEEG
%
%   Created by Alexandre Gramfort on 2009-03-30.
%   Copyright (c)  Alexandre Gramfort. All rights reserved.

if length(bndfile) == 4
    gfid = fopen(geomfile, 'w');
    fprintf(gfid,'# Domain Description 1.0\n');
    fprintf(gfid,'                        \n');
    fprintf(gfid,'Interfaces 4 Mesh       \n');
    fprintf(gfid,'                        \n');
    fprintf(gfid,'%s                      \n',bndfile{2});
    fprintf(gfid,'%s                      \n',bndfile{3});
    fprintf(gfid,'%s                      \n',bndfile{4});
    fprintf(gfid,'%s                      \n',bndfile{1});
    fprintf(gfid,'                        \n');
    fprintf(gfid,'Domains 5               \n');
    fprintf(gfid,'                        \n');
    fprintf(gfid,'Domain Skin 1 -4        \n');
    fprintf(gfid,'Domain Skull 2 -1       \n');
    fprintf(gfid,'Domain CSF 3 -2         \n');
    fprintf(gfid,'Domain Brain -3         \n');
    fprintf(gfid,'Domain Air 4            \n');
    fclose(gfid);
elseif length(bndfile) == 3
    gfid = fopen(geomfile, 'w');
    fprintf(gfid,'# Domain Description 1.0\n');
    fprintf(gfid,'                        \n');
    fprintf(gfid,'Interfaces 3 Mesh       \n');
    fprintf(gfid,'                        \n');
    fprintf(gfid,'%s                      \n',bndfile{2});
    fprintf(gfid,'%s                      \n',bndfile{3});
    fprintf(gfid,'%s                      \n',bndfile{1});
    fprintf(gfid,'                        \n');
    fprintf(gfid,'Domains 4               \n');
    fprintf(gfid,'                        \n');
    fprintf(gfid,'Domain Skin 1 -3        \n');
    fprintf(gfid,'Domain Brain -2         \n');
    fprintf(gfid,'Domain Air 3            \n');
    fprintf(gfid,'Domain Skull 2 -1       \n');
    fclose(gfid);
elseif length(bndfile) == 2
    gfid = fopen(geomfile, 'w');
    fprintf(gfid,'# Domain Description 1.0\n');
    fprintf(gfid,'                        \n');
    fprintf(gfid,'Interfaces 2 Mesh       \n');
    fprintf(gfid,'                        \n');
    fprintf(gfid,'%s                      \n',bndfile{2});
    fprintf(gfid,'%s                      \n',bndfile{1});
    fprintf(gfid,'                        \n');
    fprintf(gfid,'Domains 3               \n');
    fprintf(gfid,'                        \n');
    fprintf(gfid,'Domain Skin 1 -2        \n');
    fprintf(gfid,'Domain Air 2            \n');
    fprintf(gfid,'Domain Skull -1       \n');
    fclose(gfid);
elseif length(bndfile) == 1
    gfid = fopen(geomfile, 'w');
    fprintf(gfid,'# Domain Description 1.0\n');
    fprintf(gfid,'                        \n');
    fprintf(gfid,'Interfaces 1 Mesh       \n');
    fprintf(gfid,'                        \n');
    fprintf(gfid,'%s                      \n',bndfile{1});
    fprintf(gfid,'                        \n');
    fprintf(gfid,'Domains 2               \n');
    fprintf(gfid,'                        \n');
    fprintf(gfid,'Domain Air 1            \n');
    fprintf(gfid,'Domain Head -1          \n');
    fclose(gfid);
else
    error('OpenMEEG only supports geometry with 4 or less interfaces')
end

end %  function