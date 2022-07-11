function om_write_geom(geomfile,bndfile,names,version)
%   OM_WRITE_COND
%       [] = OM_WRITE_GEOM(GEOMFILE,BNDFILE)
%       [] = OM_WRITE_GEOM(GEOMFILE,BNDFILE,NAMES)
%       [] = OM_WRITE_GEOM(GEOMFILE,BNDFILE,NAMES,VERSION)
%
%   Write geometry file for OpenMEEG, the default Domain Description has
%   been upgraded to '1.1'.
%
%   Authors: Alexandre Gramfort alexandre.gramfort@inria.fr
%            Paul Czienskowski

% Copyright (C) 2010-2017, OpenMEEG developers
% Copyright (C) 2022, Jan-Mathijs Schoffelen
numbnd = length(bndfile);

if nargin < 4
  version = '1.1';
end

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

switch version
  case '1.0'
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
  case '1.1'
    fprintf(gfid,'# Domain Description 1.1\n');
    fprintf(gfid,'                        \n');
    fprintf(gfid,'Interfaces %d           \n', numbnd);
    fprintf(gfid,'                        \n');

    for i=1:numbnd
      fprintf(gfid,'Interface: %s         \n', bndfile{i});
    end

    fprintf(gfid,'                        \n');
    fprintf(gfid,'Domains %d              \n', numbnd+1);
    fprintf(gfid,'                        \n');

    fprintf(gfid,'Domain air: %d          \n', 1);
    for i=1:numbnd
      if i < numbnd
        fprintf(gfid,'Domain %s: +%d %d\n', names{i}, i+1, -i);
      else
        fprintf(gfid,'Domain %s: %d   \n', names{i}, -i);
      end
    end
  otherwise
    error('version should be either ''1.0'' or ''1.1''');
end
fclose(gfid);

