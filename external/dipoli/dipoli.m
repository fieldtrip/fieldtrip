function [vol] = dipoli(vol, isolated)

% DIPOLI computes a BEM system matrix
%
% Use as
%   [vol] = dipoli(vol, isolated)

% Copyright (C) 2005-2008, Robert Oostenveld
%
% $Log: dipoli.m,v $
% Revision 1.3  2008/12/24 10:25:41  roboos
% cleaned up the dipoli wrapper, replaced the binary by a better one and added a copy of the helper functions (from fileio)
%
% Revision 1.2  2008/12/24 09:20:11  roboos
% added subfunctions etc as dipoli_xxx and modified dipoli main function accordingly
%
% Revision 1.1.1.1  2008/12/24 08:52:28  roboos
% created new module that will hold the dipoli specific stuff for prepare_bemmodel
%
% Revision 1.2  2006/01/20 09:48:36  roboos
% fill remaining elements of matrix with zeros
% changed reshape and transpose
%

% find the location of the binary
str = which('dipoli.m');
[p, f, x] = fileparts(str);
dipoli = fullfile(p, f);  % without the .m extension

skin   = find_outermost_boundary(vol.bnd);
source = find_innermost_boundary(vol.bnd);

% the first compartment should be the skin, the last the source
if skin==1 && source==length(vol.bnd)
    vol.skin   = 1;
    vol.source = length(vol.bnd);
elseif skin==length(vol.bnd) && source==1
    % flip the order of the compartments
    vol.bnd    = fliplr(vol.bnd(:)');
    vol.skin   = 1;
    vol.source = length(vol.bnd);
else
    error('the first compartment should be the skin, the last  the source');
end

if isolated
    fprintf('using the isolated source approach\n');
else
    fprintf('not using isolated source approach\n');
end

% write the triangulations to file
bndfile = {};
for i=1:length(vol.bnd)
    bndfile{i} = [tempname '.tri'];
    % dipoli has another definition of the direction of the surfaces
    vol.bnd(i).tri = fliplr(vol.bnd(i).tri);
    write_tri(bndfile{i}, vol.bnd(i).pnt, vol.bnd(i).tri);
end

% these will hold the shell script and the inverted system matrix
exefile = [tempname '.sh'];
amafile = [tempname '.ama'];

fid = fopen(exefile, 'w');
fprintf(fid, '#!/bin/sh\n');
fprintf(fid, '\n');
fprintf(fid, '%s -i %s << EOF\n', dipoli, amafile);
for i=1:length(vol.bnd)
    if isolated && vol.source==i
        % the isolated potential approach should be applied using this compartment
        fprintf(fid, '!%s\n', bndfile{i});
    else
        fprintf(fid, '%s\n', bndfile{i});
    end
    fprintf(fid, '%g\n', vol.cond(i));
end
fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, 'EOF\n');
fclose(fid);
dos(sprintf('chmod +x %s', exefile));

try
    % execute dipoli and read the resulting file
    dos(exefile);
    ama = loadama(amafile);
    vol = ama2vol(ama);
catch
    warning('an error ocurred while running dipoli');
    disp(lasterr);
end

% delete the temporary files
for i=1:length(vol.bnd)
    delete(bndfile{i})
end
delete(amafile);
delete(exefile);


