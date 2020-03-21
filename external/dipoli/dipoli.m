function [vol] = dipoli(vol, isolated)

% DIPOLI computes the BEM system matrix
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

%$Id$

warning('DIPOLI is deprecated, please use FT_PREPARE_HEADMODEL with cfg.method = ''dipoli'' instead.')

% find the location of the binary
str = which('dipoli.m');
[p, f, x] = fileparts(str);
dipoli = fullfile(p, f);  % without the .m extension

dipoli = checkplatformbinary(dipoli);

if ~isempty(dipoli)

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
  bnddip = vol.bnd;
  bndfile = {};

  for i=1:length(bnddip)
    bndfile{i} = [tempname '.tri'];
    % make sure that normals on the vertices point outwards
    ok = checknormals(bnddip(i));
    if ~ok,  bnddip(i).tri = fliplr(bnddip(i).tri);end
    write_tri(bndfile{i}, bnddip(i).pos, bnddip(i).tri);
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
    vol = ama2headmodel(ama);
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
else
  error('Binary file not found or not implemented for this platform!')
end

function binname = checkplatformbinary(pathname)

binname = [];
% check for the appropriate platform dependent filename
c = computer;
is64 = ismember(lower(c),{'maci64' 'glnxa64' 'sol64' 'pcwin64'});
is32 = ismember(lower(c),{'maci' 'mac' 'pcwin'});

if ispc
  binname = [pathname '.exe'];
elseif ismac
  if is64
    allowedbinnames = {[pathname '.maci64'] [pathname '.maci'] [pathname '.mac']};
  else
    allowedbinnames = {[pathname '.maci'] [pathname '.mac']};
  end
elseif isunix
  if is64
    allowedbinnames = {[pathname '.glnxa64'] [pathname '.glnx86']};
  else
    allowedbinnames = [pathname '.glnx86'];
  end
else
  [pathname '.sol64'];
end

for i = 1:length(allowedbinnames)
  if exist(allowedbinnames{i})
    binname = allowedbinnames{i};
    return
  end
end

function ok = checknormals(bnd)
ok = 0;
pos = bnd.pos;
tri = bnd.tri;
% translate to the center
org = mean(pos,1);
pos(:,1) = pos(:,1) - org(1);
pos(:,2) = pos(:,2) - org(2);
pos(:,3) = pos(:,3) - org(3);

w = sum(solid_angle(pos, tri));

if w<0 && (abs(w)-4*pi)<1000*eps
  % FIXME: this method is rigorous only for star shaped surfaces
  warning('your normals are not oriented correctly')
  ok = 0;
elseif w>0 && abs(w-4*pi)<1000*eps
  ok = 1;
else
  error('your surface probably is irregular')
  ok = 0;
end
