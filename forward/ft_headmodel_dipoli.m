function vol = ft_headmodel_dipoli(geom, varargin)

% FT_HEADMODEL_DIPOLI creates a volume conduction model of the head
% using the boundary element method (BEM) for EEG. This function takes
% as input the triangulated surfaces that describe the boundaries and
% returns as output a volume conduction model which can be used to
% compute leadfields.
% 
% This implements
%   Oostendorp TF, van Oosterom A. "Source parameter estimation in
%   inhomogeneous volume conductors of arbitrary shape." IEEE Trans
%   Biomed Eng. 1989 Mar;36(3):382-91.
%
% The implementation of this function uses an external command-line
% executable with the name "dipoli" which is provided by Thom Oostendorp.
%
% Use as
%   vol = ft_headmodel_dipoli(geom, ...)
%
% Optional input arguments should be specified in key-value pairs and can
% include
%   isolatedsource   = string, 'yes' or 'no'
%   hdmfile          = string, filename with BEM headmodel
%   conductivity     = vector, conductivity of each compartment
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

ft_hastoolbox('dipoli', 1);

% get the optional arguments
isolatedsource  = keyval('isolatedsource', varargin);
hdmfile         = keyval('hdmfile', varargin);
conductivity    = keyval('conductivity', varargin);

% start with an empty volume conductor
vol = [];

if ~isempty(hdmfile)
  hdm = ft_read_vol(hdmfile);
  % copy the boundary of the head model file into the volume conduction model
  vol.bnd = hdm.bnd;
  if isfield(hdm, 'cond')
    % also copy the conductivities
    vol.cond = hdm.cond;
  end
else
  % copy the boundaries from the geometry into the volume conduction model
  vol.bnd = geom.bnd;
end

% determine the number of compartments
numboundaries = length(vol.bnd);

if isempty(isolatedsource)
  if numboundaries>1
    % the isolated source compartment is by default the most inner one
    isolatedsource = true;
  else
    isolatedsource = false;
  end
else
  % convert into a boolean
  isolatedsource = istrue(cfg.isolatedsource);
end

if ~isfield(vol, 'cond')
  % assign the conductivity of each compartment
  vol.cond = conductivity;
end

% determine the nesting of the compartments
nesting = zeros(numboundaries);
for i=1:numboundaries
  for j=1:numboundaries
    if i~=j
      % determine for a single vertex on each surface if it is inside or outside the other surfaces
      curpos = vol.bnd(i).pnt(1,:); % any point on the boundary is ok
      curpnt = vol.bnd(j).pnt;
      curtri = vol.bnd(j).tri;
      nesting(i,j) = bounding_mesh(curpos, curpnt, curtri);
    end
  end
end

if sum(nesting(:))~=(numboundaries*(numboundaries-1)/2)
  error('the compartment nesting cannot be determined');
end

% for a three compartment model, the nesting matrix should look like
%    0 0 0     the first is the most outside, i.e. the skin
%    0 0 1     the second is nested inside the 3rd, i.e. the outer skull
%    0 1 1     the third is nested inside the 2nd and 3rd, i.e. the inner skull
[~, order] = sort(sum(nesting,2));

fprintf('reordering the boundaries to: ');
fprintf('%d ', order);
fprintf('\n');

% update the order of the compartments
vol.bnd    = vol.bnd(order);
vol.cond   = vol.cond(order);
vol.skin_surface   = 1;
vol.source = numboundaries;

if isolatedsource
  fprintf('using compartment %d for the isolated source approach\n', vol.source);
else
  fprintf('not using the isolated source approach\n');
end

% find the location of the dipoli binary
str = which('dipoli');
[p, f, x] = fileparts(str);
dipoli = fullfile(p, f);  % without the .m extension
switch mexext
  case {'mexmaci' 'mexmaci64'}
    % apple computer
    dipoli = [dipoli '.maci'];
  case {'mexglnx86' 'mexa64'}
    % linux computer
    dipoli = [dipoli '.glnx86'];
  otherwise
    error('there is no dipoli executable for your platform');
end
fprintf('using the executable "%s"\n', dipoli);

% write the triangulations to file
bndfile = {};
for i=1:numboundaries
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
for i=1:numboundaries
  if isolatedsource && vol.source==i
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
% ensure that the temporary shell script can be executed
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
for i=1:numboundaries
  delete(bndfile{i})
end
delete(amafile);
delete(exefile);

% remember that it is a dipoli model
vol.type = 'dipoli';

