function [transform] = ft_affinecoordinates(original, target, varargin)

% FT_AFFINECOORDINATES returns the affine coordinate transformation matrix that
% converts FROM a specific head coordinate TO a specific head coordinate system.
%
% Use as
%   [transform] = ft_affinecoordinates(from, to)
%
% Note that translations are expressed in millimeters, therefore the geometrical data
% to which this coordinate transformation is applied must also be specified in
% millimeters.
%
% See also FT_CONVERT_COORDSYS, FT_CONVERT_UNITS, FT_HEADCOORDINATES, FT_WARP_APPLY

% Copyright (C) 2005-2022, Robert Oostenveld & Jan-Mathijs Schoffelen
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% these are the 48 generic axis orientation triplets, these specify the axes but no origin
%   a = anterior
%   p = posterior
%   l = left
%   r = right
%   s = superior
%   i = inferior

% get the options
requireorigin = ft_getopt(varargin, 'requireorigin', true);
computeall    = ft_getopt(varargin, 'computeall', false);

% ensure these are in lower case
original = lower(original);
target = lower(target);

generic = {
  'als'; 'ali'; 'ars'; 'ari';...
  'pls'; 'pli'; 'prs'; 'pri';...
  'las'; 'lai'; 'ras'; 'rai';...
  'lps'; 'lpi'; 'rps'; 'rpi';...
  'asl'; 'ail'; 'asr'; 'air';...
  'psl'; 'pil'; 'psr'; 'pir';...
  'sal'; 'ial'; 'sar'; 'iar';...
  'spl'; 'ipl'; 'spr'; 'ipr';...
  'sla'; 'ila'; 'sra'; 'ira';...
  'slp'; 'ilp'; 'srp'; 'irp';...
  'lsa'; 'lia'; 'rsa'; 'ria';...
  'lsp'; 'lip'; 'rsp'; 'rip'}';

% the specific ones specify the axes and also the origin
specific = {'ctf', 'bti', 'fourd', 'yokogawa', 'eeglab', 'neuromag', 'itab', 'acpc', 'spm', 'mni', 'fsaverage', 'tal', 'scanras', 'scanlps', 'dicom'};

% pre allocate all exhaustive combinations
coordsys = [specific generic];
n        = numel(coordsys);
ngeneric = numel(generic);
nspecific = numel(specific);
n1       = repmat(1:n,    n, 1);
n2       = repmat((1:n)', 1, n);

from = coordsys(n1(:))';
to   = coordsys(n2(:))';
T    = nan*zeros(4, 4, numel(to));

% numeric indices of the corresponding from/to lists into the static
% coordsys list, speeds up comparisons in the below.
from_idx = n1(:);
to_idx   = n2(:);

% do this once, and use it to check for the existence of the transform, to
% be able to break out of the exhaustive for-loop below
original_idx = strcmp(coordsys, original);
target_idx   = strcmp(coordsys, target);

% make the combined and the inverse transformations where possible
implemented = zeros(length(coordsys)); % this is only for debugging

%--------------------------------------------------------------------------
% These approximate alignments are based on transformation matrices that were
% determined by clicking on the CTF fiducial locations in the canonical T1 template
% MRI.
%
% All of the transformation matrices here are expressed in millimeter.

% this is based on the ear canals, see ALIGN_CTF2ACPC
fromthis = 'acpc';
tothis   = 'ctf';
M = [
   0.0000  0.9987  0.0517  34.7467
  -1.0000  0.0000  0.0000   0.0000
   0.0000 -0.0517  0.9987  52.2749
   0.0000  0.0000  0.0000   1.0000
  ];
[T, implemented] = transform_add(M, fromthis, tothis, T, from, to, implemented);

% this is based on the ear canals, see ALIGN_NEUROMAG2ACPC
fromthis = 'acpc';
tothis   = 'neuromag';
M = [
  1.0000  0.0000  0.0000   0.0000
  0.0000  0.9987  0.0517  34.7467
  0.0000 -0.0517  0.9987  52.2749
  0.0000  0.0000  0.0000   1.0000
  ];
[T, implemented] = transform_add(M, fromthis, tothis, T, from, to, implemented);

% see http://freesurfer.net/fswiki/CoordinateSystems
fromthis = 'fsaverage';
tothis   = 'mni';
M = [
   0.9975   -0.0073    0.0176   -0.0429
   0.0146    1.0009   -0.0024    1.5496
  -0.0130   -0.0093    0.9971    1.1840
   0.0000    0.0000    0.0000    1.0000
  ];
[T, implemented] = transform_add(M, fromthis, tothis, T, from, to, implemented);

% this is a 90 degree rotation around the z-axis
fromthis = 'ctf';
tothis   = 'neuromag';
M = [
  0.0000   -1.0000    0.0000    0.0000
  1.0000    0.0000    0.0000    0.0000
  0.0000    0.0000    1.0000    0.0000
  0.0000    0.0000    0.0000    1.0000
  ];
[T, implemented] = transform_add(M, fromthis, tothis, T, from, to, implemented);

% affine transformation from MNI to Talairach, see http://imaging.mrc-cbu.cam.ac.uk/imaging/MniTalairach
% the non-linear (i.e. piecewise linear) transform between MNI and Talairach are implemented elsewhere, see the functions MNI2TAL and TAL2MNI
fromthis = 'mni';
tothis   = 'tal';
M = [
  0.8800    0.0000    0.0000   -0.8000
  0.0000    0.9700    0.0000   -3.3200
  0.0000    0.0500    0.8800   -0.4400
  0.0000    0.0000    0.0000    1.0000
  ];
[T, implemented] = transform_add(M, fromthis, tothis, T, from, to, implemented);

% the CTF and BTI coordinate system are the same, see http://www.fieldtriptoolbox.org/faq/coordsys/
fromthis = 'ctf';
tothis   = 'bti';
[T, implemented] =  transform_add(eye(4), fromthis, tothis, T, from, to, implemented);

% the CTF and EEGLAB coordinate system are the same, see http://www.fieldtriptoolbox.org/faq/coordsys/
fromthis = 'ctf';
tothis   = 'eeglab';
[T, implemented] =  transform_add(eye(4), fromthis, tothis, T, from, to, implemented);

% the NEUROMAG and ITAB coordinate system are the same, see http://www.fieldtriptoolbox.org/faq/coordsys/
fromthis = 'neuromag';
tothis   = 'itab';
[T, implemented] =  transform_add(eye(4), fromthis, tothis, T, from, to, implemented);

% BTI and 4D are different names for exactly the same system, see http://www.fieldtriptoolbox.org/faq/coordsys/
fromthis = 'bti';
tothis   = 'fourd';
[T, implemented] =  transform_add(eye(4), fromthis, tothis, T, from, to, implemented);

% the YOKOGAWA system expresses positions relative to the dewar, not relative to the head
% see https://www.fieldtriptoolbox.org/faq/coordsys/#details-of-the-yokogawa-coordinate-system
fromthis = 'yokogawa';
tothis   = 'als';
[T, implemented] =  transform_add(eye(4), fromthis, tothis, T, from, to, implemented);

% the SPM and MNI coordinate system are the same, see http://www.fieldtriptoolbox.org/faq/acpc/
fromthis = 'spm';
tothis   = 'mni';
[T, implemented] =  transform_add(eye(4), fromthis, tothis, T, from, to, implemented);

% the SPM (aka MNI) and ACPC coordinate system are not the same, but similar enough, see http://www.fieldtriptoolbox.org/faq/acpc/
fromthis = 'spm';
tothis   = 'acpc';
[T, implemented] =  transform_add(eye(4), fromthis, tothis, T, from, to, implemented);

fromthis = 'mni';
tothis   = 'acpc';
[T, implemented] =  transform_add(eye(4), fromthis, tothis, T, from, to, implemented);

% the CTF, BTI and 4D coordinate systems are all ALS coordinate systems
% but the origin is poorly defined in ALS, hence converting from ALS to another is problematic
fromthis = 'ctf';
tothis   = 'als';
[T, implemented] =  transform_add(eye(4), fromthis, tothis, T, from, to, implemented);
fromthis = 'bti';
tothis   = 'als';
[T, implemented] =  transform_add(eye(4), fromthis, tothis, T, from, to, implemented);
fromthis = 'fourd';
tothis   = 'als';
[T, implemented] =  transform_add(eye(4), fromthis, tothis, T, from, to, implemented);

% the NEUROMAG, ITAB, ACPC, MNI, SPM and FSAVERAGE coordinate systems are all RAS coordinate systems
% but the origin is poorly defined in RAS, hence converting from RAS to another is problematic
fromthis = 'neuromag';
tothis   = 'ras';
[T, implemented] =  transform_add(eye(4), fromthis, tothis, T, from, to, implemented);
fromthis = 'itab';
tothis   = 'ras';
[T, implemented] =  transform_add(eye(4), fromthis, tothis, T, from, to, implemented);
fromthis = 'acpc';
tothis   = 'ras';
[T, implemented] =  transform_add(eye(4), fromthis, tothis, T, from, to, implemented);
fromthis = 'mni';
tothis   = 'ras';
[T, implemented] =  transform_add(eye(4), fromthis, tothis, T, from, to, implemented);
fromthis = 'spm';
tothis   = 'ras';
[T, implemented] =  transform_add(eye(4), fromthis, tothis, T, from, to, implemented);
fromthis = 'fsaverage';
tothis   = 'ras';
[T, implemented] =  transform_add(eye(4), fromthis, tothis, T, from, to, implemented);
fromthis = 'tal';
tothis   = 'ras';
[T, implemented] =  transform_add(eye(4), fromthis, tothis, T, from, to, implemented);

% the SCANRAS coordinate system is RAS with the origin at the center of the gradient coil
fromthis = 'scanras';
tothis   = 'ras';
[T, implemented] =  transform_add(eye(4), fromthis, tothis, T, from, to, implemented);

% the DICOM and SCANLPS coordinate system are the same, and rotated 180 degrees from SCANRAS
fromthis = 'dicom';
tothis   = 'scanlps';
[T, implemented] =  transform_add(eye(4), fromthis, tothis, T, from, to, implemented);
fromthis = 'dicom';
tothis   = 'lps';
[T, implemented] =  transform_add(eye(4), fromthis, tothis, T, from, to, implemented);
fromthis = 'scanlps';
tothis   = 'lps';
[T, implemented] =  transform_add(eye(4), fromthis, tothis, T, from, to, implemented);
fromthis = 'scanlps';
tothis   = 'scanras';
[T, implemented] =  transform_add(eye(4), fromthis, tothis, T, from, to, implemented);

% these are all combinations of 90 degree rotations and/or flips along one of the axes
for i=1:length(generic)
  for j=(i+1):length(generic)
    fromthis   = generic{i};
    tothis     = generic{j};

    ii = nspecific + i; %find(strcmp(coordsys, fromthis));
    jj = nspecific + j; %find(strcmp(coordsys, tothis));
    sel = ii + (jj-1)*n; %from_idx==ii & to_idx==jj;
    M = transform_generic(fromthis, tothis);
    T(:,:,sel) = M;
    sel = jj + (ii-1)*n; %from_idx==jj & to_idx==ii;
    T(:,:,sel) = inv(M);
  end
end
sel = ismember(coordsys, generic);
implemented(sel, sel) = 1;

for k=1:numel(from)
  % this is an exhaustive for-loop that covers all possible combinations,
  % but breaks once the required one has been computed
  
  fromthis = from{k};
  tothis   = to{k};
    
  % numeric indices into the (static) coordsys list, speeds up because numeric comparisons are faster than string comparisons
  i   = find(strcmp(coordsys, fromthis));
  j   = find(strcmp(coordsys, tothis));

  if implemented(original_idx, target_idx) && ~computeall
    break;
  end

  if isfinite(T(1,1,k))
    % the transform already exists already exists

  elseif isequal(fromthis, tothis)
    % construct the transformations on the diagonal
    T(:,:,k) = eye(4);
   
    implemented(i,j) = 1;

  else
    % the transformation matrix does not yet exist
    
    % check whether the inverse exists
    sel_inv = from_idx==j & to_idx==i; %strcmp(from, tothis) & strcmp(to, fromthis);
    if sum(sel_inv) && isfinite(T(1,1,sel_inv))
      T(:,:,k) = inv(T(:,:,sel_inv));
      
      implemented(i,j) = 2;
      implemented(j,i) = 2;
    else
      from_generic = ismember(fromthis, generic);
      to_generic   = ismember(tothis, generic);
        
      % a two-step approach is needed
      if ~from_generic && to_generic
        % try to make the transformation (and inverse) with a two-step approach
        % since we go from specific to generic and thereby lose the origin information anyway, it is fine to use any intermediate step
        for kk=1:numel(coordsys)
          tointermediate = coordsys{kk};
          ij = find(strcmp(coordsys, tointermediate));

          sel1 = from_idx==i  & to_idx==ij; % strcmp(from, fromthis)       & strcmp(to, tointermediate);
          sel2 = from_idx==ij & to_idx==j;  % strcmp(from, tointermediate) & strcmp(to, tothis);
          if all(isfinite(T(1,1,sel1|sel2)))
            T(:,:,k) = T(:,:,sel2) * T(:,:,sel1);
            
            % also add the inverse, if needed
            if ~isfinite(T(1,1,sel_inv))
              T(:,:,sel_inv) = inv(T(:,:,k));
            end
            implemented(i,j) = 3;
            implemented(j,i) = 3;
            break
          end
        end % for kk
    
      elseif ~from_generic && ~to_generic
        % try to make the transformation (and inverse) with a two-step approach
        % do not use the generic orientation triplets (like RAS and ALS) as intermediate steps between two specific coordinate systems
        for kk=1:numel(specific)
          tointermediate = coordsys{kk};
          ij = find(strcmp(coordsys, tointermediate));

          sel1 = from_idx==i  & to_idx==ij; % strcmp(from, fromthis)       & strcmp(to, tointermediate);
          sel2 = from_idx==ij & to_idx==j;  % strcmp(from, tointermediate) & strcmp(to, tothis);
          if all(isfinite(T(1,1,sel1|sel2)))
            T(:,:,k) = T(:,:,sel2) * T(:,:,sel1);
            
            % also add the inverse, if needed
            if all(~isfinite(T(:,:,sel_inv)))
              T(:,:,sel_inv) = inv(T(:,:,k));
            end
            implemented(i,j) = 3;
            implemented(j,i) = 3;
            break
          end
        end % for kk
      end
    end
  end 
end % for k

% these conversions should be done using FT_VOLUMENORMALISE, as they imply scaling
remove1 = {'acpc', 'spm'
           'acpc', 'mni'
           'acpc', 'fsaverage'
           'acpc', 'tal'};

% converting to/from TAL is only possible for some specific template coordinate systems
remove2 = {'bti',      'tal'
           'ctf',      'tal'
           'fourd',    'tal'
           'itab',     'tal'
           'neuromag', 'tal'
           'tal',      'bti'
           'tal',      'ctf'
           'tal',      'fourd'
           'tal',      'itab'
           'tal',      'neuromag'};

remove = cat(1, remove1, remove2);
removeidx = false(numel(to), 1);
for k=1:size(remove, 1)
  i = find(strcmp(coordsys, remove{k,1}));
  j = find(strcmp(coordsys, remove{k,2}));

  removeidx = removeidx | (from_idx==i & to_idx==j);
end

if requireorigin
  % the origin is poorly defined in generic orientation triplets (like RAS and ALS)
  % hence going from generic to specific is problematic with regards to translations
  for i=1:length(generic)
    for j=1:length(specific)
      ii = find(strcmp(coordsys, generic{i}));
      jj = find(strcmp(coordsys, specific{j}));
     
      removeidx = removeidx | (from_idx==ii & to_idx==jj);
    end
  end
end
T(:,:,removeidx) = nan;
implemented(removeidx) = 0;

if false
  % this is only for checking the coverage of all conversions
  
  figure; imagesc(implemented); caxis([0 3]);
  xticklabels(coordsys); xticks(1:numel(coordsys));
  yticklabels(coordsys); yticks(1:numel(coordsys));
end

if strcmp(original, '4d')
  fromthis = 'fourd'; % '4d' is not a valid variable name
else
  fromthis = original;
end

if strcmp(target, '4d')
  tothis = 'fourd'; % '4d' is not a valid variable name
else
  tothis = target;
end

sel = strcmp(from, fromthis) & strcmp(to, tothis);
if any(sel)
  transform = T(:,:,sel);
else
  transform = nan + zeros(4);
end

if any(isnan(transform))
  ft_error('converting from %s to %s is not supported', original, target);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to construct generic transformations such as RAS2ALS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T = transform_generic(from, to)

% ap_in  = find(to=='a'   | to=='p');
% ap_out = find(from=='a' | from=='p');
% lr_in  = find(to=='l'   | to=='r');
% lr_out = find(from=='l' | from=='r');
% si_in  = find(to=='s'   | to=='i');
% si_out = find(from=='s' | from=='i');
% 
% % index the axis according to ap,lr,si
% order_in  = [ap_in  lr_in  si_in];
% order_out = [ap_out lr_out si_out];
% 
% % check whether one of the axis needs flipping
% flip = 2.*(0.5-double(to(order_in)~=from(order_out)));

T = zeros(4);

% T(order_in+(order_out-1)*4) = flip;
% for k = 1:3 % the above avoids the for-loop
%   T(order_in(k),order_out(k)) = flip(k);
% end

in = [(to=='a')-(to=='p');
      (to=='l')-(to=='r');
      (to=='s')-(to=='i')];

out = [(from=='a')-(from=='p');
      (from=='l')-(from=='r');
      (from=='s')-(from=='i')];

%assert(isequal(in\out,T(1:3,1:3)));
T(1:3,1:3) = in\out;
T(4,4) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to fill a from/to slice in the 4x4xNtransform array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T, implemented] = transform_add(Tthis, fromthis, tothis, T, from, to, implemented)

indx = strcmp(from, fromthis) & strcmp(to, tothis);
T(:,:,indx) = Tthis;
implemented(indx) = 1;

% add the inverse as well
indx = strcmp(from, tothis) & strcmp(to, fromthis);
T(:,:,indx) = inv(Tthis);
implemented(indx) = 1;
