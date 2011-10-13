function bnd = prepare_mesh_segmentation(cfg,mri)

% PREPARE_MESH_SEGMENTATION
% Calculates the surfaces of each compartment from a segmentation or from
% an MRI volume structure to be segmented. As such the input always
% contains volumetric information.
% 
% The segmentation can be of two different types:
% - a raw segmentation: given as a 3D data matrix 
%   Requires the cfg.tissue option to be a set of numbers (intensity value of the
%   compartments)
% 
% - a mri segmentation: given as a mri volume with the segmented compartments as fields
%   Requires the cfg.tissue option to be a set of strings (the names of the
%   fields of actual input structure)
% 
% Required options:
% 
% cfg.tissue                  the tissue number/string
% cfg.numvertices             the desired number of vertices
% 
% See also PREPARE_MESH_MANUAL,PREPARE_MESH_HEADSHAPE

% Copyrights (C) 2009, Robert Oostenveld, 2011, Cristiano Micheli
%
% $Log$

% process the inputs
tissue      = ft_getopt(cfg,'tissue');
numvertices = ft_getopt(cfg,'numvertices');
smoothseg   = ft_getopt(cfg,'smoothseg',0);
threshseg   = ft_getopt(cfg,'thresholdseg',0);

% check options consistency
ntissues = numel(tissue);
if ntissues>1
  % assign one single parameter to each tissue if more than one tissue
  if numel(numvertices)==1
    numvertices(1:ntissues) = numvertices(1);
  elseif numel(numvertices)~=ntissues
    error('tissues and vertices do not match')
  end
  if numel(smoothseg)==1
    smoothseg(1:ntissues) = smoothseg(1);
  elseif numel(smoothseg)~=ntissues
    error('tissues and smooth parameter do not match')
  end  
  if numel(threshseg)==1
    threshseg(1:ntissues) = threshseg(1);
  elseif numel(threshseg)~=ntissues
    error('tissues and threshold parameter do not match')
  end
else % only one tissue, choose the first parameter
  numvertices = numvertices(1);
  smoothseg   = smoothseg(1);
  threshseg   = threshseg(1);
end

% do the mesh extrapolation
for i =1:numel(tissue)
  if ~isnumeric(tissue(i))
    comp = tissue{i};
    seg = mri.(comp);
  else
    seg = mri==tissue(i);
  end
  seg = dosmoothing(seg, smoothseg(i), num2str(i));
  seg = threshold(seg, threshseg(i), num2str(i));
  seg = fill(seg, num2str(i)); %FIXME: try fixhollow
  bnd(i) = dotriangulate(seg, numvertices(i), num2str(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bnd = dotriangulate(seg, nvert, str)
% str is just a placeholder for messages
dim = size(seg);
[mrix, mriy, mriz] = ndgrid(1:dim(1), 1:dim(2), 1:dim(3));

% construct the triangulations of the boundaries from the segmented MRI
fprintf('triangulating the boundary of compartment %s\n', str);
ori(1) = mean(mrix(seg(:)));
ori(2) = mean(mriy(seg(:)));
ori(3) = mean(mriz(seg(:)));
[pnt, tri] = triangulate_seg(seg, nvert, ori); % is tri okay? tri = projecttri(pnt);

% output
bnd.pnt = pnt;
bnd.tri = tri;
fprintf(['segmentation compartment %s completed\n'],str);

function [output] = threshold(input, thresh, str)

if thresh>0 && thresh<1
  fprintf('thresholding %s at a relative threshold of %0.3f\n', str, thresh);

  % mask by taking the negative of the brain, thus ensuring
  % that no holes are within the compartment and do a two-pass 
  % approach to eliminate potential vitamin E capsules etc.

  output   = double(input>(thresh*max(input(:))));
  [tmp, N] = spm_bwlabel(output, 6);
  for k = 1:N
    n(k,1) = sum(tmp(:)==k);
  end
  output   = double(tmp~=find(n==max(n))); clear tmp;
  [tmp, N] = spm_bwlabel(output, 6);
  for k = 1:N
    m(k,1) = sum(tmp(:)==k);
  end
  output   = double(tmp~=find(m==max(m))); clear tmp;
else
  output = input;
end

function [output] = fill(input, str)
fprintf('filling %s\n', str);
  output = input;
  dim = size(input);
  for i=1:dim(2)
    slice=squeeze(input(:,i,:));
    im = imfill(slice,8,'holes');
    output(:,i,:) = im;
  end

function [output] = dosmoothing(input, fwhm, str)
if fwhm>0
  fprintf('smoothing %s with a %d-voxel FWHM kernel\n', str, fwhm);
  spm_smooth(input, input, fwhm);
end
output = input;

% function cmprtmnt = fixhollow(cmprtmnt)
% % checks if the compartment is hollow
% [~,N] = bwlabeln(cmprtmnt);
% if N>2
% % ...and tries to fix it
%   cmprtmnt = imfill(cmprtmnt,'holes');
% end  
