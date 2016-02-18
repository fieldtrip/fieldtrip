function [source] = estimate_fwhm1(source, removecenter)

% ESTIMATE_FWHM1(SOURCE, REMOVECENTER)
%
% This function computes the fwhm of the spatial filters, according to
% Barnes et al 2003. the input source-structure should contain the filters
% The fwhm-volume is appended to the output source-structure. It is assumed
% that the dipole positions are defined on a regularly spaced 3D grid.
% 
% This function can only deal with scalar filters.


if nargin<2, removecenter = 1; end

if ~isfield(source, 'dim'),
  source.dim = pos2dim(source.pos);
end
inside  = source.inside;
ninside = length(inside);

if ~isfield(source.avg, 'filter')
  error('the input should contain spatial filters in');
end

nchan   = size(source.avg.filter{inside(1)},2);
ndir    = size(source.avg.filter{inside(1)},1);
if ndir~=1, 
  error('only scalar filters are allowed as input');
end

%create insidevol as a binary volume
insidevol         = false(source.dim);
insidevol(inside) = true;

%get filters
filter = zeros(ndir*prod(source.dim), nchan); 
index3 = reshape(1:(ndir*prod(source.dim)), [ndir prod(source.dim)]);
for k = 1:ninside
  ind                     = inside(k);
  filter(index3(:,ind),:) = source.avg.filter{ind}; 
end

%define tetraeders: the corners of each cube of voxels are
%numbered 1-8, 
% 1 = x   y   z
% 2 = x+1 y   z
% 3 = x+1 y+1 z
% 4 = x   y+1 z
% see e.g. Worsley 1999

tetra1 = [1 2 4 5;...
          3 4 2 7;...
          6 7 5 2;...
          8 5 7 4;...
          2 4 5 7];

tetra2 = [2 3 1 6;...           
          4 1 6 8;...
          5 6 8 1;...
          7 8 6 3;...
          1 3 6 8];

tetra = [tetra1;tetra2];

%index voxels
vol  = reshape(1:prod(source.dim),source.dim);
  
indx = zeros(source.dim);
f    = zeros(source.dim);
for j = 1:prod(source.dim-1)
  
  [x, y, z] = ind2sub(source.dim-1, j);
    
  %each cube of 8 is now defined as follows:
  voxindx = vol([x x+1],[y y+1],[z z+1]);
  tmpins  = insidevol([x x+1],[y y+1],[z z+1]);
  voxindx = voxindx(tmpins(:));
  
  if length(voxindx)==8,
    indf = reshape(index3(:,voxindx), [8 1]);
  else
    continue;
  end

  Cfmat = filter(indf,:)*filter(indf,:)';
  Cfmat = Cfmat./sqrt(abs(diag(Cfmat))*abs(diag(Cfmat))'); %allow for complex

  for k = 1:size(tetra,1)
    %squared correlations  
    tmp = abs(Cfmat(tetra(k,:), tetra(k,:))).^2;
          
    %Barnes&Hillebrand eq. 19
    a   = 1 + tmp(2:4, 2:4) - [tmp(1,2:4);tmp(1,2:4);tmp(1,2:4)] - ...
	                            [tmp(2:4,1) tmp(2:4,1) tmp(2:4,1)];
    %the previous line achieves the same as the next line, but much faster
	  %a   = 1 + tmp(2:4, 2:4) - repmat(tmp(1,2:4), [3 1]) - repmat(tmp(2:4,1), [1 3]);
          %if det(a)<0, keyboard;end
          %if det(a)>3, keyboard;end
          
    %Barnes&Hillebrand eq. 15 + sum tetraeders
    f(voxindx(:)) = f(voxindx(:)) + sqrt(det(a));
    %f(voxindx(:)) = f(voxindx(:)) + sqrt(abs(det(a)));
  end
  %indx(voxindx(:)) = indx(voxindx(:)) + 1;
  indx(voxindx(:)) = indx(voxindx(:)) + 2; %if using ten tetraeders you have to correct extra
end
  
%Barnes&Hillebrand eq. 15
r    = (((4*log(2)).^(-3/2))./6).*(f./indx);

% get the approximate voxel spacing, also works for non-isotropic volumes
dx = zeros(source.dim);
cnt  = zeros(source.dim);

pos  = zeros([source.dim 3]);
pos(:,:,:,1) = reshape(source.pos(:,1), source.dim);
pos(:,:,:,2) = reshape(source.pos(:,2), source.dim);
pos(:,:,:,3) = reshape(source.pos(:,3), source.dim);

% wiggle around
dx(2:end,:,:)   = dx(2:end,:,:)   + sqrt(sum((pos(1:end-1,:,:,:)-pos(2:end,:,:,:)).^2,4));
dx(1:end-1,:,:) = dx(1:end-1,:,:) + sqrt(sum((pos(1:end-1,:,:,:)-pos(2:end,:,:,:)).^2,4));
dx(:,2:end,:  ) = dx(:,2:end,:  ) + sqrt(sum((pos(:,1:end-1,:,:)-pos(:,2:end,:,:)).^2,4));
dx(:,1:end-1,:) = dx(:,1:end-1,:) + sqrt(sum((pos(:,1:end-1,:,:)-pos(:,2:end,:,:)).^2,4));
dx(:,:,2:end  ) = dx(:,:,2:end  ) + sqrt(sum((pos(:,:,1:end-1,:)-pos(:,:,2:end,:)).^2,4));
dx(:,:,1:end-1) = dx(:,:,1:end-1) + sqrt(sum((pos(:,:,1:end-1,:)-pos(:,:,2:end,:)).^2,4));
cnt(2:end,  :,:) = cnt(2:end,  :,:) + 1;
cnt(1:end-1,:,:) = cnt(1:end-1,:,:) + 1;
cnt(:,2:end,  :) = cnt(:,2:end,  :) + 1;
cnt(:,1:end-1,:) = cnt(:,1:end-1,:) + 1;
cnt(:,:,2:end  ) = cnt(:,:,2:end  ) + 1;
cnt(:,:,1:end-1) = cnt(:,:,1:end-1) + 1;

dx = dx./cnt;

%dx   = sqrt(sum((source.pos(2,:)-source.pos(1,:)).^2)); %assuming isotropic volume

%Barnes&Hillebrand eq. 16

% fwhm = dx./(3.*sqrt(r)); % this is probably a typo in the paper, spotted by Giorgos
fwhm = dx./(r.^(1/3)); % This is the correct version
  
if removecenter,
  %FIXME this is a bit strange
  %remove the centre of the head, which has artificially low fwhm
  centre = source.dim./2;
  range  = source.dim./6;
  xrange = round(centre(1)-range(1)):round(centre(1)+range(1));
  yrange = round(centre(2)-range(2)):round(centre(2)+range(2));
  zrange = round(centre(3)-range(3)):round(centre(3)+range(3));
  
  %only thought through for 3D fwhm
  tmp2 = ones(size(fwhm));
  tmp  = fwhm(xrange,yrange,zrange);
  tmp  = tmp<2*min(tmp(:));
  tmp  = convn(convn(tmp, conndef(3,'min'), 'same'), conndef(3,'min'), 'same');
  tmp2(xrange,yrange,zrange) = tmp==0;
  %fwhm(tmp2==0) = max(max(max(fwhm(xrange,yrange,zrange))));
  fwhm(tmp2==0) = nan;
  fwhm(source.outside) = nan;
  %fwhm  = fwhm.*double(tmp2);
end

source.fwhm      = fwhm(:);
source.indx      = indx(:);
source.insideold = source.inside;
source.inside    = find(isfinite(fwhm));
source.outside   = setdiff((1:size(source.pos,1))', source.inside);

%in the extreme cases a is either: [2 1 1;1 2 1;1 1 2]; det(a) = 4; no correlation      = small FWHM
%or                   a is       : [0 0 0;0 0 0;0 0 0]; det(a) = 0; perfect correlation = large FWHM
