function [rgb] = bg_rgba2rgb(bg, rgba, varargin)

% BG_RGBA2RGB overlays a transparency masked colored image on a colored background,
% and represents the result as an RGB matrix.
%
% Use as:
%   rgb = bg_rgba2rgb(bg, rgba)
%
% or
%   rgb = bg_rgba2rgb(bg, rgba, cmap, clim, alpha, amap, alim);
%
% When 2 input arguments are supplied:
%   bg   = Nx3 matrix of background rgb-coded color-values, or MxNx3
%   rgba = Nx4 matrix of rgb + alpha values, or MxNx4
%
% When 7 input arguments are supplied:
%   bg   = Nx3 matrix, Nx1 vector, 1x3 vector, MxN, or MxNx3.
%   rgba = Nx1 vector with 'functional values', or MxN.
%   cmap = Mx3 colormap, or MATLAB-supported name of colormap
%   clim = 1x2 vector denoting the color limits
%   alpha = Nx1 vector with 'alpha values', or MxN
%   amap = Mx1 alphamap, or MATLAB -supported name of alphamap ('rampup/down', 'vup/down')
%   alim = 1x2 vector denoting the opacity limits

if numel(varargin)==0
  siz1 = size(bg);
  siz2 = size(rgba);
  assert(isequal(siz1(1:end-1),siz2(1:end-1)), 'number of data points should be the same');
  assert(siz1(end)==3 && siz2(end)==4,         'inconsistent input size');
  
  bg   = reshape(bg,   [], 3);
  rgba = reshape(rgba, [], 4);
  rgb  = do_conversion(bg, rgba);
  rgb  = reshape(rgb, siz1);
else
  % this requires the data to be converted into rgb values irst, and
  % needs a cmap + clim, and an alpha, alphamap and alim
  if numel(varargin)~=5
    error('if a vector of color data is supplied, more input arguments are required');
  end
  
  % colormap handling
  cmap = varargin{1};
  if ischar(cmap)
    cmap = ft_colormap(cmap);
  elseif size(cmap,2)~=3
    error('invalid specification of colormap, should be nx3');
  end
  clim  = varargin{2}; if isempty(clim), clim(1) = min(rgba(:)); clim(2) = max(rgba(:)); end
  alpha = varargin{3};
  amap  = varargin{4};
  alim  = varargin{5}; if isempty(alim), alim(1) = min(alpha(:)); alim(2) = max(alpha(:)); end
        
  % deal with the color data
  dat     = rgba;
  siz     = size(dat);
  dat     = reshape(dat, [], 1);
  alpha   = reshape(alpha, [], 1); % assume to be same siza as input data!
  finvals = isfinite(dat);
  
  rgba              = zeros(size(dat,1),4);
  rgba(finvals,1:3) = dat2rgb(dat(finvals), cmap, clim);
    
  finvals          = isfinite(alpha);
  rgba(finvals,4)  = alpha2a(alpha(finvals), amap, alim);
  rgba(~finvals,4) = 0;
       
  % deal with the background: allow it to be 1x3, i.e. a single color per
  % pixel.
  if isequal(size(bg),[1 3])
    bg = permute(repmat(bg(:), [1 siz]), [1+(1:numel(siz)) 1]);
  else
    siz_bg = size(bg);
    if isequal(siz_bg, siz)
      % make bg rgb
      if all(bg(:)<=1) && all(bg(:)>=0)
        % don't scale
      else
        bg_min = min(bg(:));
        bg_max = max(bg(:));
        bg     = (bg-bg_min)./(bg_max-bg_min);
      end
      bg = repmat(bg, [ones(1,ndims(bg)) 3]);
    else
      %FIXME
    end
  end
  bg  = reshape(bg, [], 3);
  
  if numel(siz)==2 && siz(2)==1, siz = siz(1); end
  rgb = do_conversion(bg, rgba);
  rgb = reshape(rgb, [siz 3]);
end

function rgb = do_conversion(bg, rgba)

rgb = zeros(size(rgba,1),3);
a_  = 1-rgba(:,4);
a   = rgba(:,4);

for k = 1:3
    rgb(:,k) = bg(:,k).*a_  + rgba(:,k).*a;
end


function rgb = dat2rgb(dat, cmap, clim)

dat(end+1) = clim(1);
dat(end+1) = clim(2); % add the extremes to be sure that they are included

dat(dat<clim(1)) = clim(1);
dat(dat>clim(2)) = clim(2);


% scale between 0 and 1
dat = dat-min(dat);
dat = dat/max(dat);

ind = round(dat.*(size(cmap,1)-1))+1;
rgb = cmap(ind(1:end-2),:);

function a = alpha2a(alpha, amap, alim)

alpha(end+1) = alim(1);
alpha(end+1) = alim(2);

alpha(alpha<alim(1)) = alim(1);
alpha(alpha>alim(2)) = alim(2);

if ischar(amap)
    switch amap
        case 'rampup'
            a = alpha - min(alpha);
            a = a./max(a);
        case 'rampdown'
            a = alpha - min(alpha);
            a = a./max(a);
            a = 1-a;
        case 'vup'
            a = alpha - min(alpha);
            a = a./max(a);
            a = 1 - 2.*abs(a - 0.5);
        case 'vdown'
            a = alpha - min(alpha);
            a = a./max(a);
            a = 2.*abs(a - 0.5);
        otherwise
            error('unknown alphamap specified');
    end
else
    amap = amap(:);
    
    % assume it to be a vector that has M elements, scaled between 0 and 1
    a = alpha-min(alpha);
    a = a/max(a);
    ind = round(a.*(size(amap,1)-1))+1;
    a   = amap(ind(1:end),:);
end
a = a(1:end-2);



