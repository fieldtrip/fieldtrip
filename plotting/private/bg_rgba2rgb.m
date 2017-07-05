function [rgb] = bg_rgba2rgb(rgba, bg, varargin)

if size(bg,1) ~= size(rgba,1)
    error('the background and data should have the same number of points in the rows');
end

if size(rgba,2)==1,
    % this requires the data to be converted into rgb values first, and
    % needs a cmap + clim, and an alpha, alphamap and alim
    if numel(varargin)~=5,
        error('if a vector of data points is supplied, more input arguments are required');
    end
    cmap = varargin{1};
    if ischar(cmap)
        cmap = colormap(cmap);
    elseif size(cmap,2)~=3
        error('invalid specification of colormap, should be nx3');
    end
    
    clim = varargin{2};
    
    dat     = rgba;
    finvals = isfinite(rgba);
    
    rgba        = zeros(size(dat,1),4);
    rgba(finvals,1:3) = dat2rgb(dat(finvals), cmap, clim);
    
    alpha = varargin{3};
    amap  = varargin{4};
    alim  = varargin{5};
    
    finvals   = isfinite(alpha);
    rgba(finvals,4) = alpha2a(alpha(finvals), amap, alim);
    rgba(~finvals,4) = 0;
    % rgba(:,4) = 0.5; 
end

if size(bg,2)~=3
    % convert this into a grayscale image
    if all(bg<=1) && all(bg>=0)
        bg = repmat(bg,[1 3]);
    else
        bg = bg-min(bg);
        bg = bg/max(bg);
        bg = repmat(bg, [1 3]);
    end
end



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



