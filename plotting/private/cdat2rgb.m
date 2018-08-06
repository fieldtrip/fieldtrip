function rgb = cdat2rgb(cdat, cmap, clim, highlight)

% This funcxtion changes the color of pixels to white, regardless of colormap, without using opengl
% It does by converting by:
% 1) convert the to-be-plotted data to their respective rgb color values (determined by colormap)
% 2) convert these rgb color values to hsv values, hue-saturation-value
% 3) for to-be-masked-pixels, set saturation to 0 and value to 1 (hue is irrelevant when they are)
% 4) convert the hsv values back to rgb values

% enforce mask properties (satmask is 0 when a pixel needs to be masked, 1 if otherwise)
satmask = round(double(highlight));   % enforce binary white-masking, the hsv approach cannot be used for 'white-shading'
satmask(isnan(cdat)) = false;         % make sure NaNs are plotted as white pixels, even when using non-integer mask values

% do 1, by converting the data-values to zero-based indices of the colormap
ncolors = size(cmap, 1); % determines range of index, if a figure has been created by the caller function, gcf changes nothing, if not, a figure is created (which the below would do otherwise)
indcdat = (cdat + -clim(1)) * (ncolors / (-clim(1) + clim(2))); % transform cdat-values to have a 0-(ncolors-1) range (range depends on colormap used, and thus also on clim)
rgbcdat = ind2rgb(uint8(floor(indcdat)), cmap);
% do 2
hsvcdat = rgb2hsv(rgbcdat);
% do 3
hsvs = hsvcdat(:,:,2);
hsvs(~satmask) = 0;
hsvv = hsvcdat(:,:,3);
hsvv(~satmask) = 1;
hsvcdat(:,:,2) = hsvs;
hsvcdat(:,:,3) = hsvv;
% do 4
rgb = hsv2rgb(hsvcdat);
