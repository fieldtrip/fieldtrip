function [M] = ROIcmap(nc,opt)
% [M] = ROIcmap([nc],[opt]);
% creates a colormap with
%           no color too close to grayscale,
%           no two consecutive colors too close
%           no colors exeedingly close to another in the map
%           no colors too close to a background color (optional)
%   nc: number of colors in map. Default is 64
%   opt: optional options structure
%      .show: Figure handle in which to show the map
%             Default is 1. Set to 0 for no show.
%             In show mode, you get to pick colors with mouse
%             and read in their values. Hit enter to exit from
%             this mode.
%      .state: State of the random number generator.
%              Default is 0.
%      .write: Name of file to write colormap into
%              Default is '', no writing. Use something like
%              ROI64s0.1D.cmap, for a 64 cols, seed 0 colormap.
%      .avoid: Color to avoid getting close to.
%      .verb: verbosioty. 1 is default. 0 is for quiet
% returns
%   M: The colormap.
%
%see also readXcol, rgbdectohex, and ScaleToMap
% Ziad S. Saad SSCC/NIMH/NIH, saadz@mail.nih.gov

if (nargin == 0),
   nc = 64;
   opt.show = 1;
elseif (nargin == 1),
   opt.show = 1;
end

if (isempty(nc)), nc = 64; end

if (~isfield(opt,'show') | isempty(opt.show)), opt.show = 1; end
if (~isfield(opt,'state') | isempty(opt.state)), opt.state = 0; end
if (~isfield(opt,'write') | isempty(opt.write)), opt.write = ''; end
if (~isfield(opt,'verb') | isempty(opt.verb)), opt.verb = 1; end
if (~isfield(opt,'avoid') | isempty(opt.avoid)), opt.avoid = []; end

%initialize rng
rand('state',opt.state);

M = zeros(nc,3);
alldiff_lim = 0.5; %between 0 and 1, controls how different all colors in map are.
               %The first few colors can be quite different, high alldiff_lim
               %The difference is adjusted as more colors are demanded.
g_lim = 0.2; %limit for too gray (0-1)
d_lim = 0.40; %limit for too dim (0-3)
b_lim = 2.2; %limit for too bright  (0-3)
for (i=1:1:nc),
   M(i,:) = rand(1,3);
   cnt = 0;
   %reject if too gray or too close to previous color
   while (  toogray(M(i,:), g_lim, d_lim, b_lim) |  tooclose(M,i, 0.6, alldiff_lim) | (~isempty(opt.avoid) & (sum(abs(M(i,:)-opt.avoid)) < 0.6))),
      M(i,:) = rand(1,3);
      cnt = cnt + 1;
      if (cnt > 2000), % too tight, relax
         alldiff_lim = max([0.95.*alldiff_lim 0.02]) ;
         d_lim = max([0.95.*d_lim 0.01]);
         b_lim = min([b_lim*1.05, 8.0]);
         if (opt.verb) fprintf(1,'Reduced alldiff_lim to %g, d_lim to %g, b_lim to %g\n', alldiff_lim, d_lim, b_lim); end
         cnt = 0;
      end
   end
   if (opt.verb) fprintf(1,'Color %d OK\n', i); end
end
if (opt.verb) fprintf(1,'alldiff_lim final was %g, d_lim final was %g, b_lim final was %g\n', alldiff_lim, d_lim, b_lim); end

if (~isempty(opt.write)),
   optw.OverWrite = 'p';
   wryte3(M, opt.write, optw);
end

if (opt.show),
   ShowCmap(M, opt.show);
end

return;


function [a] = toogray(c, g_lim, d_lim, b_lim)

   a = 0;
   dc = abs(c - mean(c));
   cs = sum(c);
   if (dc(1) < g_lim & dc(2) < g_lim & dc(3) < g_lim), a = 1; return; end
   if (cs < d_lim | cs > b_lim), a = 1; return; end
   return;

function [a] = tooclose(M,i,prev_lim, alldiff_lim)

   if (i==1), a = 0; return; end

   a = 1;

   %too close to previous ?
   dc = abs(M(i,:)-M(i-1,:));
   if (sum(dc) < prev_lim), return; end

   %too close to one before?
   if (i > 2),
      for (j=1:1:i-2),
         dc = abs(M(i,:)-M(j,:));
         if (dc(1) < alldiff_lim & dc(2) < alldiff_lim & dc(3) < alldiff_lim), return; end
      end
   end
   %OK if you get here
   a = 0;
   return;
