function [c] = ROIcol(i,nc, verb)
% [c] = ROIcol(i,nc);
% returns a color from an ROI colormap.
% useful for autocoloring of multiplots
%   i: index of color, starts at 1
%   nc: number of colors in map. Default is 32.
%       larger numbers cause longer delays at map generation.
% returns
%   c: The color RGB triplets
%
%example:
%  plot (sin([0:0.1:3]),'Color', ROIcol); hold on
%  plot (cos([0:0.1:3]),'Color', ROIcol); hold on
%see also ROIcmap
% Ziad S. Saad SSCC/NIMH/NIH, saadz@mail.nih.gov

persistent M;
persistent ic;
persistent ncl;
if (nargin == 2),
   verb = 1;
elseif (nargin == 1),
   nc = 32;
   verb = 1;
elseif (nargin == 0),
   nc = 32;
   verb = 0;
   i = [];
end

if (isempty(nc)), nc = 32; end

if (isempty(ncl)) ncl = nc; end
if (ncl ~= nc),
   M = [];
end

if (isempty(M)),
   if (verb),
       fprintf(1,'Generating colormap of %d colors.\n', nc);
   end
   opt.write = '';opt.show = 0; opt.state = 0; opt.verb = 0; opt.avoid = [1 1 1];
   M = ROIcmap(nc, opt);
end

if (isempty(i)),
   if (isempty(ic)),
      ic = 1;
   else
      ic = ic + 1;
   end
else
   ic = i;
end

ic = rem(ic,nc)+1;

c = M(ic,:);
return;
