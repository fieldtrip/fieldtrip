function cfg = topoplotCC(cfg, freq)

% TOPOPLOTCC plots the connections between significantly coherent
% sensor pairs
%
% Use as
%  topoplotCC(cfg, freq)
%
% The configuration should contain:
%   cfg.feedback    = string (default = 'textbar')
%   cfg.layout      = specification of the layout, see PREPARE_LAYOUT
%   cfg.foi         = the frequency of interest which is to be plotted (default is the first frequency bin)
%   cfg.widthparam  = string, parameter to be used to control the line width
%   cfg.alphaparam  = string, parameter to be used to control the opacity
%   cfg.colorparam  = string, parameter to be used to control the line color
%
% The widthparam should be indicated in pixels, e.g. usefull numbers are 1
% and larger.
%
% The alphaparam should be indicated as opacity between 0 (fully transparent) 
% and 1 (fully opaque).
%
% See also: PREPARE_LAYOUT, MULTIPLOTCC

% $Log: topoplotCC.m,v $
% Revision 1.11  2009/06/15 12:57:39  roboos
% added checkdata and checkconfig calls
% changed default colorparam into cohspctrm
% changed default width into 1
% updated documentation
%
% Revision 1.10  2009/04/15 13:07:14  roboos
% largely rewritten the whole function to work with dimord=chan_chan_freq
% now also with specification of the width, opacity and color for each line
%
% Revision 1.9  2008/09/22 20:17:44  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.8  2007/03/21 15:44:05  chrhes
% updated documentation regarding the fact that cfg.layout can also contain a
% layout structure obtained using the function prepare_layout.m
%
% Revision 1.7  2007/03/14 08:43:12  roboos
% replaced call to createlayout to prepare_layout, made some small changes to the lay structure
%
% Revision 1.6  2006/06/26 10:13:15  ingnie
% now works also if no stat.freq present
%
% Revision 1.5  2006/04/20 09:58:34  roboos
% updated documentation
%
% Revision 1.4  2006/03/09 08:20:54  roboos
% changed freq.foi into freq.freq
%
% Revision 1.3  2005/10/01 16:19:04  jansch
% added some help-information. have fun!
%
% Revision 1.2  2005/09/05 06:40:44  jansch
% added defaults in configuration
%
% Revision 1.1  2005/08/25 10:59:37  jansch
% New implementation
%

fieldtripdefs

% check if the input data is valid for this function
freq = checkdata(freq, 'cmbrepresentation', 'full');

% check if the input configuration is valid for this function
cfg = checkconfig(cfg, 'trackconfig', 'on');
cfg = checkconfig(cfg, 'required', {'foi', 'layout'});

% set the defaults
if ~isfield(cfg, 'feedback'), cfg.feedback = 'textbar'; end
if ~isfield(cfg, 'alphaparam'), cfg.alphaparam = []; end
if ~isfield(cfg, 'widthparam'), cfg.widthparam = []; end
if ~isfield(cfg, 'colorparam'), cfg.colorparam = 'cohspctrm'; end

lay = prepare_layout(cfg, freq);

reflabel  = freq.label;
chanlabel = freq.label;
refindx   = match_str(lay.label, reflabel);
chanindx  = match_str(lay.label, chanlabel);
nref      = length(refindx);
nchan     = length(chanindx);

% select the data to be used in the figure
fbin = nearest(freq.freq, cfg.foi);

if isfield(freq, cfg.widthparam)
  widthparam = freq.(cfg.widthparam)(:,:,fbin);
else
  widthparam = ones(nref,nchan);
end

if isfield(freq, cfg.alphaparam)
  alphaparam = freq.(cfg.alphaparam)(:,:,fbin);
else
  alphaparam = [];
end

if isfield(freq, cfg.colorparam)
  colorparam = freq.(cfg.colorparam)(:,:,fbin);
else
  colorparam = [];
end

k = 0;
progress('init', cfg.feedback, 'plotting connections...');

figure
hold on

rgb  = colormap;
cmin = min(colorparam(:));
cmax = max(colorparam(:));
colorparam = (colorparam - cmin)./(cmax-cmin);
colorparam = round(colorparam * (size(rgb,1)-1) + 1);

for i=1:nref
  for j=1:nchan
    if i==j
      % don't plot the connection with itself
      continue
    end
    k = k+1;

    progress(k/(nref*nchan), 'plotting connection %d from %d (%d -> %d)\n',k, nref*nchan, i, j);

    if widthparam(i,j)>0
      xbeg = lay.pos(refindx(i),1);
      ybeg = lay.pos(refindx(i),2);
      xend = lay.pos(chanindx(j),1);
      yend = lay.pos(chanindx(j),2);
      x = [xbeg xend]';
      y = [ybeg yend]';
      % h = line(x, y);
      h = patch(x, y, 1);

      if ~isempty(widthparam)
        set(h, 'LineWidth', widthparam(i,j));
      end

      if ~isempty(alphaparam)
        set(h, 'EdgeAlpha', alphaparam(i,j));
      end

      if ~isempty(colorparam)
        set(h, 'EdgeColor', rgb(colorparam(i,j),:));
      end

    end
  end
end
progress('close');

% also plot the position of the electrodes
plot(lay.pos(:,1), lay.pos(:,2), 'k.');

% also plot the outline, i.e. head shape or sulci
if isfield(lay, 'outline')
  fprintf('solid lines indicate the outline, e.g. head shape or sulci\n');
  for i=1:length(lay.outline)
    if ~isempty(lay.outline{i})
      X = lay.outline{i}(:,1);
      Y = lay.outline{i}(:,2);
      h = line(X, Y);
      set(h, 'color', 'k');
      set(h, 'linewidth', 1.5);
      set(h, 'linestyle', '-');
    end
  end
end

% also plot the mask, i.e. global outline for masking the topoplot
if isfield(lay, 'mask')
  fprintf('dashed lines indicate the mask for topograpic interpolation\n');
  for i=1:length(lay.mask)
    if ~isempty(lay.mask{i})
      X = lay.mask{i}(:,1);
      Y = lay.mask{i}(:,2);
      % the polygon representing the mask should be closed
      X(end+1) = X(1);
      Y(end+1) = Y(1);
      h = line(X, Y);
      set(h, 'color', 'k');
      set(h, 'linewidth', 1.5);
      set(h, 'linestyle', '-');
    end
  end
end

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');

if nargout<1
  clear cfg
end
