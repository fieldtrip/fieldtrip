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

% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

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
