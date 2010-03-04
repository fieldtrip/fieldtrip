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
freq = checkdata(freq, 'cmbrepresentation', 'sparse');

% check if the input configuration is valid for this function
cfg = checkconfig(cfg, 'trackconfig', 'on');
cfg = checkconfig(cfg, 'required', {'foi', 'layout'});

% set the defaults
if ~isfield(cfg, 'feedback'),   cfg.feedback = 'text';        end
if ~isfield(cfg, 'alphaparam'), cfg.alphaparam = [];          end
if ~isfield(cfg, 'widthparam'), cfg.widthparam = [];          end
if ~isfield(cfg, 'colorparam'), cfg.colorparam = 'cohspctrm'; end
if ~isfield(cfg, 'newfigure'),  cfg.newfigure = 'yes';        end

if ~isfield(cfg, 'arrowhead'),   cfg.arrowhead = 'no';         end % no, stop, start, both
if ~isfield(cfg, 'arrowlength'), cfg.arrowlength = 1;          end % relative to the complete line
if ~isfield(cfg, 'arrowoffset'), cfg.arrowoffset = 0;          end % absolute, should be in the same units as the layout

lay = prepare_layout(cfg, freq);

beglabel = freq.labelcmb(:,1);
endlabel = freq.labelcmb(:,2);
ncmb     = size(freq.labelcmb,1);

% select the data to be used in the figure
fbin = nearest(freq.freq, cfg.foi);

if isfield(freq, cfg.widthparam)
  widthparam = freq.(cfg.widthparam)(:,fbin);
else
  widthparam = ones(ncmb,1);
end

if isfield(freq, cfg.alphaparam)
  alphaparam = freq.(cfg.alphaparam)(:,fbin);
else
  alphaparam = [];
end

if isfield(freq, cfg.colorparam)
  colorparam = freq.(cfg.colorparam)(:,:,fbin);
else
  colorparam = [];
end

if strcmp(cfg.newfigure, 'yes')
  figure
end

hold on

rgb  = colormap;
cmin = min(colorparam(:));
cmax = max(colorparam(:));
colorparam = (colorparam - cmin)./(cmax-cmin);
colorparam = round(colorparam * (size(rgb,1)-1) + 1);

if strcmp(cfg.newfigure, 'yes')
  % also plot the position of the electrodes
  plot_vector(lay.pos(:,1), lay.pos(:,2), 'style','k.');

  % also plot the outline, i.e. head shape or sulci
  if isfield(lay, 'outline')
    fprintf('solid lines indicate the outline, e.g. head shape or sulci\n');
    for i=1:length(lay.outline)
      if ~isempty(lay.outline{i})
        X = lay.outline{i}(:,1);
        Y = lay.outline{i}(:,2);
        plot_line(X, Y, 'color', 'k', 'linewidth', 1.5, 'linestyle', '-');
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
        plot_line(X, Y, 'color', 'k', 'linewidth', 1.5, 'linestyle', '-');
      end
    end
  end
end % if newfigure

% fix the limits for the axis
axis(axis);

progress('init', cfg.feedback, 'plotting connections...');

for i=1:ncmb

  if strcmp(beglabel{i}, endlabel{i})
    % skip autocombinations
    continue
  end

  progress(i/ncmb, 'plotting connection %d from %d (%s -> %s)\n', i, ncmb, beglabel{i}, endlabel{i});

  if widthparam(i)>0
    begindx = strmatch(beglabel{i}, lay.label);
    endindx = strmatch(endlabel{i}, lay.label);
    xbeg = lay.pos(begindx,1);
    ybeg = lay.pos(begindx,2);
    xend = lay.pos(endindx,1);
    yend = lay.pos(endindx,2);
    
    if strcmp(cfg.arrowhead, 'no')
      x = [xbeg xend]';
      y = [ybeg yend]';
      % h = line(x, y);
      h = patch(x, y, 1);
    else
      arrowbeg  = [xbeg ybeg];
      arrowend  = [xend yend];
      center    = (arrowbeg+arrowend)/2;
      direction = (arrowend - arrowbeg);
      direction = direction/norm(direction);
      offset    = [direction(2) -direction(1)];
      arrowbeg  = cfg.arrowlength * (arrowbeg-center) + center + cfg.arrowoffset * offset;
      arrowend  = cfg.arrowlength * (arrowend-center) + center + cfg.arrowoffset * offset;

      switch cfg.arrowhead
        case {'yes' 'stop'}
          h = arrow(arrowbeg, arrowend, 'Ends', 'stop');
        case 'start'
          h = arrow(arrowbeg, arrowend, 'Ends', 'start');
        case 'both'
          h = arrow(arrowbeg, arrowend, 'Ends', 'both');
        case 'none'
          h = arrow(arrowbeg, arrowend, 'Ends', 'none');
        otherwise
          error('unsupported value for cfg.arrowhead')
      end % switch arrowhead
    end % if arrow

    if ~isempty(widthparam)
      set(h, 'LineWidth', widthparam(i));
    end

    if ~isempty(alphaparam)
      set(h, 'EdgeAlpha', alphaparam(i));
    end

    if ~isempty(colorparam)
      set(h, 'EdgeColor', rgb(colorparam(i),:));
    end

  end
end
progress('close');

% improve the fit in the axis
axis tight

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');

if nargout<1
  clear cfg
end
