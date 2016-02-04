function [cfg] = ft_intracranplotER(cfg, data)

% FT_INTRACRANPLOTER plots results from surface/depth
% recordings of a single subject, on a surface/slices. 
% Only single values per electrode can be visualized.
% The input data should contain an elec-field, describing electrode location.
% The input data should also contain anatomy, being either single subject anatomy, 
% template anatomy, pial surface, or an MRI.
%
% Use as
%   ft_intracranplotER(cfg, data)
%
% The input data can be output of FT_TIMELOCKANALYSIS, or FT_FREQANALYSIS, and must be 
% a single value or a 1D vector per electrode (which can be averaged selectively.
%
% Two visualization methods are supported, 'point' or 'area'
%    'point' will visualize single values per electrode using a circle, which will be colored according to value.
%       The size of the circle can be set by 'size masking'
%    'area' will visualize single values using using a coloring of the surface/matter surrounding the electrode.
%       Colors will be interpolated in between electrodes. Transparency of the area can be set using masking.
% 
% The configuration can have the following parameters:
%   cfg.method        = string, visualization method, 'point', 'area', 'both'
%   cfg.channel       = nx1 cell-array with selection of channels (default = 'all') (see FT_CHANNELSELECTION)
%   cfg.parameter     = string, field in data to be visualized (default = [])
%   cfg.maskparameter = string, field in the data to be used for size/transparency masking of cfg.parameter.
%   cfg.paramsel      = scalar/1x2 vector, index or range of indices of data to select (default is to use all data, averaged)
%                       (if range, data will be averaged over the range)
%   cfg.paramlim      = 1x2 vector or 'minmax' (default), 'maxabs', 'zeromax', 'minzero', range to be used for visualization 
%   cfg.anatomypial   = pial surface (for use with SURFACE electrodes) (see FT_READ_HEADSHAPE)
%   cfg.anatomymri    = MRI (for use with DEPTH electrodes) (see FT_READ_MRI)
%   cfg.anatomyalpha  = opaqueness of the surface/MRI, between 0 (invisible) to 1 (fully opaque) (default = 1) 
%   cfg.pialcolor     = 1xN vector, color of the Pial surface (default = [.75 .75 .75])
%   cfg.colormap      = string, colormap for to-be-vizualized data, see COLORMAP (default = 'parula')
%   cfg.colorbar      = 'yes' or 'no' (default = 'yes')
%   cfg.renderer      = 'painters', 'zbuffer',' opengl' or 'none' (default = [])
%
%
% to facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
% if you specify this option the input data will be read from a *.mat
% file on disk. this mat files should contain only a single variable named 'data',
% corresponding to the input structure.
%
% See also FT_SINGLEPLOTER, FT_MULTIPLOTER, FT_SOURCEPLOT, FT_ELECTRODEPLACEMENT


% Copyright (C) 2016, Roemer van der Meij
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%



%%% TO DO LIST (non-exhaustive)
% add check for data being intracranial
% expand possible input data types
% add viewpoint manipulation
% add multiple viewpoint possibility
% add viewpoint detection (!)
% improve area interpolation
% add area interpolation controls
% add depth plotting (!)
% implement maskparameter for method = area
% add controls for mask sizing for method = point


warning('under development, appropriateness not likeley')

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar data
ft_preamble provenance data
ft_preamble trackconfig

% the abort variable is set to true or false in ft_preamble_init
if abort
  return
end

% set the defaults 
cfg.method           = ft_getopt(cfg, 'method',           'point');
cfg.parameter        = ft_getopt(cfg, 'parameter',        []);
cfg.maskparameter    = ft_getopt(cfg, 'maskparameter',    []);
cfg.paramsel         = ft_getopt(cfg, 'paramsel',         []);
cfg.paramlim         = ft_getopt(cfg, 'paramlim',         'minmax');
cfg.colormap         = ft_getopt(cfg, 'colormap',         'parula');
cfg.colorbar         = ft_getopt(cfg, 'colorbar',         'yes');
cfg.anatomypial      = ft_getopt(cfg, 'anatomypial',      []);
cfg.anatomymri       = ft_getopt(cfg, 'anatomymri',       []);
cfg.channel          = ft_getopt(cfg, 'channel',          'all');
cfg.renderer         = ft_getopt(cfg, 'renderer',         []);
cfg.pialcolor        = ft_getopt(cfg, 'pialcolor',        [.75 .75 .75]);
cfg.anatomyalpha     = ft_getopt(cfg, 'anatomyalpha',     1);
 
% ensure anatomy is given
if isempty(cfg.anatomypial) && isempty(cfg.anatomymri)
  error('need to provide either pial surface, MRI, or both')
end

% ensure data has elec
if ~isfield(data,'elec')
  error('electrode positions need to be provided in elec field')
end

% ensure that the input type is correct
data  = ft_checkdata(data, 'datatype', {'timelock', 'freq'},'dimord',{'chan_time','chan_freq'});
dtype = ft_datatype(data);

% set parameter defaults according to datatype
switch dtype
  case 'timelock'
    cfg.parameter = ft_getopt(cfg,  'parameter', 'avg');
  case 'freq'
      cfg.parameter = ft_getopt(cfg,  'parameter', 'powspctrm');
  otherwise
    % not supported
end
if ~isfield(cfg,'parameter')
  error('cfg.parameter not present as a field in the data')
end

% ensure cfg.maskparameter is valid
hasmask = false;
if ~isempty(cfg.maskparameter)
  hasmask = true;
  if numel(data.label) ~= numel(data.(cfg.maskparameter))
    error('field to be used for masking should have the same size as the number of channels in the data, and should follow that order')
  end
end

% fetch anatomy and ensure units are the same as elec
elec = data.elec;
hasanatpial = false;
hasanatmri  = false;
if ~isempty(cfg.anatomypial)
  hasanatpial = true;
  pial        = cfg.anatomypial;
  pial        = ft_convert_units(pial,elec.unit);
end
if ~isempty(cfg.anatomymri)
  hasanatmri = true;
  mri        = cfg.anatomymri;
  mri        = ft_convert_units(mri,elec.unit);
end


% select channels
if isfield(cfg, 'channel') && isfield(data, 'label')
  cfg.channel = ft_channelselection(cfg.channel, data.label);
end
if isempty(cfg.channel)
  error('no channels selected');
end
selchan = cfg.channel;


% determine electrode type, surf vs depth
selind   = match_str(elec.label,selchan);
electype = elec.chantype(selind);
if any(strcmp(electype,'surf'))
  hassurfelec = true;
else
  hassurfelec = false;
end
if any(strcmp(electype,'depth'))
  hasdepthelec = true;
else
  hasdepthelec = false;
end
if hassurfelec && ~hasanatpial
  error('cannot plot surface electrodes without surface (pial) anatomy')
end
if ~hassurfelec && hasanatpial
  warning('surface (pial) anatomy given but no surface electrodes will be plotted')
end
if hasdepthelec && ~hasanatmri
  error('cannot plot depth electrodes without MRI')
end
if ~hasdepthelec && hasanatmri
  warning('MRI given but no depth electrodes will be plotted')
end
if ~any(strcmp(electype,'surf')) && ~any(strcmp(electype,'depth'))
  error('electrodes of unknown type detected, remove electrodes from cfg.channel and/or confirm elec field')
end


% fetch electrode locations
selind  = match_str(elec.label,selchan);
if numel(selind)~=numel(selchan)
  error('not all selected electrodes have a position defined in the elec field')
end
elecpos = elec.elecpos(selind,:);


% fetch data and select/average when applicable (dimensionality checked using ft_checkdata above)
dat = data.(cfg.parameter);
% select channels and extract mask
selind  = match_str(data.label, selchan);
dat     = dat(selind,:);
if hasmask
  mask   = data.(cfg.maskparameter);
  mask   = mask(selind); % size checked above
end
% apply cfg.paramsel if needed
if size(dat,2)>1
  if isempty(cfg.paramsel)
    dat = nanmean(dat,2);
  else
    if numel(cfg.paramsel)==2
      dat = nanmean(dat(:,cfg.paramsel(1):cfg.paramsel(2)),2);
    elseif numel(cfg.paramsel)==1
      dat = dat(:,cfg.paramsel);
    else
      error('cfg.paramsel needs to be either a scaler or a 1x2 vector')
    end
  end
end
% check for appropriateness of mask in case of method = area
if hasmask && any(strcmp(cfg.method,{'area','both'})) && (min(mask)<0 || max(mask)>1)
  error('masking values should be between 0 and 1 when cfg.method = ''area'' or ''both''')
end

% determine data limits and convert dat to color index values
switch cfg.paramlim
  case 'minmax'
    cmin = min(dat(:));
    cmax = max(dat(:));
  case 'maxabs'
    cmin = -max(abs(dat(:)));
    cmax = max(abs(dat(:)));
  case 'zeromax'
    cmin = 0;
    cmax = max(dat(:));
  case 'minzero'
    cmin = min(dat(:));
    cmax = 0;
  otherwise
    if isnumeric(cfg.paramlim)
      cmin = cfg.paramlim(1);
      cmax = cfg.paramlim(2);
    else
      error('unsupported option for cfg.paramlim')
    end
end
if (cmax-cmin)<0
  error('inappropriate cfg.paramlim given the data')
end
% create colormat based on colormap (make it huge for indexing precision)
colormat = eval([cfg.colormap '(1e5);']); % damn evil eval
% transform dat to go from 0 1 and multiply with 10e4 to act as index into colormat
datcind = round(((dat + -cmin) .* (1 / (-cmin + cmax))) .* 10e4);
datcind(datcind==0) = 1; % if a color index is zero, zet it to the first color-value (which has index 1)
datcind(datcind<0) = 1;
datcind(datcind>1e5) = 1e5;
% % convert values at the extremes to gray
% colormat(1,:) = [.5 .5 .5];
% colormat(1e5,:) = [.5 .5 .5];


% set the figure window title suffix
funcname = mfilename;
if isfield(cfg, 'dataname')
  dataname = cfg.dataname;
elseif nargin > 1
  dataname = inputname(2);
else % data provided through cfg.inputfile
  dataname = cfg.inputfile;
end
fignamesuffix = [funcname '_' dataname];




if hassurfelec
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Plot SURFACE electrodes
  
  % select only surface electrodes
  selind   = match_str(elec.label,selchan);
  electype = elec.chantype(selind);
  selind   = find(strcmp(electype,'surf'));
  nselec   = numel(selind);
  sdatcind = datcind(selind);
  selecpos = elecpos(selind,:);
  if hasmask
    smask  = mask(selind);
  end
  
  % create figure
  figure('name',['surf_' fignamesuffix],'numbertitle','off')
  hold on
  
  % switch between display methods
  switch cfg.method
       
    
    case 'point'
      
      % first, plot surface
      ft_plot_mesh(pial, 'facecolor',cfg.pialcolor,'edgecolor','none','edgealpha',0);
      lighting phong
      material dull
      set(gca,'DataAspectRatio',[1 1 1], 'PlotBoxAspectRatio', [1 1 1]);
      alpha(cfg.anatomyalpha)
      % set viewpoint
      %view([0 0])
      % set camlight (this is always done w.r.t. to the current camera position, hence the name...)
      camlight(0,0);
      
      % then, plot electrodes loop over channels, as plot3 cannot take a color matrix as input
      for ichan = 1:nselec
        if hasmask
          plot3(selecpos(ichan,1),selecpos(ichan,2),selecpos(ichan,3),'o','markerfacecolor',colormat(sdatcind(ichan),:),'markeredgecolor',[0 0 0],'markersize',smask(ichan));
        else
          plot3(selecpos(ichan,1),selecpos(ichan,2),selecpos(ichan,3),'o','markerfacecolor',colormat(sdatcind(ichan),:),'markeredgecolor',[0 0 0]);
        end
      end
      axis off 
      
      
    case {'area','both'}
    
      % first, calculate electrode spacing
      distmat = NaN(nselec,nselec);
      for iseed = 1:nselec
        xs = selecpos(iseed,1);
        ys = selecpos(iseed,2);
        zs = selecpos(iseed,3);
        for itarg = 1:nselec
          xt = selecpos(itarg,1);
          yt = selecpos(itarg,2);
          zt = selecpos(itarg,3);
          % calculate euclidian distance
          distmat(iseed,itarg) = ((xs-xt)^2 + (ys-yt)^2 + (zs-zt) ^2 ) ^ 0.5;
        end
      end
      distmat(logical(eye(nselec))) = inf;
      elecspacing = mean(min(distmat));
      
      %       % then, calculate the distance between each electrode, and all vertices
      %       vertelecdist = NaN(nselec,size(pial.pos,1)); %
      %       for ichan = 1:nselec
      %         xs = selecpos(ichan,1);
      %         ys = selecpos(ichan,2);
      %         zs = selecpos(ichan,3);
      %         vertelecdist(ichan,:) = sqrt((pial.pos(:,1) - xs).^2 + (pial.pos(:,2) - ys).^2 + (pial.pos(:,3) - zs).^2);
      %       end
      %       % then, compute the color of each vertex as a function of the distance to all electrodes
      %       col = datcind.' * (1./vertelecdist);
      %       % the min of col should be 0, the max (max(datcind)-1)
      %       col = col - min(col);
      %       col = (col ./ max(col)) .* (max(datcind)-1);
      %       col = col + 1; % col is an index
      
      % % then, compute the color of each vertex as a function using a sphere surrounding the electrode
      %col = interp_ungridded(selecpos,pial.pos,'data',sdatcind,'projmethod','sphere_avg','sphereradius',2*elecspacing);
      col = interp_ungridded(selecpos,pial.pos,'data',sdatcind,'projmethod','nearest');
      % the min of col should be 0, the max (max(datcind)-1)
      col = col - min(col);
      col = (col ./ max(col)) .* (max(datcind)-1);
      col = col + 1; % col is an index
      
      % then, index the colors based on a distance limit
      vertcol = repmat(cfg.pialcolor,size(pial.pos,1),1);
      for ivert = 1:size(vertcol,1)
        % detect whether it falls within patchdistlim
        if all(sum(abs(repmat(pial.pos(ivert,:),[nselec,1]) - selecpos) < elecspacing))
          % index the color with col
          currcol = round(col(ivert));
          if currcol==0
            currcol = 1; % if color index is zero, set it to the first color-value (which has index 1)
          end
          vertcol(ivert,:) = colormat(currcol,:);
        end
      end
      
      % then, finally, plot the surface
      ft_plot_mesh(pial, 'facecolor',cfg.pialcolor,'vertexcolor',vertcol,'edgecolor','none','edgealpha',0);
      lighting phong
      material dull
      set(gca,'DataAspectRatio',[1 1 1], 'PlotBoxAspectRatio', [1 1 1]);
      alpha(cfg.anatomyalpha)
      % set viewpoint
      %view([0 0])
      % set camlight (this is always done w.r.t. to the current camera position, hence the name...)
      camlight(0,0);
      
      if strcmp(cfg.method,'both')
        for ichan = 1:nselec
          if hasmask
            plot3(selecpos(ichan,1),selecpos(ichan,2),selecpos(ichan,3),'o','markerfacecolor',colormat(sdatcind(ichan),:),'markeredgecolor',[0 0 0],'markersize',smask(ichan));
          else
            plot3(selecpos(ichan,1),selecpos(ichan,2),selecpos(ichan,3),'o','markerfacecolor',colormat(sdatcind(ichan),:),'markeredgecolor',[0 0 0]);
          end
        end
        axis off
      end
      
      
    otherwise
      error('unsupported cfg.method')
  end
  
  % colorbar
  if strcmp(cfg.colorbar,'yes')
    colorbar
    caxis([cmin cmax])
  end
  
  % set renderer if specified
  if ~isempty(cfg.renderer)
    set(gcf, 'renderer', cfg.renderer)
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end




if hasdepthelec
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Plot DEPTH electrodes
  
  error('TO BE IMPLEMENTED')
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end



% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous data
ft_postamble provenance
