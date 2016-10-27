function ft_realtime_headlocalizer(cfg)

% FT_REALTIME_HEADLOCALIZER is a realtime application for online visualization of the
% head position indicator (HPI) coils in CTF275 and Elekta/Neuromag systems.
%
% Repositioning within a recording session can be achieved by marking the HPI
% coil positions at an arbitrary point, i.e. by clicking the 'Update' button.
% Black unfilled markers should appear which indicate the positions of the coils
% at the moment of buttonpress. Distance to these marked positions then
% become colorcoded, i.e. green, orange, or red.
%
% Repositioning between a recording session, i.e. to a previous recording session,
% can be achieved by specifying a template; e.g. by pointing to another dataset;
% e.g. cfg.template = 'subject01xxx.ds' (CTF only), or by pointing to a textfile
% created during a previous recording; e.g. cfg.template = '29-Apr-2013-xxx.txt'.
% The latter textfile is created automatically with each 'Update' buttonpress.
%
% Use as
%   ft_realtime_headlocalizer(cfg)
% with the following configuration options
%   cfg.dataset         = string, name or location of a dataset/buffer (default = 'buffer://odin:1972')
%   cfg.template        = string, name of a template dataset for between-session repositioning (default = [])
%   cfg.bufferdata      = whether to start on the 'first or 'last' data that is available (default = 'last')
%   cfg.coilfreq        = single number in Hz or list of numbers (Neuromag default = [293, 307, 314, 321, 328])
%   cfg.blocksize       = number, size of the blocks/chuncks that are processed (default = 1 second)
%   cfg.accuracy_green  = distance from fiducial coordinate; green when within limits (default = 0.15 cm)
%   cfg.accuracy_orange = orange when within limits, red when out (default = 0.3 cm)
%
% This method is described in Stolk A, Todorovic A, Schoffelen JM, Oostenveld R.
% "Online and offline tools for head movement compensation in MEG." 
% Neuroimage. 2013 Mar;68:39-48. doi: 10.1016/j.neuroimage.2012.11.047.

% Copyright (C) 2008-2013,  Arjen Stolk & Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
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
% $Id$

agreement = {
  'By using this realtime headlocalizer tool in your research, you agree to citing the publication below.'
  ''
  'Stolk A, Todorovic A, Schoffelen JM, Oostenveld R.'
  '"Online and offline tools for head movement compensation in MEG."'
  'Neuroimage. 2013 Mar;68:39-48.'
  };

if ~strcmp(questdlg(agreement,'User agreement', 'Yes', 'Cancel', 'Cancel'), 'Yes')
  return
end

% do the general setup of the function
ft_defaults

% set the defaults
cfg.dataset         = ft_getopt(cfg, 'dataset', 'buffer://odin:1972'); % location of the buffer/dataset
cfg.accuracy_green  = ft_getopt(cfg, 'accuracy_green',  .15); % green when within this distance from reference
cfg.accuracy_orange = ft_getopt(cfg, 'accuracy_orange',  .3); % orange when within this distance from reference
cfg.template        = ft_getopt(cfg, 'template',         []); % template dataset containing the references
cfg.blocksize       = ft_getopt(cfg, 'blocksize',         1); % in seconds
cfg.bufferdata      = ft_getopt(cfg, 'bufferdata',   'last'); % first (replay) or last (real-time)
cfg.coilfreq        = ft_getopt(cfg, 'coilfreq',   [293, 307, 314, 321, 328]); % Hz, Neuromag

% ensure pesistent variables are cleared
clear ft_read_header

% start by reading the header from the realtime buffer
cfg = ft_checkconfig(cfg, 'dataset2files', 'yes'); % translate dataset into datafile+headerfile
hdr = ft_read_header(cfg.headerfile, 'cache', true, 'coordsys', 'dewar');

% determine the size of blocks to process
blocksize   = round(cfg.blocksize * hdr.Fs);
prevSample  = 0;
count       = 0;

% determine MEG system type
isneuromag = ft_senstype(hdr.grad, 'neuromag');
isctf      = ft_senstype(hdr.grad, 'ctf275');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read template head position, to reposition to, if template file is specified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isctf
  if ~isempty(cfg.template)
    [p, f, x]=fileparts(cfg.template);
    if strcmp(x, '.ds')
      shape = ft_read_headshape(cfg.template, 'coordsys', 'dewar', 'format', 'ctf_ds');
      template(1,:) = [shape.fid.pos(1,1), shape.fid.pos(1,2), shape.fid.pos(1,3)]; % chan X pos
      template(2,:) = [shape.fid.pos(2,1), shape.fid.pos(2,2), shape.fid.pos(2,3)];
      template(3,:) = [shape.fid.pos(3,1), shape.fid.pos(3,2), shape.fid.pos(3,3)];
    elseif strcmp(x, '.txt')
      template = dlmread(cfg.template);
    else
      error('incorrect template file specified');
    end
  else
    template = [];
  end
  
  % remove CTF REF sensors, for plotting purposes
  chansel = match_str(hdr.grad.chantype,'meggrad');
  hdr.grad.chanpos  = hdr.grad.chanpos(chansel,:);
  hdr.grad.chanori  = hdr.grad.chanori(chansel,:);
  hdr.grad.chantype = hdr.grad.chantype(chansel,:);
  hdr.grad.label    = hdr.grad.label(chansel,:);
  hdr.grad.tra      = hdr.grad.tra(chansel,:);
  
elseif isneuromag
  if ~isempty(cfg.template)
    template = dlmread(cfg.template);
  else
    template = [];
  end
  
else
  error('the data does not resemble ctf, nor neuromag')
end % if ctf or neuromag


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read digitized head position (for dipole fitting)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isctf
  sens = hdr.grad;
  % not needed for CTF275 systems
  dip  = [];
  vol  = [];
  coilsignal = [];
elseif isneuromag
    shape = ft_read_headshape(cfg.headerfile, 'coordsys', 'dewar', 'format', 'neuromag_fif','unit','m'); % ensure SI units
  for i = 1:min(size(shape.pos,1),length(cfg.coilfreq)) % for as many digitized or specified coils
    if ~isempty(strfind(shape.label{i},'hpi'))
      dip(i).pos = shape.pos(i,:); % chan X pos, initial guess for each of the dipole/coil positions
      dip(i).mom = [0 0 0]';
    end
  end
  if ~exist('dip', 'var')
    error('head localization requires digitized positions for Neuromag systems')
  end
  
  % prepare the forward model and the sensor array for subsequent fitting
  % note that the forward model is a magnetic dipole in an infinite vacuum
  %cfg.channel = ft_channelselection('MEG', hdr.label); % because we want to planars as well (previously only magnetometers)
  cfg.channel = ft_channelselection('MEGMAG', hdr.label); % old
  %cfg.channel = setdiff(ft_channelselection('MEG', hdr.label),ft_channelselection('MEGMAG', hdr.label)); % just trying out (planar mags)
  %cfg.channel = ft_channelselection('IAS*',hdr.label); % internal active shielding
  [vol, sens] = ft_prepare_vol_sens([], hdr.grad, 'channel', cfg.channel);
  sens = ft_datatype_sens(sens, 'version', 'upcoming', 'scaling', 'amplitude/distance', 'distance', 'm'); % ensure SI units
  coilsignal = [];
  
  % update distances, given that sensor units are m an not cm
  cfg.accuracy_green = cfg.accuracy_green/100;
  cfg.accuracy_orange = cfg.accuracy_orange/100;
else
  error('the data does not resemble ctf, nor neuromag')
end % if ctf or neuromag


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define a subset of channels for reading
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isctf
  [dum, chanindx] = match_str('headloc', hdr.chantype);
elseif isneuromag
  % 102 magnetometers
  [dum, chanindx] = match_str('megmag', hdr.chantype);
  
  % all 306
  %[dum, chanindx1] = match_str('megmag', hdr.chantype);
  %[dum, chanindx2] = match_str('megplanar', hdr.chantype); % because we want to planars as well
  %chanindx = sort([chanindx1; chanindx2],1); % because we want to planars as well
  
  % 204 planar gradiometers
  %[dum, chanindx] = match_str('megplanar', hdr.chantype);
  
  % 11 IAS
  % chanindx = 1:11;  
end

if isempty(chanindx)
  error('the data does not seem to have head localization channels');
end

% this information is passed between the GUI callback functions
info                    = [];
info.hdr                = hdr;
info.blocksize          = blocksize;
info.isctf              = isctf;
info.isneuromag         = isneuromag;
info.cfg                = cfg;
info.template           = template;
info.sens               = sens;
info.vol                = vol;
info.dip                = dip;
info.continue           = true;
clear hdr blocksize isctf isneuromag cfg template sens vol dip

% initiate main figure
hMainFig = figure;

% attach the info in the figure
guidata(hMainFig, info);

% initiate gui controls
uicontrol_sub(hMainFig);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the general BCI loop where realtime incoming data is handled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while ishandle(hMainFig) && info.continue % while the flag is one, the loop continues
  
  % get the potentially updated information from the main window
  info = guidata(hMainFig);
  
  % determine number of samples available in buffer
  info.hdr = ft_read_header(info.cfg.headerfile, 'cache', true, 'coordsys', 'dewar');
  
  % see whether new samples are available
  newsamples = (info.hdr.nSamples*info.hdr.nTrials-prevSample);
  
  if newsamples>=info.blocksize
    
    if strcmp(info.cfg.bufferdata, 'last')
      begsample  = info.hdr.nSamples*info.hdr.nTrials - info.blocksize + 1;
      endsample  = info.hdr.nSamples*info.hdr.nTrials;
    elseif strcmp(info.cfg.bufferdata, 'first')
      begsample  = prevSample + 1;
      endsample  = prevSample + info.blocksize ;
    else
      error('unsupported value for cfg.bufferdata');
    end
    
    % remember up to where the data was read
    prevSample  = endsample;
    count       = count + 1;
    fprintf('processing segment %d from sample %d to %d\n', count, begsample, endsample);
    
    % read data segment from buffer
    dat = ft_read_data(info.cfg.datafile, 'header', info.hdr, 'begsample', begsample, 'endsample', endsample, 'chanindx', chanindx, 'checkboundary', false);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % from here onward it is specific to the head localization
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % put the data in a fieldtrip-like raw structure
    data.trial{1} = double(dat);
    data.time{1}  = offset2time(begsample, info.hdr.Fs, endsample-begsample+1);
    data.label    = info.hdr.label(chanindx);
    data.hdr      = info.hdr;
    data.fsample  = info.hdr.Fs;
    
    if info.isneuromag && size(coilsignal,2)~=info.blocksize
      % construct the reference signal for each of the coils
      % this needs to be updated if the blocksize changes
      ncoil = length(info.cfg.coilfreq);
      if ncoil==0
        error('no coil frequencies were specified');
      else
        time = (1:info.blocksize)./info.hdr.Fs;
        coilsignal = zeros(ncoil, info.blocksize);
        for i=1:ncoil
          coilsignal(i,:) = exp(time*info.cfg.coilfreq(i)*1i*2*pi);
          coilsignal(i,:) = coilsignal(i,:) / norm(coilsignal(i,:));
        end
      end
    end
    
    % compute the HPI coil positions, this takes some time
    [hpi, info.dip] = data2hpi(data, info.dip, info.vol, info.sens, coilsignal, info.isctf, info.isneuromag); % for neuromag datasets this is relatively slow
    
    if ~ishandle(hMainFig)
      % the figure has been closed
      break
    end
    
    % get the potentially updated information from the main window
    info = guidata(hMainFig);
    % update the info
    info.hpi = hpi;
    % store the updated gui variables
    guidata(hMainFig, info);
    
    % DRAW LEFT PANEL - TOP VIEW
    a = subplot(1,2,1);
    h = get(a, 'children');
    hold on;
    
    if ~isempty(h)
      % done on every iteration
      delete(h);
    end
    
    % draw the color-coded head and distances from the templates
    draw_sub(hMainFig);
    
    % show current timesample
    title(sprintf('top view, runtime = %d s\n', round(mean(data.time{1}))));
    
    % not needed any more
    clear data;
    
    % viewing angle
    if info.isctf
      view(-45, 90)
    elseif info.isneuromag
      view(0, 90)
    end
    
    % DRAW RIGHT PANEL - FRONT/REAR VIEW
    b = subplot(1,2,2);
    i = get(b, 'children');
    hold on;
    
    if ~isempty(i)
      % done on every iteration
      delete(i);
    end
    
    % draw the color-coded head and distances from the templates
    draw_sub(hMainFig);
    
    % viewing angle
    if get(info.hViewRadioButton1,'Value') == 1
      if info.isctf
        view(135, 0)
      elseif info.isneuromag
        view(180, 0)
      end
      title(sprintf('anterior view, clock time %s', datestr(now))); % show current data & time
      
    elseif get(info.hViewRadioButton2,'Value') == 1
      if info.isctf
        view(-45, 0)
      elseif info.isneuromag
        view(0, 0)
      end
      title(sprintf('posterior view, clock time %s', datestr(now))); % show current data & time
    end
    
    % force Matlab to update the figure
    drawnow
    
  end % if enough new samples
  
end % while true
close(hMainFig); % close the figure


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that initiates the figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function uicontrol_sub(handle, eventdata)

% get the info
info = guidata(handle);

% initiate figure
set(handle, 'KeyPressFcn', {@key_sub});

hUpdateButton = uicontrol(...
  'Parent', handle,...
  'Style', 'pushbutton',...
  'String', 'Update',...
  'Units', 'normalized',...
  'Position', [.65 .0875 .15 .075],...
  'FontSize', 12,...
  'Callback', {@update_ButtonDownFcn});

hQuitButton = uicontrol(...
  'Parent', handle,...
  'Style', 'pushbutton',...
  'String', 'Quit',...
  'Units', 'normalized',...
  'Position', [.8 .0875 .15 .075],...
  'FontSize', 12,...
  'Callback', {@quit_ButtonDownFcn});

hCoilCheckBox = uicontrol(...
  'Parent', handle,...
  'Style', 'checkbox',...
  'String', 'Coils',...
  'Units', 'normalized',...
  'Position', [.05 .1 .075 .05],...
  'FontSize', 8,...
  'BackgroundColor', [.8 .8 .8],...
  'Value', 1,...
  'Callback', {@coil_CheckBox});

hHeadCheckBox = uicontrol(...
  'Parent', handle,...
  'Style', 'checkbox',...
  'String', 'Head',...
  'Units', 'normalized',...
  'Position', [.125 .1 .075 .05],...
  'FontSize', 8,...
  'BackgroundColor', [.8 .8 .8],...
  'Value', 1,...
  'Callback', {@head_CheckBox});

hSensorCheckBox = uicontrol(...
  'Parent', handle,...
  'Style', 'checkbox',...
  'String', 'Sensors',...
  'Units', 'normalized',...
  'Position', [.2 .1 .075 .05],...
  'FontSize', 8,...
  'BackgroundColor', [.8 .8 .8],...
  'Value', 0,...
  'Callback', {@sensor_CheckBox});

hViewRadioButton1 = uicontrol(...
  'Parent', handle,...
  'Style', 'radiobutton',...
  'String', 'Anterior view',...
  'Units', 'normalized',...
  'Position', [.275 .1 .1 .05],...
  'FontSize', 8,...
  'BackgroundColor', [.8 .8 .8],...
  'Value', 0,... % by default switched off
  'Callback', {@view_RadioButton1});

hViewRadioButton2 = uicontrol(...
  'Parent', handle,...
  'Style', 'radiobutton',...
  'String', 'Posterior view',...
  'Units', 'normalized',...
  'Position', [.375 .1 .1 .05],...
  'FontSize', 8,...
  'BackgroundColor', [.8 .8 .8],...
  'Value', 1,... % by default switched on
  'Callback', {@view_RadioButton2});

hBlocksizeMenu = uicontrol(...
  'Parent', handle,...
  'Style', 'popupmenu',...
  'String', {'.1 second','.2 second','.5 second','1 second','1.5 second','2 seconds','5 seconds','10 seconds','30 seconds'},...
  'Units', 'normalized',...
  'Position', [.475 .0925 .1 .05],...
  'FontSize', 8,...
  'BackgroundColor', [.8 .8 .8],...
  'Value', 4,... % default
  'Callback', {@blocksize_Menu});

info.hQuitButton        = hQuitButton;
info.hCoilCheckBox      = hCoilCheckBox;
info.hHeadCheckBox      = hHeadCheckBox;
info.hSensorCheckBox    = hSensorCheckBox;
info.hViewRadioButton1  = hViewRadioButton1;
info.hViewRadioButton2  = hViewRadioButton2;
info.hBlocksizeMenu     = hBlocksizeMenu;

% put the info back
guidata(handle, info);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that computes the HPI coil positions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hpi, dip] = data2hpi(data, dip, vol, sens, coilsignal, isctf, isneuromag)

% The CTF275 system localizes the HPI coil positions online, and writes them
% to the dataset. For the Neuromag systems the signals evoked by the HPI coils
% are superimposed on the other signals. This requires additional online
% dipolefitting of those HPI coils.

if isctf
  % assign the channels to the resp. coil coordinates
  [dum, x1] = match_str('HLC0011', data.label);
  [dum, y1] = match_str('HLC0012', data.label);
  [dum, z1] = match_str('HLC0013', data.label);
  [dum, x2] = match_str('HLC0021', data.label);
  [dum, y2] = match_str('HLC0022', data.label);
  [dum, z2] = match_str('HLC0023', data.label);
  [dum, x3] = match_str('HLC0031', data.label);
  [dum, y3] = match_str('HLC0032', data.label);
  [dum, z3] = match_str('HLC0033', data.label);
  
  % convert from meter to cm and assign to the resp. coil
  hpi{1} = data.trial{1}([x1 y1 z1],end) * 100;
  hpi{2} = data.trial{1}([x2 y2 z2],end) * 100;
  hpi{3} = data.trial{1}([x3 y3 z3],end) * 100;
  
elseif isneuromag
  % estimate the complex-valued MEG topography for each coil
  % this implements a discrete Fourier transform (DFT)
  topo = [];
  %[x, ut] = svdfft( data.trial{1} );
  %data.trial{1} = x;
  topo = ft_preproc_detrend(data.trial{1}) * ctranspose(coilsignal);
  
  % ignore the out-of-phase spectral component in the topography
  topo = real(topo); % THIS SEEMS TO BE CRUCIAL
  
  % fit a magnetic dipole to each of the topographies
  constr.sequential = true; % for BTI systems this would be 'false' as all coils have the same frequency
  constr.rigidbody = true;
  
  % fit the coils together
  dipall = [];
  ncoil = numel(dip);
  for i=1:ncoil
    dipall.pos(i,:) = dip(i).pos;
  end
  dipall = dipole_fit(dipall, sens, vol, topo, 'constr', constr, 'display', 'off');
  for i=1:ncoil
    sel = (1:3) + 3*(i-1);
    dip(i).pos = dipall.pos(i,:);
    dip(i).mom = real(dipall.mom(sel,i)); % ignore the complex phase information
    hpi{i} = dip(i).pos;
  end
  
else
  error('the data does not resemble ctf, nor neuromag')
end % if ctf or neuromag


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that does the timing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [time] = offset2time(offset, fsample, nsamples)

offset   = double(offset);
nsamples = double(nsamples);
time = (offset + (0:(nsamples-1)))/fsample;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which computes the circumcenter(x,y,z) of the 3D triangle (3 coils)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cc] = circumcenter(hpi)

% use coordinates relative to point `a' of the triangle
xba = hpi{2}(1) - hpi{1}(1);
yba = hpi{2}(2) - hpi{1}(2);
zba = hpi{2}(3) - hpi{1}(3);
xca = hpi{3}(1) - hpi{1}(1);
yca = hpi{3}(2) - hpi{1}(2);
zca = hpi{3}(3) - hpi{1}(3);

% squares of lengths of the edges incident to `a'
balength = xba * xba + yba * yba + zba * zba;
calength = xca * xca + yca * yca + zca * zca;

% cross product of these edges
xcrossbc = yba * zca - yca * zba;
ycrossbc = zba * xca - zca * xba;
zcrossbc = xba * yca - xca * yba;

% calculate the denominator of the formulae
denominator = 0.5 / (xcrossbc * xcrossbc + ycrossbc * ycrossbc + zcrossbc * zcrossbc);

% calculate offset (from `a') of circumcenter
xcirca = ((balength * yca - calength * yba) * zcrossbc - (balength * zca - calength * zba) * ycrossbc) * denominator;
ycirca = ((balength * zca - calength * zba) * xcrossbc - (balength * xca - calength * xba) * zcrossbc) * denominator;
zcirca = ((balength * xca - calength * xba) * ycrossbc - (balength * yca - calength * yba) * xcrossbc) * denominator;

cc(1) = xcirca + hpi{1}(1,end);
cc(2) = ycirca + hpi{1}(2,end);
cc(3) = zcirca + hpi{1}(3,end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which draws the color-coded head and distances to the template
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function draw_sub(handle)

% get the info
info = guidata(handle);

if get(info.hSensorCheckBox, 'Value') && ~isempty(info.sens)
  % plot the sensors
  hold on; ft_plot_sens(info.sens, 'style', 'k.');
end

% plot the template fiducial positions
if ~isempty(info.template)
  if info.isctf
    plot3(info.template(1,1), info.template(1,2), info.template(1,3), 'k^', 'MarkerSize', 27, 'LineWidth', 2); % chan X pos
    plot3(info.template(2,1), info.template(2,2), info.template(2,3), 'ko', 'MarkerSize', 27, 'LineWidth', 2);
    plot3(info.template(3,1), info.template(3,2), info.template(3,3), 'ko', 'MarkerSize', 27, 'LineWidth', 2);
    
    text(-8,8, info.template(2,3), 'Left', 'FontSize', 15);
    text(6,-6, info.template(3,3), 'Right', 'FontSize', 15);
  elseif info.isneuromag
    for j = 1:size(info.template,1)
      plot3(info.template(j,1), info.template(j,2), info.template(j,3), 'ko', 'MarkerSize', 27, 'LineWidth', 2); % chan X pos
    end
  end
end

% plot the HPI coil positions
for j = 1:numel(info.hpi)
  plot3(info.hpi{j}(1), info.hpi{j}(2), info.hpi{j}(3), 'ko', 'LineWidth', 1,'MarkerSize', 5)
end

if get(info.hCoilCheckBox, 'Value')
  if info.isctf
    
    % draw nasion position
    if ~isempty(info.template)
      if abs(info.template(1,1))-info.cfg.accuracy_green < abs(info.hpi{1}(1)) && abs(info.hpi{1}(1)) < abs(info.template(1,1))+info.cfg.accuracy_green ...
          && abs(info.template(1,2))-info.cfg.accuracy_green < abs(info.hpi{1}(2)) && abs(info.hpi{1}(2)) < abs(info.template(1,2))+info.cfg.accuracy_green ...
          && abs(info.template(1,3))-info.cfg.accuracy_green < abs(info.hpi{1}(3)) && abs(info.hpi{1}(3)) < abs(info.template(1,3))+info.cfg.accuracy_green
        plot3(info.hpi{1}(1),info.hpi{1}(2),info.hpi{1}(3),'g^', 'MarkerFaceColor',[.5 1 .5],'MarkerSize',25)
        head1 = true;
      elseif abs(info.template(1,1))-info.cfg.accuracy_orange < abs(info.hpi{1}(1)) && abs(info.hpi{1}(1)) < abs(info.template(1,1))+info.cfg.accuracy_orange ...
          && abs(info.template(1,2))-info.cfg.accuracy_orange < abs(info.hpi{1}(2)) && abs(info.hpi{1}(2)) < abs(info.template(1,2))+info.cfg.accuracy_orange ...
          && abs(info.template(1,3))-info.cfg.accuracy_orange < abs(info.hpi{1}(3)) && abs(info.hpi{1}(3)) < abs(info.template(1,3))+info.cfg.accuracy_orange
        plot3(info.hpi{1}(1),info.hpi{1}(2),info.hpi{1}(3),'y^', 'MarkerFaceColor',[1 .5 0],'MarkerEdgeColor',[1 .5 0],'MarkerSize',25)
        head1 = false;
      else % when not in correct position
        plot3(info.hpi{1}(1),info.hpi{1}(2), info.hpi{1}(3),'r^', 'MarkerFaceColor',[1 0 0],'MarkerSize',25)
        head1 = false;
      end
    else
      plot3(info.hpi{1}(1),info.hpi{1}(2), info.hpi{1}(3),'r^', 'MarkerFaceColor',[1 0 0],'MarkerSize',25)
      head1 = false;
    end
    
    % draw left ear position
    if ~isempty(info.template)
      if abs(info.template(2,1))-info.cfg.accuracy_green < abs(info.hpi{2}(1)) && abs(info.hpi{2}(1)) < abs(info.template(2,1))+info.cfg.accuracy_green ...
          && abs(info.template(2,2))-info.cfg.accuracy_green < abs(info.hpi{2}(2)) && abs(info.hpi{2}(2)) < abs(info.template(2,2))+info.cfg.accuracy_green ...
          && abs(info.template(2,3))-info.cfg.accuracy_green < abs(info.hpi{2}(3)) && abs(info.hpi{2}(3)) < abs(info.template(2,3))+info.cfg.accuracy_green
        plot3(info.hpi{2}(1),info.hpi{2}(2),info.hpi{2}(3),'go', 'MarkerFaceColor',[.5 1 .5],'MarkerSize',25)
        head2 = true;
      elseif abs(info.template(2,1))-info.cfg.accuracy_orange < abs(info.hpi{2}(1)) && abs(info.hpi{2}(1)) < abs(info.template(2,1))+info.cfg.accuracy_orange ...
          && abs(info.template(2,2))-info.cfg.accuracy_orange < abs(info.hpi{2}(2)) && abs(info.hpi{2}(2)) < abs(info.template(2,2))+info.cfg.accuracy_orange ...
          && abs(info.template(2,3))-info.cfg.accuracy_orange < abs(info.hpi{2}(3)) && abs(info.hpi{2}(3)) < abs(info.template(2,3))+info.cfg.accuracy_orange
        plot3(info.hpi{2}(1),info.hpi{2}(2),info.hpi{2}(3),'yo', 'MarkerFaceColor',[1 .5 0],'MarkerEdgeColor',[1 .5 0],'MarkerSize',25)
        head2 = false;
      else % when not in correct position
        plot3(info.hpi{2}(1),info.hpi{2}(2), info.hpi{2}(3),'ro', 'MarkerFaceColor',[1 0 0],'MarkerSize',25)
        head2 = false;
      end
    else
      plot3(info.hpi{2}(1),info.hpi{2}(2), info.hpi{2}(3),'ro', 'MarkerFaceColor',[1 0 0],'MarkerSize',25)
      head2 = false;
    end
    
    % draw right ear position
    if ~isempty(info.template)
      if abs(info.template(3,1))-info.cfg.accuracy_green < abs(info.hpi{3}(1)) && abs(info.hpi{3}(1)) < abs(info.template(3,1))+info.cfg.accuracy_green  ...
          && abs(info.template(3,2))-info.cfg.accuracy_green  < abs(info.hpi{3}(2)) && abs(info.hpi{3}(2)) < abs(info.template(3,2))+info.cfg.accuracy_green  ...
          && abs(info.template(3,3))-info.cfg.accuracy_green  < abs(info.hpi{3}(3)) && abs(info.hpi{3}(3)) < abs(info.template(3,3))+info.cfg.accuracy_green
        plot3(info.hpi{3}(1),info.hpi{3}(2),info.hpi{3}(3),'go', 'MarkerFaceColor',[.5 1 .5],'MarkerSize',25)
        head3 = true;
      elseif abs(info.template(3,1))-info.cfg.accuracy_orange < abs(info.hpi{3}(1)) && abs(info.hpi{3}(1)) < abs(info.template(3,1))+info.cfg.accuracy_orange ...
          && abs(info.template(3,2))-info.cfg.accuracy_orange < abs(info.hpi{3}(2)) && abs(info.hpi{3}(2)) < abs(info.template(3,2))+info.cfg.accuracy_orange ...
          && abs(info.template(3,3))-info.cfg.accuracy_orange < abs(info.hpi{3}(3)) && abs(info.hpi{3}(3)) < abs(info.template(3,3))+info.cfg.accuracy_orange
        plot3(info.hpi{3}(1),info.hpi{3}(2),info.hpi{3}(3),'yo', 'MarkerFaceColor',[1 .5 0],'MarkerEdgeColor',[1 .5 0],'MarkerSize',25)
        head3 = false;
      else % when not in correct position
        plot3(info.hpi{3}(1),info.hpi{3}(2), info.hpi{3}(3),'ro', 'MarkerFaceColor',[1 0 0],'MarkerSize',25)
        head3 = false;
      end
    else
      plot3(info.hpi{3}(1),info.hpi{3}(2), info.hpi{3}(3),'ro', 'MarkerFaceColor',[1 0 0],'MarkerSize',25)
      head3 = false;
    end
    
    if get(info.hHeadCheckBox, 'Value')
      % draw 3d head
      cc = circumcenter(info.hpi);
      x_radius = sqrt((info.hpi{2}(1) - cc(1))^2 + (info.hpi{2}(2) - cc(2))^2);
      y_radius = sqrt((info.hpi{3}(1) - cc(1))^2 + (info.hpi{3}(2) - cc(2))^2);
      [xe, ye, ze] = ellipsoid(cc(1),cc(2),cc(3),x_radius,y_radius,11);
      hh = surfl(xe, ye, ze);
      shading interp
      if get(info.hCoilCheckBox, 'Value') % this only works if 'coils' are updated
        if head1 == true && head2 == true && head3 == true
          colormap cool
        else
          colormap hot
        end
      end
      alpha(.15)
    end
    
  elseif info.isneuromag
    
    % plot fitted positions of each coil
    if ~isempty(info.template)
      for j = 1:size(info.template,1)
        if abs(info.template(j,1))-info.cfg.accuracy_green < abs(info.hpi{j}(1)) && abs(info.hpi{j}(1)) < abs(info.template(j,1))+info.cfg.accuracy_green ...
            && abs(info.template(j,2))-info.cfg.accuracy_green < abs(info.hpi{j}(2)) && abs(info.hpi{j}(2)) < abs(info.template(j,2))+info.cfg.accuracy_green ...
            && abs(info.template(j,3))-info.cfg.accuracy_green < abs(info.hpi{j}(3)) && abs(info.hpi{j}(3)) < abs(info.template(j,3))+info.cfg.accuracy_green
          plot3(info.hpi{j}(1),info.hpi{j}(2),info.hpi{j}(3),'go', 'MarkerFaceColor',[.5 1 .5],'MarkerSize',25)
        elseif abs(info.template(j,1))-info.cfg.accuracy_orange < abs(info.hpi{j}(1,end)) && abs(info.hpi{j}(1)) < abs(info.template(j,1))+info.cfg.accuracy_orange ...
            && abs(info.template(j,2))-info.cfg.accuracy_orange < abs(info.hpi{j}(2)) && abs(info.hpi{j}(2)) < abs(info.template(j,2))+info.cfg.accuracy_orange ...
            && abs(info.template(j,3))-info.cfg.accuracy_orange < abs(info.hpi{j}(3)) && abs(info.hpi{j}(3)) < abs(info.template(j,3))+info.cfg.accuracy_orange
          plot3(info.hpi{j}(1),info.hpi{j}(2),info.hpi{j}(3),'yo', 'MarkerFaceColor',[1 .5 0],'MarkerEdgeColor',[1 .5 0],'MarkerSize',25)
        else % when not in correct position
          plot3(info.hpi{j}(1,end),info.hpi{j}(2), info.hpi{j}(3),'ro', 'MarkerFaceColor',[1 0 0],'MarkerSize',25);
        end
      end
    else
      for j = 1:numel(info.hpi)
        plot3(info.hpi{j}(1),info.hpi{j}(2), info.hpi{j}(3),'ro', 'MarkerFaceColor',[1 0 0],'MarkerSize',25);
      end
    end
  end
end

% axis
grid on
xlabel('x (cm)');
ylabel('y (cm)');
zlabel('z (cm)');
set(gca, 'xtick', -10:2:10)
set(gca, 'ytick', -10:2:10)
set(gca, 'ztick', -40:2:-10) % note the different scaling
axis square

% put the info back
guidata(handle, info);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION which handles hot keys in the current plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key_sub(handle, eventdata)

% get the info
info = guidata(handle);

switch eventdata.Key
  case 'u'
    % update the template positions
    fprintf('updating template coordinates \n')
    for j = 1:numel(info.hpi)
      info.template(j,:) = info.hpi{j}(:); % chan X pos
    end
    
    % write template position to text file for later re-positioning
    template_time = [date datestr(now,'-HH-MM-SS')];
    fprintf('writing to %s.txt \n', template_time);
    dlmwrite([template_time '.txt'], info.template, ' ');
    
  case 'q'
    % stop the application
    fprintf('stopping the application \n')
    info.continue = false;
  case 'c'
    % display the sensors/dewar
    if get(info.hCoilCheckBox,'Value') == 0;
      fprintf('displaying coils \n')
      set(info.hCoilCheckBox, 'Value', 1); % toggle on
    elseif get(info.hCoilCheckBox,'Value') == 1;
      set(info.hCoilCheckBox, 'Value', 0); % toggle off
    end
  case 'h'
    % display the sensors/dewar
    if get(info.hHeadCheckBox,'Value') == 0;
      fprintf('displaying head \n')
      set(info.hHeadCheckBox, 'Value', 1); % toggle on
    elseif get(info.hHeadCheckBox,'Value') == 1;
      set(info.hHeadCheckBox, 'Value', 0); % toggle off
    end
  case 's'
    % display the sensors/dewar
    if get(info.hSensorCheckBox,'Value') == 0;
      fprintf('displaying sensors/dewar \n')
      set(info.hSensorCheckBox, 'Value', 1); % toggle on
    elseif get(info.hSensorCheckBox,'Value') == 1;
      set(info.hSensorCheckBox, 'Value', 0); % toggle off
    end
  otherwise
    fprintf('no command executed \n')
end

% put the info back
guidata(handle, info);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONs which handle button presses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function update_ButtonDownFcn(handle, eventdata)

% get the info
info = guidata(handle);

% update the template positions
fprintf('updating template coordinates \n')
for j = 1:numel(info.hpi)
  info.template(j,:) = info.hpi{j}(:); % chan X pos
end

% write template position to text file for later re-positioning
template_time = [date datestr(now,'-HH-MM-SS')];
fprintf('writing to %s.txt \n', template_time);
dlmwrite([template_time '.txt'], info.template, ' ');

% put the info back
guidata(handle, info);

function quit_ButtonDownFcn(handle, eventdata)

% get the info
info = guidata(handle);

% stop the application
fprintf('stopping the application \n')
info.continue = false;

% put the info back
guidata(handle, info);

function coil_CheckBox(hObject, eventdata)
% toggle coils display

function head_CheckBox(hObject, eventdata)
% toggle head display

function sensor_CheckBox(hObject, eventdata)
% toggle sensors display

function view_RadioButton1(handle, eventdata)

% get the info
info = guidata(handle);

% toggle front view - in combination with view_RadioButton2
if get(info.hViewRadioButton1,'Value') == 1;
  set(info.hViewRadioButton2, 'Value', 0); % toggle off radiobutton 2
  set(info.hViewRadioButton1, 'Value', 1); % toggle on radiobutton 1
elseif get(info.hViewRadioButton1,'Value') == 0;
  set(info.hViewRadioButton2, 'Value', 1); % toggle on radiobutton 2
  set(info.hViewRadioButton1, 'Value', 0); % toggle off radiobutton 1
end

% put the info back
guidata(handle, info);

function view_RadioButton2(handle, eventdata)

% get the info
info = guidata(handle);

% toggle back view - in combination with view_RadioButton1
if get(info.hViewRadioButton2,'Value') == 1;
  set(info.hViewRadioButton1, 'Value', 0); % toggle off radiobutton 1
  set(info.hViewRadioButton2, 'Value', 1); % toggle on radiobutton 2
elseif get(info.hViewRadioButton2,'Value') == 0;
  set(info.hViewRadioButton1, 'Value', 1); % toggle on radiobutton 1
  set(info.hViewRadioButton2, 'Value', 0); % toggle off radiobutton 2
end

% put the info back
guidata(handle, info);

function blocksize_Menu(handle, eventdata)

% get the info
info = guidata(handle);

val = get(info.hBlocksizeMenu, 'Value');
switch val
  case 1
    info.blocksize = round(0.1 * info.hdr.Fs); % 0.1 s
  case 2
    info.blocksize = round(0.2 * info.hdr.Fs);
  case 3
    info.blocksize = round(0.5 * info.hdr.Fs);
  case 4
    info.blocksize = round(1 * info.hdr.Fs);
  case 5
    info.blocksize = round(1.5 * info.hdr.Fs);
  case 6
    info.blocksize = round(2 * info.hdr.Fs);
  case 7
    info.blocksize = round(5 * info.hdr.Fs);
  case 8
    info.blocksize = round(10 * info.hdr.Fs);
  case 9
    info.blocksize = round(30 * info.hdr.Fs);
end
fprintf('changing blocksize to %d samples\n', info.blocksize);

% put the info back
guidata(handle, info);
