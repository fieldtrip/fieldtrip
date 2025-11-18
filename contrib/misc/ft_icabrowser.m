function [rej_comp, artifact] = ft_icabrowser(cfg, comp)

% FT_ICABROWSER takes as an input comp structure from FieldTrip ft_componentanalysis
% and presents a GUI interface showing the power spectrum, variance over
% time and the topography of the components, as well as the possibility to
% save a PDF, inspect the timecourse and toggle components to be rejected vs
% kept. Also, it allows for the specification of artifactual segments,
% through the functionality of ft_databrowser. This may be handy in order
% to identify data segments that cause infrequent data features loading on
% a limited set of channels (e.g. isolated SQUID jumps in MEG). This type
% of artitfact defendibly would require the extreme segments to be rejected
% from the data a priori, i.e. before application of the independent
% component analysis.
%
% Use as
%    [rej_comp, artifact] = ft_icabrowser(cfg, comp)
%
% where the input comp structure should be obtained from FT_COMPONENTANALYSIS.
%
% The configuration must contain:
%   cfg.layout     = filename of the layout, see FT_PREPARE_LAYOUT
%
% further optional configuration parameters are
%   cfg.rejcomp       = list of components which shall be initially marked for rejection, e.g. [1 4 7]
%   cfg.blocksize     = blocksize of time course length for visualization (default = [])
%   cfg.powscale      = scaling of y axis in power plot, 'lin' or 'log10', (default = 'log10')
%   cfg.freqscale     = scaling of x axis in power plot, 'lin' or 'log', (default = 'lin')
%   cfg.foilim
%   cfg.zlim          = plotting limits for color dimension of topoplot, 'maxmin', 'maxabs', 'zeromax', 'minzero', or [zmin zmax] (default = 'maxmin')
%   cfg.outputfolder  = where pdfs will be saved (default = pwd)
%   cfg.prefix        = prefix of the pdf files (default = 'ICA')
%   cfg.colormap      = any sized colormap, see FT_COLORMAP
%   cfg.outputfile    = MAT file which contains indices of all components to reject
%   cfg.showcallinfo  = show call info, 'yes' or 'no' (default: 'no')
%   cfg.chunklength
%
% original written by Thomas Pfeffer
% adapted by Jonathan Daume and Anne Urai
% University Medical Center Hamburg-Eppendorf, 2015
%
% modified by Daniel Matthes
% Max Planck Institute for Human Cognitive and Brain Sciences, 2019
%
% Jan-Mathijs did a big overhaul of the functionality, allowing for more
% interaction. Future versions may lose some functionality, e.g. w.r.t
% saving pdfs.
%
% See also FT_COMPONENTANALYSIS, FT_TOPOPLOTIC, FT_PREPARE_LAYOUT, FT_DATABROWSER

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if the input data is valid for this function
comp = ft_checkdata(comp, 'datatype', 'raw+comp', 'hassampleinfo', 'yes');

% set the defaults
layout        = ft_getopt(cfg, 'layout');
rejcomp       = ft_getopt(cfg, 'rejcomp',   []);
blocksize     = ft_getopt(cfg, 'blocksize', []);
powscale      = ft_getopt(cfg, 'powscale',  'log10');
freqscale     = ft_getopt(cfg, 'freqscale', 'lin');
cfg.zlim      = ft_getopt(cfg, 'zlim',      'maxmin');
outputfolder  = ft_getopt(cfg, 'outputfolder',      pwd);
prefix        = ft_getopt(cfg, 'prefix',    'ICA');
colormap      = ft_getopt(cfg, 'colormap',  'default');
outputfile    = ft_getopt(cfg, 'outputfile',   []);
showcallinfo  = ft_getopt(cfg, 'showcallinfo', 'no');
cfg.chunklength = ft_getopt(cfg, 'chunklength', 'trial');
cfg.foilim   = ft_getopt(cfg, 'foilim', []);
cfg.colormap = ft_getopt(cfg, 'colormap', []);

% cfg.rejcomp can be indicated by number or by label
rejcomp = ft_channelselection(rejcomp, comp.label);
reject = match_str(comp.label, rejcomp);

% create folder if not exist
if ~exist(outputfolder, 'dir')
  mkdir outputfolder;
end

% test config
if ~ismember(powscale, {'lin', 'log10'})
  powscale = 'log10';
end

if ~ismember(freqscale, {'lin', 'log'})
  freqscale = 'lin';
end

if numel(comp.trial)==1 && isequal(cfg.chunklength, 'trial')
  fprintf('Using a default chunk length = 2 seconds\n');
  cfg.chunklength = 2;
end

if ~isequal(cfg.chunklength, 'trial')
  % NEW: precompute a ncomponent x nchunk matrix of variance values, for the data cut into 2 second chunks
  fprintf('Computing the variance in %d-second chunks\n', cfg.chunklength);
  cfgredef         = [];
  cfgredef.length  = cfg.chunklength;
  cfgredef.keeppartial = 'yes';
  cfgredef.updatetrialinfo = 'yes';
  newcomp  = ft_redefinetrial(cfgredef, comp);
else
  newcomp = comp;
end
trialinfo = (1:numel(newcomp.trial))';

if istable(trialinfo)
  trialinfo = table2array(trialinfo);
end
comp_var = zeros(numel(comp.label), numel(newcomp.trial));
for s = 1:numel(newcomp.trial)
  comp_var(:,s) = var(newcomp.trial{s}, [], 2, 'omitnan');
end
comp_time = 1:size(comp_var,2);

% NEW: precompute the spectra
fprintf('Computing power spectra\n');
hasnan = false(numel(newcomp.trial),1);
for s = 1:numel(newcomp.trial)
  hasnan(s) = any(~isfinite(newcomp.trial{s}(:)));
end
cfgfreq = [];
cfgfreq.method = 'mtmfft';
cfgfreq.taper  = 'hanning';
cfgfreq.pad    = 10;
cfgfreq.trials = find(~hasnan);
if ~isempty(cfg.foilim)
  cfgfreq.foilim = cfg.foilim;
end
freq = ft_freqanalysis(cfgfreq, newcomp);

% preallocate rejected components
rej_comp = false(size(comp.label,1),1);
rej_comp(reject) = true;

numOfPlots = 4;                            % number of subplots per page
page = 1;                                  % number of current page
numOfPages = ceil(size(comp.label, 1)/4);  % maximum number of pages

% to save time redoing this for each topo
cfglay              = [];
cfglay.layout       = layout;
cfglay.showcallinfo = showcallinfo;
lay                 = ft_prepare_layout(cfglay);

% match the labels in the layout with the labels in the data
[i1, i2] = match_str(lay.label, comp.topolabel);
xpos = lay.pos(i1,1);
ypos = lay.pos(i1,2);

cfgtopo               = [];
cfgtopo.layout        = lay;                                              
cfgtopo.comment       = 'no';
cfgtopo.zlim          = cfg.zlim;
cfgtopo.highlight     = 'off';
cfgtopo.marker        = 'off';
cfgtopo.style         = 'straight';
cfgtopo.showcallinfo  = showcallinfo;
if ~isempty(colormap)
  cfgtopo.colormap    = ft_colormap(cfg.colormap);
  close;
end

err      = 0;
manpos   = [0.1 0.1 0.8 0.8]; % figure position, can be updated later
artifact = zeros(0,2);

rej_str = {'Keep' 'Reject'};
rej_col = [0 0.5 0; 0.5 0 0];

f = figure('units','normalized','outerposition', manpos, 'CloseRequestFcn', @(h, evt)quitme);
set(f, 'Name', sprintf('%d: %s', double(f), 'ft_icabrowser'));

% ------------------------------------------------
% FIRST/PREVIOUS/NEXT/LAST PAGE BUTTONS
% ------------------------------------------------
first = uicontrol('Units','normalized','Position',[0.05 0.01 0.075 0.05],'Style','pushbutton','String','First','Callback',@(h, evt)firstpage);
first.Enable = 'off';

prev = uicontrol('Units','normalized','Position',[0.15 0.01 0.075 0.05],'Style','pushbutton','String','Prev','Callback',@(h, evt)prevpage);
prev.Enable = 'off';

next = uicontrol('Units','normalized','Position',[0.25 0.01 0.075 0.05],'Style','pushbutton','String','Next','Callback',@(h, evt)nextpage);
next.Enable = 'off';

last = uicontrol('Units','normalized','Position',[0.35 0.01 0.075 0.05],'Style','pushbutton','String','Last','Callback',@(h, evt)lastpage);
last.Enable = 'off';

% ------------------------------------------------
% SWITCH POWER AXIS LOG/LINEAR
% ------------------------------------------------
powlog = uicontrol('Units','normalized','Position',[0.45 0.01 0.075 0.05],'Style','pushbutton', 'Callback',@(h, evt)powplotlog);
powlog.Enable = 'off';
if strcmp(powscale, 'log10')
  powlog.String = 'Linear Power Scale';
else
  powlog.String = 'Log Power Scale';
end

% ------------------------------------------------
% SWITCH FREQUENCY AXIS LOG/LINEAR
% ------------------------------------------------
freqlog = uicontrol('Units','normalized','Position',[0.55 0.01 0.090 0.05],'Style','pushbutton', 'Callback',@(h, evt)freqplotlog);
freqlog.Enable = 'off';
if strcmp(freqscale, 'log')
  freqlog.String = 'Linear Frequency Scale';
else
  freqlog.String = 'Log Frequency Scale';
end

% ------------------------------------------------
% SAVE AND QUIT
% ------------------------------------------------
save_it = uicontrol('Units','normalized','Position',[0.80 0.01 0.075 0.05], ...
  'Style','pushbutton','String','Save','Callback',@(h, evt)save_callback);
save_it.Enable = 'off';
quit_it = uicontrol('Units','normalized','Position',[0.90 0.01 0.075 0.05],...
  'Style','pushbutton','String','Quit','Callback',@(h, evt)quitme);
quit_it.Enable = 'off';

orient landscape
drawnow;

trial_for_disp = 1;
subcomp = cell(4,numOfPlots);
while err == 0 % KEEP GOING UNTIL THERE IS AN ERROR

  for row = 1:numOfPlots % il is the subplot count
    freqlog.Enable = 'off';
    powlog.Enable  = 'off';
      
    compNum = (page-1) * numOfPlots + row; % number of current component
   
    % ------------------------------------------------
    % PLOT FIRST TRIAL
    % ------------------------------------------------
    if isempty(subcomp{1,row})
      subcomp{1,row} = subplot(numOfPlots,4,row*4-3, 'tag', sprintf('data_%s', comp.label{compNum}));
      xlabel('Time (s)');
      plot(comp.time{trial_for_disp}, comp.trial{trial_for_disp}(compNum,:));
      ylabel(comp.label{compNum});
      if row==1
        set(gca, 'XaxisLocation', 'top');
        xlabel(sprintf('trial %d', trial_for_disp));
      end
      axis tight;
    else
      l = get(subcomp{1,row}, 'children');
      set(l, 'XData', comp.time{trial_for_disp}, 'YData', comp.trial{trial_for_disp}(compNum,:));
      set(subcomp{1,row}, 'tag', sprintf('data_%s', comp.label{compNum}));
      set(get(subcomp{1,row}, 'ylabel'), 'string', comp.label{compNum});
      if row==1
        set(get(subcomp{1,row}, 'xlabel'), 'string', sprintf('trial %d', trial_for_disp));
      end
    end

    % ------------------------------------------------
    % PLOT POWER SPECTRUM
    % ------------------------------------------------
    strt = find(freq.freq > 1,1,'first');
    stp  = find(freq.freq < 200,1,'last');
    subcomp{2,row} = subplot(numOfPlots,4,row*4-2, 'tag', sprintf('freq_%s', comp.label{compNum}));
    xlabel('Frequency (Hz)'); grid on;
    
    if strcmp(powscale, 'log10')
      dat = log10(freq.powspctrm(compNum, strt:stp));
      plot(freq.freq(strt:stp), dat);
      set(gca, 'ylim', max(dat)+[-2 0], 'xlim', freq.freq([strt stp]));
      ylabel('PSD (dB/Hz)');
    else
      plot(freq.freq(strt:stp),freq.powspctrm(compNum, strt:stp));
      ylabel('PSD (T^2/Hz)');
    end

    if strcmp(freqscale, 'log')
      set(gca, 'XScale', 'log');
      set(gca,'TickDir','out','XTick', [1 10 100]);
    else
      set(gca, 'XScale', 'linear');
      set(gca,'TickDir','out','XTick',0:25:200);
    end
    
    % ------------------------------------------------
    % PLOT VARIANCE OVER TIME
    % ------------------------------------------------
    if isempty(subcomp{3,row})
      subcomp{3,row} = subplot(numOfPlots,4,row*4-1);
      plot(comp_time,comp_var(compNum,:),'k.');
      xlabel('Chunk'); ylabel('Variance');
      axis tight; set(gca, 'tickdir', 'out');
      set(gca, 'ButtonDownFcn', @(h,evt)toggle_visual, 'tag', sprintf('var_%s', comp.label{compNum}));
    else
      l = get(subcomp{3,row}, 'children');
      set(l, 'YData', comp_var(compNum,:));
      set(subcomp{3,row}, 'tag',  sprintf('var_%s', comp.label{compNum}));
    end

    % ------------------------------------------------
    % PLOT COMPONENT TOPOGRAPHY
    % ------------------------------------------------
    subcomp{4,row} = subplot(numOfPlots,4,row*4, 'Tag', sprintf('topo_%s', comp.label{compNum}));
    ft_plot_layout(lay, 'label', 'off', 'box', 'no', 'point', 'no');
    dat = comp.topo(i2, compNum);
    % Get physical min/max range of z:
    if strcmp(cfg.zlim, 'maxmin')
      zmin = min(dat);
      zmax = max(dat);
    elseif strcmp(cfg.zlim, 'maxabs')
      zmin = -max(max(abs(dat)));
      zmax = max(max(abs(dat)));
    elseif strcmp(cfg.zlim, 'zeromax')
      zmin = 0;
      zmax = max(dat);
    elseif strcmp(cfg.zlim, 'minzero')
      zmin = min(dat);
      zmax = 0;
    else
      zmin = cfg.zlim(1);
      zmax = cfg.zlim(2);
    end
    ft_plot_topo(xpos, ypos, dat, 'outline', lay.outline, 'mask', lay.mask, 'clim', [zmin zmax], 'gridscale', 95);
    
    % decorate with buttons
    ypos_button = [0.73 0.51 0.29 0.07];

    % ------------------------------------------------
    % SHOW TIMECOURSE OF THIS COMPONENT BUTTON
    % ------------------------------------------------
    button{row,1} = uicontrol('Units','normalized','Position',[0.9 ypos_button(row) 0.075 0.035],...
      'Style','pushbutton','String','Timecourse','Callback',@(h, evt)tc_cb(compNum));

    % ------------------------------------------------
    % REJECT COMPONENT BUTTON
    % ------------------------------------------------
    button{row,2} = uicontrol('Units','normalized', 'Tag', sprintf('rej%d',row), 'Position',[0.75 ypos_button(row) 0.11 0.035],...
      'Style','checkbox','String','select for rejection','Callback',@(h, evt)rej_callback(row));

    % ------------------------------------------------
    % SAVE COMPONENT PDF BUTTON
    % ------------------------------------------------
    button{row,3} = uicontrol('Units','normalized','Position',[0.9 ypos_button(row)+0.05 0.075 0.035],'Style','pushbutton','String','Save PDF','Callback',@(h, evt)sc_cb(1));
    if isempty(outputfolder)
      button{row,3}.Enable = 'off';
    end

    if mod(compNum, numOfPlots)==0 || compNum == numel(comp.label)
      % ------------------------------------------------
      % ENABLE/DISABLE NEXT/LAST PAGE BUTTONS
      % ------------------------------------------------
      if compNum > 4
        prev.Enable = 'on';
        first.Enable = 'on';
      else
        prev.Enable = 'off';
        first.Enable = 'off';
      end

      if compNum < size(comp.label,1)-3
        next.Enable = 'on';
        last.Enable = 'on';
      else
        next.Enable = 'off';
        last.Enable = 'off';
      end

      if isempty(outputfile)
        save_it.Enable = 'Off';
      else
        save_it.Enable = 'On';
      end
      freqlog.Enable = 'on';
      powlog.Enable  = 'on';
      quit_it.Enable = 'on';
      uiwait
    end
  end % for
end % while
rej_comp = find(rej_comp);

% ------------------------------------------------
% DEFINE NESTED CALLBACK FUNCTIONS
% ------------------------------------------------
  function rej_callback(k)
    rejk = findobj('Tag', sprintf('rej%d',k));
    if (rej_comp(compNum+k-4) == 0)
      %set(rejk,'Backgroundcolor',[0.5 0 0]);
      %set(rejk,'String', 'Reject');
      rej_comp(compNum+k-4) = true;
    else
      %set(rejk,'Backgroundcolor',[0 0.5 0]);
      %set(rejk,'String', 'Keep');
      rej_comp(compNum+k-4) = false;
    end
  end

  function toggle_visual()
    % copied from FT_SELECT_BOX, but without the waitforbuttonpress command since here it is triggered by the ButtonDown event
    point1 = get(gca, 'CurrentPoint');    % button down detected
    rbbox;                                % this draws the rubber-band box
    point2 = get(gca, 'CurrentPoint');    % button up detected
    point1 = point1(1, 1:2);              % extract x and y
    point2 = point2(1, 1:2);
    x = sort([point1(1) point2(1)]);
    y = sort([point1(2) point2(2)]);

    thistag = get(gca, 'tag');
    if ~contains(thistag, 'var')
      return;
    else
      % identify the intended trial to be displayed
      l = get(gca,'children');
      xsel = find(l.XData>x(1) & l.XData<x(2));
      ysel = l.YData(xsel);
      iy = xsel(ysel>y(1) & ysel<y(2));
      if numel(iy)>1
        iy = iy(1);
      end

      if ~isempty(iy)
        trial_for_disp = trialinfo(iy);
        for k = 1:4
          p = get(subcomp{1, k}, 'children');
          set(p, 'XData', comp.time{trial_for_disp}, 'YData', comp.trial{trial_for_disp}(compNum-row+k, :));
        end
        set(get(subcomp{1,1}, 'xlabel'), 'string', sprintf('trial %d', trial_for_disp));
      end
    end
    uiresume;
    
  end

% timecourse funcs
  function tc_cb(compNum)
    thistag = get(gca, 'tag');
    if ~startsWith(thistag, 'topo')
      return;
    else
      thistag = thistag(6:end);
    end

    cfgtc = [];
    cfgtc.layout        = lay;
    cfgtc.viewmode      = 'butterfly';
    cfgtc.channel       = compNum;
    cfgtc.blocksize     = blocksize;
    cfgtc.showcallinfo  = showcallinfo;
    cfgtc.linecolor     = 'k';

    ft_info off;
    cfgout = ft_databrowser(cfgtc, comp);
    if isfield(cfgout, 'artfctdef') && isfield(cfgout.artfctdef, 'visual')
      artifact = cat(1, artifact, cfgout.artfctdef.visual.artifact);
    end
    ft_info on;
    uiresume;
  end

% save to figure
  function sc_cb(whichcomp)
    h = figure;
    set(h,'Position',[200 200 1000 300]);
    set(h,'Units','inches');
    screenposition = get(h,'Position');
    set(h, 'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));
    new = copyobj(subcomp{1}{whichcomp},h);
    set(new,'Position',[.05 .1 0.25 0.85]);
    new = copyobj(subcomp{2}{whichcomp},h);
    set(new,'Position',[.35 .1 0.25 0.85]);
    set(new,'LineWidth',2);
    new = copyobj(subcomp{3}{whichcomp},h); set(new,'Position',[.55 .05 0.5 0.95]);

    % save under the correct comp nr
    compnrs = compNum-3:compNum;
    print(h,'-dpdf',sprintf('%s/%s_comp%d.pdf', outputfolder, prefix, compnrs(whichcomp)));
    fprintf('saved pdf to %s/%s_comp%d.pdf\n', outputfolder, prefix, compnrs(whichcomp));
    close(h)
  end

% gui
  function firstpage()
    manpos = get(f,'OuterPosition');
    page = 1;
    uiresume;
  end

  function prevpage()
    manpos = get(f,'OuterPosition');
    page = page - 1;
    uiresume;
  end

  function nextpage()
    manpos = get(f,'OuterPosition');
    page = page + 1;
    uiresume;
  end

  function lastpage()
    manpos = get(f,'OuterPosition');
    page = numOfPages;
    uiresume;
  end

  function powplotlog()
    manpos = get(f,'OuterPosition');
    if strcmp(powscale, 'log10')
      powscale = 'lin';
    else
      powscale = 'log10';
    end
    uiresume;
  end

  function freqplotlog()
    manpos = get(f,'OuterPosition');
    if strcmp(freqscale, 'log')
      freqscale = 'lin';
    else
      freqscale = 'log';
    end
    clf;
    uiresume;
  end

  function save_callback()
    if ~isempty(outputfile)
      idx = find(rej_comp == 1);
      save(outputfile, 'idx', 'rej_comp');
    end
    delete(f);
    disp('saved');
    err = 1;
  end

  function quitme()
    delete(f);
    err = 1;
  end

  function trl = sampleinfo2trl(data)

    % SAMPLEINFO2TRL constructs the trial definition from the sampleinfo, the time axes
    % and optionally from the trialinfo
    %
    % Use as
    %   trl = sampleinfo2trl(data)
    %
    % See also ARTIFACT2BOOLVEC, ARTIFACT2EVENT, ARTIFACT2TRL, BOOLVEC2ARTIFACT, BOOLVEC2EVENT, BOOLVEC2TRL, EVENT2ARTIFACT, EVENT2BOOLVEC, EVENT2TRL, TRL2ARTIFACT, TRL2BOOLVEC, TRL2EVENT

    % get the begin and end sample of each trial
    begsample = data.sampleinfo(:,1);
    endsample = data.sampleinfo(:,2);

    % recreate the offset
    offset = zeros(numel(data.trial), 1);
    for i=1:numel(data.trial)
      offset(i) = time2offset(data.time{i}, data.fsample);
    end

    if isfield(data, 'trialinfo') && istable(data.trialinfo)
      trl = table(begsample, endsample, offset);
      trl = horzcat(trl, data.trialinfo);
    elseif isfield(data, 'trialinfo') && isnumeric(data.trialinfo)
      trl = [begsample endsample offset data.trialinfo];
    else
      trl = [begsample endsample offset];
    end
  end

  function offset = time2offset(time, fsample)

    % TIME2OFFSET converts a time-axis of a trial into the offset in samples
    % according to the definition from DEFINETRIAL
    %
    % Use as
    %   [offset] = time2offset(time, fsample)
    %
    % The trialdefinition "trl" is an Nx3 matrix. The first column contains
    % the sample-indices of the begin of the trial relative to the begin
    % of the raw data , the second column contains the sample_indices of
    % the end of the trials, and the third column contains the offset of
    % the trigger with respect to the trial. An offset of 0 means that
    % the first sample of the trial corresponds to the trigger. A positive
    % offset indicates that the first sample is later than the trigger, a
    % negative offset indicates a trial beginning before the trigger.

    if ~isempty(time)
      offset = round(time(1)*fsample);
    else
      offset = nan;
    end
  end

end
