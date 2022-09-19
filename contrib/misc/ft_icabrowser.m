function [rej_comp] = ft_icabrowser(cfg, comp)

% FT_ICABROWSER loads in comp structure from FieldTrip ft_componentanalysis
% and presents a GUI interface showing the power spectrum, variance over
% time and the topography of the components, as well as the possibility to
% save a PDF, view the timecourse and toggle components to be rejected vs
% kept.
%
% Use as
%    [rej_comp] = ft_icabrowser(cfg, comp)
%
% where the input comp structure should be obtained from FT_COMPONENTANALYSIS.
%
% The configuration must contain:
%   cfg.layout     = filename of the layout, see FT_PREPARE_LAYOUT
%
% further optional configuration parameters are
%   cfg.rejcomp       = list of components which shall be initially marked for rejection, e.g. [1 4 7]
%   cfg.blocksize     = blocksize of time course (default = 1 sec)
%   cfg.powscale      = scaling of y axis in power plot, 'lin' or 'log10', (default = 'log10')
%   cfg.zlim          = plotting limits for color dimension of topoplot, 'maxmin', 'maxabs', 'zeromax', 'minzero', or [zmin zmax] (default = 'maxmin')
%   cfg.path          = where pdfs will be saves (default = pwd)
%   cfg.prefix        = prefix of the pdf files (default = 'ICA')
%   cfg.colormap      = any sized colormap, see COLORMAP
%   cfg.outputfile    = MAT file which contains indices of all components to reject
%   cfg.showcallinfo  = show call info, 'yes' or 'no' (default: 'no')
%
% original written by Thomas Pfeffer
% adapted by Jonathan Daume and Anne Urai
% University Medical Center Hamburg-Eppendorf, 2015
%
% modified by Daniel Matthes
% Max Planck Institute for Human Cognitive and Brain Sciences, 2019
%
% See also FT_COMPONENTANALYSIS, FT_TOPOPLOTIC, FT_PREPARE_LAYOUT

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if the input data is valid for this function
comp = ft_checkdata(comp, 'datatype', 'raw+comp');

% set the defaults
layout        = ft_getopt(cfg, 'layout');
rejcomp       = ft_getopt(cfg, 'rejcomp', []);
blocksize     = ft_getopt(cfg, 'blocksize', 1);
powscale      = ft_getopt(cfg, 'powscale', 'log10');
zlim          = ft_getopt(cfg, 'zlim', 'maxmin');
path          = ft_getopt(cfg, 'path', pwd);
prefix        = ft_getopt(cfg, 'prefix', 'ICA');
colormap      = ft_getopt(cfg, 'colormap', 'default');
outputfile    = ft_getopt(cfg, 'outputfile', []);
showcallinfo  = ft_getopt(cfg, 'showcallinfo', 'no');

% cfg.rejcomp can be indicated by number or by label
rejcomp = ft_channelselection(rejcomp, comp.label);
reject = match_str(comp.label, rejcomp);

% create folder if not exist
if ~exist(path, 'dir')
	mkdir path;
end

% test config
if ~ismember(powscale, {'lin', 'log10'})
  powscale = 'log10';
end

% setup
var_data = cat(2,comp.trial{:});
var_time = (1:(size(var_data,2)))/comp.fsample;
% only do the fft on a subset of trials, saves time
fft_data = cat(2,comp.trial{1:5:end});
% preallocate rejected components
rej_comp = false(size(comp.label,1),1);
rej_comp(reject) = true;

numOfPlots = 4;                                                             % number of subplots per page
page = 1;                                                                   % number of current page
numOfPages = ceil(size(comp.label, 1)/4);                                   % maximum number of pages
row = 0;                                                                    % number of current row in page

% to save time redoing this for each topo
cfglay              = [];
cfglay.layout       = layout;
cfglay.showcallinfo = showcallinfo;
lay                 = ft_prepare_layout(cfglay);

cfgtopo = [];
cfgtopo.layout        = lay;                                                % specify the layout file that should be used for plotting
cfgtopo.comment       = 'no';
cfgtopo.zlim          = zlim;
cfgtopo.highlight     = 'off';
cfgtopo.marker        = 'off';
cfgtopo.style         = 'straight';
cfgtopo.showcallinfo  = showcallinfo;
if ~isempty(colormap)
  cfgtopo.colormap    = colormap;
end

err = 0;
manpos = [0.1 0.1 0.8 0.8];                                                 % figure position, can be updated later

% ------------------------------------------------
% COMPUTE LATENCY FOR 2s-WINDOWS
% ------------------------------------------------

slen = floor(2*comp.fsample);
smax = floor(size(var_data,2)/slen);
comp_var  = nan(numOfPlots, smax); % preallocate
comp_time = nan(1, smax);     % preallocate
for s = 1 : smax
  comp_time(s) = mean(var_time(1,(s-1)*slen+1:s*slen));
end

f = figure('units','normalized','outerposition', manpos, 'CloseRequestFcn', @(h, evt)quitme);

while err == 0 % KEEP GOING UNTIL THERE IS AN ERROR

  while row < numOfPlots % il is the subplot count
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
    % SWITCH FREQUENCY AXIS LOG/LINEAR
    % ------------------------------------------------
    log = uicontrol('Units','normalized','Position',[0.45 0.01 0.075 0.05],'Style','pushbutton', 'Callback',@(h, evt)plotlog);
    log.Enable = 'off';
    if strcmp(powscale, 'log10')
      log.String = 'Linear';
    else
      log.String = 'Log10';
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

    row = row + 1;
    compNum = (page-1) * numOfPlots + row;                              % estimate number of component

    % ------------------------------------------------
    % COMPUTE VARIANCE FOR 2s-WINDOWS
    % ------------------------------------------------

    for s = 1 : smax
        comp_var(compNum,s)=var(var_data(compNum,(s-1)*slen+1:s*slen));
    end

    % ------------------------------------------------
    % COMPUTE POWER SPECTRUM
    % ------------------------------------------------

    smo = 50;
    steps = 10;
    Fs = comp.fsample;
    N = floor(size(fft_data,2));
    xdft = fft(fft_data(compNum,:));
    xdft = xdft(1:floor(N/2)+1);
    psdx = (1/(Fs*N)).*abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);

    smoothed = conv(psdx, ones(smo,1), 'valid');
    smoothed = smoothed(1:steps:fix(length(psdx) - smo))/smo;

    freq = linspace(0,Fs/2,size(smoothed,2));
    strt = find(freq > 2,1,'first');
    stp  = find(freq < 200,1,'last');

    % ------------------------------------------------
    % PLOT POWER SPECTRUM
    % ------------------------------------------------
    subcomp{1}{row} = subplot(numOfPlots,3,row*3-2);
    if strcmp(powscale, 'log10')
      plot(freq(strt:stp),log10(smoothed(strt:stp)));
      ylabel('PSD (dB/Hz)');
    else
      plot(freq(strt:stp),smoothed(strt:stp));
      ylabel('PSD (T^2/Hz)');
    end
    set(gca,'TickDir','out','XTick',0:25:200)
    xlabel('Frequency (Hz)'); grid on;
    axis tight;

    % ------------------------------------------------
    % PLOT VARIANCE OVER TIME
    % ------------------------------------------------
    subcomp{2}{row} = subplot(numOfPlots,3,row*3-1);
    scatter(comp_time,comp_var(compNum,:),'k.');
    xlabel('Time (s)'); ylabel('Variance');
    axis tight; set(gca, 'tickdir', 'out');

    % ------------------------------------------------
    % PLOT COMPONENT TOPOGRAPHY
    % ------------------------------------------------
    subcomp{3}{row} = subplot(numOfPlots,3,row*3);
    cfgtopo.component = compNum;       % specify the component(s) that should be plotted
    ft_topoplotIC(cfgtopo, comp);

    if mod(compNum,numOfPlots)==0 || compNum == 80
      set(f, 'Name', sprintf('%d: %s', double(f), 'ft_icabrowser'));

      pos = [0.76 0.73 0.075 0.035; ...
          0.76 0.51 0.075 0.035; ...
          0.76 0.29 0.075 0.035; ...
          0.76 0.07 0.075 0.035];

      bgc     = cell(1,4);
      rej_str = cell(1,4);
      for ibgc = 1 : 4
        if rej_comp(compNum+(ibgc-4)) == 1
          bgc{ibgc} = 'r';
          rej_str{ibgc} = 'Reject';
        else
          bgc{ibgc} = 'g';
          rej_str{ibgc} = 'Keep';
        end
      end

      % ------------------------------------------------
      % SHOW TIMECOURSE OF THIS COMPONENT
      % ------------------------------------------------
      uicontrol('Units','normalized','Position',[0.86 0.73 0.075 0.035],...
          'Style','pushbutton','String','Timecourse','Callback',@(h, evt)tc_cb(1));
      uicontrol('Units','normalized','Position',[0.86 0.51 0.075 0.035],...
          'Style','pushbutton','String','Timecourse','Callback',@(h, evt)tc_cb(2));
      uicontrol('Units','normalized','Position',[0.86 0.29 0.075 0.035],...
          'Style','pushbutton','String','Timecourse','Callback',@(h, evt)tc_cb(3));
      uicontrol('Units','normalized','Position',[0.86 0.07 0.075 0.035],...
          'Style','pushbutton','String','Timecourse','Callback',@(h, evt)tc_cb(4));

      % ------------------------------------------------
      % REJECT COMPONENT
      % ------------------------------------------------
      rej1 = uicontrol('Units','normalized', 'Tag', 'rej1', 'Position',pos(1,:),'Style','pushbutton','String', rej_str{1}, ...
          'Backgroundcolor',bgc{1},'Callback',@(h, evt)rej_callback1);
      rej2 = uicontrol('Units','normalized','Tag', 'rej2','Position',pos(2,:),'Style','pushbutton','String', rej_str{2}, ...
          'Backgroundcolor',bgc{2},'Callback',@(h, evt)rej_callback2);
      rej3 = uicontrol('Units','normalized','Tag', 'rej3','Position',pos(3,:),'Style','pushbutton','String', rej_str{3}, ...
          'Backgroundcolor',bgc{3},'Callback',@(h, evt)rej_callback3);
      rej4 = uicontrol('Units','normalized','Tag', 'rej4','Position',pos(4,:),'Style','pushbutton','String', rej_str{4}, ...
          'Backgroundcolor',bgc{4},'Callback',@(h, evt)rej_callback4);

      % ------------------------------------------------
      % SAVE COMPONENT PDF
      % ------------------------------------------------
      savecomp1 = uicontrol('Units','normalized','Position',[0.86 0.78 0.075 0.035],'Style','pushbutton','String','Save PDF','Callback',@(h, evt)sc_cb(1));
      savecomp2 = uicontrol('Units','normalized','Position',[0.86 0.56 0.075 0.035],'Style','pushbutton','String','Save PDF','Callback',@(h, evt)sc_cb(2));
      savecomp3 = uicontrol('Units','normalized','Position',[0.86 0.34 0.075 0.035],'Style','pushbutton','String','Save PDF','Callback',@(h, evt)sc_cb(3));
      savecomp4 = uicontrol('Units','normalized','Position',[0.86 0.12 0.075 0.035],'Style','pushbutton','String','Save PDF','Callback',@(h, evt)sc_cb(4));
      if isempty(path)
        savecomp1.Enable = 'off';
        savecomp2.Enable = 'off';
        savecomp3.Enable = 'off';
        savecomp4.Enable = 'off';
      end

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
      log.Enable = 'on';
      quit_it.Enable = 'on';
      uiwait
    end
  end
end

% ------------------------------------------------
% DEFINE NESTED CALLBACK FUNCTIONS
% ------------------------------------------------
    function rej_callback1()
      rej1 = findobj('Tag','rej1');
      if (rej_comp(compNum-3) == 0)
        set(rej1,'Backgroundcolor','r');
        set(rej1,'String', 'Reject');
        rej_comp(compNum-3) = true;
      else
        set(rej1,'Backgroundcolor','g');
        set(rej1,'String', 'Keep');
        rej_comp(compNum-3) = false;
      end
    end

    function rej_callback2()
      rej2 = findobj('Tag','rej2');
      if (rej_comp(compNum-2) == 0)
        set(rej2,'Backgroundcolor','r');
        set(rej2,'String', 'Reject');
        rej_comp(compNum-2) = true;
      else
        set(rej2,'Backgroundcolor','g');
        set(rej2,'String', 'Keep');
        rej_comp(compNum-2) = false;
      end
    end
    function rej_callback3()
      rej3 = findobj('Tag','rej3');
      if (rej_comp(compNum-1) == 0)
        set(rej3,'Backgroundcolor','r');
        set(rej3,'String', 'Reject');
        rej_comp(compNum-1) = true;
      else
        set(rej3,'Backgroundcolor','g');
        set(rej3,'String', 'Keep');
        rej_comp(compNum-1) = false;
      end
    end
    function rej_callback4()
      rej4 = findobj('Tag','rej4');
      if (rej_comp(compNum-0) == 0)
        set(rej4,'Backgroundcolor','r');
        set(rej4,'String', 'Reject');
        rej_comp(compNum-0) = true;
      else
        set(rej4,'Backgroundcolor','g');
        set(rej4,'String', 'Keep');
        rej_comp(compNum-0) = false;
      end
    end

% timecourse funcs
    function tc_cb(whichcomp)
      cfgtc = [];
      cfgtc.layout        = lay;
      cfgtc.viewmode      = 'butterfly';
      cfgtc.channel       = compNum + (whichcomp-4);
      cfgtc.blocksize     = blocksize;
      cfgtc.showcallinfo  = showcallinfo;


      ft_info off;
      ft_databrowser(cfgtc, comp);
      ft_info on;
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
      print(h,'-dpdf',sprintf('%s/%s_comp%d.pdf', path, prefix, compnrs(whichcomp)));
      fprintf('saved pdf to %s/%s_comp%d.pdf\n', path, prefix, compnrs(whichcomp));
      close(h)
    end

% gui
    function firstpage()
      manpos = get(f,'OuterPosition');
      page = 1;
      row = 0;
      clf;
      uiresume;
    end

    function prevpage()
      manpos = get(f,'OuterPosition');
      page = page - 1;
      row = 0;
      clf;
      uiresume;
    end

    function nextpage()
      manpos = get(f,'OuterPosition');
      page = page + 1;
      row = 0;
      clf;
      uiresume;
    end

    function lastpage()
      manpos = get(f,'OuterPosition');
      page = numOfPages;
      row = 0;
      clf;
      uiresume;
    end

    function plotlog()
      manpos = get(f,'OuterPosition');
      row = 0;
      if strcmp(powscale, 'log10')
        powscale = 'lin';
      else
        powscale = 'log10';
      end
      clf;
      uiresume;
    end

    function save_callback()
      if ~isempty(outputfile)
          idx = find(rej_comp == 1); %#ok<NASGU>
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

end
