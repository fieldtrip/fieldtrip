function browse_simpleFFT(cfg, data)

% BROWSE_SIMPLEFFT is a helper function for FT_DATABROWSER that shows a
% simple FFT of the data. 
%
% Included are a button to switch between log and non-log space, and a selection button to deselect channels,
% for the purpose of zooming in on bad channels.
%
%
%
% See also BROWSE_MOVIEPLOTER, BROWSE_TOPOPLOTER, BROWSE_MULTIPLOTER, BROWSE_TOPOPLOTVAR, BROWSE_SIMPLEFFT

% Copyright (C) 2011, Roemer van der Meij

% call ft_freqanalysis on data
cfgfreq = [];
cfgfreq.method = 'mtmfft';
cfgfreq.taper  = 'boxcar';
freqdata = ft_freqanalysis(cfgfreq,data);

% remove zero-bin for plotting
if freqdata.freq(1) == 0
  freqdata.freq = freqdata.freq(2:end);
  freqdata.powspctrm = freqdata.powspctrm(:,2:end);
end
  
%%% FIXME: call ft_singleplotER for everything below 

% make figure window for fft
ffth = figure('name',cfg.figurename,'numbertitle','off','units','normalized');
% set button
butth = uicontrol('tag', 'simplefft_l2', 'parent', ffth, 'units', 'normalized', 'style', 'checkbox', 'string','log10 of power','position', [0.87, 0.6 , 0.12, 0.05],  'value',1,'callback',{@draw_simple_fft_cb},'backgroundcolor',[.8 .8 .8]);
uicontrol('tag', 'simplefft_l2_chansel', 'parent', ffth, 'units', 'normalized', 'style', 'pushbutton', 'string','select channels','position', [0.87, 0.45 , 0.12, 0.10],'callback',{@selectchan_fft_cb});

% put dat in fig (sparse)
fftopt = [];
fftopt.freqdata   = freqdata;
fftopt.chancolors = cfg.chancolors;
fftopt.chansel    = 1:numel(data.label);
fftopt.butth      = butth;
setappdata(ffth, 'fftopt', fftopt);

% draw fig
draw_simple_fft_cb(butth)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function selectchan_fft_cb(h,eventdata)
ffth = get(h,'parent');
fftopt = getappdata(ffth, 'fftopt');

% open chansel dialog
chansel = select_channel_list(fftopt.opt.freq.label, fftopt.chansel, 'select channels for viewing power');

% output data
fftopt.chansel = chansel;
setappdata(ffth, 'fftopt', fftopt);
draw_simple_fft_cb(fftopt.butth)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function draw_simple_fft_cb(h,eventdata)
togglestate = get(h,'value');
ffth        = get(h,'parent');
fftopt = getappdata(ffth, 'fftopt');
% clear axis (for switching)
cla

% switch log10 or nonlog10
if togglestate == 0
  dat = fftopt.freqdata.powspctrm;
 elseif togglestate == 1
  dat = log10(fftopt.freqdata.powspctrm);
end

% select data and chanel colors
chancolors = fftopt.chancolors(fftopt.chansel,:);
dat        = dat(fftopt.chansel,:);

% plot using specified colors
set(0,'currentFigure',ffth)
for ichan = 1:size(dat,1)
  color = chancolors(ichan,:);
  ft_plot_vector(fftopt.freqdata.freq, dat(ichan,:), 'box', false, 'color', color)
end
ylabel('log10(power)')
xlabel('frequency (hz)')
yrange = abs(max(max(dat)) - min(min(dat)));
axis([fftopt.freqdata.freq(1) fftopt.freqdata.freq(end) (min(min(dat)) - yrange.*.1) (max(max(dat)) + yrange*.1)]) 

% switch log10 or nonlog10
if togglestate == 0
  ylabel('power')
  axis([fftopt.freqdata.freq(1) fftopt.freqdata.freq(end) 0 (max(max(dat)) + yrange*.1)]) 
elseif togglestate == 1
  ylabel('log10(power)')
  axis([fftopt.freqdata.freq(1) fftopt.freqdata.freq(end) (min(min(dat)) - yrange.*.1) (max(max(dat)) + yrange*.1)]) 
end
set(gca,'Position', [0.13 0.11 0.725 0.815])

