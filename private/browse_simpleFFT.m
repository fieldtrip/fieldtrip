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



% select  the requested data segment and do fft
dat     = data.trial{1};
fsample = data.fsample;
nsample = size(dat,2);
fftdat = transpose(fft(transpose(dat))); % double
% scale the same way as mtmfft
fftdat = fftdat .* sqrt(2 ./ size(dat,2));
% compute powerspectrum (my pref would be to actually don't square, to make line spectra easier to visualize in case of low freq stuff, but doing it so people dont get confused)
fftdat = abs(fftdat) .^ 2;

% create proper xaxis
freqboilim = round([0 fsample/2] ./ (fsample ./ nsample)) + 1;
freqboi    = freqboilim(1):1:freqboilim(2);
freqoi     = (freqboi-1) ./ (nsample * (1/fsample));
xaxis = freqoi;

% make figure window for fft
ffth = figure('name',cfg.figurename,'numbertitle','off','units','normalized');
% set button
butth = uicontrol('tag', 'simplefft_l2', 'parent', ffth, 'units', 'normalized', 'style', 'togglebutton', 'string','log of power','position', [0.87, 0.6 , 0.12, 0.10],  'value',1,'callback',{@draw_simple_fft_cb});
uicontrol('tag', 'simplefft_l2_chansel', 'parent', ffth, 'units', 'normalized', 'style', 'pushbutton', 'string','select channels','position', [0.87, 0.45 , 0.12, 0.10],'callback',{@selectchan_fft_cb});

% put dat in fig (sparse)
fftopt = [];
fftopt.chancolors = cfg.chancolors;
fftopt.xaxis      = xaxis;
fftopt.xaxis      = xaxis;
fftopt.chanlabel  = data.label;
fftopt.chansel    = 1:size(fftdat,1);
fftopt.fftdat     = fftdat(:,freqboi);
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
chansel = select_channel_list(fftopt.chanlabel, fftopt.chansel, 'select channels for viewing power');

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

% switch log or nonlog
if togglestate == 0
  dat = fftopt.fftdat;
  set(h,'backgroundColor',[0.8 0.8 0.8])
elseif togglestate == 1
  dat = log(fftopt.fftdat);
  set(h,'backgroundColor','g')
end

% select data and chanel colors
chancolors = fftopt.chancolors(fftopt.chansel,:);
dat        = dat(fftopt.chansel,:);

% plot using specified colors
set(0,'currentFigure',ffth)
for ichan = 1:size(dat,1)
  color = chancolors(ichan,:);
  ft_plot_vector(fftopt.xaxis, dat(ichan,:), 'box', false, 'color', color)
end
ylabel('log(power)')
xlabel('frequency (hz)')
yrange = abs(max(max(dat)) - min(min(dat)));
axis([fftopt.xaxis(1) fftopt.xaxis(end) (min(min(dat)) - yrange.*.1) (max(max(dat)) + yrange*.1)]) 

% switch log or nonlog
if togglestate == 0
  ylabel('power')
  axis([fftopt.xaxis(1) fftopt.xaxis(end) 0 (max(max(dat)) + yrange*.1)]) 
elseif togglestate == 1
  ylabel('log(power)')
  axis([fftopt.xaxis(1) fftopt.xaxis(end) (min(min(dat)) - yrange.*.1) (max(max(dat)) + yrange*.1)]) 
end
set(gca,'Position', [0.13 0.11 0.725 0.815])

