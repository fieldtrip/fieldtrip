function browse_simpleFFT(cfg, data)

% BROWSE_SIMPLEFFT is a helper function for FT_DATABROWSER that shows a
% simple FFT of the data.
%
% Included are a button to switch between log and non-log space, and a
% selection button to deselect channels, for the purpose of zooming in on
% bad channels.
%
% See also BROWSE_MOVIEPLOTER, BROWSE_TOPOPLOTER, BROWSE_MULTIPLOTER, BROWSE_TOPOPLOTVAR, BROWSE_SIMPLEFFT

% Copyright (C) 2011-2014, Roemer van der Meij

% call ft_freqanalysis on data
cfgfreq = [];
cfgfreq.method = 'mtmfft';
cfgfreq.taper = 'boxcar';
freqdata = ft_freqanalysis(cfgfreq, data);

% remove zero-bin for plotting
if freqdata.freq(1) == 0
  freqdata.freq = freqdata.freq(2:end);
  freqdata.powspctrm = freqdata.powspctrm(:, 2:end);
end

%%% FIXME: it would be nice to use ft_singleplotER for everything below,
%%% this would mean that ft_singleplotER needs to get lin/log buttons

% make figure window for fft
ffth = figure('name', cfg.figurename, 'numbertitle', 'off', 'units', 'normalized', 'toolbar', 'figure');
% set button
button_x       = uicontrol('tag', 'semilogx', 'parent', ffth, 'units', 'normalized', 'style', 'checkbox', 'string', 'log10(frequency)', 'position', [0.87, 0.75 , 0.12, 0.05], 'value', true, 'callback', {@togglex}, 'backgroundcolor', [.8 .8 .8]);
button_y       = uicontrol('tag', 'semilogy', 'parent', ffth, 'units', 'normalized', 'style', 'checkbox', 'string', 'log10(power)', 'position', [0.87, 0.60 , 0.12, 0.05], 'value', true, 'callback', {@toggley}, 'backgroundcolor', [.8 .8 .8]);
button_chansel = uicontrol('tag', 'chansel', 'parent', ffth, 'units', 'normalized', 'style', 'pushbutton', 'string', 'select channels', 'position', [0.87, 0.45 , 0.12, 0.10], 'callback', {@selectchan_fft_cb});

% put dat in fig (sparse)
fftopt = [];
fftopt.freqdata = freqdata;
fftopt.chancolors = cfg.chancolors;
fftopt.chansel = 1:numel(data.label);
fftopt.semilogx = true;
fftopt.semilogy = true;
setappdata(ffth, 'fftopt', fftopt);

% draw fig
draw_simple_fft_cb(button_chansel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function togglex(h, eventdata)
ffth = get(h, 'parent');
fftopt = getappdata(ffth, 'fftopt');
fftopt.semilogx = ~fftopt.semilogx;
setappdata(ffth, 'fftopt', fftopt);
draw_simple_fft_cb(h)

function toggley(h, eventdata)
ffth = get(h, 'parent');
fftopt = getappdata(ffth, 'fftopt');
fftopt.semilogy = ~fftopt.semilogy;
setappdata(ffth, 'fftopt', fftopt);
draw_simple_fft_cb(h)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function selectchan_fft_cb(h, eventdata)
ffth = get(h, 'parent');
fftopt = getappdata(ffth, 'fftopt');

% open chansel dialog
chansel = select_channel_list(fftopt.freqdata.label, fftopt.chansel, 'select channels for viewing power');

% output data
fftopt.chansel = chansel;
setappdata(ffth, 'fftopt', fftopt);
draw_simple_fft_cb(h)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function draw_simple_fft_cb(h, eventdata)
ffth = get(h, 'parent');
fftopt = getappdata(ffth, 'fftopt');
% clear axis (for switching)
cla

% toggle log10 or not
if fftopt.semilogx
  freq = log10(fftopt.freqdata.freq);
else
  freq = fftopt.freqdata.freq;
end

if fftopt.semilogy
  dat = log10(fftopt.freqdata.powspctrm);
else
  dat = fftopt.freqdata.powspctrm;
end


% select data and chanel colors
chancolors = fftopt.chancolors(fftopt.chansel, :);
dat = dat(fftopt.chansel, :);

% plot using specified colors
set(0, 'currentFigure', ffth)
for ichan = 1:size(dat, 1)
  color = chancolors(ichan, :);
  ft_plot_vector(freq, dat(ichan, :), 'box', false, 'color', color)
end

minx = min(freq);
maxx = max(freq);
miny = min(dat(:));
maxy = max(dat(:));
yrange = maxy-miny;
axis([minx maxx miny-yrange*0.1 maxy+yrange*.01])

if fftopt.semilogx
  xlabel('log10(frequency)')
else
  xlabel('frequency')
end

if fftopt.semilogy
  ylabel('log10(power)')
else
  ylabel('power')
end

set(gca, 'Position', [0.13 0.11 0.725 0.815])
