function [spectrum,freqoi,timeoi] = ft_specest_wavelet(dat, time, varargin)

% SPECEST_WAVELET performs time-frequency analysis on any time series trial
% data using the 'wavelet method' based on Morlet wavelets, doing
% convolution in the time domain by multiplaction in the frequency domain
%
% Use as
%   [spectrum,freqoi,timeoi] = specest_wavelet(dat,time...)
% where
%   dat      = matrix of chan*sample
%   time     = vector, containing time in seconds for each sample
%   spectrum = array of chan*freqoi*timeoi of fourier coefficients
%   freqoi   = vector of frequencies in spectrum
%   timeoi   = vector of timebins in spectrum
%
% Optional arguments should be specified in key-value pairs and can include:
%   pad        = number, total length of data after zero padding (in seconds)
%   freqoi     = vector, containing frequencies of interest
%   timeoi     = vector, containing time points of interest (in seconds)
%   width      = number or vector, width of the wavelet, determines the temporal and spectral resolution
%   gwidth     = number, determines the length of the used wavelets in standard deviations of the implicit Gaussian kernel
%
% See also SPECEST_MTMCONVOL, SPECEST_CONVOL, SPECEST_HILBERT, SPECEST_MTMFFT

% Copyright (C) 2010, Donders Institute for Brain, Cognition and Behaviour
%
% $Log$

% get the optional input arguments
keyvalcheck(varargin, 'optional', {'pad','width','gwidth','freqoi','timeoi','polyremoval'});
freqoi    = keyval('freqoi',      varargin);  if isempty(freqoi),   freqoi  = 'all';   end
timeoi    = keyval('timeoi',      varargin);  if isempty(timeoi),   timeoi  = 'all';   end
width     = keyval('width',       varargin);  if isempty(width),    width    = 7;      end
gwidth    = keyval('gwidth',      varargin);  if isempty(gwidth),   gwidth   = 3;      end
pad       = keyval('pad',         varargin);
polyorder = keyval('polyremoval', varargin); if isempty(polyorder), polyorder = 1; end


% Set n's
[nchan,ndatsample] = size(dat);

% Remove polynomial fit from the data -> default is demeaning
if polyorder >= 0
  dat = ft_preproc_polyremoval(dat, polyorder, 1, ndatsample);
end

% Determine fsample and set total time-length of data
fsample = 1/(time(2)-time(1));
dattime = ndatsample / fsample; % total time in seconds of input data

% Zero padding
if round(pad * fsample) < ndatsample
  error('the padding that you specified is shorter than the data');
end
if isempty(pad) % if no padding is specified padding is equal to current data length
  pad = dattime;
end
postpad = zeros(1,round((pad - dattime) * fsample));
endnsample = round(pad * fsample);  % total number of samples of padded data
endtime    = pad;            % total time in seconds of padded data

% Set freqboi and freqoi
if isnumeric(freqoi) % if input is a vector
  freqboi   = round(freqoi ./ (fsample ./ endnsample)) + 1;
  freqboi   = unique(freqboi);
  freqoi    = (freqboi-1) ./ endtime; % boi - 1 because 0 Hz is included in fourier output
elseif strcmp(freqoi,'all') % if input was 'all'
  freqboilim = round([0 fsample/2] ./ (fsample ./ endnsample)) + 1;
  freqboi    = freqboilim(1):1:freqboilim(2);
  freqoi     = (freqboi-1) ./ endtime;
end
% check for freqoi = 0 and remove it, there is no wavelet for freqoi = 0
if freqoi(1)==0
  freqoi(1)  = [];
  freqboi(1) = [];
end
nfreqboi = length(freqboi);
nfreqoi  = length(freqoi);

% Set timeboi and timeoi
offset = round(time(1)*fsample);
if isnumeric(timeoi) % if input is a vector
  timeboi  = round(timeoi .* fsample - offset) + 1;
  ntimeboi = length(timeboi);
  timeoi   = round(timeoi .* fsample) ./ fsample;
elseif strcmp(timeoi,'all') % if input was 'all'
  timeboi  = 1:length(time);
  ntimeboi = length(timeboi);
  timeoi   = time;
end

% Creating wavelets
% expand width to array if constant width
if numel(width) == 1
  width = ones(1,nfreqoi) * width;
end
wltspctrm = cell(nfreqoi,1);
for ifreqoi = 1:nfreqoi
  dt = 1/fsample;
  sf = freqoi(ifreqoi) / width(ifreqoi);
  st = 1/(2*pi*sf);
  toi2 = -gwidth*st:dt:gwidth*st;
  A = 1/sqrt(st*sqrt(pi));
  tap = (A*exp(-toi2.^2/(2*st^2)))';
  acttapnumsmp = size(tap,1);
  taplen(ifreqoi) = acttapnumsmp;
  ins = ceil(endnsample./2) - floor(acttapnumsmp./2);
  prezer = zeros(ins,1);
  pstzer = zeros(endnsample - ((ins-1) + acttapnumsmp)-1,1);
  
  % produce angle with convention: cos must always be 1  and sin must always be centered in upgoing flank, so the centre of the wavelet (untapered) has angle = 0
  ind  = (-(acttapnumsmp-1)/2 : (acttapnumsmp-1)/2)'   .*  ((2.*pi./fsample) .* freqoi(ifreqoi));
  
  % create wavelet and fft it
  wavelet = complex(vertcat(prezer,tap.*cos(ind),pstzer), vertcat(prezer,tap.*sin(ind),pstzer));
  wltspctrm{ifreqoi} = complex(zeros(1,endnsample));
  wltspctrm{ifreqoi} = fft(wavelet,[],1)';
  
  
  %%%% debug plotting
%   figure('name',['wavelet @ ' num2str(freqoi(ifreqoi)) 'Hz' ],'NumberTitle','off');
%   subplot(2,1,1);
%   hold on;
%   plot(real(wavelet));
%   plot(imag(wavelet),'color','r');
%   legend('real','imag');
%   tline = length(wavelet)/2;
%   if mod(tline,2)==0
%     line([tline tline],[-max(abs(wavelet)) max(abs(wavelet))],'color','g','linestyle','--')
%   else
%     line([ceil(tline) ceil(tline)],[-max(abs(wavelet)) max(abs(wavelet))],'color','g','linestyle','--');
%     line([floor(tline) floor(tline)],[-max(abs(wavelet)) max(abs(wavelet))],'color','g','linestyle','--');
%   end;
%   subplot(2,1,2);
%   plot(angle(wavelet),'color','g');
%   if mod(tline,2)==0,
%     line([tline tline],[-pi pi],'color','r','linestyle','--')
%   else
%     line([ceil(tline) ceil(tline)],[-pi pi],'color','r','linestyle','--')
%     line([floor(tline) floor(tline)],[-pi pi],'color','r','linestyle','--')
%   end
  %%%% debug plotting
  
end

% Compute fft
spectrum = complex(nan(nchan,nfreqoi,ntimeboi),nan(nchan,nfreqoi,ntimeboi));
datspectrum = transpose(fft(transpose([dat repmat(postpad,[nchan, 1])]))); % double explicit transpose to speedup fft
for ifreqoi = 1:nfreqoi
  fprintf('processing frequency %d (%.2f Hz)\n', ifreqoi,freqoi(ifreqoi));
  % compute indices that will be used to extracted the requested fft output
  nsamplefreqoi    = taplen(ifreqoi);
  reqtimeboiind    = find((timeboi >=  (nsamplefreqoi ./ 2)) & (timeboi < (ndatsample - (nsamplefreqoi ./2))));
  reqtimeboi       = timeboi(reqtimeboiind);
  
  % compute datspectrum*wavelet, if there are reqtimeboi's that have data
  if ~isempty(reqtimeboi)
    dum = fftshift(transpose(ifft(transpose(datspectrum .* repmat(wltspctrm{ifreqoi},[nchan 1])))),2); % double explicit transpose to speedup fft
    dum = dum .* sqrt(2 ./ fsample);
    spectrum(:,ifreqoi,reqtimeboiind) = dum(:,reqtimeboi);
  end
end
