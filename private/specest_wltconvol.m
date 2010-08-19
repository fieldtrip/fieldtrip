function [spectrum,freqoi,timeoi] = specest_wltconvol(dat, time, varargin) 

% SPECEST_WLTCONVOL performs time-frequency analysis on any time series trial data using the 'wavelet method' based on Morlet wavelets,
% doing convolution in the time domain by multiplaction in the frequency domain
%
%
% Use as
%   [spectrum,freqoi,timoei] = specest_wltconvol(dat,time...)   
%
%   dat      = matrix of chan*sample 
%   time     = vector, containing time in seconds for each sample
%   spectrum = matrix of chan*freqoi*timeoi of fourier coefficients
%   freqoi   = vector of frequencies in spectrum
%   timeoi   = vector of timebins in spectrum
%
%
%
% Optional arguments should be specified in key-value pairs and can include:
%   pad        = number, total length of data after zero padding (in seconds)
%   freqoi     = vector, containing frequencies of interest                                           
%   timeoi     = vector, containing time points of interest (in seconds)
%   width      = 
%   gwidth     = 
%
%
%
%
%
% FFT SPEED NOT YET OPTIMIZED (e.g. matlab version, transpose or not)
% SHOULD FREQOI = 'ALL' BE REMOVED OR NOT?
%
%
% See also SPECEST_MTMCONVOL, SPECEST_TFR, SPECEST_HILBERT, SPECEST_MTMWELCH, SPECEST_NANFFT, SPECEST_MVAR, SPECEST_MTMCONVOL


% get the optional input arguments
keyvalcheck(varargin, 'optional', {'pad','width','gwidth','freqoi','timeoi'});
freqoi    = keyval('freqoi',      varargin);  if isempty(freqoi),   freqoi  = 'all';   end
timeoi    = keyval('timeoi',      varargin);  if isempty(timeoi),   timeoi  = 'all';   end
width     = keyval('width',       varargin);  if isempty(width),    width    = 7;      end
gwidth    = keyval('gwidth',      varargin);  if isempty(gwidth),   gwidth   = 3;      end
pad       = keyval('pad',         varargin);



% Set n's
[nchan,ndatsample] = size(dat);


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
prepad  = zeros(1,floor((pad - dattime) * fsample ./ 2));
postpad = zeros(1,ceil((pad - dattime) * fsample ./ 2));
endnsample = pad * fsample;  % total number of samples of padded data
endtime    = pad;            % total time in seconds of padded data



% Set freqboi and freqoi
if isnumeric(freqoi) % if input is a vector
  freqboi   = round(freqoi ./ (fsample ./ endnsample)) + 1;
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
% expand width to array if constant width
if numel(width) == 1
  width = ones(1,nfreqoi) * width;
end


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


% minoffset = min(data.offset);
% timboi = round(cfg.toi .* data.fsample - minoffset);
% toi = round(cfg.toi .* data.fsample) ./ data.fsample;
% numtoi = length(toi);
% numfoi = length(cfg.foi);



% Creating wavelets
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
  ind = (0:acttapnumsmp-1)' .* ((2.*pi./fsample) .* freqoi(ifreqoi));
  wltspctrm{ifreqoi} = complex(zeros(1,endnsample));
  wltspctrm{ifreqoi} = fft(complex(vertcat(prezer,tap.*cos(ind),pstzer), vertcat(prezer,tap.*sin(ind),pstzer)),[],1)';
end


% Compute fft
spectrum = complex(nchan,nfreqoi,ntimeboi);
datspectrum = fft([repmat(prepad,[nchan, 1]) dat repmat(postpad,[nchan, 1])],[],2);
for ifreqoi = 1:nfreqoi
  fprintf('processing frequency %d (%.2f Hz)\n', ifreqoi,freqoi(ifreqoi));
  for ichan = 1:nchan
    % compute indices that will be used to extracted the requested fft output
    nsamplefreqoi    = taplen(ifreqoi) .* fsample;
    reqtimeboiind    = find((timeboi >=  (nsamplefreqoi ./ 2)) & (timeboi <    ndatsample - (nsamplefreqoi ./2)));
    reqtimeboi       = timeboi(reqtimeboiind);
    
    % compute datspectrum*wavelet, if there are reqtimeboi's that have data
    if ~isempty(reqtimeboi)
      dum = fftshift(ifft(datspectrum(ichan,:) .* wltspctrm{ifreqoi},[],2)); % why is this fftshift necessary?
      spectrum(ichan,ifreqoi,reqtimeboiind) = dum(reqtimeboi);
    end
  end
end
















