function [freq, state] = specest_tfr(dat, fsample, foi, varargin)

state         = keyval('state',         varargin);
waveletwidth  = keyval('waveletwidth',  varargin); if isempty(waveletwidth), waveletwidth = 7; end
downsample    = keyval('downsample',    varargin);
time          = keyval('time',          varargin);
toi           = keyval('toi',           varargin);

% FIXME add width

if ~isempty(downsample) && ~isempty(toi)
  error('the downsample and toi options are mutually exclusive');
end

nchans   = size(dat,1);
nsamples = size(dat,2);

if isempty(time)
  if ~isempty(toi)
    error('specification of toi without time is not permitted');
  end
  time = (0:(nsamples-1))/fsample;
end

if ~isempty(state) && isequal(state.fsample, fsample)  && isequal(state.foi, foi) && isequal(state.waveletwidth, waveletwidth)
  % reuse the wavelet family from the previous state
  disp('using previous state');
  M = state.M;
else
  if ~isempty(state)
    state = [];
    warning('recomputing the state');
  end
  % recompute the wavelet family
  M = waveletfam(foi,fsample,waveletwidth);
end

% compute freq by convolving the wavelets with the data
nfreq = length(foi);
freq = zeros(nchans, nsamples, nfreq);

for f=1:nfreq
  wavelet = M{f};
  for c=1:nchans
    freq(c,:,f) = conv(dat(c,:),  wavelet, 'same');
  end
  % pad the edges with nans to indicate that the wavelet was not fully immersed in the data
  pad = ceil(length(M{f})/2);
  % the padding should not be longer than the actual data
  pad = min(pad, nsamples);
  begpad = 1:pad;
  endpad = (nsamples-pad+1):nsamples;
  freq(:,begpad,f) = nan;
  freq(:,endpad,f) = nan;
end

if ~isempty(downsample)
  % this would be done in case of downsample
  tbin = 1:downsample:nsamples;
elseif ~isempty(toi)
  % this would be done in case of toi/time specification
  tbin = zeros(size(toi));
  for i=1:length(tbin)
    tbin = nearest(time, toi(i));
  end
else
  tbin = 1:nsamples;
end

% select the samples for the output
freq = freq(:,tbin,:);

% remember the wavelets so that they can be reused in a subsequent call
state.fsample = fsample;
state.foi     = foi;
state.waveletwidth = waveletwidth;
state.M       = M;

% %convolves data in chan by time matrix (from one trial) with 1 wavelet for
% %each specified frequency of interest
% freq = zeros(size(data,1),length(foi), size(data,2));
% for k=1:size(data,1) %nchans
%   for j=1:length(foi)
%     cTmp = conv(data(k,:),M{j});
%     cTmp = 2*(abs(cTmp).^2)/fsample;
%     cTmp = cTmp(ceil(length(M{j})/2):length(cTmp)-floor(length(M{j})/2));
%     freq(k,j,:) = cTmp;
%   end
%
% end

%output should be chan by freq by time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for waveletanalysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = waveletfam(foi,fsample,waveletwidth)
dt = 1/fsample;
for k=1:length(foi)
  sf  = foi(k)/waveletwidth;
  st  = 1/(2*pi*sf);
  toi = -3.5*st:dt:3.5*st;
  A   = 1/sqrt(st*sqrt(pi));
  M{k}= A*exp(-toi.^2/(2*st^2)).*exp(i*2*pi*foi(k).*toi);
end