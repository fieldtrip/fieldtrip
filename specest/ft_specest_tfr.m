function [spectrum, freqoi, timeoi] = ft_specest_tfr(dat, time, varargin)

% FT_SPECEST_TFR performs time-frequency analysis on any time series trial data using
% the 'wavelet method' based on Morlet wavelets, doing convolution in the time
% domain.
%
% Use as
%   [spectrum,freqoi,timeoi] = ft_specest_convol(dat,time,...)
% where
%   dat       = matrix of chan*sample
%   time      = vector, containing time in seconds for each sample
%   spectrum  = array of chan*freqoi*timeoi of fourier coefficients
%   freqoi    = vector of frequencies in spectrum
%   timeoi    = vector of timebins in spectrum
%
% Optional arguments should be specified in key-value pairs and can include
%   timeoi    = vector, containing time points of interest (in seconds, analysis window will be centered around these time points)
%   freqoi    = vector, containing frequencies (in Hz)
%   width     = number or vector, width of the wavelet, determines the temporal and spectral resolution (default = 7)
%   gwidth    = number, determines the length of the used wavelets in standard deviations of the implicit Gaussian kernel
%   verbose   = output progress to console (0 or 1, default 1)
%   polyorder = number, the order of the polynomial to fitted to and removed from the data prior to the fourier transform (default = 0 -> remove DC-component)
%
% See also FT_FREQANALYSIS, FT_SPECEST_MTMFFT, FT_SPECEST_MTMCONVOL, FT_SPECEST_HILBERT, FT_SPECEST_NANFFT, FT_SPECEST_WAVELET

% Copyright (C) 2010, Donders Institute for Brain, Cognition and Behaviour
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

% get the optional input arguments
freqoi    = ft_getopt(varargin, 'freqoi', 'all');
timeoi    = ft_getopt(varargin, 'timeoi', 'all');
width     = ft_getopt(varargin, 'width', 7);
gwidth    = ft_getopt(varargin, 'gwidth', 3);
pad       = ft_getopt(varargin, 'pad');
polyorder = ft_getopt(varargin, 'polyorder', 0);
fbopt     = ft_getopt(varargin, 'feedback');
verbose   = ft_getopt(varargin, 'verbose', true);

if isempty(fbopt),
  fbopt.i = 1;
  fbopt.n = 1;
end

% Set n's
[nchan,ndatsample] = size(dat);

% This does not work on integer data
typ = class(dat);
if ~strcmp(typ, 'double') && ~strcmp(typ, 'single')
  dat = cast(dat, 'double');
end

% Remove polynomial fit from the data -> default is demeaning
if polyorder >= 0
  dat = ft_preproc_polyremoval(dat, polyorder, 1, ndatsample);
end

% Determine fsample and set total time-length of data
fsample = 1./mean(diff(time));
dattime = ndatsample / fsample; % total time in seconds of input data

% Zero padding
if round(pad * fsample) < ndatsample
  error('the padding that you specified is shorter than the data');
end
if isempty(pad) % if no padding is specified padding is equal to current data length
  pad = dattime;
end
endnsample = round(pad * fsample);  % total number of samples of padded data
endtime    = pad;            % total time in seconds of padded data
prepad     = zeros(1,floor(((pad - dattime) * fsample) ./ 2));
postpad    = zeros(1,ceil(((pad - dattime) * fsample) ./ 2));

% Set freqboi and freqoi
freqoiinput = freqoi;
if isnumeric(freqoi) % if input is a vector
  freqboi   = round(freqoi ./ (fsample ./ endnsample)) + 1; % is equivalent to: round(freqoi .* endtime) + 1;
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

% throw a warning if input freqoi is different from output freqoi
if isnumeric(freqoiinput)
  % check whether padding is appropriate for the requested frequency resolution
  rayl = 1/endtime;
  if any(rem(freqoiinput,rayl)) % not always the case when they mismatch
    ft_warning('padding not sufficient for requested frequency resolution, for more information please see the FAQs on www.ru.nl/neuroimaging/fieldtrip');
  end
  if numel(freqoiinput) ~= numel(freqoi) % freqoi will not contain double frequency bins when requested
    ft_warning('output frequencies are different from input frequencies, multiples of the same bin were requested but not given');
  else
    if any(abs(freqoiinput-freqoi) >= eps*1e6)
      ft_warning('output frequencies are different from input frequencies');
    end
  end
end


% Set timeboi and timeoi
timeoiinput = timeoi;
offset = round(time(1)*fsample);
if isnumeric(timeoi) % if input is a vector
  timeoi   = unique(round(timeoi .* fsample) ./ fsample);
  timeboi  = round(timeoi .* fsample - offset) + 1;
  ntimeboi = length(timeboi);
elseif strcmp(timeoi,'all') % if input was 'all'
  timeboi  = 1:length(time);
  ntimeboi = length(timeboi);
  timeoi   = time;
end

% throw a warning if input timeoi is different from output timeoi
if isnumeric(timeoiinput)
  if numel(timeoiinput) ~= numel(timeoi) % timeoi will not contain double time-bins when requested
    ft_warning('output time-bins are different from input time-bins, multiples of the same bin were requested but not given');
  else
    if any(abs(timeoiinput-timeoi) >= eps*1e6) 
      ft_warning('output time-bins are different from input time-bins');
    end
  end
end


% Creating wavelets
% expand width to array if constant width
if numel(width) == 1
  width = ones(1,nfreqoi) * width;
end
wavelet = cell(nfreqoi,1);
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
  %prezer = zeros(ins,1);
  %pstzer = zeros(endnsample - ((ins-1) + acttapnumsmp)-1,1);
  
  % produce angle with convention: cos must always be 1  and sin must always be centered in upgoing flank, so the centre of the wavelet (untapered) has angle = 0
  ind  = (-(acttapnumsmp-1)/2 : (acttapnumsmp-1)/2)'   .*  ((2.*pi./fsample) .* freqoi(ifreqoi));
  
  % create wavelet and fft it
  %wavelet{ifreqoi} = complex(vertcat(prezer,tap.*cos(ind),pstzer), vertcat(prezer,tap.*sin(ind),pstzer));
  wavelet{ifreqoi} = complex(vertcat(tap.*cos(ind)), vertcat(tap.*sin(ind)));
  
  
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


% compute spectrum by convolving the wavelets with the data
spectrum = complex(nan(nchan,nfreqoi,ntimeboi),nan(nchan,nfreqoi,ntimeboi));
for ifreqoi = 1:nfreqoi
  str = sprintf('frequency %d (%.2f Hz)', ifreqoi,freqoi(ifreqoi));
  [st, cws] = dbstack;
  if length(st)>1 && strcmp(st(2).name, 'ft_freqanalysis') && verbose
    % specest_convol has been called by ft_freqanalysis, meaning that ft_progress has been initialised
    ft_progress(fbopt.i./fbopt.n, ['trial %d, ',str,'\n'], fbopt.i);
  elseif verbose
    fprintf([str, '\n']);
  end
  
  % compute indices that will be used to extracted the requested output (this keeps nans when the wavelet is not fully immersed in the data)
  nsamplefreqoi    = taplen(ifreqoi);
  reqtimeboiind    = find((timeboi >=  (nsamplefreqoi ./ 2)) & (timeboi < (ndatsample - (nsamplefreqoi ./2))));
  reqtimeboi       = timeboi(reqtimeboiind);
  
  % do convolution, if there are reqtimeboi's that have data
  if ~isempty(reqtimeboi)
    dum = complex(zeros(nchan,numel(reqtimeboi)));
    for ichan = 1:nchan
      dumconv = conv(dat(ichan,:),  wavelet{ifreqoi}, 'same');
      dum(ichan,:) = dumconv(reqtimeboi); % keeping nans nans when the wavelet is not fully immersed in the data
    end
    spectrum(:,ifreqoi,reqtimeboiind) = dum;
  end
end



% old code used to create wavelets
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % SUBFUNCTION for waveletanalysis
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function M = waveletfam(foi,fsample,waveletwidth)
% dt = 1/fsample;
% for k=1:length(foi)
%   sf  = foi(k)/waveletwidth;
%   st  = 1/(2*pi*sf);
%   toi = -3.5*st:dt:3.5*st;
%   A   = 1/sqrt(st*sqrt(pi));
%   M{k}= A*exp(-toi.^2/(2*st^2)).*exp(i*2*pi*foi(k).*toi);
% end

