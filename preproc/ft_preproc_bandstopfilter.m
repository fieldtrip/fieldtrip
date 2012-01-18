function [filt] = ft_preproc_bandstopfilter(dat,Fs,Fbp,N,type,dir)

% FT_PREPROC_BANDSTOPFILTER applies a band-stop filter to the data and thereby
% removes the spectral components in the specified frequency band 
% 
% Use as
%   [filt] = ft_preproc_bandstopfilter(dat, Fsample, Fbp, N, type, dir)
% where
%   dat        data matrix (Nchans X Ntime)
%   Fsample    sampling frequency in Hz
%   Fbp        frequency band, specified as [Fhp Flp]
%   N          optional filter order, default is 4 (but) or dependent upon
%              frequency band and data length (fir/firls)
%   type       optional filter type, can be
%                'but' Butterworth IIR filter (default)
%                'fir' FIR filter using Matlab fir1 function 
%                'firls' FIR filter using Matlab firls function (requires Matlab Signal Processing Toolbox)
%   dir        optional filter direction, can be
%                'onepass'         forward filter only
%                'onepass-reverse' reverse filter only, i.e. backward in time
%                'twopass'         zero-phase forward and reverse filter (default)
%
% Note that a one- or two-pass filter has consequences for the
% strength of the filter, i.e. a two-pass filter with the same filter
% order will attenuate the signal twice as strong.
%
% See also PREPROC

% Copyright (c) 2007-2008, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

% set the default filter order later
if nargin<4 || isempty(N)
    N = [];
end

% set the default filter type
if nargin<5 || isempty(type)
  type = 'but';
end

% set the default filter direction
if nargin<6|| isempty(dir)
  dir = 'twopass';
end

% Nyquist frequency
Fn = Fs/2;

% compute filter coefficients 
switch type
  case 'but'
    if isempty(N)
      N = 4;
    end
    [B, A] = butter(N, [min(Fbp)/Fn max(Fbp)/Fn], 'stop');
  case 'fir'
    if isempty(N)
      N = 3*fix(Fs / Fbp(1));
      if rem(N,2)==1,   N=N+1;    end
    end
    if N > floor( (size(dat,2) - 1) / 3)
      N=floor(size(dat,2)/3) - 2;
      if rem(N,2)==1,   N=N+1;    end
    end
    [B, A] = fir1(N, [min(Fbp)/Fn max(Fbp)/Fn], 'stop');
  case 'firls' % from NUTMEG's implementation
    if isempty(N)
      N = 3*fix(Fs / Fbp(1));
    end
    if N > floor( (size(dat,2) - 1) / 3)
      N=floor(size(dat,2)/3) - 1;
    end
    
    f = 0:0.001:1;
    if rem(length(f),2)~=0
      f(end)=[];
    end
    z = ones(1,length(f));
    if(isfinite(min(Fbp)))
      [val,pos1] = min(abs(Fs*f/2 - min(Fbp)));
    else 
      [val,pos2] = min(abs(Fs*f/2 - max(Fbp)));
      pos1=pos2;
    end
    if(isfinite(max(Fbp)))
      [val,pos2] = min(abs(Fs*f/2 - max(Fbp)));
    else
      pos2 = length(f);
    end
    z(pos1:pos2) = 0;
    A = 1;
    B = firls(N,f,z); % requires Matlab signal processing toolbox
end

meandat=mean(dat,2);
dat=dat-repmat(meandat,[1 size(dat,2)]);
filt = filter_with_correction(B,A,dat,dir);
filt=filt+repmat(meandat,[1 size(dat,2)]);

%SK: I think the following is non-sense. Approximating a high-order
% bandstop filter by a succession of low-order bandstop filters
% will most likely give you very bad accuracy.

% check for filter instabilities and try to solve them
rangedat  = max(dat,[],2)  - min(dat,[],2);
rangefilt = max(filt,[],2) - min(filt,[],2);
result_instable = any(isnan(filt(:))) || (max(rangefilt)/max(rangedat)>2);
if result_instable && N>1
  warning('instable filter detected, applying two sequential filters');
  step1 = floor(N/2);  
  step2 = N - step1;
  % apply the filter in two steps, note that this is recursive
  filt = dat;
  filt = ft_preproc_bandstopfilter(filt,Fs,Fbp,step1,type,dir);
  filt = ft_preproc_bandstopfilter(filt,Fs,Fbp,step2,type,dir);
end
