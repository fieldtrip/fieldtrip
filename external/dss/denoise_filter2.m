function [params, s_new] = denoise_filter2(params, s, state)
% Generic filter function. Can be used as DSS denoising function.
%   [params, s_new] = denoise_filter2(params, s, state)
%     params  Defines the used filter
%       params.filter_conv  Inpulse response for convolution filter
%       params.filter_dct   Mask for DCT filter
%       params.filter_fft   Mask for FFT filter
%       none of above       Pass the signal unfiltered
%     state   DSS algorithm state (not used)
%     s       Source signal, matrix of row vector signals
%     s_new   Filtered signal
%
% This function has been adjusted from the original denoise_filter
% function in the DSS package, allowing for filtfilt filtering

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id: denoise_filter.m,v 1.15 2005/11/30 08:29:40 jaakkos Exp $

if nargin<2
    params.name = 'Generic filter';
    params.description = 'Generic filter function for convolutive, DCT, FFT, and filtfilt filtering';
    params.param = {'filter_conv', 'filter_dct', 'filter_fft' 'filter_filtfilt'};
    params.param_value ={[], [], []};
    params.param_type = {'vector', 'vector', 'vector'};
    params.param_desc = {'Convolution coefficients', 'DCT filter coefficients', 'FFT filter coefficients', 'filtfilt coefficients'};
    params.approach = {'pca', 'defl', 'symm'};
    params.beta = {'beta_global'};
    return;
end

if isfield(params, 'filter_conv')
  % -- Filter with impulse response
  %s_new = conv(s, repmat(params.filter_conv, size(s, 1), 1));
  %s_new = s_new(round(length(params.filter_conv)/2)+[1:length(s)]);
  s_new = convolution(s, repmat(params.filter_conv, size(s, 1), 1));
elseif isfield(params, 'filter_dct')
  % DCT filtering
  s_new = idct(repmat(params.filter_dct, size(s,1), 1)' .* dct(s'))';
elseif isfield(params, 'filter_fft')
  % FFT filtering
  s_new = real(ifft(repmat(params.filter_fft, size(s,1), 1)' .* fft(s'))');
elseif isfield(params, 'filter_filtfilt')
  if isfield(params.filter_filtfilt, 'function')
    fhandle = str2func(params.filter_filtfilt.function);
  else
    fhandle = @filtfilt;
  end
  if iscell(s)
    s_new = s;
    for k = 1:numel(s)
      tmp = bsxfun(@minus,s{k},mean(s{k},2));
      s_new{k} = fhandle(params.filter_filtfilt.B, params.filter_filtfilt.A, tmp')';
    end
  else
    tr_begin = params.tr_begin(:);
    tr_end   = params.tr_end(:);
    s_new    = zeros(size(s));
    for k = 1:numel(tr_begin)
      indx = tr_begin(k):tr_end(k);
      tmp  = s(:,indx);
      tmp  = tmp - mean(tmp,2)*ones(1,numel(indx));
      s_new(:,indx) = fhandle(params.filter_filtfilt.B, params.filter_filtfilt.A, tmp')';
    end
  end
elseif isfield(params, 'filter_bandpass')
  if iscell(s)
    s_new = s;
    dobp = size(params.filter_bandpass.bpfreq,1)==1;
    for k = 1:numel(s)
      tmp = bsxfun(@minus,s{k},mean(s{k},2));
      if dobp
        s_new{k} = ft_preproc_bandpassfilter(tmp, params.filter_bandpass.fsample, params.filter_bandpass.bpfreq,[],'firws');
      else
        % do a series of bandstops
        tmporig = tmp;
        for m = 1:size(params.filter_bandpass.bpfreq,1)
          tmp = ft_preproc_bandstopfilter(tmp, params.filter_bandpass.fsample, params.filter_bandpass.bpfreq(m,:),[],'firws');
        end
        s_new{k} = tmporig - tmp;
      end
    end
  else
    tr_begin = params.tr_begin(:);
    tr_end   = params.tr_end(:);
    s_new    = zeros(size(s));
    for k = 1:numel(tr_begin)
      indx = tr_begin(k):tr_end(k);
      tmp  = s(:,indx);
      tmp  = tmp - mean(tmp,2)*ones(1,numel(indx));
      s_new(:,indx) = ft_preproc_bandpassfilter(tmp, params.filter_bandpass.fsample, params.filter_bandpass.bpfreq);
    end
  end
elseif isfield(params, 'filter_highpass')
  tr_begin = params.tr_begin(:);
  tr_end   = params.tr_end(:);
  s_new    = zeros(size(s));
  for k = 1:numel(tr_begin)
    indx = tr_begin(k):tr_end(k);
    tmp  = s(:,indx);
    tmp  = tmp - mean(tmp,2)*ones(1,numel(indx));
    s_new(:,indx) = ft_preproc_highpassfilter(tmp, params.filter_highpass.fsample, params.filter_highpass.hpfreq,[],[],[],false);
  end
elseif isfield(params, 'filter_lowpass')
  tr_begin = params.tr_begin(:);
  tr_end   = params.tr_end(:);
  s_new    = zeros(size(s));
  for k = 1:numel(tr_begin)
    indx = tr_begin(k):tr_end(k);
    tmp  = s(:,indx);
    mtmp = mean(tmp,2);
    tmp  = tmp - mtmp*ones(1,numel(indx));
    tmp  = ft_preproc_lowpassfilter(tmp, params.filter_lowpass.fsample, params.filter_lowpass.lpfreq) + mtmp*ones(1,numel(indx));
    tmp(tmp<1e12*eps) = 1e12*eps;
    s_new(:,indx) = tmp;
  end
elseif isfield(params, 'filter_boxcar')
  if iscell(s)
    s_new = cell(size(s));
    for k = 1:numel(s)
      s_new{k} = ft_preproc_smooth(s{k}, params.filter_boxcar);
    end
  else
    tr_begin = params.tr_begin(:);
    tr_end   = params.tr_end(:);
    s_new    = zeros(size(s));
    nbox     = params.filter_boxcar;
    krn      = ones(1,nbox)./nbox;
    for k = 1:numel(tr_begin)
      indx = tr_begin(k):tr_end(k);
      tmp  = s(:,indx);
      mtmp = mean(tmp,2);
      tmp  = tmp - mtmp*ones(1,numel(indx));
      tmp  = convn(tmp,krn,'same') + mtmp*ones(1,numel(indx));
      tmp(:,[1:round(nbox/2) (end-round(nbox/2)):end]) = mtmp;
      s_new(:,indx) = tmp;
    end
  end
else 
  % -- No filtering
  s_new = s;
end

% --------
% Calculate convolution for multiple row vector pairs and
% cut tails based on length of B.
function R = convolution(A, B)

R = zeros(size(A, 1), size(A,2)+size(B,2)-1);
for i=1:size(A,1)
  R(i,:) = conv(A(i,:), B(i,:));
end
R = R(:,round(size(B,2)/2)+[1:size(A,2)]);

% fir_filterdcpadded() - Pad data with DC constant and filter
%
% Usage:
%   >> data = fir_filterdcpadded(b, a, data, causal);
%
% Inputs:
%   b             - vector of filter coefficients
%   a             - 1
%   data          - raw data (times x chans)
%   causal        - boolean perform causal filtering {default 0}
%   usefftfilt    - boolean use fftfilt instead of filter
%
% Outputs:
%   data      - smoothed data
%
% Note:
%   firfiltdcpadded always operates (pads, filters) along first dimension.
%   Not memory optimized.
%
% Author: Andreas Widmann, University of Leipzig, 2014

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2013 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
% $Id$

function [ data ] = fir_filterdcpadded(b, a, data, causal, usefftfilt)

% Defaults
if nargin <= 4 || isempty(usefftfilt)
    usefftfilt = 0;
end
if nargin <= 3 || isempty(causal)
    causal = 0;
end

% Check arguments
if nargin < 3
    ft_error('Not enough input arguments.');
end

% Is FIR?
if ~isscalar(a) || a ~= 1
    ft_error('Not a FIR filter. onepass-zerophase and onepass-minphase filtering is available for FIR filters only.')
end

% Group delay
if mod(length(b), 2) ~= 1
    ft_error('Filter order is not even.');
end
groupDelay = (length(b) - 1) / 2;

% Filter symmetry
isSym = all(b(1:groupDelay) == b(end:-1:groupDelay + 2));
isAntisym = all([b(1:groupDelay) == -b(end:-1:groupDelay + 2) b(groupDelay + 1) == 0]);
if causal == 0 && ~(isSym || isAntisym)
    ft_error('Filter is not anti-/symmetric. For onepass-zerophase filtering the filter must be anti-/symmetric.')
end

% Padding
if causal
    startPad = repmat(data(1, :), [2 * groupDelay 1]);
    endPad = [];
else
    startPad = repmat(data(1, :), [groupDelay 1]);
    endPad = repmat(data(end, :), [groupDelay 1]);
end

% Filter data (with double precision)
isSingle = isa(data, 'single');

if usefftfilt
    data = fftfilt(double(b), double([startPad; data; endPad]));
else
    data = filter(double(b), 1, double([startPad; data; endPad])); % Pad and filter with double precision
end

% Convert to single
if isSingle
    data = single(data);
end

% Remove padded data
data = data(2 * groupDelay + 1:end, :);
