function [params, s_new] = denoise_filter(params, s, state)
% Generic filter function. Can be used as DSS denoising function.
%   [params, s_new] = denoise_filter(params, s, state)
%     params  Defines the used filter
%       params.filter_conv  Inpulse response for convolution filter
%       params.filter_dct   Mask for DCT filter
%       params.filter_fft   Mask for FFT filter
%       none of above       Pass the signal unfiltered
%     state   DSS algorithm state (not used)
%     s       Source signal, matrix of row vector signals
%     s_new   Filtered signal

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

if nargin<3
    params.name = 'Generic filter';
    params.description = 'Generic filter function for convolutive, DCT and FFT filtering';
    params.param = {'filter_conv', 'filter_dct', 'filter_fft'};
    params.param_value ={[], [], []};
    params.param_type = {'vector', 'vector', 'vector'};
    params.param_desc = {'Convolution coefficients', 'DCT filter coefficients', 'FFT filter coefficients.'};
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
