function [y, opt] = specest_nanfft(x, varargin)

% SPECEST_NANFFT computes a fast Fourier transform in the presence of NaNs
% in the data
%
% Use as
%   [y] = specest_nanfft(x, ...)
%
% Optional arguments should be specified in key-value pairs and can include
%   basis      = precomputes set of basis functions (sines/cosines)
%   datataype  = 0, 1, 2

% Copyright (C) 2008, Robert Oostenveld
%
% $Log: specest_nanfft.m,v $
% Revision 1.4  2008/10/01 13:17:36  roboos
% correct number of channels for zeros, deal with case when all data is nan (thanks to Thilo)
%
% Revision 1.3  2008/10/01 11:24:13  roboos
% implemented spectral estimate for data with odd number of samples (sofar only an even number of samples would work)
%
% Revision 1.2  2008/10/01 08:19:45  roboos
% complete update, variable nan-locations are supported, using pseudo-inverse
%

% get the optional arguments
basis     = keyval('basis',     varargin);
datatype  = keyval('datatype',  varargin);

[m, n] = size(x);

if mod(n,2)==0
  % the number of samples is even
  k = n/2+1;
else
  % the number of samples is odd
  k = floor(n/2+1);
end

% determine the type of data and thereby the most suitable algorithm to use
nancount = sum(isnan(x), 1);
if isempty(datatype)
  if all(nancount==0)
    % there is no missing data
    datatype = 0;
  elseif all(nancount==0 | nancount==m)
    % the missing data is at the same location for all channels
    datatype = 1;
  else
    % the missing data is at different timepoints for different channels
    datatype = 2;
  end
end

if datatype==0
  % no basis functions are needed, because the standard FFT routine will be used

elseif datatype~=0 && isempty(basis)
  % create a seperate set of basis functions for the cosine and sine
  basis_c = zeros(k, n);
  basis_s = zeros(k, n);

  % create the time axis
  t = linspace(0, 2*pi, n+1);
  t = t(1:end-1);

  for w=1:k
    c = cos((w-1)*t);
    s = sin((w-1)*t);
    if w==1 || (w==(k) && mod(n,2)==0)
      % the normalization for the lowest (DC) and the highest frequency component is different
      s = s/(n);
      c = c/(n);
    else
      s = s/(n/2);
      c = c/(n/2);
    end
    basis_c(w,:) = c;
    basis_s(w,:) = s;
  end
  % concatenate the sine and cosine basis functions
  % leaving the first and last sine functions out, since those are all zero
  if mod(n,2)==0
    % the number of samples is even -> the last sine wave basis function is zero
    basis = cat(1, basis_c, basis_s(2:end-1, :));
  else
    % the number of samples is odd -> also include the last sine wave basis function
    basis = cat(1, basis_c, basis_s(2:end, :));
  end
end


switch datatype
  case 0
    % there is no missing data
    % use the standard FFT implementation
    y = fft(x, [], 2);

  case 1
    % the missing data is at the same location for all channels
    % remove that piece from the data and from the basis functions and use linear estimation
    
    keep = ~isnan(x(1,:));

    if all(~keep)
      % the data is all NaN, no reason to try to estimate the basis
      % functions
      y = nan(size(x));
      return
    end

    basis = basis(:,keep);
    x = x(:,keep);

    % do the linear estimation based on x=y*basis
    % y = x / basis;
    y = x * pinv(basis);

    % disentagle the estimated components

    if mod(n,2)==0
      % the number of samples is even -> the last sine wave basis function is zero
      sel1 = 1;       % lowest cosine, i.e. DC
      sel2 = 2:(k-1); % all cosines in between
      sel3 = k;       % highest cosine
      sel4 = (k+1):n; % all sines

      est1 = y(:,sel1);
      est2 = y(:,sel2);
      est3 = y(:,sel3);
      est4 = y(:,sel4);

      % combine the various estimates into a complex representation compatible with standard FFT
      y_real = cat(2, est1,        est2, est3,          fliplr(est2));
      y_imag = cat(2, zeros(m,1), -est4, zeros(m,1),    fliplr(est4));
      y = y_real + i*y_imag;

    else
      % the number of samples is odd -> also include the last sine wave basis function
      sel1 = 1;       % lowest cosine, i.e. DC
      sel2 = 2:k;     % all other cosines
      sel3 = (k+1):n; % all sines

      est1 = y(:,sel1);
      est2 = y(:,sel2);
      est3 = y(:,sel3);

      % combine the various estimates into a complex representation compatible with standard FFT
      y_real = cat(2, est1,        est2, fliplr(est2));
      y_imag = cat(2, zeros(m,1), -est3, fliplr(est3));
      y = y_real + i*y_imag;
    end

  case 2
    % the missing data is at different timepoints for different channels
    % use recursion to compute the nanfft for each channel
    y = zeros(size(x));
    for k=1:m
      y(k,:) = specest_nanfft(x(k,:), 'basis', basis);
    end

  otherwise
    error('unsupported configuration of NaNs in the data');
end
