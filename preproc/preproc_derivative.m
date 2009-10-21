function [dat] = preproc_derivative(dat, order, padding)

% PREPROC_DERIVATIVE computes the temporal Nth order derivative of the
% data
%
% Use as
%   [dat] = preproc_derivative(dat, order, padding)
% where
%   dat        data matrix (Nchans X Ntime)
%   order      number representing the Nth derivative (default = 1)
%   padding    string that determines whether and how the data will be
%              padded to keep the number of samples the same, can be
%                'none'  do not apply padding, the output will be N samples shorter
%                'both'  apply padding to both sides
%                'beg'   apply padding at the beginning of the data
%                'end'   apply padding at the end of the data (default)
%
% See also PREPROC

% Copyright (C) 2008, Robert Oostenveld
%
% $Log: preproc_derivative.m,v $
% Revision 1.2  2008/05/23 09:13:58  roboos
% cleaned up code and documentation, ensure that all functions are consistent, added proper implementation to the scratch functions
%

% determine the size of the data
[Nchans, Nsamples] = size(dat);

% set the defaults if options are not specified
if nargin<2 || isempty(order)
  order = 1;
end
if nargin<3 || isempty(padding)
  padding = 'end';
end

% compute the derivative
dat = diff(dat, order, 2);

% pad the resulting data to keep the number of samples the same
switch padding
  case 'beg'
    dat = cat(2, zeros(Nchans, order), dat);
  case 'end'
    dat = cat(2, dat, zeros(Nchans, order));
  case 'both'
    if rem(order,2)
      error('padding can only be applied to both sides if the order is an even number');
    end
    dat = cat(2, zeros(Nchans, order/2), dat, zeros(Nchans, order/2));
end

