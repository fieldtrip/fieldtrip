function [dat] = preproc_hilbert(dat, option)

% PREPROC_HILBERT computes the Hilbert transpose of the data and optionally
% performs post-processing on the complex representation, e.g. the absolute
% value of the Hilbert transform of a band-pass filtered signal corresponds
% with the amplitude envelope.
%
% Use as
%   [dat] = preproc_hilbert(dat, option)
% where
%   dat        data matrix (Nchans X Ntime)
%   option     string that determines whether and how the Hilbert transform
%              should be post-processed, can be
%                'abs'
%                'complex'
%                'real'
%                'imag'
%                'absreal'
%                'absimag'
%                'angle'
%
% The default is to return the absolute value of the Hilbert transform.
%
% See also PREPROC

% Copyright (C) 2008, Robert Oostenveld
%
% $Log: preproc_hilbert.m,v $
% Revision 1.2  2008/05/23 09:13:58  roboos
% cleaned up code and documentation, ensure that all functions are consistent, added proper implementation to the scratch functions
%

% set the defaults if option is not specified
if nargin<2 || isempty(option)
  option = 'abs';
end

% use the non-conjugate transpose to be sure
dat = transpose(hilbert(transpose(dat)));

% do postprocessing of the complex representation
switch option
  case {'yes' 'abs'}
    dat = abs(dat);   % this is the default if 'yes' is specified
  case {'no' 'complex'}
    dat = dat;        % this is the default if 'no' is specified
  case 'real'
    dat = real(dat);
  case 'imag'
    dat = imag(dat);
  case 'absreal'
    dat = abs(real(dat));
  case 'absimag'
    dat = abs(imag(dat));
  case 'angle'
    dat = unwrap(angle(dat));
  otherwise
    error('incorrect specification of the optional input argument');
end
