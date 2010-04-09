function [dat] = ft_preproc_hilbert(dat, option)

% FT_PREPROC_HILBERT computes the Hilbert transpose of the data and optionally
% performs post-processing on the complex representation, e.g. the absolute
% value of the Hilbert transform of a band-pass filtered signal corresponds
% with the amplitude envelope.
%
% Use as
%   [dat] = ft_preproc_hilbert(dat, option)
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
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

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
        dat = (angle(dat./abs(dat)));
    case 'unwrap_angle'
        dat = unwrap(angle(dat./abs(dat)),[],2);
    otherwise
        error('incorrect specification of the optional input argument');
end
