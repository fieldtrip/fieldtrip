function [spectrum, freqoi] = ft_specest_neuvar(dat, time, varargin)

% FT_SPECEST_NEUVAR computes a time-domain estimation of overall signal 
% power, having compensated for the 1/f distribution of spectral content.
%
% Use as
%   [spectrum,ntaper,freqoi] = ft_specest_neuvar(dat,time...)
% where
%   dat        = matrix of chan*sample
%   time       = vector, containing time in seconds for each sample
%   neuvar     = matrix of chan*neuvar
%
% Optional arguments should be specified in key-value pairs and can include
%   order      = number, the order of differentation for compensating for the 1/f (default: 1)
%   pad        = number, total length of data after zero padding (in seconds)
%   padtype    = string, indicating type of padding to be used (see ft_preproc_padding, default: 0)
%   verbose    = output progress to console (0 or 1, default 1)
%
% See also FT_FREQANALYSIS, FT_SPECEST_MTMFFT, FT_SPECEST_MTMCONVOL, FT_SPECEST_TFR, FT_SPECEST_HILBERT, FT_SPECEST_WAVELET

% Copyright (C) 2016, Arjen Stolk
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
order     = ft_getopt(varargin, 'order', 1);
pad       = ft_getopt(varargin, 'pad');
padtype   = ft_getopt(varargin, 'padtype', 'zero');
fbopt     = ft_getopt(varargin, 'feedback');
verbose   = ft_getopt(varargin, 'verbose', true);

if isempty(fbopt)
  fbopt.i = 1;
  fbopt.n = 1;
end

% this does not work on integer data
dat = cast(dat, 'double');

% Set n's
[nchan,ndatsample] = size(dat);

% This does not work on integer data
if ~isa(dat, 'double') && ~isa(dat, 'single')
  dat = cast(dat, 'double');
end

% Determine fsample and set total time-length of data
fsample = 1./mean(diff(time));
dattime = ndatsample / fsample; % total time in seconds of input data

% Zero padding
if round(pad * fsample) < ndatsample
  ft_error('the padding that you specified is shorter than the data');
end
if isempty(pad) % if no padding is specified padding is equal to current data length
  pad = dattime;
end
postpad    = ceil((pad - dattime) * fsample);
endnsample = round(pad * fsample);  % total number of samples of padded data
endtime    = pad;                   % total time in seconds of padded data

% Filter throught a differentiatior to compensate for the 1/f
dat = ft_preproc_derivative(dat, order);

% Calculate overall power using variance
str = sprintf('neuronal variance: %d samples, datalength: %d samples',endnsample,ndatsample);
[st, cws] = dbstack;
if length(st)>1 && strcmp(st(2).name, 'ft_freqanalysis')
  % specest_mtmfft has been called by ft_freqanalysis, meaning that ft_progress has been initialised
  ft_progress(fbopt.i./fbopt.n, ['processing trial %d/%d ',str,'\n'], fbopt.i, fbopt.n);
elseif verbose
  fprintf([str, '\n']);
end
spectrum = var(ft_preproc_padding(dat, padtype, 0, postpad), [], 2)'; % freq X chan
spectrum = sqrt(spectrum); % take square root since cfg.output = 'pow' assumes the output is amplitude and squares it 
freqoi = 0; % assign to DC frequency
