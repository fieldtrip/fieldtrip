function [dat] = ft_preproc_derivative(dat, order, padding)

% FT_PREPROC_DERIVATIVE computes the temporal Nth order derivative of the
% data
%
% Use as
%   [dat] = ft_preproc_derivative(dat, order, padding)
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

