function [p] = ft_connectivity_plm(input, varargin)

% FT_CONNECTIVITY_PLM computes the phase linearity measurement from a cell
% array of time-domain data, where each cell is an epoch. This function implements
% the metric described in Baselice et al. "Phase Linearity Measurement:
% a novel index for brain functional connectivity", IEEE Transactions
% on Medical Imaging, 2018. Please reference the paper in case of use.
%
% Use as
%   [p] = ft_connectivity_plm(input, ...)
%
% The input data input should be organized as a cell-array, one element for each epoch.
% Each cell element should be a matrix of of nchan x nsamples values.
%
% Additional optional input arguments come as key-value pairs:
%   bandwidth	=	scalar, half-bandwidth parameter: the frequency range
%			across which to integrate
%   fsample     =       sampling frequency, needed to convert bandwidth to number of bins
%
% The output p contains the phase linearity measurement in the [0, 1] interval.
% The output p is organized as a 3D matrix of nepoch x nchan x  nchan dimensions.
%
% See also FT_CONNECTIVITYANALYSIS

% Copyright (C) 2018, Fabio Baselice, Pierpaolo Sorrentino, Jan-Mathijs Schoffelen
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

% the sequence of steps is as follows:
%  - Hilbert transformation
%  - multiply with complex conjugate
%  - fft
%  - remove volume conduction component
%  - integrate over bandwidth


% NOTE BY JM: if the user inputs data with different length trials, the fft per trial is going
% to have different frequency resolutions, which is not good. Better to throw an error in that
% case.
fs = ft_getopt(varargin, 'fsample');
B = ft_getopt(varargin, 'bandwidth');
if isempty(fs)
  error('sampling rate is not defined');
end
if isempty(B)
  warning('bandwidth parameter is not defined, assumed 1Hz');
  B=1;
end

nsmp = cellfun('size', input, 2);
assert(all(nsmp==nsmp(1)), 'currently there is no support for input, where the trials are of different length');

nrpt=numel(input);
for k = 1:numel(input)
  input{k} = hilbert(input{k}')';
end
% NOTE by JM: Is it expected that the data has been bandpassfiltered at
% this point? How would this be checked? 

nchan=size(input{1},1);
trial_length=size(input{1},2);
ph_min=0.1;        % Eps of Eq.(17) of the manuscript
f=(fs/trial_length)*(0:(trial_length-1));
f_integr=(abs(f)<B) | (abs(f-fs)<B);
p=zeros(nchan, nchan, nrpt);

for ktime=1:nrpt
  for kchan1=1:(nchan-1)
    for kchan2=(kchan1+1):nchan
      temp=fft(input{ktime}(kchan1,:).*conj(input{ktime}(kchan2,:)));    % NOTE BY FB: The inner cycle can be vectorized
      temp(1)=temp(1).*(abs(angle(temp(1)))>ph_min);  % Volume conduction suppression
      temp=(abs(temp)).^2;
      p_temp=sum(temp(f_integr))./sum(temp);
      p(kchan1, kchan2, ktime)=p_temp;
      p(kchan2, kchan1, ktime)=p_temp;
    end
  end
end

p = permute(p, [3 1 2]); % permute to adhere to the conventional matrix shape of the ft_connectivity_* codebase
