function [filt] = ft_preproc_dftfilter(dat, Fs, Fl, varargin)

% FT_PREPROC_DFTFILTER reduces power line noise (50 or 60Hz) via two 
% alternative methods:
% A) DFT filter (Flreplace = 'zero') or
% B) Spectrum Interpolation (Flreplace = 'neighbour').
%
% A) The DFT filter applies a notch filter to the data to remove the 50Hz
% or 60Hz line noise components ('zeroing'). This is done by fitting a sine 
% and cosine at the specified frequency to the data and subsequently 
% subtracting the estimated components. The longer the data is, the sharper 
% the spectral notch will be that is removed from the data.
% Preferably the data should have a length that is a multiple of the
% oscillation period of the line noise (i.e. 20ms for 50Hz noise). If the
% data is of different lenght, then only the first N complete periods are
% used to estimate the line noise. The estimate is subtracted from the
% complete data.
%
% B) Alternatively line noise is reduced via spectrum interpolation 
% (introduced by Mewett et al., 2004, Med. Biol. Eng. Comput. 42, 
% doi:10.1007/BF02350994). 
% The signal is:
% I)   transformed into the frequency domain via a discrete Fourier 
%       transform (DFT), 
% II)  the line noise component (e.g. 50Hz, Flwidth = 1 (±1Hz): 49-51Hz) is 
%       interpolated in the amplitude spectrum by replacing the amplitude 
%       of this frequency bin by the mean of the adjacent frequency bins 
%       ('neighbours', e.g. 49Hz and 51Hz). 
%       Neighwidth defines frequencies considered for the mean (e.g. 
%       Neighwidth = 2 (±2Hz) implies 47-49 Hz and 51-53 Hz). 
%       The original phase information of the noise frequency bin is
%       retained.
% III) the signal is transformed back into the time domain via inverse DFT
%       (iDFT).
% If Fline is a vector (e.g. [50 100 150]), harmonics are also considered. 
% Preferably the data should be continuous or consist of long data segments
% (several seconds) to avoid edge effects. If the sampling rate and the
% data length are such, that a full cycle of the line noise and the harmonics
% fit in the data and if the line noise is stationary (e.g. no variations
% in amplitude or frequency), then spectrum interpolation can also be 
% applied to short trials. But it should be used with caution and checked 
% for edge effects.
%
% Use as
%   [filt] = ft_preproc_dftfilter(dat, Fsample, Fline, varargin)
% where
%   dat             data matrix (Nchans X Ntime)
%   Fsample         sampling frequency in Hz
%   Fline           line noise frequency (and harmonics)
%
% Additional input arguments come as key-value pairs:
%
%   Flreplace       'zero' or 'neighbour', method used to reduce line noise, 'zero' implies DFT filter, 'neighbour' implies spectrum interpolation  
%   Flwidth         bandwidth of line noise frequencies, applies to spectrum interpolation, in Hz
%   Neighwidth      width of frequencies neighbouring line noise frequencies, applies to spectrum interpolation (Flreplace = 'neighbour'), in Hz 
%
% The line frequency should be specified as a single number for the DFT filter.
% If omitted, a European default of 50Hz will be assumed
%
% See also PREPROC

% Undocumented option:
%   Fline can be a vector, in which case the regression is done for all
%   frequencies in a single shot. Prerequisite is that the requested
%   frequencies all fit with an integer number of cycles in the data.
%
% Copyright (C) 2003, Pascal Fries
% Copyright (C) 2003-2015, Robert Oostenveld
% Copyright (C) 2016, Sabine Leske 
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

% defaults
Flreplace = ft_getopt(varargin,'dftreplace','zero');
Flwidth = ft_getopt(varargin, 'dftbandwidth', [1 2 3]);
Neighwidth = ft_getopt(varargin, 'dftneighbourwidth', [2 2 2]);

% determine the size of the data
[nchans, nsamples] = size(dat);

% set the default filter frequency
if nargin<3 || isempty(Fl)
  Fl = 50;
end

% ensure to be a column  vector
Fl = Fl(:);

% preprocessing fails on channels that contain NaN
if any(isnan(dat(:)))
  ft_warning('FieldTrip:dataContainsNaN', 'data contains NaN values');
end

% Method A): DFT filter 
if strcmp(Flreplace,'zero')
    % determine the largest integer number of line-noise cycles that fits in the data
    n = round(floor(nsamples .* Fl./Fs) * Fs./Fl);
    if all(n==n(1))
        % make a selection of samples such that the line-noise fits the data
        sel = 1:n(1);
    else
        % the different frequencies require different numbers of samples, apply the filters sequentially
        filt = dat;
        for i=1:numel(Fl)
            optarg = varargin;
            filt = ft_preproc_dftfilter(filt, Fs, Fl(i), optarg{:});
        end
        return
    end
    
    % temporarily remove mean to avoid leakage
    meandat = mean(dat(:,sel),2);
    for i=1:nsamples
        % demean the data
        dat(:,i) = dat(:,i) - meandat;
    end
    
    % fit a sine and cosine to each channel in the data and subtract them
    time = (0:nsamples-1)/Fs;
    tmp  = exp(1i*2*pi*Fl*time);                   % complex sin and cos
    ampl = 2*dat(:,sel)/tmp(:,sel);                % estimated amplitude of complex sin and cos on integer number of cycles
    est  = ampl*tmp;                               % estimated signal at this frequency
    filt = dat - est;                              % subtract estimated signal
    filt = real(filt);
    
    for i=1:nsamples
        % add the mean back to the filtered data
        filt(:,i) = filt(:,i) + meandat;
    end

% Method B): Spectrum Interpolation
elseif strcmp(Flreplace,'neighbour')
   Flwidth = Flwidth(:);
   Neighwidth = Neighwidth(:);
   
   if (length(Fl) ~= length(Flwidth)) || (length(Fl) ~= length(Neighwidth))
       error('The number of frequencies to interpolate (cfg.dftfreq) should be the same as the number of bandwidths (cfg.dftbandwidth) and bandwidths of neighbours (cfg.neighbourwidth)');
   end

    % frequencies to interpolate
    for i = 1:length(Fl)
        f2int(i,:) = [Fl(i)-Flwidth(i) Fl(i)+Flwidth(i)];
    end
    % frequencies used for interpolation
    for i = 1:length(Neighwidth)
        f4int(i,:) = [f2int(i,1)-Neighwidth(i) f2int(i,:) f2int(i,2)+Neighwidth(i)];
    end
    
    data_fft = fft(dat,nsamples,2); % calculate fft to obtain spectrum that will be interpolated
    frq = Fs*linspace(0,1,nsamples+1);
    
    % interpolate 50Hz (and harmonics) amplitude in spectrum
    for i = 1:length(Fl)
        smpl2int = nearest(frq,f2int(i,1)):nearest(frq,f2int(i,2));                                                   % samples of frequencies that will be interpolated
        smpl4int = [(nearest(frq,f4int(i,1)):nearest(frq,f4int(i,2))-1),(nearest(frq,f4int(i,3))+1:nearest(frq,f4int(i,4)))]; % samples of neighbouring frequencies used to calculate the mean
        
        % new amplitude is calculated as the mean of the neighbouring frequencies
        mns4int= bsxfun(@times, ones(size(data_fft(:,smpl2int))), mean(abs(data_fft(:,smpl4int)),2));
        
        % Eulers formula: replace noise components with new mean amplitude combined with phase, that is retained from the original data
        data_fft(:,smpl2int) = bsxfun(@times, exp(bsxfun(@times,angle(data_fft(:,smpl2int)),1i)), mns4int); 
    end
    
    % complex fourier coefficients are transformed back into time domin, fourier coefficients are treated as conjugate 'symmetric'
    % to ensure a real valued signal after iFFT
    filt = ifft(data_fft,[],2,'symmetric');
   
end

