function [all_IAC, v, n] = ft_connectivity_iac(input, varargin)

% FT_CONNECTIVITY_IAC computes instantaneous amplitude coupling (IAC)
% between pairs of channels from Fourier-domain MEG/EEG data.
%
%
% This implementation follows amplitude-based
% connectivity (IAC) based on:
%
% Wei, H.T., Francois-Nienaber, A., Deschamps, T., Bellana, B.,
% Hebscher, M., Sivaratnam, G., Zadeh, M., & Meltzer, J.A. (2021).
% Sensitivity of amplitude and phase based MEG measures of
% interhemispheric connectivity during unilateral finger movements.
% NeuroImage, 242, 118457.
% https://doi.org/10.1016/j.neuroimage.2021.118457
%
% FieldTrip integration by:
% monaragala, K. (2026) 
%
% IAC is motivated by amplitude envelope correlation:
% when two signals exhibit high amplitudes at the same time, their
% amplitude envelopes covary. A time-resolved (instantaneous) measure
% can be formed from the pointwise product (Hadamard product) of the
% two amplitude-envelope time series, followed by averaging
% over time/trials to obtain a scalar coupling estimate.
%
% Data Input:
%
% The input data must be complex Fourier spectra organized as:
%
% Repetitions (rpttap) x Channel x Frequency x Time
%
% with dimord = 'rpttap_chan_freq_time'.
%
% The first dimension is treated as repetitions (e.g. trials).
%
% Notes:
% - This function operates directly on Fourier-domain input and does not
%   compute cross-spectra internally.
% - If input contains series of zeros in |X|, divisions in the orthogonalization step
%   may produce NaNs/Infs.
%
% See also FT_CONNECTIVITYANALYSIS 
%
% Copyright (C) 2009-2010 Donders Institute, Jan-Mathijs Schoffelen
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


hasjack   = ft_getopt(varargin, 'hasjack', 0);
feedback  = ft_getopt(varargin, 'feedback', 'none');
dimord    = ft_getopt(varargin, 'dimord');
powindx   = ft_getopt(varargin, 'powindx');
normalize = ft_getopt(varargin, 'normalize', 'no');
nbin      = ft_getopt(varargin, 'nbin');

disp('dimord:')
disp(dimord)

%%

if isempty(dimord)
  ft_error('input parameters should contain a dimord');
end

if ~strcmp(dimord, 'rpttap_chan_freq_time')
ft_error([ ...
  'Expected dimord ''rpttap_chan_freq_time'', got ''%s''.\n' ...
  'Please ensure input data follows rpttap_chan_freq_time for FieldTrip to process it correctly.' ...
], dimord);
end

siz = size(input);
n = siz(1);
v = [];

outsum = zeros(siz(2:end));
outssq = zeros(siz(2:end));
pvec   = [2 setdiff(1:numel(siz),2)];

ft_progress('init', feedback, 'computing metric...');

figure(111)
imagesc(squeeze(abs(squeeze(mean(input(:,4,:,:),1)))))

for ch1 = 1:size(input, 2)
    for ch2 = 1:size(input, 2)
        %Repetitions x Channel (x Frequency) (x Time)
        HEM1 = squeeze(input(:,ch1,:,:) - mean(input(:,ch1,:,:),4));
        
        HEM2 = squeeze(input(:,ch2,:,:) - mean(input(:,ch2,:,:),4));
        
        ortho_HEM1HEM2 = [];
        ortho_HEM2HEM1 = [];
        
        % perform orthoganalization between M1s
        for tr = 1:size(HEM1,1) %trials/epochs
            for f = 1:size(HEM1,2) %frequency
                for t = 1:size(HEM1,3) %time
        
                 ortho_HEM1HEM2(tr,f,t) = imag(HEM1(tr,f,t).*conj(HEM2(tr,f,t))/abs(HEM2(tr,f,t)));
                 ortho_HEM2HEM1(tr,f,t) = imag(HEM2(tr,f,t).*conj(HEM1(tr,f,t))/abs(HEM1(tr,f,t)));
        
                end
            end
        end
        
        %% compute envelopes for M1s anterior
        env_ortho_HEM1HEM2 = log(abs(ortho_HEM1HEM2).^2);  % [trial x freqs x time]
        env_HEM1 = log(abs(HEM1).^2);
        env_HEM2 = log(abs(HEM2).^2);
        env_ortho_HEM2HEM1 = log(abs(ortho_HEM2HEM1).^2);  % [trial x freqs x time]
        
        %% permute to put time and trial as the last two dimensions
        perm_env_ortho_HEM1HEM2 = permute(env_ortho_HEM1HEM2,[2 1 3]);
        perm_env_HEM1 = permute(env_HEM1,[2 1 3]);
        perm_env_HEM2 = permute(env_HEM2,[2 1 3]);
        perm_env_ortho_HEM2HEM1 = permute(env_ortho_HEM2HEM1,[2 1 3]);
        
        %%  pool the last two dimensions as one
        pool_env_ortho_HEM1HEM2 = perm_env_ortho_HEM1HEM2(:,:); 
        pool_env_HEM1 = perm_env_HEM1(:,:);
        pool_env_HEM2 = perm_env_HEM2(:,:);
        pool_env_ortho_HEM2HEM1 = perm_env_ortho_HEM2HEM1(:,:); 
        
        %% take the mean and std of the last dimension (across time and trial)
        env_ortho_HEM1HEM2_mean = mean(pool_env_ortho_HEM1HEM2,2);
        env_HEM1_mean = mean(pool_env_HEM1,2);
        env_HEM2_mean = mean(pool_env_HEM2,2);
        env_ortho_HEM2HEM1_mean = mean(pool_env_ortho_HEM2HEM1,2);
        
        env_ortho_HEM1HEM2_std = std(pool_env_ortho_HEM1HEM2,0,2);
        env_HEM1_std = std(pool_env_HEM1,0,2);
        env_HEM2_std = std(pool_env_HEM2,0,2);
        env_ortho_HEM2HEM1_std = std(pool_env_ortho_HEM2HEM1,0,2);
        
        
        %% take z-score
        for tr = 1:size(env_ortho_HEM1HEM2,1)
            for t = 1:size(env_ortho_HEM1HEM2,3)
                env_ortho_HEM1HEM2(tr,:,t) = (squeeze(env_ortho_HEM1HEM2(tr,:,t)) - permute(env_ortho_HEM1HEM2_mean,[2 1]))./permute(env_ortho_HEM1HEM2_std,[2 1]);
                env_HEM1z(tr,:,t) = (squeeze(env_HEM1(tr,:,t)) - permute(env_HEM1_mean,[2 1]))./permute(env_HEM1_std,[2 1]);
                env_HEM2z(tr,:,t) = (squeeze(env_HEM2(tr,:,t)) - permute(env_HEM2_mean,[2 1]))./permute(env_HEM2_std,[2 1]);
                env_ortho_HEM2HEM1(tr,:,t) = (squeeze(env_ortho_HEM2HEM1(tr,:,t)) - permute(env_ortho_HEM2HEM1_mean,[2 1]))./permute(env_ortho_HEM2HEM1_std,[2 1]);
            end
        end
        
        
        %% now do element-wise multiplication between left and right M1 envelopes
        IAC_ortho_HEM1HEM2 = env_ortho_HEM1HEM2.*env_HEM1z;  % [trial x freqs x time]
        IAC_ortho_HEM2HEM1 = env_ortho_HEM2HEM1.*env_HEM2z;  % [trial x freqs x time]
        
        %% average IAC over all trials
        IAC_ortho_HEM1HEM2 = squeeze(mean(IAC_ortho_HEM1HEM2,1));
        IAC_ortho_HEM2HEM1 = squeeze(mean(IAC_ortho_HEM2HEM1,1));
        
        
        % average across orthogonalization directions
        IAC_HEM_1221 = (IAC_ortho_HEM1HEM2 + IAC_ortho_HEM2HEM1)./2;
        
        all_IAC(ch1, ch2, :, :) = IAC_HEM_1221;

    end
end

end