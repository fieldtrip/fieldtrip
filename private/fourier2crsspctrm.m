function output = fourier2crsspctrm(cfg, freq)

% FOURIER2CRSSPCTRM transforms a fourier-containing freq-structure 
% into a crsspctrm-containing freq-structure, in which the
% powerspectra are also contained in the cross-spectra, being a 
% channelcombination of a channel with itself.
%
% Use as
%   [freq] = fourier2crsspctrm(cfg, freq)
%
% where you have the following configuration options:
%   cfg.channel    = cell-array with selection of channels,
%                    see CHANNELSELECTION for details
%   cfg.channelcmb = cell-array with selection of combinations between
%                    channels, see CHANNELCOMBINATION for details
%   cfg.keeptrials = 'yes' or 'no' (default)
%   cfg.foilim     = 2-element vector defining your frequency limits of 
%                    interest. By default the whole frequency range of the 
%                    input is taken.
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

if ~isfield(cfg, 'channel'),     cfg.channel     = {'all'};                       end
if ~isfield(cfg, 'channelcmb'),  cfg.channelcmb  = {};                            end
if ~isfield(cfg, 'foilim'),      cfg.foilim      = [freq.freq(1) freq.freq(end)]; end
if ~isfield(cfg, 'keepfourier'), cfg.keepfourier = 'no';                          end
if ~isfield(cfg, 'feedback'),    cfg.feedback    = 'text';                        end

%select the channels on which the power-spectra will be computed
chn = ft_channelselection(cfg.channel,freq.label);
for j = 1:length(chn)
  chnindx(j,1) = find(strcmp(chn(j), freq.label));
  %chnindx(j,1) = find(strcmp(chn{j}, freq.label));
end

%convert the channelcombinations to indices
chncmb  = ft_channelcombination(cfg.channelcmb, freq.label);
cmbindx = zeros(size(chncmb,1),2);
for j = 1:size(chncmb,1)
  cmbindx(j,1) = find(strcmp(chncmb(j,1), freq.label));
  cmbindx(j,2) = find(strcmp(chncmb(j,2), freq.label));
  %cmbindx(j,1) = find(strcmp(chncmb{j,1}, freq.label));
  %cmbindx(j,2) = find(strcmp(chncmb{j,2}, freq.label));
end

%dimensionality of the input data
Nrpt   = length(freq.cumtapcnt);
Nfrq   = size(freq.fourierspctrm,3);
Ntim   = size(freq.fourierspctrm,4);
Nchn   = length(chnindx);
Ncmb   = size(cmbindx,1);

%%FIXME
%if Ntim>1, ft_error('correct handling of time-frequency data is not yet implemented, no information about tapers is available'); end

%keeping track of the tapers
%in the case of tfr fourier-data cumtapcnt is highly redundant; for each frequency
%the number of tapers is equal, as well as for each trial, thus it is sufficient to
%reduce the original cumtapcnt to a vector of Ntrlx1 containing the number of tapers
cumtapcnt = freq.cumtapcnt;
if ~isempty(strfind(freq.dimord, 'time')),
  %cumtapcnt is NtrlxNfrqxNtim; create just one column-vector
  cumtapcnt = ones(size(cumtapcnt,1),1).*unique(cumtapcnt(~isnan(cumtapcnt(:))));
end
sumtapcnt = cumsum([0; cumtapcnt(:)]);

powspctrm = zeros(Nrpt, Nchn, Nfrq, Ntim);
progress('init', cfg.feedback, 'computing single-trial power-spectral densities');
for j = 1:Nrpt
  progress(j/Nrpt, 'trial %d/%d\n', j, Nrpt);
  tmp1     = freq.fourierspctrm([1+sumtapcnt(j):sumtapcnt(j+1)], chnindx, :, :);
  tmp2     = freq.fourierspctrm([1+sumtapcnt(j):sumtapcnt(j+1)], chnindx, :, :);
  powspctrm(j, :, :, :) = squeeze(sum( tmp1 .* conj(tmp2), 1) ./ size(tmp1,1));
end
progress('close');

crsspctrm = complex(zeros(Nrpt, Ncmb, Nfrq, Ntim), zeros(Nrpt, Ncmb, Nfrq, Ntim));
progress('init', cfg.feedback, 'computing single-trial cross-spectral densities');
for j = 1:Nrpt
  progress(j/Nrpt, 'trial %d/%d\n', j, Nrpt);
  tmp1     = freq.fourierspctrm([1+sumtapcnt(j):sumtapcnt(j+1)], cmbindx(:,1), :, :);
  tmp2     = freq.fourierspctrm([1+sumtapcnt(j):sumtapcnt(j+1)], cmbindx(:,2), :, :);
  crsspctrm(j, :, :, :) = squeeze(sum( tmp1 .* conj(tmp2), 1) ./ size(tmp1,1));
end
progress('close');

output.dimord = freq.dimord;
output.freq   = freq.freq;
output.label  = chn;
output.labelcmb(:,1) = freq.label(cmbindx(:,1));
output.labelcmb(:,2) = freq.label(cmbindx(:,2));
output.cumtapcnt = freq.cumtapcnt;
try, output.grad = freq.grad; end
try, output.time = freq.time; end
output.powspctrm = powspctrm;
output.crsspctrm = crsspctrm;
if strcmp(cfg.keepfourier, 'yes'), output.fourierspctrm = freq.fourierspctrm; end 

if isempty(output.crsspctrm), output = rmfield(output, 'crsspctrm'); end
if isempty(output.labelcmb ), output = rmfield(output, 'labelcmb' ); end

% add information about the version of this function to the configuration
cfg.version.name = mfilename('fullpath');
cfg.version.id = '$Id$';

% remember the configuration details of the input data
try, cfg.previous = freq.cfg; end

% remember the exact configuration details in the output 
output.cfg = cfg;

