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

%$Log: fourier2crsspctrm.m,v $
%Revision 1.14  2008/11/27 09:04:48  kaigoe
%added default cfg.feedback=text
%
%Revision 1.13  2007/08/31 07:12:57  jansch
%added feedback to be specified by cfg.feedback instead of hardcoded textbar
%
%Revision 1.12  2007/07/31 08:28:36  jansch
%some cosmetic changes
%
%Revision 1.11  2006/06/23 10:49:33  jansch
%added cfg.keepfourier as an option
%
%Revision 1.10  2006/03/22 13:59:18  jansch
%implemented support for tfr-structures containing fourierspectra
%
%Revision 1.9  2006/03/20 11:25:12  jansch
%made changes, to be called from the newly implemented freqdescriptives
%
%Revision 1.8  2006/02/28 12:43:55  roboos
%changed a foi into freq
%
%Revision 1.7  2006/02/23 10:28:17  roboos
%changed dimord strings for consistency, changed toi and foi into time and freq, added fixdimord where neccessary
%
%Revision 1.6  2006/02/01 12:26:04  roboos
%made all uses of dimord consistent with the common definition of data dimensions, see the fixdimord() function
%
%Revision 1.5  2005/08/18 12:17:47  jansch
%added a try-catch in the assignment of the gradiometer-description to the
%output.
%
%Revision 1.4  2005/08/16 07:49:55  jansch
%included the powerspectra in the crsspctrm, removed powspctrm from output
%
%Revision 1.3  2005/08/15 14:46:57  jansch
%fixed small bug in creation of output-structure
%
%Revision 1.2  2005/08/15 10:36:08  jansch
%added version information to the configuration
%
%Revision 1.1  2005/08/15 10:33:19  jansch
%new implementation, to be used as a subfunction for freqdescriptives and
%other stuff
%

if ~isfield(cfg, 'channel'),     cfg.channel     = {'all'};                       end
if ~isfield(cfg, 'channelcmb'),  cfg.channelcmb  = {};                            end
if ~isfield(cfg, 'foilim'),      cfg.foilim      = [freq.freq(1) freq.freq(end)]; end
if ~isfield(cfg, 'keepfourier'), cfg.keepfourier = 'no';                          end
if ~isfield(cfg, 'feedback'),    cfg.feedback    = 'text';                        end

%select the channels on which the power-spectra will be computed
chn     = channelselection(cfg.channel,freq.label);
for j = 1:length(chn)
  chnindx(j,1) = find(strcmp(chn(j), freq.label));
  %chnindx(j,1) = strmatch(chn{j}, freq.label, 'exact');
end

%convert the channelcombinations to indices
chncmb  = channelcombination(cfg.channelcmb, freq.label);
cmbindx = zeros(size(chncmb,1),2);
for j = 1:size(chncmb,1)
  cmbindx(j,1) = find(strcmp(chncmb(j,1), freq.label));
  cmbindx(j,2) = find(strcmp(chncmb(j,2), freq.label));
  %cmbindx(j,1) = strmatch(chncmb{j,1}, freq.label, 'exact');
  %cmbindx(j,2) = strmatch(chncmb{j,2}, freq.label, 'exact');
end

%dimensionality of the input data
Nrpt   = length(freq.cumtapcnt);
Nfrq   = size(freq.fourierspctrm,3);
Ntim   = size(freq.fourierspctrm,4);
Nchn   = length(chnindx);
Ncmb   = size(cmbindx,1);

%%FIXME
%if Ntim>1, error('correct handling of time-frequency data is not yet implemented, no information about tapers is available'); end

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
try, output.grad = freq.grad; end;
try, output.time = freq.time; end;
output.powspctrm = powspctrm;
output.crsspctrm = crsspctrm;
if strcmp(cfg.keepfourier, 'yes'), output.fourierspctrm = freq.fourierspctrm; end 

if isempty(output.crsspctrm), output = rmfield(output, 'crsspctrm'); end;
if isempty(output.labelcmb ), output = rmfield(output, 'labelcmb' ); end;

% add information about the version of this function to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i1] = dbstack;
  cfg.version.name = st(i1);
end
cfg.version.id = '$Id: fourier2crsspctrm.m,v 1.14 2008/11/27 09:04:48 kaigoe Exp $';
% remember the configuration details of the input data
try, cfg.previous = freq.cfg; end
% remember the exact configuration details in the output 
output.cfg = cfg;

