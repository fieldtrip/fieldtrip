function data = ft_removetemplateartifact(cfg, data, template)

% FT_REMOVETEMPLATEARTIFACT removes an artifact from preprocessed data by template
% subtraction. The template can for example be formed by averaging an ECG-triggered
% MEG timecourse.
%
% Use as
%   dataclean = ft_removetemplateartifact(cfg, data, template)
% where data is raw data as obtained from FT_PREPROCESSING and template is a averaged
% timelock structure as obtained from FT_TIMELOCKANALYSIS. The configuration should
% be according to
%
%   cfg.channel  = Nx1 cell-array with selection of channels (default = 'all'), see FT_CHANNELSELECTION for details
%   cfg.artifact = Mx2 matrix with sample numbers of the artifact segments, e.g. obtained from FT_ARTIFACT_EOG
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_ARTIFACT_ECG, FT_PREPROCESSING, FT_TIMELOCKANALYSIS, FT_REJECTCOMPONENT

% Copyright (C) 2014, Robert Oostenveld
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the initial part deals with parsing the input options and data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar data template
ft_preamble provenance data template
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  % do not continue function execution in case the outputfile is present and the user indicated to keep it
  return
end

% ensure that the input data is valid for this function, this will also do
% backward-compatibility conversions of old data that for example was
% read from an old *.mat file
data     = ft_checkdata(data, 'datatype', 'raw', 'feedback', 'yes', 'hassampleinfo', 'yes');
template = ft_checkdata(template, 'datatype', 'timelock');

% get the options
cfg.channel = ft_getopt(cfg, 'method', data.label);
cfg.feedback = ft_getopt(cfg, 'method', 'text');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the actual computation is done in the middle part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ensure that the same channels are in the data and template
cfg.channel = ft_channelselection(cfg.channel, data.label);
cfg.channel = ft_channelselection(cfg.channel, template.label);

tmpcfg   = keepfields(cfg, {'channel'});
data     = ft_selectdata(tmpcfg, data);
template = ft_selectdata(tmpcfg, template);
% restore the provenance information
[cfg, data]     = rollback_provenance(cfg, data);
[cfg, template] = rollback_provenance(cfg, template);

ntrial    = length(data.trial);
nchan     = length(cfg.channel);
nartifact = size(cfg.artifact,1);

ft_progress('init', cfg.feedback, 'removing artifacts');

for i=1:ntrial
  datbegsample = data.sampleinfo(i,1);
  datendsample = data.sampleinfo(i,2);
  dattrllength = datendsample-datbegsample+1;

  model = zeros(nchan, dattrllength);
  count = 0;
  for j=1:nartifact
    artbegsample = cfg.artifact(j,1);
    artendsample = cfg.artifact(j,2);
    if artendsample>=datbegsample && artbegsample<=datendsample
      % one of the artifacts overlaps with this trial
      count = count + 1;

      % express the artifact relative to the trial
      artbegsample = artbegsample-datbegsample+1;
      artendsample = artendsample-datbegsample+1;
      % the artifact might partially fall outside the trial, in that case it needs to be trimmed
      artbegtrim = 0;
      artendtrim = 0;

      if artbegsample<1
        artbegtrim   = 1 - artbegsample;
        artbegsample = 1;
      end
      if artendsample>dattrllength
        artendtrim   = artendsample - dattrllength;
        artendsample = dattrllength;
      end

      % insert the trimmed template in the model for this trial
      model(:, artbegsample:artendsample) = template.avg(:, (1+artbegtrim):(end-artendtrim));

    end
  end % for each artifact
  ft_progress(i/ntrial, 'removing %d artifacts from trial %d of %d\n', count, i, ntrial);

  % remove the artifact model from the actual data
  data.trial{i} = data.trial{i} - model;

end % for each trial
ft_progress('close');

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   data template
ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data
