function [data] = ft_combineplanar(cfg, data)

% FT_COMBINEPLANAR computes the planar gradient magnitude over both directions
% combining the two gradients at each sensor to a single positive-valued number.
% This can be done for averaged ERFs or TFRs (i.e. powerspectra).
%
% Use as
%   [data] = ft_combineplanar(cfg, data)
% where data contains an averaged planar gradient (either ERF or TFR).
%
% In the case of ERFs, the configuration can contain
%   cfg.blc       = 'yes' or 'no' (default)
%   cfg.blcwindow = [begin end]
%
% After combining the planar data, the planar gradiometer definition does not
% match the data any more and therefore it is removed from the data. With
%   cfg.combinegrad  = 'yes'
% the function will try to reconstruct the axial gradiometer definition.
%
% See also FT_MEGPLANAR

% Undocumented local options:
% cfg.baseline
% cfg.combinemethod
% cfg.foilim
% cfg.trials
% cfg.inputfile     = one can specifiy preanalysed saved data as input
% cfg.outputfile    = one can specify output as file to save to disk

% Copyright (C) 2004, Ole Jensen, Robert Oostenveld
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

fieldtripdefs

% set the defaults
if ~isfield(cfg, 'blc'),           cfg.blc            = 'no';       end
if ~isfield(cfg, 'foilim'),        cfg.foilim         = [-inf inf]; end
if ~isfield(cfg, 'blcwindow'),     cfg.blcwindow      = [-inf inf]; end
if ~isfield(cfg, 'combinemethod'), cfg.combinemethod  = 'sum';      end
if ~isfield(cfg, 'trials'),        cfg.trials         = 'all';      end
if ~isfield(cfg, 'feedback'),      cfg.feedback       = 'none';     end
if ~isfield(cfg, 'inputfile'),     cfg.inputfile      = [];         end
if ~isfield(cfg, 'outputfile'),    cfg.outputfile     = [];         end

hasdata = (nargin>1);
if ~isempty(cfg.inputfile)
  % the input data should be read from file
  if hasdata
    error('cfg.inputfile should not be used in conjunction with giving input data to this function');
  else
    data = loadvar(cfg.inputfile, 'data');
  end
end

% check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', {'raw', 'freq', 'timelock'}, 'feedback', 'yes', 'senstype', {'ctf151_planar', 'ctf275_planar', 'neuromag122', 'neuromag306', 'bti248_planar', 'bti148_planar', 'itab153_planar'});

cfg = ft_checkconfig(cfg, 'trackconfig', 'on');
cfg = ft_checkconfig(cfg, 'forbidden', {'combinegrad'});
cfg = ft_checkconfig(cfg, 'deprecated', {'baseline'});

if isfield(cfg, 'baseline')
  warning('only supporting cfg.baseline for backwards compatibility, please update your cfg');
  cfg.blc       = 'yes';
  cfg.blcwindow = cfg.baseline;
end

israw      = datatype(data, 'raw');
isfreq     = datatype(data, 'freq');
istimelock = datatype(data, 'timelock');
try, dimord = data.dimord; end

% select trials of interest
if ~strcmp(cfg.trials, 'all')
  error('trial selection has not been implemented yet') % first fix ft_checkdata (see above)
end

% find the combination of horizontal and vertical channels that should be combined
planar    = planarchannelset(data);
sel_dH    = match_str(data.label, planar(:,1));  % indices of the horizontal channels
sel_dV    = match_str(data.label, planar(:,2));  % indices of the vertical   channels
lab_dH    = data.label(sel_dH);
lab_dV    = data.label(sel_dV);

if length(sel_dH)~=length(sel_dV)
  error('not all planar channel combinations are complete')
end

% find the other channels that are present in the data
sel_other = setdiff(1:length(data.label), [sel_dH(:)' sel_dV(:)']);
lab_other = data.label(sel_other);

% define the channel names after combining the planar combinations
% they should be sorted according to the order of the planar channels in the data
[dum, sel_planar] = match_str(data.label, planar(:,1));
lab_comb          = planar(sel_planar,3);

% perform baseline correction
if strcmp(cfg.blc, 'yes')
  if ~(istimelock || israw)
    error('baseline correction is only supported for timelocked or raw input data')
  end
  if ischar(cfg.blcwindow) && strcmp(cfg.blcwindow, 'all')
    cfg.blcwindow = [-inf inf];
  end
  % find the timebins corresponding to the baseline interval
  tbeg = nearest(data.time, cfg.blcwindow(1));
  tend = nearest(data.time, cfg.blcwindow(2));
  cfg.blcwindow(1) = data.time(tbeg);
  cfg.blcwindow(2) = data.time(tend);
  data.avg = blc(data.avg, tbeg, tend);
end

if isfreq
  
  switch cfg.combinemethod
    case 'sum'
      if isfield(data, 'powspctrm'),
        % compute the power of each planar channel, by summing the horizontal and vertical gradients
        dimtok = tokenize(dimord,'_');
        catdim = strmatch('chan',dimtok);
        if catdim==1,
          combined = data.powspctrm(sel_dH,:,:,:) + data.powspctrm(sel_dV,:,:,:);
          other    = data.powspctrm(sel_other,:,:,:);
        elseif catdim==2,
          combined = data.powspctrm(:,sel_dH,:,:,:) + data.powspctrm(:,sel_dV,:,:,:);
          other    = data.powspctrm(:,sel_other,:,:,:);
        else
          error('unsupported dimension order of frequency data');
        end
        data.powspctrm = cat(catdim, combined, other);
        data.label     = cat(1, lab_comb(:), lab_other(:));
      else
        error('cfg.combinemethod = ''%s'' only works for frequency data with powspctrm', cfg.combinemethod);
      end
    case 'svd'
      if isfield(data, 'fourierspctrm'),
        fbin = nearest(data.freq, cfg.foilim(1)):nearest(data.freq, cfg.foilim(2));
        Nrpt   = size(data.fourierspctrm,1);
        Nsgn   = length(sel_dH);
        Nfrq   = length(fbin);
        Ntim   = size(data.fourierspctrm,4);
        %fourier= complex(zeros(Nrpt,Nsgn,Nfrq,Ntim),zeros(Nrpt,Nsgn,Nfrq,Ntim));
        fourier= zeros(Nrpt,Nsgn,Nfrq,Ntim)+nan;
        progress('init', cfg.feedback, 'computing the svd');
        for j = 1:Nsgn
          progress(j/Nsgn, 'computing the svd of signal %d/%d\n', j, Nsgn);
          for k = 1:Nfrq
            dum = reshape(data.fourierspctrm(:,[sel_dH(j) sel_dV(j)],fbin(k),:), [Nrpt 2 Ntim]);
            dum = permute(dum, [2 3 1]);
            dum = reshape(dum, [2 Ntim*Nrpt]);
            timbin = ~isnan(dum(1,:));
            dum2   = svdfft(dum(:,timbin),1,data.cumtapcnt);
            dum(1,timbin) = dum2;
            dum = reshape(dum(1,:),[Ntim Nrpt]);
            fourier(:,j,k,:) = transpose(dum);
            
            %for m = 1:Ntim
            %  dum                     = data.fourierspctrm(:,[sel_dH(j) sel_dV(j)],fbin(k),m);
            %  timbin                  = find(~isnan(dum(:,1)));
            %  [fourier(timbin,j,k,m)] = svdfft(transpose(dum(timbin,:)),1);
            %end
          end
        end
        progress('close');
        other              = data.fourierspctrm(:,sel_other,fbin,:);
        data               = rmfield(data,'fourierspctrm');
        data.fourierspctrm = cat(2, fourier, other);
        data.label         = cat(1, lab_comb(:), lab_other(:));
        data.freq          = data.freq(fbin);
      else
        error('cfg.combinemethod = ''%s'' only works for frequency data with fourierspctrm', cfg.combinemethod);
      end
    otherwise
      error('cfg.combinemethod = ''%s'' is not supported for frequency data', cfg.combinemethod);
  end
  
elseif (israw || istimelock)
  if istimelock,
    % convert timelock to raw
    data = ft_checkdata(data, 'datatype', 'raw', 'feedback', 'yes');
  end
  
  switch cfg.combinemethod
    case 'sum'
      Nrpt = length(data.trial);
      for k = 1:Nrpt
        combined = sqrt(data.trial{k}(sel_dH,:).^2 + data.trial{k}(sel_dV,:).^2);
        other    = data.trial{k}(sel_other,:);
        data.trial{k} = [combined; other];
      end
      data.label = cat(1, lab_comb(:), lab_other(:));
    case 'svd'
      Nrpt = length(data.trial);
      Nsgn = length(sel_dH);
      Nsmp = cellfun('size', data.trial, 2);
      Csmp = cumsum([0 Nsmp]);
      % do a 'fixed orientation' across all trials approach here
      % this is different from the frequency case FIXME
      tmpdat = zeros(2, sum(Nsmp));
      for k = 1:Nsgn
        for m = 1:Nrpt
          tmpdat(:, (Csmp(m)+1):Csmp(m+1)) = data.trial{m}([sel_dH(k) sel_dV(k)],:);
        end
        tmpdat2 = abs(svdfft(tmpdat,1));
        tmpdat2 = mat2cell(tmpdat2, 1, Nsmp);
        for m = 1:Nrpt
          if k==1, trial{m} = zeros(Nsgn, Nsmp(m)); end
          trial{m}(k,:) = tmpdat2{m};
        end
      end
      for m = 1:Nrpt
        other = data.trial{m}(sel_other,:);
        trial{m} = [trial{m}; other];
      end
      data.trial = trial;
      data.label = cat(1, lab_comb(:), lab_other(:));
    otherwise
      error('cfg.combinemethod = ''%s'' is not supported for timelocked or raw data', cfg.combinemethod);
  end
  
  if istimelock,
    % convert raw to timelock
    data = ft_checkdata(data, 'datatype', 'timelock', 'feedback', 'yes');
  end
  
else
  error('unsupported input data');
end % which datatype

% remove the fields for which the planar gradient could not be combined
try, data = rmfield(data, 'crsspctrm');   end
try, data = rmfield(data, 'labelcmb');    end

% accessing this field here is needed for the configuration tracking
% by accessing it once, it will not be removed from the output cfg
cfg.outputfile;

% get the output cfg
cfg = ft_checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');

% store the configuration of this function call, including that of the previous function call
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id  = '$Id$';
% remember the configuration details of the input data
try, cfg.previous = data.cfg; end

% remember the exact configuration details in the output
data.cfg = cfg;

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'data', data); % use the variable name "data" in the output file
end
