function [data] = ft_combineplanar(cfg, data)

% FT_COMBINEPLANAR computes the planar gradient magnitude over both directions
% combining the two gradients at each sensor to a single positive-valued number. This
% can be done for single-trial/averaged planar gradient ERFs or single-trial/averaged
% TFRs.
%
% Use as
%   [data] = ft_combineplanar(cfg, data)
% where data contains an averaged planar gradient ERF or single-trial/averaged TFR.
%
% The configuration can contain
%   cfg.method         = 'sum', 'svd', 'abssvd', or 'complex' (default = 'sum')
%   cfg.updatesens     = 'no' or 'yes' (default = 'yes')
% and for timelocked input data (i.e. ERFs), the configuration can also contain
%   cfg.demean         = 'yes' or 'no' (default = 'no')
%   cfg.baselinewindow = [begin end]
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_MEGPLANAR

% Undocumented local options:
% cfg.foilim
% cfg.trials

% Copyright (C) 2004, Ole Jensen
% Copyright (C) 2004-2013, Robert Oostenveld
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar data
ft_preamble provenance data
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', {'raw', 'freq', 'timelock'}, 'feedback', 'yes', 'senstype', {'ctf151_planar', 'ctf275_planar', 'neuromag122', 'neuromag306', 'bti248_planar', 'bti148_planar', 'itab153_planar', 'yokogawa160_planar', 'yokogawa64_planar', 'yokogawa440_planar'});

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden',   {'combinegrad'});
cfg = ft_checkconfig(cfg, 'deprecated',  {'baseline'});
cfg = ft_checkconfig(cfg, 'renamed',     {'blc', 'demean'});
cfg = ft_checkconfig(cfg, 'renamed',     {'blcwindow', 'baselinewindow'});
cfg = ft_checkconfig(cfg, 'renamed',     {'combinemethod', 'method'});

% set the defaults
cfg.demean         = ft_getopt(cfg, 'demean',         'no');
cfg.foilim         = ft_getopt(cfg, 'foilim',         [-inf inf]);
cfg.baselinewindow = ft_getopt(cfg, 'baselinewindow', [-inf inf]);
cfg.trials         = ft_getopt(cfg, 'trials',         'all', 1);
cfg.feedback       = ft_getopt(cfg, 'feedback',       'none');
cfg.method         = ft_getopt(cfg, 'method',         'sum');
cfg.updatesens     = ft_getopt(cfg, 'updatesens',     'yes');

if isfield(cfg, 'baseline')
  ft_warning('only supporting cfg.baseline for backwards compatibility, please update your cfg');
  cfg.demean         = 'yes';
  cfg.baselinewindow = cfg.baseline;
end

israw      = ft_datatype(data, 'raw');
istimelock = ft_datatype(data, 'timelock');
isfreq     = ft_datatype(data, 'freq');
if isfield(data, 'dimord')
  dimord = data.dimord;
end

% select trials of interest
if ~strcmp(cfg.trials, 'all')
  ft_error('trial selection has not been implemented yet') % first fix ft_checkdata (see above)
end

% find the combination of horizontal and vertical channels that should be combined
planar = ft_senslabel(ft_senstype(data), 'output', 'planarcombined');
[sel_pH, sel_dH]    = match_str(planar(:,1), data.label);  % indices of the horizontal channels
[sel_pV, sel_dV]    = match_str(planar(:,2), data.label);  % indices of the vertical   channels

% identify and remove unnpaired channels
[dum,iH,iV]     = intersect(sel_pH,sel_pV);
sel_dH=sel_dH(iH);
sel_dV=sel_dV(iV);

% find the other channels that are present in the data
sel_other = setdiff(1:length(data.label), [sel_dH(:)' sel_dV(:)']);
lab_other = data.label(sel_other);

% define the channel names after combining the planar combinations
% they should be sorted according to the order of the planar channels in the data
[dum, sel_planar] = match_str(data.label(sel_dH),planar(:,1));
lab_comb          = planar(sel_planar,end);

% perform baseline correction
if strcmp(cfg.demean, 'yes')
  if ~(istimelock || israw)
    ft_error('baseline correction is only supported for timelocked or raw input data')
  end
  if ischar(cfg.baselinewindow) && strcmp(cfg.baselinewindow, 'all')
    cfg.baselinewindow = [-inf inf];
  end
  % find the timebins corresponding to the baseline interval
  tbeg = nearest(data.time, cfg.baselinewindow(1));
  tend = nearest(data.time, cfg.baselinewindow(2));
  cfg.baselinewindow(1) = data.time(tbeg);
  cfg.baselinewindow(2) = data.time(tend);
  data.avg = ft_preproc_baselinecorrect(data.avg, tbeg, tend);
end

if isfreq
  
  switch cfg.method
    case 'sum'
      if isfield(data, 'powspctrm')
        % compute the power of each planar channel, by summing the horizontal and vertical gradients
        dimtok = tokenize(dimord, '_');
        catdim = strmatch('chan',dimtok);
        if catdim==1
          combined = data.powspctrm(sel_dH,:,:,:) + data.powspctrm(sel_dV,:,:,:);
          other    = data.powspctrm(sel_other,:,:,:);
        elseif catdim==2
          combined = data.powspctrm(:,sel_dH,:,:,:) + data.powspctrm(:,sel_dV,:,:,:);
          other    = data.powspctrm(:,sel_other,:,:,:);
        else
          ft_error('unsupported dimension order of frequency data');
        end
        data.powspctrm = cat(catdim, combined, other);
        data.label     = cat(1, lab_comb(:), lab_other(:));
      else
        ft_error('cfg.method = ''%s'' only works for frequency data with powspctrm', cfg.method);
      end
    case 'svd'
      if isfield(data, 'fourierspctrm')
        fbin = nearest(data.freq, cfg.foilim(1)):nearest(data.freq, cfg.foilim(2));
        Nrpt   = size(data.fourierspctrm,1);
        Nsgn   = length(sel_dH);
        Nfrq   = length(fbin);
        Ntim   = size(data.fourierspctrm,4);
        %fourier= complex(zeros(Nrpt,Nsgn,Nfrq,Ntim),zeros(Nrpt,Nsgn,Nfrq,Ntim));
        fourier= nan(Nrpt,Nsgn,Nfrq,Ntim);
        ft_progress('init', cfg.feedback, 'computing the svd');
        for j = 1:Nsgn
          ft_progress(j/Nsgn, 'computing the svd of signal %d/%d\n', j, Nsgn);
          for k = 1:Nfrq
            dum = reshape(data.fourierspctrm(:,[sel_dH(j) sel_dV(j)],fbin(k),:), [Nrpt 2 Ntim]);
            dum = permute(dum, [2 3 1]);
            dum = reshape(dum, [2 Ntim*Nrpt]);
            timbin = ~isnan(dum(1,:));
            [loading, ~,  ori, sin_val] = svdfft(dum(:,timbin),2,data.cumtapcnt);
            dum2   = loading(1,:);
            dum(1,timbin) = dum2;
            dum = reshape(dum(1,:),[Ntim Nrpt]);
            fourier(:,j,k,:) = transpose(dum);
            data.ori{k} = ori; % to change into a cell
            data.eta{k} = sin_val(1)/sum(sin_val(2:end)); % to change into a cell
            
            %for m = 1:Ntim
            %  dum                     = data.fourierspctrm(:,[sel_dH(j) sel_dV(j)],fbin(k),m);
            %  timbin                  = find(~isnan(dum(:,1)));
            %  [fourier(timbin,j,k,m)] = svdfft(transpose(dum(timbin,:)),1);
            %end
          end
        end
        ft_progress('close');
        other              = data.fourierspctrm(:,sel_other,fbin,:);
        data               = rmfield(data, 'fourierspctrm');
        data.fourierspctrm = cat(2, fourier, other);
        data.label         = cat(1, lab_comb(:), lab_other(:));
        data.freq          = data.freq(fbin);
      else
        ft_error('cfg.method = ''%s'' only works for frequency data with fourierspctrm', cfg.method);
      end
    otherwise
      ft_error('cfg.method = ''%s'' is not supported for frequency data', cfg.method);
  end % switch method
  
elseif (israw || istimelock)
  if istimelock
    % convert timelock to raw
    data = ft_checkdata(data, 'datatype', 'raw', 'feedback', 'yes');
  end
  
  switch cfg.method
    case 'sum'
      Nrpt = length(data.trial);
      for k = 1:Nrpt
        combined = sqrt(data.trial{k}(sel_dH,:).^2 + data.trial{k}(sel_dV,:).^2);
        other    = data.trial{k}(sel_other,:);
        data.trial{k} = [combined; other];
      end
      data.label = cat(1, lab_comb(:), lab_other(:));
      
    case 'complex'
      Nrpt = length(data.trial);
      for k = 1:Nrpt
        combined = data.trial{k}(sel_dH,:)*1i + data.trial{k}(sel_dV,:);
        other    = data.trial{k}(sel_other,:);
        data.trial{k} = [combined; other];
      end
      data.label = cat(1, lab_comb(:), lab_other(:));
      
    case {'svd' 'abssvd'}
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
        if strcmp(cfg.method, 'abssvd')||strcmp(cfg.method, 'svd')
          [loading, ~,  ori, sin_val] = svdfft(tmpdat,2);
          data.ori{k} = ori; % to change into a cell
          data.eta{k} = sin_val(1)/sum(sin_val(2:end)); % to change into a cell
          if strcmp(cfg.method, 'abssvd')
            tmpdat2 = abs(loading(1,:));
          else
            tmpdat2 = loading(1,:);
          end
        end
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
      ft_error('cfg.method = ''%s'' is not supported for timelocked or raw data', cfg.method);
  end % switch method
  
  if istimelock
    % convert raw to timelock
    data = ft_checkdata(data, 'datatype', 'timelock', 'feedback', 'yes');
  end
  
else
  ft_error('unsupported input data');
end % which ft_datatype

% remove the fields for which the planar gradient could not be combined
data = removefields(data, {'crsspctrm', 'labelcmb'});

if strcmp(cfg.updatesens, 'yes') && isfield(data, 'grad')
  % update the grad and only retain the channel related info
  [sel_dH, sel_comb] = match_str(data.grad.label, planar(:,1));  % indices of the horizontal channels
  [sel_dV          ] = match_str(data.grad.label, planar(:,2));  % indices of the vertical   channels
  
  % find the other channels that are present in the data
  sel_other = setdiff(1:length(data.grad.label), [sel_dH(:)' sel_dV(:)']);
  lab_other = data.grad.label(sel_other);
  lab_comb  = planar(sel_comb,end);
  
  % compute the average position
  newpos   = [
    (data.grad.chanpos(sel_dH,:)+data.grad.chanpos(sel_dV,:))/2
    data.grad.chanpos(sel_other,:)
    ];
  % compute the average orientation
  newori   = [
    (data.grad.chanori(sel_dH,:)+data.grad.chanori(sel_dV,:))/2
    data.grad.chanori(sel_other,:)
    ];
  newlabel = [
    lab_comb
    lab_other
    ];
  newtype = [
    repmat({'unknown'}, numel(sel_comb), 1) % combined planar
    data.grad.chantype(sel_other(:))        % keep the known channel details
    ];
  newunit = [
    repmat({'unknown'}, numel(sel_comb), 1) % combined planar
    data.grad.chanunit(sel_other(:))        % keep the known channel details
    ];
  
  newgrad.chanpos  = newpos;
  newgrad.chanori  = newori;
  newgrad.label    = newlabel;
  newgrad.chantype = newtype;
  newgrad.chanunit = newunit;
  newgrad.unit     = data.grad.unit;
  newgrad.type     = [data.grad.type '_combined'];
  
  % remember the original channel position details
  if isfield(data.grad, 'chanposold')
    newgrad = copyfields(data.grad, newgrad, {'chanposold', 'chanoriold', 'labelold', 'chantypeold', 'chanunitold'});
  else
    newgrad.labelold     = data.grad.label;
    newgrad.chanposold   = data.grad.chanpos;
    newgrad.chanoriold   = data.grad.chanori;
    newgrad.chantypeold  = data.grad.chantype;
    newgrad.chanunitold  = data.grad.chanunit;
  end
  
  % replace it with the updated gradiometer description
  data.grad = newgrad;
end

% convert back to input type if necessary
if istimelock
  data = ft_checkdata(data, 'datatype', 'timelock');
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   data
ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data
