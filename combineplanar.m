function [data] = combineplanar(cfg, data)

% COMBINEPLANAR computes the planar gradient magnitude over both directions
% combining the two gradients at each sensor to a single positive-valued number.
% This can be done for averaged ERFs or TFRs (i.e. powerspectra). 
% 
% Use as
%   [data] = combineplanar(cfg, data)
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
% See also MEGPLANAR

% Undocumented local options:
% cfg.baseline
% cfg.combinemethod
% cfg.foilim
% cfg.trials

% Copyright (C) 2004, Ole Jensen, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

fieldtripdefs

cfg = checkconfig(cfg, 'trackconfig', 'on');

% check if the input data is valid for this function
data = checkdata(data, 'datatype', {'raw', 'freq', 'timelock'}, 'feedback', 'yes', 'senstype', {'ctf151_planar', 'ctf275_planar', 'neuromag122', 'neuromag306', 'bti248_planar', 'bti148_planar'});

israw      = datatype(data, 'raw');
isfreq     = datatype(data, 'freq');
istimelock = datatype(data, 'timelock');
try, dimord = data.dimord; end

% set the defaults
if ~isfield(cfg, 'blc'),           cfg.blc           = 'no';  end
if ~isfield(cfg, 'blcwindow'),     cfg.blcwindow     = 'all'; end
if ~isfield(cfg, 'combinegrad'),   cfg.combinegrad   = 'no';  end
if ~isfield(cfg, 'combinemethod'), cfg.combinemethod = 'sum'; end 
if ~isfield(cfg, 'foilim'),        cfg.foilim        = [];    end
if ~isfield(cfg, 'trials'),        cfg.trials        = 'all'; end
if ~isfield(cfg, 'feedback'),      cfg.feedback      = 'none'; end
if isfield(cfg, 'baseline')
  warning('only supporting cfg.baseline for backwards compatibility, please update your cfg');
  cfg.blc = 'yes';
  cfg.blcwindow = cfg.baseline;
end
if strcmp(cfg.blc, 'yes') && isempty(cfg.blcwindow)
  cfg.blcwindow = [min(data.time) max(data.time)];
end

% select trials of interest
if ~strcmp(cfg.trials, 'all')
  error('trial selection has not been implemented yet') % first fix checkdata (see above)
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
lab_comb  = planar(sel_planar,3);

% perform baseline correction
if strcmp(cfg.blc, 'yes') 
  if ~istimelock,
    error('baseline correction is only supported for ERFs')
  else
    if ischar(cfg.blcwindow) && strcmp(cfg.blcwindow, 'all')
      cfg.blcwindow = [min(data.time) max(data.time)];
    end
    % find the timebins corresponding to the baseline interval
    tbeg = nearest(data.time, cfg.blcwindow(1));
    tend = nearest(data.time, cfg.blcwindow(2));
    cfg.blcwindow(1) = data.time(tbeg);
    cfg.blcwindow(2) = data.time(tend);
    data.avg = blc(data.avg, tbeg, tend); 
  end
end

if isfreq
  switch cfg.combinemethod
    case 'sum'
      if isfield(data, 'powspctrm'),
        % compute the power of each planar channel, by summing the horizontal and vertical gradients
    dimtok = tokenize(dimord,'_');
    catdim = strmatch('chan',dimtok);
    if catdim==1,
      tmp1 = data.powspctrm(sel_dH,:,:,:) + data.powspctrm(sel_dV,:,:,:);
      tmp2 = data.powspctrm(sel_other,:,:,:);
    elseif catdim==2,
      tmp1 = data.powspctrm(:,sel_dH,:,:,:) + data.powspctrm(:,sel_dV,:,:,:);
      tmp2 = data.powspctrm(:,sel_other,:,:,:);
        else
          error('unsupported dimension order of frequency data');
        end
        data.powspctrm = cat(catdim, tmp1, tmp2);
      else
        error('cfg.combinemethod = ''sum'' only works for frequency data with powspctrm');
      end
    case 'svd'
      if isfield(data, 'fourierspctrm'), 
        if isempty(cfg.foilim), cfg.foilim = [data.freq(1) data.freq(end)]; end;
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
        data.freq          = data.freq(fbin);
      else
        error('cfg.combinemethod = ''svd'' only works for frequency data with fourierspctrm');
      end
    otherwise
  end
else
  if istimelock,
    data = checkdata(data, 'datatype', 'raw', 'feedback', 'yes');
  end
  
  switch cfg.combinemethod
    case 'sum'
      Nrpt = length(data.trial);
      for k = 1:Nrpt
        tmp1 = sqrt(data.trial{k}(sel_dH,:).^2 + data.trial{k}(sel_dV,:).^2);
    tmp2 = data.trial{k}(sel_other,:);
    data.trial{k} = [tmp1;tmp2];
      end
    case 'svd'
      Nrpt = length(data.trial);
      Nsgn = length(sel_dH);
      Nsmp = cellfun('size', data.trial, 2);
      Csmp = cumsum([0 Nsmp]);
      %do a 'fixed orientation' across all trials approach here
      %this is different from the frequency case FIXME
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

    otherwise
  end

  if istimelock,
    data = checkdata(data, 'datatype', 'timelock', 'feedback', 'yes');
  end
end

if strcmp(cfg.combinegrad, 'no') && ~isfield(data, 'grad')
  % the planar gradiometer definition was already removed
  % nothing needs to be done here
elseif strcmp(cfg.combinegrad, 'no') && isfield(data, 'grad')
  % remove the planar gradiometer definition since it does not match the data any more
  data = rmfield(data, 'grad');
elseif strcmp(cfg.combinegrad, 'yes') && ~isfield(data, 'grad')
  % there is no gradiometer definition, impossible to reconstruct it
  error('the planar gradiometer definition is missing, cannot convert it back to axial');
elseif strcmp(cfg.combinegrad, 'yes') && isfield(data, 'grad')
  warning('trying to convert planar to axial gradiometers, this is experimental');
  % try to reconstruct the original axial gradiometer array from the planar gradiometer definition
  orig = data.grad;
  if all(size(orig.pnt)==[302 3]) && ...
    all(size(orig.pnt)==[302 3]) && ...
    all(size(orig.tra)==[302 302]) && ...
    length(orig.label)==302 && ...
    all(sum(orig.tra~=0,1)>2)
    % This looks as if it was made using the MEGPLANAR nearest neighbour approach
    % which means that the coil position and orientation still correspond
    % with those of the original axial gradiometer. Only the label and tra
    % have been modified and have to be restored to their original values.
    axial.pnt = orig.pnt;
    axial.ori = orig.ori;
    for i=1:151
      axial.label{i} = orig.label{i}(1:(end-3));
    end
    if all(orig.ori(1,:)==orig.ori(152,:))
      % orientation is the same, the subtraction should be in "tra"
      axial.tra = [eye(151) -eye(151)];
    else
      % orientation is opposite, the subtraction should not be in "tra"
      axial.tra = [eye(151) eye(151)];
    end
    try
      axial.unit = orig.unit;
    end
  else
    error('cannot convert gradiometer definition back to axial, please contact Robert');
  end
  data.grad = axial;
end

% reconstruct the original channel labels
data.label = [lab_comb(:); lab_other(:)];

% remove the fields for which the planar gradient could not be combined
try, data = rmfield(data, 'crsspctrm');   end
try, data = rmfield(data, 'labelcmb');      end

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes'); 

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

