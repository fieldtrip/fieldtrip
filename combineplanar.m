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
% $Log: combineplanar.m,v $
% Revision 1.44  2009/09/30 12:48:14  jansch
% added cumtapcnt as input to svdfft
%
% Revision 1.43  2009/07/23 08:11:29  crimic
% fixed tiny bug
%
% Revision 1.42  2009/07/17 08:17:24  jansch
% rewriting of big parts of the code; incorporating checkdata etc. implementation
% of 'svd' combinemethod also for time domain data
%
% Revision 1.41  2009/01/20 13:01:31  sashae
% changed configtracking such that it is only enabled when BOTH explicitly allowed at start
% of the fieldtrip function AND requested by the user
% in all other cases configtracking is disabled
%
% Revision 1.40  2008/11/21 12:48:17  sashae
% added call to checkconfig at start and end of function
%
% Revision 1.39  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.38  2008/07/16 10:20:42  jansch
% added support for bti248_planar data
%
% Revision 1.37  2008/01/31 17:20:17  sashae
% added option for trial selection
%
% Revision 1.36  2007/05/29 14:26:44  ingnie
% changed sensortype to senstype (in calling checkdata) to avoid overlap with function name sensortype
%
% Revision 1.35  2007/05/29 12:51:31  roboos
% added new options for checkdata
%
% Revision 1.34  2007/05/07 09:38:09  chrhes
% added support for single trial time-domain data
%
% Revision 1.33  2007/05/02 15:56:41  roboos
% added some comments to the code
%
% Revision 1.32  2007/04/03 15:37:07  roboos
% renamed the checkinput function to checkdata
%
% Revision 1.31  2007/03/30 17:05:40  ingnie
% checkinput; only proceed when input data is allowed datatype
%
% Revision 1.30  2007/03/27 11:05:19  ingnie
% changed call to fixdimord in call to checkinput
%
% Revision 1.29  2006/04/20 09:58:33  roboos
% updated documentation
%
% Revision 1.28  2006/02/23 10:28:16  roboos
% changed dimord strings for consistency, changed toi and foi into time and freq, added fixdimord where neccessary
%
% Revision 1.27  2006/02/07 22:18:31  roboos
% changed the dimord chancmb (used to be sgncmb) into chan
%
% Revision 1.26  2006/02/01 12:26:00  roboos
% made all uses of dimord consistent with the common definition of data dimensions, see the fixdimord() function
%
% Revision 1.25  2006/01/30 14:12:13  jansch
% included svdfft instead of svd, for the combination of planar fourier-components
%
% Revision 1.24  2005/08/16 09:05:38  jansch
% implemented svd for fourier-input
%
% Revision 1.23  2005/06/29 12:46:29  roboos
% the following changes were done on a large number of files at the same time, not all of them may apply to this specific file
% - added try-catch around the inclusion of input data configuration
% - moved cfg.version, cfg.previous and the assignment of the output cfg to the end
% - changed help comments around the configuration handling
% - some changes in whitespace
%
% Revision 1.22  2005/06/28 15:56:34  roboos
% removed subfunction planarchannelset from the code, it is now a separate file in private/planarchannelset.m
%
% Revision 1.21  2005/06/02 16:00:20  roboos
% fixed bug: instead of combining the average many times, it should combine each trial
%
% Revision 1.20  2005/06/02 12:27:07  roboos
% changed from processing only the average ERF to processing only raw trials
% added call to data2raw and raw2data, which converts various types of averages to raw trials
% for keeptrial and keepsubject, the average and variance are recomputed after the combination of planar data for each trial/subject
%
% Revision 1.19  2005/05/17 17:50:36  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.18  2005/05/02 10:15:39  roboos
% corrected output labels for planar combinations
% fixed bug in removal of fields that are not combined (such as var and crsspctrm)
%
% Revision 1.17  2005/04/29 10:28:56  roboos
% changed selection of planar combinations into a lookup table
% added support for Neuromag 306
%
% Revision 1.16  2005/02/18 12:57:50  roboos
% added support for singletrial power spectra
%
% Revision 1.15  2004/11/12 09:21:01  roboos
% added detection for 306 channel Neuromag data, no actual implementation yet
%
% Revision 1.14  2004/06/23 20:28:54  roberto
% added support for 122 channels, indicating Neuromag data
%
% Revision 1.13  2004/06/03 15:46:27  roberto
% added warning message for planar gradiometers
%
% Revision 1.12  2004/04/28 09:33:44  roberto
% added check for missing gradiometer definition
%
% Revision 1.11  2004/04/22 08:48:17  roberto
% added cfg.combinegrad, wrapped around planar to axial conversion
%
% Revision 1.10  2004/04/20 20:09:29  roberto
% convert input (planar) gradiometer definition into axial gradiometer definition,
% some changes in splaces and tabs
%
% Revision 1.9  2004/04/13 16:31:09  roberto
% fixed bug in dbstack selection of function filename for Matlab 6.1
%
% Revision 1.8  2004/04/13 14:25:24  roberto
% wrapped code-snippet around mfilename to make it compatible with Matlab 6.1
%
% Revision 1.7  2004/04/06 20:00:40  roberto
% minor changes in documentation
%
% Revision 1.6  2004/03/29 15:13:29  roberto
% added version and history to output configuration
%
% Revision 1.5  2004/03/23 11:01:09  roberto
% implemented support to keep non-MEG channels in the data after combining the planar gradients
%
% Revision 1.4  2004/02/10 15:45:04  roberto
% fixed bug in cfg.blcwindow and changed its default into 'all'
%
% Revision 1.3  2004/01/27 17:02:49  roberto
% added support for frequency and time-frequency data (only singletrial)
% removed custom baselinecorrection and replaced with blc function
%
% Revision 1.2  2004/01/27 09:49:26  roberto
% fixed small bug in length(data.label)~=302
%
% Revision 1.1  2004/01/27 09:42:06  roberto
% took over this function from Ole, improved compatibility with framework,
% added help and built in some checks for invalid input data
%

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
cfg.version.id  = '$Id: combineplanar.m,v 1.44 2009/09/30 12:48:14 jansch Exp $';
% remember the configuration details of the input data
try, cfg.previous = data.cfg; end
% remember the exact configuration details in the output 
data.cfg = cfg;

