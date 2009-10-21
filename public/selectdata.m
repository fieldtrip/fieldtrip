function [data] = selectdata(varargin)

% this function serves to concatenate the input data-structures along the
% compatible dimensions and thus is a more general implementation of
% appenddata, which only raw data. Moreover, it can be used to equate the
% data of different conditions to match e.g. in channels time-axis etc
% Moreover, it can be used as a generalization to ...average with 'keepindividual'
%
% Finally, this function serves to subselect regions-of-interest from the input data,
% either or not averaging across the specified dimensions.
%
% Supported input data:
%   freq
%   timelock
%   source   (not yet)
%   volume   (not yet)
%
% supported options:
%   foilim
%   toilim
%   roi
%   channel (FIXME this is also done by preprocessing?)
%   avgoverchan
%   avgoverfreq
%   avgovertime
%   avgoverroi
%   avgoverrpt

% Copyright (C) 2009, Jan-Mathijs Schoffelen
%
% $Log: selectdata.m,v $
% Revision 1.13  2009/10/17 18:06:43  jansch
% built-in support to compute jackknife samples
%
% Revision 1.12  2009/10/17 17:50:25  jansch
% allowed conversion back to raw datatype
%
% Revision 1.11  2009/10/01 12:14:05  jansch
% some changes
%
% Revision 1.10  2009/08/18 09:55:46  jansch
% included possibility to concatenate over grid positions, allowing for cutting
% the dipole grid and glueing it together later on
%
% Revision 1.9  2009/08/17 08:41:19  jansch
% multiple changes
%
% Revision 1.8  2009/07/15 12:11:57  jansch
% fixed small bug
%
% Revision 1.7  2009/07/06 09:41:18  jansch
% multiple changes. allowing for selection of rpt in frequency data when input
% data has rpttap. allowing for grandaveraging functionality in the case of
% multiple inputs with the same dimensionalities. this is equivalent to the
% XXXgrandaverage functions with keepindividual = 'yes'.
%
% Revision 1.6  2009/04/14 18:29:32  roboos
% deleted the subfunction istrue, since it now is a seperate function
%
% Revision 1.5  2009/03/18 19:49:54  roboos
% use the smart xxxdim functions
%
% Revision 1.4  2009/01/26 20:53:11  roboos
% ensure that some options are true|false
% changes some whitespace
%
% Revision 1.3  2009/01/12 17:05:58  roboos
% fixed some whitespace
%

% check the input data and options
isdata  = find(cellfun(@isstruct,varargin));
keyvals = setdiff(1:length(varargin),isdata);

data   = varargin(isdata);
kvp    = varargin(keyvals);
dtype  = cell(1,length(data));
dimord = cell(1,length(data));

for k = 1:length(data)
  data{k} = checkdata(data{k}, 'datatype', {'freq' 'timelock' 'source', 'volume', 'freqmvar', 'raw'});
  [dtype{k}, dimord{k}]  = datatype(data{k});
  if strcmp(dtype{k}, 'raw'),
    %convert to timelock and keep track of this
    israw = 1;
    data{k} = checkdata(data{k}, 'datatype', 'timelock');
    [dtype{k}, dimord{k}] = datatype(data{k});
  else
    israw = 0;
  end
end

if any(~strmatch(dtype{1},dtype))
  error('different types of input data is not supported');
end

% check consistency of input data
if any(~strmatch(dimord{1},dimord))
  error('a different dimord in the input data is not supported');
end

isfreq   = datatype(data{1},'freq');
istlck   = datatype(data{1},'timelock');
issource = datatype(data{1},'source');
isvolume = datatype(data{1},'volume');
isfreqmvar = datatype(data{1},'freqmvar');

selchan  = keyval('channel', kvp); selectchan = ~isempty(selchan);
selfoi   = keyval('foilim',  kvp); selectfoi  = ~isempty(selfoi);
seltoi   = keyval('toilim',  kvp); selecttoi  = ~isempty(seltoi);
selroi   = keyval('roi',     kvp); selectroi  = ~isempty(selroi);
selrpt   = keyval('rpt',     kvp); selectrpt  = ~isempty(selrpt);
param    = keyval('param',   kvp); if isempty(param), param = 'all'; end

avgoverchan  = keyval('avgoverchan',  kvp); if isempty(avgoverchan), avgoverchan = false; end
avgoverfreq  = keyval('avgoverfreq',  kvp); if isempty(avgoverfreq), avgoverfreq = false; end
avgovertime  = keyval('avgovertime',  kvp); if isempty(avgovertime), avgovertime = false; end
avgoverroi   = keyval('avgoverroi',   kvp); if isempty(avgoverroi),  avgoverroi  = false; end
avgoverrpt   = keyval('avgoverrpt',   kvp); if isempty(avgoverrpt),  avgoverrpt  = false; end
dojack       = keyval('jackknife',    kvp); if isempty(dojack),      dojack      = false; end

% create anonymous function and apply it to the boolean input arguments
istrue = @(x)(ischar(x) && (strcmpi(x, 'yes') || strcmpi(x, 'true')) || (~isempty(x) && numel(x)==1 && x==1));

% ensure that these are boolean arguments, optionally convert from "yes"/"no" to true/false
avgoverchan = istrue(avgoverchan);
avgoverfreq = istrue(avgoverfreq);
avgovertime = istrue(avgovertime);
avgoverroi  = istrue(avgoverroi);
avgoverrpt  = istrue(avgoverrpt);
dojack      = istrue(dojack);

if dojack && avgoverrpt,
  error('it is not possible to do both a jackknife and to average across replicates');
end

if length(data)>1 && selectrpt,
  error('multiple data structures as input is not supported in combination with subselection of trials');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% concatenate the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(data)>1,
  % determine the way to concatenate

  if issource || isvolume,
    param = parameterselection(param, data{1}); % FIXME check consistency across input data of presence of specific parameters
  else
    param = {param};
  end

  dimtok                           = tokenize(dimord{1}, '_');
  dimtok(strmatch('chan', dimtok)) = {'label'}; % data.chan does not exist

  dimmat      = zeros(length(dimtok), length(data));
  dimmat(:,1) = 1;
  for k = 1:length(dimtok)
    if isempty(strfind(dimtok{k},'rpt')),
      dimdat = getfield(data{1}, dimtok{k});
    else
      % dimtok is 'rpt' or 'rpttap'
      dimdat = size(getsubfield(data{1}, param{1}),1);
    end
    for m = 2:length(data)
      if isempty(strfind(dimtok{k},'rpt')),
        dimdat2 = getfield(data{m},dimtok{k});
      else
        % dimtok is 'rpt' or 'rpttap'
        dimdat2 = size(getsubfield(data{m}, param{1}),1);
      end
      try, dimmat(k,m) = all(dimdat(:)==dimdat2(:));            catch end;
      try, dimmat(k,m) = all(cellfun(@isequal,dimdat,dimdat2)); catch end;
    end
  end
  catdim = find(sum(dimmat,2)<length(data));
  
  if length(catdim)>1,
    error('ambiguous dimensions for concatenation');
  elseif isempty(catdim) && isempty(strmatch('rpt',dimtok)) && isempty(strmatch('rpttap',dimtok)),
    %treat as individual observations: prepend a first dimension 'rpt'
    %(so this part should be able to cover the functionality of ...grandaverage)
    catdim = 0;
  elseif isempty(catdim) && (~isempty(strmatch('rpt',dimtok)) || ~isempty(strmatch('rpttap',dimtok)))
    %append observations
    catdim = 1;
  elseif isempty(catdim)
    error('don''t know how to concatenate the data');
  end

  % concatenate the data
  for k = 1:length(param)
    tmp = cell(1,length(data));
    for m = 1:length(tmp)
      tmp{m} = getsubfield(data{m},param{k});
    end
    if catdim==0,
      ndim    = length(size(tmp{1}));
      data{1} = setsubfield(data{1}, param{k}, permute(cat(ndim+1,tmp{:}),[ndim+1 1:ndim]));
    else
      data{1} = setsubfield(data{1}, param{k}, cat(catdim,tmp{:}));
    end
  end

  if catdim==0,
    %a dimension has been prepended
    dimtok    = ['rpt' dimtok];
    catdim    = 1;
    dimord{1} = ['rpt_',dimord{1}];
    if issubfield(data{1}, 'dim'),
      dim       = [length(data) data{1}.dim];
    end
  else 
    if issubfield(data{1}, 'dim'),
      dim       = data{1}.dim;
    end
  end
  
  % concatenate the relevant descriptive fields in the data-structure
  if ~strcmp(dimtok{catdim},'rpt') && ~strcmp(dimtok{catdim},'rpttap'),
    for k = 1:length(data)
      if k==1,
        tmp = getsubfield(data{k}, dimtok{catdim});
	if strcmp(dimtok{catdim},'pos') && isfield(data{k},'inside'),
	  tmpinside  = getfield(data{k}, 'inside');
	  tmpoutside = getfield(data{k}, 'outside'); 
	  tmpnvox    = numel(tmpinside)+numel(tmpoutside);
	end
      else
        if strcmp(dimtok{catdim},'pos'),
          tmp = [tmp; getsubfield(data{k}, dimtok{catdim})]; sortflag = 0;
	  
	  %FIXME make this robust, now inside as vector is assumed
	  if exist('tmpinside', 'var')
	    tmpx       = getfield(data{k}, 'inside');
	    tmpx2      = getfield(data{k}, 'outside');
	    tmpnvox    = numel(tmpinside)+numel(tmpoutside);
	    tmpinside  = [tmpinside(:)'  tmpnvox(end)+tmpx(:)'];
	    tmpoutside = [tmpoutside(:)' tmpnvox(end)+tmpx2(:)'];
	  end
	else
          tmp = [tmp  getsubfield(data{k}, dimtok{catdim})]; sortflag = 1;
        end
      end
    end
    data{1} = setsubfield(data{1}, dimtok{catdim}, tmp);
    if exist('tmpinside', 'var')
      data{1} = setfield(data{1}, 'inside',  tmpinside);
      data{1} = setfield(data{1}, 'outside', tmpoutside);
    end
    %FIXME think about this
    tryfields = {'dof'};
  else
    % no such field as {'label','time','freq','pos'} has to be concatenated
    sortflag  = 0;
    tryfields = {'cumsumcnt','cumtapcnt'}; 
  end
  
  % concatenate the relevant descriptive fields in the data-structure (continued)
  for k = 1:length(tryfields)
    try,
      for m = 1:length(data)
        if m==1,
          tmpfield = getfield(data{m}, tryfields{k});
        else
          tmpfield = [tmpfield; getfield(data{m}, tryfields{k})];
        end
      end
      data{1} = setfield(data{1}, tryfields{k}, tmpfield);
    catch
    end
  end

  % FIXME this is ugly: solve it
  %if issource || isvolume,
  %  data{1}.dim(catdim) = max(size(tmp));
  %end

  % sort concatenated data FIXME this is also ugly and depends on tmp
  if sortflag && ~iscell(tmp),
    [srt, ind] = sort(tmp, 2);
    data{1} = setsubfield(data{1}, dimtok{catdim}, tmp(ind));
    for k = 1:length(param)
      tmp     = getsubfield(data{1}, param{k});
      tmp     = permute(tmp, [catdim setdiff(1:length(size(tmp)), catdim)]);
      tmp     = ipermute(tmp(ind,:,:,:,:), [catdim setdiff(1:length(size(tmp)), catdim)]);
      data{1} = setsubfield(data{1}, param{k}, tmp);
    end
  elseif exist('tmp', 'var') && iscell(tmp)
    %in this case (ugly!) tmp is probably a cell-array containing functional data
  end
  
  % remove unspecified parameters
  rmparam = setdiff(parameterselection('all',data{1}),[param 'pos' 'inside' 'outside']);
  for k = 1:length(rmparam)
    data{1} = rmsubfield(data{1}, rmparam{k});
  end
  
  % keep the first structure only
  data        = data{1};
  dimord      = dimord{1};
  data.dimord = dimord;
  if isfield(data, 'dim'),
    %data.dim    = dim;
    data.dim = size(data.(param{1}));
  end

else
  % nothing to do
  data = data{1};
  dimord = dimord{1};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% from here on the data is concatenated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine the subselection in the data
if selectrpt,
  dimtok = tokenize(data.dimord, '_');
  if strcmp(dimtok{1}, 'rpttap'),
    %account for the tapers
    sumtapcnt = [0;cumsum(data.cumtapcnt(:))];
    begtapcnt = sumtapcnt(1:end-1)+1;
    endtapcnt = sumtapcnt(2:end);
    begtapcnt = begtapcnt(selrpt);
    endtapcnt = endtapcnt(selrpt);
    tapers = zeros(1,sumtapcnt(end));
    for k = 1:length(begtapcnt)
      tapers(begtapcnt(k):endtapcnt(k)) = 1;
    end
    selrpt = find(tapers);
  else
    % do nothing
  end
end

if selectchan,
  selchan = match_str(data.label, channelselection(selchan, data.label));
end

if selectfoi,
  if length(selfoi)==1, selfoi(2) = selfoi; end;
  if length(selfoi)==2,
    %treat selfoi as lower limit and upper limit
    selfoi = nearest(data.freq, selfoi(1)):nearest(data.freq, selfoi(2));
  else
    %treat selfoi as a list of frequencies
    for k=1:length(selfoi)
      tmpfoi(k) = nearest(data.freq, selfoi(k));
    end
    selfoi = tmpfoi;
  end
end

if selecttoi,
  if length(seltoi)==1, seltoi(2) = seltoi; end;
  seltoi = nearest(data.time, seltoi(1)):nearest(data.time, seltoi(2));
end

if selectroi,
  error('not yet implemented');
end

if isfreq,
  if isfield(data, 'labelcmb'),
    %there is a crsspctrm field, this will only be selectdimmed
    %if we apply a trick
    tmpdata = data;
    tmpdata.label = data.labelcmb;
    if selectrpt,  tmpdata = seloverdim(tmpdata, 'rpt',  selrpt);  end
    if selectchan, tmpdata = seloverdim(tmpdata, 'chan', selchan); end
    if selectfoi,  tmpdata = seloverdim(tmpdata, 'freq', selfoi);  end
    if selecttoi,  tmpdata = seloverdim(tmpdata, 'time', seltoi);  end
    % average over dimensions
    if avgoverrpt,  tmpdata = avgoverdim(tmpdata, 'rpt');   end
    if avgoverfreq, tmpdata = avgoverdim(tmpdata, 'freq');  end
    if avgovertime, tmpdata = avgoverdim(tmpdata, 'time');  end
    if dojack,      tmpdata = leaveoneout(tmpdata);         end    
    crsspctrm = tmpdata.crsspctrm; clear tmpdata;
  else
    crsspctrm = [];
  end
  % make the subselection
  if selectrpt,  data = seloverdim(data, 'rpt',  selrpt);  end
  if selectchan, data = seloverdim(data, 'chan', selchan); end
  if selectfoi,  data = seloverdim(data, 'freq', selfoi);  end
  if selecttoi,  data = seloverdim(data, 'time', seltoi);  end
  % average over dimensions
  if avgoverrpt,  data = avgoverdim(data, 'rpt');   end
  if avgoverchan, data = avgoverdim(data, 'chan');  end
  if avgoverfreq, data = avgoverdim(data, 'freq');  end
  if avgovertime, data = avgoverdim(data, 'time');  end
  if dojack,      data = leaveoneout(data);         end    
  if ~isempty(crsspctrm), data.crsspctrm = crsspctrm; end

elseif istlck,
  % make the subselection
  if selectrpt,  data = seloverdim(data, 'rpt',  selrpt);  end
  if selectchan, data = seloverdim(data, 'chan', selchan); end
  if selectfoi,  data = seloverdim(data, 'freq', selfoi);  end
  if selecttoi,  data = seloverdim(data, 'time', seltoi);  end
  % average over dimensions
  if avgoverrpt,  data = avgoverdim(data, 'rpt');   end
  if avgoverchan, data = avgoverdim(data, 'chan');  end
  if avgoverfreq, data = avgoverdim(data, 'freq');  end
  if avgovertime, data = avgoverdim(data, 'time');  end
  if dojack,      data = leaveoneout(data);         end    

elseif issource,
  %FIXME fill in everything
  if selectrpt,  data = seloverdim(data, 'rpt',  selrpt);  end
  if selectfoi,  data = seloverdim(data, 'freq', selfoi);  end
  if avgoverrpt,  data = avgoverdim(data, 'rpt');  end
  if avgoverfreq, data = avgoverdim(data, 'freq'); end

elseif isvolume,
  error('this is not yet implemented');
elseif isfreqmvar,
  % make the subselection
  if selectrpt,  data = seloverdim(data, 'rpt',  selrpt);  end
  if selectchan, data = seloverdim(data, 'chan', selchan); end
  if selectfoi,  data = seloverdim(data, 'freq', selfoi);  end
  if selecttoi,  data = seloverdim(data, 'time', seltoi);  end
  % average over dimensions
  if avgoverrpt,  data = avgoverdim(data, 'rpt');   end
  if avgoverchan, data = avgoverdim(data, 'chan');  end
  if avgoverfreq, data = avgoverdim(data, 'freq');  end
  if avgovertime, data = avgoverdim(data, 'time');  end
  if dojack,      data = leaveoneout(data);         end    
end

%convert back to raw
if israw
  data = checkdata(data, 'datatype', 'raw');
end
