function [data] = selectdata(varargin)

% this function serves to concatenate the input data-structures along the
% compatible dimensions and thus is a more general implementation of
% appenddata, which deals only with raw data. Moreover, it can be used to equate the
% data of different conditions to match e.g. in channels time-axis etc
% Moreover, it can be used as a generalization to ...average with 'keepindividual'
%
% Finally, this function serves to subselect regions-of-interest from the input data,
% either or not averaging across the specified dimensions.
%
% Supported input data:
%   freq
%   timelock
%   source  
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
    %ensure it to have an offset
    data{k} = checkdata(data{k}, 'datatype', 'raw', 'hasoffset', 'yes');
  end
  if strcmp(dtype{k}, 'source'),
    data{k} = checkdata(data{k}, 'sourcerepresentation', 'new');
  end
end

if any(~strmatch(dtype{1},dtype))
  error('different types of input data is not supported');
end

% check consistency of input data
if any(~strmatch(dimord{1},dimord))
  error('a different dimord in the input data is not supported');
end

israw    = strcmp(dtype{1},'raw');
isfreq   = strcmp(dtype{1},'freq');
istlck   = strcmp(dtype{1},'timelock');
issource = strcmp(dtype{1},'source');
isvolume = strcmp(dtype{1},'volume');
isfreqmvar = strcmp(dtype{1},'freqmvar');

selchan  = keyval('channel', kvp); selectchan = ~isempty(strmatch('channel', kvp(cellfun(@ischar, kvp))));
selfoi   = keyval('foilim',  kvp); selectfoi  = ~isempty(strmatch('foilim',  kvp(cellfun(@ischar, kvp))));
seltoi   = keyval('toilim',  kvp); selecttoi  = ~isempty(strmatch('toilim',  kvp(cellfun(@ischar, kvp))));
selroi   = keyval('roi',     kvp); selectroi  = ~isempty(strmatch('roi', kvp(cellfun(@ischar, kvp))));
selrpt   = keyval('rpt',     kvp); selectrpt  = ~isempty(strmatch('rpt', kvp(cellfun(@ischar, kvp))));
selpos   = keyval('pos',     kvp); selectpos  = ~isempty(strmatch('pos', kvp(cellfun(@ischar, kvp))));
param    = keyval('param',   kvp); if isempty(param), param = 'all'; end % FIXME think about this

avgoverchan  = keyval('avgoverchan',  kvp); if isempty(avgoverchan), avgoverchan = false; end
avgoverfreq  = keyval('avgoverfreq',  kvp); if isempty(avgoverfreq), avgoverfreq = false; end
avgovertime  = keyval('avgovertime',  kvp); if isempty(avgovertime), avgovertime = false; end
avgoverroi   = keyval('avgoverroi',   kvp); if isempty(avgoverroi),  avgoverroi  = false; end
avgoverrpt   = keyval('avgoverrpt',   kvp); if isempty(avgoverrpt),  avgoverrpt  = false; end
dojack       = keyval('jackknife',    kvp); if isempty(dojack),      dojack      = false; end

fb       = keyval('feedback', kvp); if isempty(fb), fb = 'yes'; end
if isstr(fb) && strcmp(fb, 'yes'), 
  fb = 1;
elseif isstr(fb) && strcmp(fb, 'no'),
  fb = 0;
end

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

if length(data)>1 && ~israw,
  % determine the way to concatenate
  if issource,
    if isfield(data{1}, 'inside')
      isboolean = islogical(data{1}.inside); 
      if isboolean,
        for k = 1:numel(data)
          data{k} = fixinside(data{k}, 'index');
        end
      end
    end
  end
   
  %if issource || isvolume,
  %  param = parameterselection(param, data{1}); % FIXME check consistency across input data of presence of specific parameters
  %else
  if ~iscell(param), param = {param}; end
  %end
  if issource || isvolume
    if numel(param)>1,
      error('selectdata for source inputs works only for one parameter at a time');
    end
    dimord(:) = {getfield(data{1}, [param{1},'dimord'])};
  end
  dimtok                           = tokenize(dimord{1}, '_');
  dimtok(strmatch('chan', dimtok)) = {'label'}; % data.chan does not exist
  
  dimmat      = zeros(length(dimtok), length(data));
  dimmat(:,1) = 1;
  for k = 1:length(dimtok)
    if isempty(strfind(dimtok{k},'rpt')) && isempty(strfind(dimtok{k},'{pos}')) && isempty(strfind(dimtok{k},'ori')),
      dimdat = getfield(data{1}, dimtok{k});
    elseif ~isempty(strfind(dimtok{k},'{pos}')),
      dimdat = getfield(data{1}, dimtok{k}(2:end-1)); 
    elseif isempty(strfind(dimtok{k},'ori')),
      % dimtok is 'rpt' or 'rpttap'
      dimdat = size(getfield(data{1}, param{1}),1);
    end
    for m = 2:length(data)
      if isempty(strfind(dimtok{k},'rpt')) && isempty(strfind(dimtok{k}, '{pos}')) && isempty(strfind(dimtok{k},'ori')),
        dimdat2 = getfield(data{m},dimtok{k});
      elseif ~isempty(strfind(dimtok{k},'{pos}')),
        dimdat2 = getfield(data{m},dimtok{k}(2:end-1));
      elseif isempty(strfind(dimtok{k},'ori')),
        % dimtok is 'rpt' or 'rpttap'
        dimdat2 = size(getfield(data{m}, param{1}),1);
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
  elseif ~isempty(strfind(dimtok{catdim},'pos'))
    dimtok{catdim} = 'pos';
  elseif isempty(catdim)
    error('don''t know how to concatenate the data');
  end

  % concatenate the data
  % FIXME this works for source data, does this also work for volume data?
  for k = 1:length(param)
    tmp = cell(1,length(data));
    for m = 1:length(tmp)
      tmp{m} = getfield(data{m},param{k});
    end
    if ~iscell(tmp{1}),
      % this is for the 'normal' case
      if catdim==0,
        ndim    = length(size(tmp{1}));
        data{1} = setfield(data{1}, param{k}, permute(cat(ndim+1,tmp{:}),[ndim+1 1:ndim]));
      else
        data{1} = setfield(data{1}, param{k}, cat(catdim,tmp{:}));
      end
    else
      % this is for source data with the positions in a cell-array
      npos = numel(tmp{1});
      if catdim==0,
        error('not implemented yet');
      elseif catdim==1,
        data{1}.(param{k}) = cat(1, tmp{:});
      else
        for kk = 1:npos
          tmpsiz = size(tmp{1}{kk});
          if ~all(tmpsiz==0)
            for kkk = 1:numel(data)
              tmp2{kkk} = tmp{kkk}{kk};
            end          
            data{1}.(param{k}){kk} = cat(catdim-1, tmp2{:});
          else
            %keep empty
          end
        end %for kk = 1:npos
      end %if catdim==0
    end %if ~iscell(tmp{1}) 
    paramdimord{k} = [param{k},'dimord'];
  end %for k = 1:numel(param)
  
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
        tmp       = getsubfield(data{k}, dimtok{catdim})';
	if isfield(data{k}, 'inside'),
	  tmpnvox   = numel(data{k}.inside)+numel(data{k}.outside);
	  tmpinside = data{k}.inside(:);
	end
      else
        if strcmp(dimtok{catdim},'pos')
	  tmp       = [tmp;       getsubfield(data{k}, dimtok{catdim})];
	  tmpinside = [tmpinside; data{k}.inside(:)+tmpnvox];
	  tmpnvox   = tmpnvox+numel(data{k}.inside)+numel(data{k}.outside);
	  sortflag  = 0;
	else
          tmp       = [tmp       getsubfield(data{k}, dimtok{catdim})'];
	  sortflag  = 1;
	end
      end
    end
    data{1} = setsubfield(data{1}, dimtok{catdim}, tmp);
    if isfield(data{1}, 'inside'),
      data{1} = setsubfield(data{1}, 'inside',  tmpinside);
      data{1} = setsubfield(data{1}, 'outside', setdiff(1:size(data{1}.pos,1)', tmpinside));  
    end

    %FIXME think about this
    tryfields = {'dof'};
  else
    % no such field as {'label','time','freq','pos'} has to be concatenated
    sortflag  = 0;
    tryfields = {'cumsumcnt','cumtapcnt','trialinfo'}; 
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
  % FIXME handle inside in previous loop

  % FIXME this is ugly: solve it
  %if issource || isvolume,
  %  data{1}.dim(catdim) = max(size(tmp));
  %end
  
  % sort concatenated data FIXME this is also ugly and depends on tmp
  % FIXME if functional data in cell-array no sorting takes place
  if sortflag && ~iscell(tmp) && ~iscell(data{1}.(param{1})),
    [srt, ind] = sort(tmp, 2);
    data{1} = setfield(data{1}, dimtok{catdim}, tmp(ind));
    for k = 1:length(param)
      tmp     = getsubfield(data{1}, param{k});
      tmp     = permute(tmp, [catdim setdiff(1:length(size(tmp)), catdim)]);
      tmp     = ipermute(tmp(ind,:,:,:,:), [catdim setdiff(1:length(size(tmp)), catdim)]);
      data{1} = setfield(data{1}, param{k}, tmp);
    end
  elseif exist('tmp', 'var') && iscell(tmp)
    %in this case (ugly!) tmp is probably a cell-array containing functional data
  end
  % remove unspecified parameters
  if ~issource,
    rmparam = setdiff(parameterselection('all',data{1}),[param 'pos' 'inside' 'outside']);
  else
    rmparam = setdiff(fieldnames(data{1}), [param(:)' paramdimord(:)' 'pos' 'inside' 'outside' 'dim' 'cfg' 'vol' 'cumtapcnt' 'orilabel']);
  end
  for k = 1:length(rmparam)
    data{1} = rmfield(data{1}, rmparam{k});
  end
  
  % keep the first structure only
  data        = data{1};
  dimord      = dimord{1};
  if ~issource,
    data.dimord = dimord;
  else
    data.([param{1},'dimord']) = dimord;
  end
  if isfield(data, 'dim') & ~issource,
    %data.dim    = dim;
    data.dim = size(data.(param{1}));
  elseif isfield(data, 'dim')
    data     = rmfield(data, 'dim'); %source data should not contain a dim
    %FIXME this should be handled by checkdata once the source structure is
    %unequivocally defined
  end

elseif length(data)>1 && israw
  error('concatenation of several raw data-structures is done by ''ft_appenddata''');
else
  % nothing to do
  data   = data{1};
  dimord = dimord{1};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% from here on the data is concatenated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine the subselection in the data
if selectrpt && ~israw
  if islogical(selrpt),
    selrpt = find(selrpt);
  elseif isempty(selrpt),
    warning('you request all repetitions to be thrown away');
  end
  
  if ~issource
    rpttapflag = ~isempty(strfind(data.dimord, 'rpttap')); 
    rptflag    = ~isempty(strfind(data.dimord, 'rpt')) && ~rpttapflag; 
  else
    fn    = fieldnames(data);
    selfn = find(~cellfun('isempty', strfind(fn, 'dimord')));
    fn    = fn(selfn);
    for k = 1:numel(fn)
      rpttapflag(k,1) = ~isempty(strfind(data.(fn{k}), 'rpttap'));
      rptflag(k,1)    = ~isempty(strfind(data.(fn{k}), 'rpt')) && ~rpttapflag(k,1);
    end
  end    
  
  if any(rpttapflag),
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
elseif selectrpt && israw
  if islogical(selrpt),
    selrpt = find(selrpt);
  elseif isempty(selrpt),
    warning('you request all repetitions to be thrown away');
  end
end

if selectchan,
  %FIXME give selchan according to the order requested in selchan
  %this does not work 
  tmp            = ft_channelselection(selchan, data.label);
  [dum, selchan] = match_str(tmp, data.label);
end

if selectfoi,
  if numel(selfoi)==1, selfoi(2) = selfoi; end;
  if numel(selfoi)==2,
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

if israw,
  if selectrpt,
    data = selfromraw(data, 'rpt', selrpt);
  end

  if selectchan,
    data = selfromraw(data, 'chan', selchan); 
  end

elseif isfreq,
  if isfield(data, 'labelcmb') && isfield(data, 'label') && (selectchan || avgoverchan)
    error('selection of or averaging across channels in the presence of both label and labelcmb is ambiguous');
  end
  
  if isfield(data, 'labelcmb'),
    %there is a crsspctrm or powcovspctrm field, 
    %this will only be selectdimmed
    %if we apply a trick
    tmpdata = data;
    tmpdata.label = data.labelcmb;
    if selectrpt,  tmpdata = seloverdim(tmpdata, 'rpt',  selrpt,  fb); end
    if selectchan, tmpdata = seloverdim(tmpdata, 'chan', selchan, fb); end
    if selectfoi,  tmpdata = seloverdim(tmpdata, 'freq', selfoi,  fb); end
    if selecttoi,  tmpdata = seloverdim(tmpdata, 'time', seltoi,  fb); end
    % average over dimensions
    if avgoverrpt,  tmpdata = avgoverdim(tmpdata, 'rpt',  fb);  end
    if avgoverfreq, tmpdata = avgoverdim(tmpdata, 'freq', fb);  end
    if avgovertime, tmpdata = avgoverdim(tmpdata, 'time', fb);  end
    if avgoverchan, error('avgoverchan not yet implemented in this case'); end
    if dojack,      tmpdata = leaveoneout(tmpdata);         end    

    if isfield(tmpdata, 'crsspctrm'),    crsspctrm = tmpdata.crsspctrm;    end
    if isfield(tmpdata, 'powcovspctrm'), crsspctrm = tmpdata.powcovspctrm; end
    if isfield(tmpdata, 'crsspctrm') || isfield(tmpdata, 'powcovspctrm'), clear tmpdata; end 
  else
    crsspctrm = [];
  end
  % make the subselection
  if selectrpt,  data = seloverdim(data, 'rpt',  selrpt,  fb); end
  if selectchan, data = seloverdim(data, 'chan', selchan, fb); end
  if selectfoi,  data = seloverdim(data, 'freq', selfoi,  fb); end
  if selecttoi,  data = seloverdim(data, 'time', seltoi,  fb); end
  % average over dimensions
  if avgoverrpt,  data = avgoverdim(data, 'rpt',  fb);  end
  if avgoverchan, data = avgoverdim(data, 'chan', fb);  end
  if avgoverfreq, data = avgoverdim(data, 'freq', fb);  end
  if avgovertime, data = avgoverdim(data, 'time', fb);  end
  if dojack,      data = leaveoneout(data);             end    
  
  if ~isempty(crsspctrm),
    if isfield(data, 'crsspctrm'),    data.crsspctrm    = crsspctrm; end
    if isfield(data, 'powcovspctrm'), data.powcovspctrm = crsspctrm; end
  end

elseif istlck,
  % make the subselection
  if selectrpt,  data = seloverdim(data, 'rpt',  selrpt,  fb); end
  if selectchan, data = seloverdim(data, 'chan', selchan, fb); end
  if selectfoi,  data = seloverdim(data, 'freq', selfoi,  fb); end
  if selecttoi,  data = seloverdim(data, 'time', seltoi,  fb); end
  % average over dimensions
  if avgoverrpt,  data = avgoverdim(data, 'rpt',  fb); end
  if avgoverchan, data = avgoverdim(data, 'chan', fb); end
  if avgoverfreq, data = avgoverdim(data, 'freq', fb); end
  if avgovertime, data = avgoverdim(data, 'time', fb); end
  if dojack,      data = leaveoneout(data);         end    

elseif issource,
  %FIXME fill in everything
  if selectrpt,  data = seloverdim(data, 'rpt',  selrpt, fb); end
  if selectfoi,  data = seloverdim(data, 'freq', selfoi, fb); end
  if selectpos,  data = seloverdim(data, 'pos',  selpos, fb); end
  if avgoverrpt,  data = avgoverdim(data, 'rpt'); end
  if avgoverfreq, data = avgoverdim(data, 'freq'); end

elseif isvolume,
  error('this is not yet implemented');
elseif isfreqmvar,
  % make the subselection
  if selectrpt,  data = seloverdim(data, 'rpt',  selrpt,  fb); end
  if selectchan, data = seloverdim(data, 'chan', selchan, fb); end
  if selectfoi,  data = seloverdim(data, 'freq', selfoi,  fb); end
  if selecttoi,  data = seloverdim(data, 'time', seltoi,  fb); end
  % average over dimensions
  if avgoverrpt,  data = avgoverdim(data, 'rpt',  fb); end
  if avgoverchan, data = avgoverdim(data, 'chan', fb); end
  if avgoverfreq, data = avgoverdim(data, 'freq', fb); end
  if avgovertime, data = avgoverdim(data, 'time', fb); end
  if dojack,      data = leaveoneout(data);            end    
end
