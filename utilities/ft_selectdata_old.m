function [data] = ft_selectdata_old(varargin)

% FT_SELECTDATA_OLD is deprecated, please use FT_SELECTDATA instead.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Old documentation for reference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function serves to subselect regions-of-interest from the input data,
% with or without averaging across the specified dimensions. It also
% concatenates multiple input data structures along the compatible
% dimension and is thus a more general implementation of ft_appenddata,
% which deals only with raw data.
%
% Use as
%  [data] = ft_selectdata(data1, data2, ..., key1, value1, key2, value2, ...)
%
% Supported input data:
%   freq
%   timelock
%   source
%   volume   (not yet)
%   raw      only for subselection of channels and replicates. this is the
%              same functionality as ft_preprocessing
%
% Optional input arguments should be specified as key-value pairs and may include
%   param         string
%   foilim        [begin end]  edges of frequency band to be retained
%   toilim        [begin end]  edges of time window to be retained
%   roi           [Nx1]        indices of voxels of region-of-interest
%   rpt           [Nx1]        indices of replicates to be retained
%   channel       {Nx1}        list of channels labels to be retained
%   avgoverchan   'no' ('yes') average across channels
%   avgoverfreq   'no' ('yes') average across frequency bins
%   avgovertime   'no' ('yes') average across time points
%   avgoverroi    'no' ('yes') average across voxels in ROI
%   avgoverrpt    'no' ('yes') average across replicates
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2009-2011, Jan-Mathijs Schoffelen
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

% FIXME ROI selection is not yet implemented

% the input consists of one or multiple data structures, followed by the optional key-value pairs
isdata   = find(cellfun(@isstruct,varargin));  % FIXME this fails in case one of the key-value pairs contains a structure
iskeyval = setdiff(1:length(varargin),isdata);
data     = varargin(isdata);
varargin = varargin(iskeyval);

% go over all input data structures
dtype   = cell(1,length(data));
dimord  = cell(1,length(data));
hassubj = false(1, length(data));
for k = 1:length(data)
  data{k} = ft_checkdata(data{k}, 'datatype', {'freq', 'freq+comp', 'timelock', 'timelock+comp', 'raw', 'raw+comp', 'source', 'volume', 'freqmvar', 'chan'});
  if isfield(data{k}, 'dimord') && ~isempty(strfind(data{k}.dimord, 'subj'))
    hassubj(k)     = true;
    data{k}.dimord = strrep(data{k}.dimord, 'subj', 'rpt');
  end
  [dtype{k}, dimord{k}]  = ft_datatype(data{k});
  if strcmp(dtype{k}, 'raw'),
    data{k} = ft_checkdata(data{k}, 'datatype', 'raw');
  end
  if strcmp(dtype{k}, 'source'),
    data{k} = ft_checkdata(data{k}, 'sourcerepresentation', 'new');
  end
end

if ~all(strcmp(dtype{1},dtype))
  error('the data type is not consistent for all inputs');
end

% check consistency of input data
if ~all(strcmp(dimord{1},dimord))
  error('the dimord is not consistent for all inputs');
end

israw      = strcmp(dtype{1},'raw') || strcmp(dtype{1},'comp'); % comp can be treated as raw
isfreq     = strcmp(dtype{1},'freq');
istlck     = strcmp(dtype{1},'timelock');
issource   = strcmp(dtype{1},'source');
isvolume   = strcmp(dtype{1},'volume');
isfreqmvar = strcmp(dtype{1},'freqmvar');

% get the optional arguments
selchan      = ft_getopt(varargin, 'channel',      'all', 1); selectchan = ~(ischar(selchan) && strcmp(selchan, 'all'));
selfoi       = ft_getopt(varargin, 'foilim',       'all', 1); selectfoi  = ~(ischar(selfoi)  && strcmp(selfoi,  'all'));
seltoi       = ft_getopt(varargin, 'toilim',       'all', 1); selecttoi  = ~(ischar(seltoi)  && strcmp(seltoi,  'all'));
selroi       = ft_getopt(varargin, 'roi',          []); selectroi  = ~isempty(selroi);
selrpt       = ft_getopt(varargin, 'rpt',          'all', 1); selectrpt  = ~(ischar(selrpt)  && strcmp(selrpt,  'all'));
selpos       = ft_getopt(varargin, 'pos',          []); selectpos  = ~isempty(selpos);
param        = ft_getopt(varargin, 'param',        'all'); % FIXME think about this
avgoverchan  = ft_getopt(varargin, 'avgoverchan',  false);
avgoverfreq  = ft_getopt(varargin, 'avgoverfreq',  false);
avgovertime  = ft_getopt(varargin, 'avgovertime',  false);
avgoverroi   = ft_getopt(varargin, 'avgoverroi',   false);
avgoverrpt   = ft_getopt(varargin, 'avgoverrpt',   false);
dojack       = ft_getopt(varargin, 'jackknife',    false);
fb           = ft_getopt(varargin, 'feedback',     true);
% FIXME implement toi and foi

% ensure that these are boolean arguments, optionally convert from "yes"/"no" to true/false
fb          = istrue(fb);
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

% a quick check to ensure the user does not use this function in cases
% where it is known to contain a bug
%if selectrpt && (isempty(selrpt) || ~any(selrpt))
%  error('ft_selectdata_old does not work when selecting 0 trials; please use ft_selectdata_new instead (use a cfg input, instead of key-value pairs, to ft_selectdata)');
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% concatenate the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(data)>1 && ~israw,
  % determine the way to concatenate
  
  % force inside field for source data to be indices
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
      error('ft_selectdata for source inputs works only for one parameter at a time');
    end
    dimord(:) = {data{1}.([param{1},'dimord'])};
  end
  dimtok                           = tokenize(dimord{1}, '_');
  dimtok(strmatch('chan', dimtok)) = {'label'}; % data.chan does not exist
  
     dimmat      = zeros(length(dimtok), length(data));
     dimmat(:,1) = 1;
     for k = 1:length(dimtok)
       if isempty(strfind(dimtok{k},'rpt')) && isempty(strfind(dimtok{k},'{pos}')) && isempty(strfind(dimtok{k},'ori')),
         dimdat = data{1}.(dimtok{k});
       elseif ~isempty(strfind(dimtok{k},'{pos}')),
         dimdat = data{1}.(dimtok{k}(2:end-1));
       elseif isempty(strfind(dimtok{k},'ori')),
         % dimtok is 'rpt' or 'rpttap'
         dimdat = size(data{1}.(param{1}),1);
       end
       for m = 2:length(data)
         if isempty(strfind(dimtok{k},'rpt')) && isempty(strfind(dimtok{k}, '{pos}')) && isempty(strfind(dimtok{k},'ori')),
           dimdat2 = data{m}.(dimtok{k});
         elseif ~isempty(strfind(dimtok{k},'{pos}')),
           dimdat2 = data{m}.(dimtok{k}(2:end-1));
         elseif isempty(strfind(dimtok{k},'ori')),
           % dimtok is 'rpt' or 'rpttap'
           dimdat2 = size(data{m}.(param{1}),1);
         end
         try, dimmat(k,m) = all(dimdat(:)==dimdat2(:));            catch end;
         try, dimmat(k,m) = all(cellfun(@isequal,dimdat,dimdat2)); catch end;
       end
     end
    catdim = find(sum(dimmat,2)<length(data));

%   if any(strcmp(dimtok, 'rpt'))
%     catdim = find(strcmp(dimtok, 'rpt'));
%   elseif any(strcmp(dimtok, 'rpttap'))
%     catdim = find(strcmp(dimtok, 'rpttap'));
%   elseif any(strcmp(dimtok, 'subj'))
%     catdim = find(strcmp(dimtok, 'subj'));
%   else
%     catdim = [];
%   end
%   
  if length(catdim)>1,
    error('ambiguous dimensions for concatenation');
  elseif isempty(catdim) && isempty(intersect(dimtok, {'rpt', 'rpttap', 'subj'}))
    % treat as individual observations: prepend a first dimension 'rpt'
    % (so this part should be able to cover the functionality of ...grandaverage)
    catdim = 0;
  elseif isempty(catdim) && (~isempty(strmatch('rpt',dimtok)) || ~isempty(strmatch('rpttap',dimtok)))
    %append observations
    catdim = find(~cellfun('isempty',strfind(dimtok, 'rpt')));
  elseif ~isempty(strfind(dimtok{catdim},'pos'))
    dimtok{catdim} = 'pos';
  elseif isempty(catdim)
    error('don''t know how to concatenate the data');
  end
  
  % concatenate the data
  % FIXME this works for source data, does this also work for volume data?
  for k = 1:length(param)
    tmp = cell(1,length(data));
    % try to get the numeric data 'param{k}' if present
    try
      for m = 1:length(tmp)
        tmp{m} = data{m}.(param{k});
      end
    catch
      continue;
    end
    if ~iscell(tmp{1}),
      % this is for the 'normal' case
      if catdim==0,
        ndim    = length(size(tmp{1}));
        datacat.(param{k}) = permute(cat(ndim+1,tmp{:}),[ndim+1 1:ndim]);
      else
        datacat.(param{k}) = cat(catdim,tmp{:});
      end
    else
      % this is for source data with the positions in a cell-array
      npos = numel(tmp{1});
      if catdim==0,
        error('not implemented yet');
      elseif catdim==1,
        datacat.(param{k}) = cat(1, tmp{:});
      else
        for kk = 1:npos
          tmpsiz = size(tmp{1}{kk});
          if ~all(tmpsiz==0)
            for kkk = 1:numel(data)
              tmp2{kkk} = tmp{kkk}{kk};
            end
            datacat.(param{k}){kk} = cat(catdim-1, tmp2{:});
          else
            %keep empty
          end
        end %for kk = 1:npos
      end %if catdim==0
    end %if ~iscell(tmp{1})
    paramdimord{k} = [param{k},'dimord'];
  end %for k = 1:numel(param)
  
  if catdim==0,
    % a dimension has been prepended
    dimtok    = ['rpt' dimtok];
    catdim    = 1;
    dimord{1} = ['rpt_',dimord{1}];
    
    % adjust dim-field
    if issubfield(data{1}, 'dim'),
      dim       = [length(data) data{1}.dim];
    end
    
    % adjust inside-field according to intersection
    if isfield(data{1}, 'inside')
      isboolean = islogical(data{1}.inside);
      for k = 1:numel(data)
        if k==1,
          if isboolean
            inside = double(data{k}.inside);
          else
            inside = zeros(numel(data{k}.inside)+numel(data{k}.outside),1);
            inside(data{k}.inside) = inside(data{k}.inside) + 1;
          end
        else
          if isboolean
            inside = double(data{k}.inside) + inside;
          else
            inside(data{k}.inside) = inside(data{k}.inside) + 1;
          end
        end
      end
      
      % determine which sources were inside or outside the brain in all subjects
      nalloutside = sum(inside(:)==0);
      nsomeinside = sum(inside(:)>0 & inside(:)~=numel(data));
      inside      = inside==numel(data);
      nallinside  = sum(inside(:));
      
      fprintf('%d voxels are inside the brain of all subjects\n',               nallinside);
      fprintf('%d voxels are inside the brain of some, but not all subjects\n', nsomeinside);
      fprintf('%d voxels are outside the brain of all subjects\n',              nalloutside);
      warning('marking only voxels inside the brain of all subjects as ''inside''');
      
      if isboolean
        data{1}.inside = inside;
      else
        data{1}.inside  = find(inside);
        data{1}.outside = setdiff((1:numel(inside))', data{1}.inside);
      end
      
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
        tmp = getsubfield(data{k}, dimtok{catdim})';
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
        elseif strcmp(dimtok{catdim}, 'time') || strcmp(dimtok{catdim}, 'freq')
          tmp       = [tmp(:)' data{k}.(dimtok{catdim})];
          sortflag  = 1;
        else
          tmp       = [tmp(:); data{k}.(dimtok{catdim})];
          sortflag  = 0;
        end
      end
    end
    datacat = setsubfield(datacat, dimtok{catdim}, tmp);
    if isfield(datacat, 'inside'),
      datacat = setsubfield(datacat, 'inside',  tmpinside);
      datacat = setsubfield(datacat, 'outside', setdiff(1:size(data{1}.pos,1)', tmpinside));
      datacat = setsubfield(datacat, 'pos', data{1}.pos);
    end
    
    %FIXME think about this
    tryfields = {'dof'};
  else
    % no such field as {'label','time','freq','pos'} has to be concatenated
    sortflag  = 0;
    tryfields = {'cumsumcnt','cumtapcnt','trialinfo','sampleinfo'};
  end
  % add additional descriptive fields
  if isfield(data{1}, 'label'), datacat.label = data{1}.label; end
  if isfield(data{1}, 'freq'),  datacat.freq  = data{1}.freq;  end
  if isfield(data{1}, 'time'),  datacat.time  = data{1}.time;  end
  if isfield(data{1}, 'cumtapcnt'), datacat.cumtapcnt = data{1}.cumtapcnt; end
  if isfield(data{1}, 'cumsumcnt'), datacat.cumsumcnt = data{1}.cumsumcnt; end
  if isfield(data{1}, 'trialinfo'), datacat.trialinfo = data{1}.trialinfo; end
  if isfield(data{1}, 'sampleinfo'), datacat.sampleinfo = data{1}.sampleinfo; end
  if isfield(data{1}, 'labelcmb'),  datacat.labelcmb  = data{1}.labelcmb; end
  
  
  % concatenate the relevant descriptive fields in the data-structure (continued)
  for k = 1:length(tryfields)
    try
      for m = 1:length(data)
        if m==1,
          tmpfield = data{m}.(tryfields{k});
        else
          tmpfield = [tmpfield; data{m}.(tryfields{k})];
        end
      end
      datacat.(tryfields{k}) = tmpfield;
    catch
    end
  end
  % FIXME handle inside in previous loop
  
  % FIXME this is ugly: solve it
  %if issource || isvolume,
  %  datacat.dim(catdim) = max(size(tmp));
  %end
  
  % sort concatenated data FIXME this is also ugly and depends on tmp
  % FIXME if functional data in cell-array no sorting takes place
  if sortflag && ~iscell(tmp) && ~iscell(data{1}.(param{1})),
    [srt, ind] = sort(tmp, 2);
    data{1}.(dimtok{catdim}) = tmp(ind);
    for k = 1:length(param)
      try
        tmp     = data{1}.(param{k});
      catch
        continue;
      end
      tmp     = permute(tmp, [catdim setdiff(1:length(size(tmp)), catdim)]);
      tmp     = ipermute(tmp(ind,:,:,:,:), [catdim setdiff(1:length(size(tmp)), catdim)]);
      datacat.(param{k}) = tmp;
    end
  elseif exist('tmp', 'var') && iscell(tmp)
    %in this case (ugly!) tmp is probably a cell-array containing functional data
  end
  % remove unspecified parameters
  if ~issource,
    %rmparam = setdiff(parameterselection('all',datacat),[param 'pos' 'inside' 'outside' 'freq' 'time']);
    rmparam = {};
  else
    rmparam = setdiff(fieldnames(datacat), [param(:)' paramdimord(:)' 'pos' 'inside' 'outside' 'dim' 'cfg' 'vol' 'cumtapcnt' 'orilabel' 'time' 'freq']);
  end
  for k = 1:length(rmparam)
    datacat = rmfield(datacat, rmparam{k});
  end
  
  
  % keep the first structure only
  data        = datacat;
  dimord      = dimord{1};
  if ~issource,
    data.dimord = dimord;
  else
    data.([param{1},'dimord']) = dimord;
  end
  if isfield(data, 'dim'),
    data.dim    = dim;
    %data.dim = size(data.(param{1}));
  elseif isfield(data, 'dim')
    data     = rmfield(data, 'dim'); %source data should not contain a dim
    %FIXME this should be handled by ft_checkdata once the source structure is
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
      tapers(begtapcnt(k):endtapcnt(k)) = k;
    end
    selrpt   = find(tapers);
    [srt,ix] = sort(tapers(tapers~=0));
    selrpt   = selrpt(ix);
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
  
  if isfield(data, 'label')
    tmp            = ft_channelselection(selchan, data.label);
    [dum, selchan] = match_str(tmp, data.label);
  elseif isfield(data, 'labelcmb')
    tmp            = ft_channelselection(selchan, unique(data.labelcmb(:)));
    [dum, selchan1] = match_str(tmp, data.labelcmb(:,1));
    [dum, selchan2] = match_str(tmp, data.labelcmb(:,2));
    selchan         = intersect(selchan1, selchan2);
  end
  
end

if selectfoi,
  if numel(selfoi)==1, selfoi(2) = selfoi; end;
  if numel(selfoi)==2,
    % treat selfoi as lower limit and upper limit
    selfoi = nearest(data.freq, selfoi(1)):nearest(data.freq, selfoi(2));
    % selfoi = find(data.freq>=selfoi(1) & data.freq<=selfoi(2));
  else
    % treat selfoi as a list of frequencies
    tmpfoi = zeros(1,numel(selfoi));
    for k=1:length(selfoi)
      tmpfoi(k) = nearest(data.freq, selfoi(k));
    end
    selfoi = tmpfoi;
  end % numel
end

if selecttoi && ~israw,
  if length(seltoi)==1, seltoi(2) = seltoi; end;
  if numel(seltoi)==2,
    % treat seltoi as lower limit and upper limit
    toitmp=nearest(data.time,[seltoi(1) seltoi(2)]);
    seltoi=toitmp(1):toitmp(2);
%     seltoi = nearest(data.time, seltoi(1)):nearest(data.time, seltoi(2));
%     seltoi = find(data.time>=seltoi(1) & data.time<=seltoi(2));
  else
    % treat seltoi as a list of timepoints
    tmptoi = zeros(1,numel(seltoi));
    for k=1:length(seltoi)
      tmptoi(k) = nearest(data.time, seltoi(k));
    end
    seltoi = tmptoi;
  end % numel
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
  
  if selecttoi,
    data = selfromraw(data, 'latency', seltoi);
  end
  
elseif isfreq,
  %   if isfield(data, 'labelcmb') && isfield(data, 'label') && (selectchan || avgoverchan)
  %     error('selection of or averaging across channels in the presence of both label and labelcmb is not possible');
  %   end
  tmpdata = [];
  if isfield(data, 'labelcmb'),
    % there is a crsspctrm or powcovspctrm field,
    % this will only be selectdimmed
    % if we apply a trick
    tmpdata = data;
    tmpdata.label = data.labelcmb;
    try, tmpdata = rmfield(tmpdata, 'powspctrm'); end
    
    % selection of combinations
    if selectchan && isfield(data, 'label')
      selcmb = find(sum(ismember(tmpdata.label, data.label(selchan)),2)==2);
    elseif selectchan
      error('this is not yet implemented');
    end
    
    % make the subselection
    if selectrpt,  tmpdata = seloverdim(tmpdata, 'rpt',  selrpt,  fb); end
    if selectchan, tmpdata = seloverdim(tmpdata, 'chan', selcmb,  fb); end
    if selectfoi,  tmpdata = seloverdim(tmpdata, 'freq', selfoi,  fb); end
    if selecttoi,  tmpdata = seloverdim(tmpdata, 'time', seltoi,  fb); end
    % average over dimensions
    if avgoverrpt,  tmpdata = avgoverdim(tmpdata, 'rpt',  fb);  end
    if avgoverfreq, tmpdata = avgoverdim(tmpdata, 'freq', fb);  end
    if avgovertime, tmpdata = avgoverdim(tmpdata, 'time', fb);  end
    if avgoverchan, tmpdata = avgoverdim(tmpdata, 'chan', fb);  end
    if dojack,      tmpdata = leaveoneout(tmpdata);         end
  end
  
  if isfield(data, 'label'),
    % make the subselection
    if selectrpt,  data = seloverdim(data, 'rpt',  selrpt,  fb); end
    if selectchan, data = seloverdim(data, 'chan', selchan, fb); end
    if selectfoi,  data = seloverdim(data, 'freq', selfoi,  fb); end
    if selecttoi,  data = seloverdim(data, 'time', seltoi,  fb); end
    % average over dimensions
    if avgoverrpt,  data = avgoverdim(data, 'rpt',  fb);  end
    if avgoverfreq, data = avgoverdim(data, 'freq', fb);  end
    if avgovertime, data = avgoverdim(data, 'time', fb);  end
    if avgoverchan, data = avgoverdim(data, 'chan', fb); end
    if dojack,      data = leaveoneout(data);            end
  end
  
  if isstruct(tmpdata)
    param = selparam(tmpdata);
    for k = 1:numel(param)
      data.(param{k}) = tmpdata.(param{k});
    end
  end
  
elseif istlck,
  % make the subselection
  if selectrpt,  data = seloverdim(data, 'rpt',  selrpt,  fb); end
  if selectchan, data = seloverdim(data, 'chan', selchan, fb); end
  if selectfoi,  data = seloverdim(data, 'freq', selfoi,  fb); end
  if selecttoi,  data = seloverdim(data, 'time', seltoi,  fb); end
  if isfield(data,'trial') && isfield(data,'avg') %&& size(data.trial,3)~=size(data.avg,2)
    warning('Warning: .avg, .var, .dof and .cov not updated.');
    if isfield(data, 'avg'), data = rmfield(data, 'avg'); end
    if isfield(data, 'cov'), data = rmfield(data, 'cov'); end
    if isfield(data, 'dof'), data = rmfield(data, 'dof'); end
    if isfield(data, 'var'), data = rmfield(data, 'var'); end
  end
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

if hassubj(1)
  data.dimord = strrep(data.dimord, 'rpt', 'subj');
end
