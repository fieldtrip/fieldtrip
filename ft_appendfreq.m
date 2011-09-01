function [freq] = ft_appendfreq(cfg, varargin)

% FT_APPENDFREQ concatenates multiple frequency or time-frequency data
% structures that have been processed seperately. If the input data
% structures contain different channels, it will be concatenated along the
% channel direction. If the channels are identical in the input data
% structures, the data will be concatenated along the repetition dimension.
%
% Use as
%   combined = ft_appendfreq(cfg, freq1, freq2, ...)
%
% See also FT_FREQANALYSIS, FT_APPENDDATA, FT_APPENDTIMELOCK, FT_APPENDSOURCE

% Copyright (C) 2011, Robert Oostenveld
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

ft_defaults

% record start time and total processing time
ftFuncTimer = tic();
ftFuncClock = clock();

cfg = ft_checkconfig(cfg, 'required', 'parameter');

% set the defaults
cfg.inputfile  = ft_getopt(cfg, 'inputfile',  []);
cfg.outputfile = ft_getopt(cfg, 'outputfile', []);
cfg.appenddim  = ft_getopt(cfg, 'appenddim',  'auto');


hasdata = nargin>1;
if ~isempty(cfg.inputfile) % the input data should be read from file
  if hasdata
    error('cfg.inputfile should not be used in conjunction with giving input data to this function');
  elseif ~iscell(cfg.inputfile)
    error('you should specify cfg.inpoutfile as cell-array with multiple file names');
  else
    for i=1:numel(cfg.inputfile)
      varargin{i} = loadvar(cfg.inputfile{i}, 'freq'); % read datasets from array inputfile
    end
  end
end

Ndata = numel(varargin);
if Ndata==1
  % nothing to do
  freq = varargin{1};
  return;
end

% check if the input data is valid for this function
for i=1:Ndata
  varargin{i} = ft_checkdata(varargin{i}, 'datatype', 'freq', 'feedback', 'yes');
end

% do a basic check to see whether the dimords match
dimordmatch = true;
dimord      = cell(1,Ndata);
for i=1:Ndata
  dimord{i} = varargin{i}.dimord;
  if i>1
    dimordmatch = strcmp(dimord{1}, dimord{i}) && dimordmatch;
  end
end

if ~dimordmatch
  error('the dimords of the input data structures are not equal');
end

% create the output structure from scratch
freq   = [];

dimtok = tokenize(dimord{1}, '_');
switch cfg.appenddim
  case 'auto'
    % to be implemented
  case {'rpt' 'rpttap' 'subj'}
    catdim = find(ismember(dimtok, {'rpt' 'rpttap' 'subj'}));
    if numel(catdim)==0
      catdim = 0;
    elseif numel(catdim)==1
      % this is OK
    elseif numel(catdim)>1
      error('ambiguous dimord for concatenation');
    end
    
    % if any of these are present, concatenate
    % if not prepend the dimord with rpt (and thus shift the dimensions)
    
    % here we need to check whether the other dimensions are the same. if
    % not, consider some tolerance.
    boolval1 = ~checkchan(varargin{:});
    boolval2 = ~checkfreq(varargin{:});
    if isfield(varargin{1}, 'time'),
      boolval3 = ~checktime(varargin{:});
    else 
      boolval3 = true;
    end
    
    if any([boolval1 boolval2 boolval3]==false)
      error('appending across observations is not possible, because the dimensions are incompatible');
    end
    
    % update the dimord
    if catdim==0
      freq.dimord = ['rpt_',varargin{1}.dimord];
    else
      freq.dimord = varargin{1}.dimord;
    end
    
    % fill in the rest of the descriptive fields
    freq.label = varargin{1}.label;
    freq.freq  = varargin{1}.freq;
    if isfield(freq, 'time'), freq.time = varargin{1}.time; end
    
  case 'chan'
    catdim = strmatch('chan', dimtok);
    if isempty(catdim)
      % try chancmb
      catdim = strmatch('chancmb', dimtok);
    end
    
    % check whether all channels are unique and throw an error if not
    [boolval, list] = checkchan(varargin{:});
    if ~boolval
      error('the input data structures have non-unique channels, concatenation across channel is not possible');
    end
    
    % update the channel description
    freq.label = list;
    
    % fill in the rest of the descriptive fields
    freq.freq  = varargin{1}.freq;
    if isfield(freq, 'time'), freq.time = varargin{1}.time; end
    freq.dimord = varargin{1}.dimord;
    
  case 'freq'
    catdim = strmatch('freq', dimtok);
    
    % check whether all frequencies are unique and throw an error if not
    [boolval, list] = checkfreq(varargin{:});
    if ~boolval
      error('the input data structures have non-unique frequency bins, concatenation across frequency is not possible');
    end
    
    % update the frequency description
    varargin{1}.freq = list;
    
    % fill in the rest of the descriptive fields
    freq.label  = varargin{1}.label;
    freq.dimord = varargin{1}.dimord;
    if isfield(freq, 'time'), freq.time = varargin{1}.time; end
    
  case 'time'
    catdim = strmatch('time', dimtok);
    
    % check whether all time points are unique and throw an error if not 
    [boolval, list] = checktime(varargin{:});
    if ~boolval
      error('the input data structures have non-unique time bins, concatenation across time is not possible');
    end
    
    % update the time description
    varargin{1}.time = list;
    
    % fill in the rest of the descriptive fields
    freq.label  = varargin{1}.label;
    freq.freq   = varargin{1}.freq;
    freq.dimord = varargin{1}.dimord;
    
  otherwise
end

%param = selparam(varargin{2}); % use the second one because the dimord of the first has been adjusted
% FIXME do a check on whether the parameters are present in all datasets
param = cfg.parameter;
if ~iscell(param), param = {param}; end

% concatenate the data
for k = 1:numel(param)
  tmp = cell(1,Ndata);
  % get the numeric data 'param{k}' if present
  for m = 1:Ndata
    tmp{m} = varargin{m}.(param{k});
  end
  
  if catdim==0,
    ndim    = length(size(tmp{1}));
    freq.(param{k}) = permute(cat(ndim+1,tmp{:}),[ndim+1 1:ndim]);
  else
    freq.(param{k}) = cat(catdim,tmp{:});
  end
end %for k = 1:numel(param)





% % use a helper function to select the consistent parts of the data and to concatenate it
% freq = ft_selectdata(varargin{:}, 'param', {'powspctrm' 'crsspctrm' 'fourierspctrm'});
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % concatenate the data
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% if length(data)>1 && ~israw,
%   % determine the way to concatenate
%   
%   if ~iscell(param), param = {param}; end
%   %end
%   
%   dimtok                           = tokenize(dimord{1}, '_');
%   dimtok(strmatch('chan', dimtok)) = {'label'}; % data.chan does not exist
%   
%   dimmat      = zeros(length(dimtok), length(data));
%   dimmat(:,1) = 1;
%   for k = 1:length(dimtok)
%     if isempty(strfind(dimtok{k},'rpt')) && isempty(strfind(dimtok{k},'{pos}')) && isempty(strfind(dimtok{k},'ori')),
%       dimdat = data{1}.(dimtok{k});
%     elseif ~isempty(strfind(dimtok{k},'{pos}')),
%       dimdat = data{1}.(dimtok{k}(2:end-1));
%     elseif isempty(strfind(dimtok{k},'ori')),
%       % dimtok is 'rpt' or 'rpttap'
%       dimdat = size(data{1}.(param{1}),1);
%     end
%     for m = 2:length(data)
%       if isempty(strfind(dimtok{k},'rpt')) && isempty(strfind(dimtok{k}, '{pos}')) && isempty(strfind(dimtok{k},'ori')),
%         dimdat2 = data{m}.(dimtok{k});
%       elseif ~isempty(strfind(dimtok{k},'{pos}')),
%         dimdat2 = data{m}.(dimtok{k}(2:end-1));
%       elseif isempty(strfind(dimtok{k},'ori')),
%         % dimtok is 'rpt' or 'rpttap'
%         dimdat2 = size(data{m}.(param{1}),1);
%       end
%       try, dimmat(k,m) = all(dimdat(:)==dimdat2(:));            catch end;
%       try, dimmat(k,m) = all(cellfun(@isequal,dimdat,dimdat2)); catch end;
%     end
%   end
%   catdim = find(sum(dimmat,2)<length(data));
%   
%   if length(catdim)>1,
%     error('ambiguous dimensions for concatenation');
%   elseif isempty(catdim) && isempty(strmatch('rpt',dimtok)) && isempty(strmatch('rpttap',dimtok)),
%     %treat as individual observations: prepend a first dimension 'rpt'
%     %(so this part should be able to cover the functionality of ...grandaverage)
%     catdim = 0;
%   elseif isempty(catdim) && (~isempty(strmatch('rpt',dimtok)) || ~isempty(strmatch('rpttap',dimtok)))
%     %append observations
%     catdim = find(~cellfun('isempty',strfind(dimtok, 'rpt')));
%   elseif ~isempty(strfind(dimtok{catdim},'pos'))
%     dimtok{catdim} = 'pos';
%   elseif isempty(catdim)
%     error('don''t know how to concatenate the data');
%   end
%   
%   % concatenate the data
%   % FIXME this works for source data, does this also work for volume data?
%   for k = 1:length(param)
%     tmp = cell(1,length(data));
%     % try to get the numeric data 'param{k}' if present
%     try
%       for m = 1:length(tmp)
%         tmp{m} = data{m}.(param{k});
%       end
%     catch
%       continue;
%     end
%     if ~iscell(tmp{1}),
%       % this is for the 'normal' case
%       if catdim==0,
%         ndim    = length(size(tmp{1}));
%         data{1}.(param{k}) = permute(cat(ndim+1,tmp{:}),[ndim+1 1:ndim]);
%       else
%         data{1}.(param{k}) = cat(catdim,tmp{:});
%       end
%     else
%       % this is for source data with the positions in a cell-array
%       npos = numel(tmp{1});
%       if catdim==0,
%         error('not implemented yet');
%       elseif catdim==1,
%         data{1}.(param{k}) = cat(1, tmp{:});
%       else
%         for kk = 1:npos
%           tmpsiz = size(tmp{1}{kk});
%           if ~all(tmpsiz==0)
%             for kkk = 1:numel(data)
%               tmp2{kkk} = tmp{kkk}{kk};
%             end
%             data{1}.(param{k}){kk} = cat(catdim-1, tmp2{:});
%           else
%             %keep empty
%           end
%         end %for kk = 1:npos
%       end %if catdim==0
%     end %if ~iscell(tmp{1})
%     paramdimord{k} = [param{k},'dimord'];
%   end %for k = 1:numel(param)
%   
%   if catdim==0,
%     % a dimension has been prepended
%     dimtok    = ['rpt' dimtok];
%     catdim    = 1;
%     dimord{1} = ['rpt_',dimord{1}];
%   end
%   
%   % concatenate the relevant descriptive fields in the data-structure
%   if ~strcmp(dimtok{catdim},'rpt') && ~strcmp(dimtok{catdim},'rpttap'),
%     for k = 1:length(data)
%       if k==1,
%         tmp = getsubfield(data{k}, dimtok{catdim})';
%         if isfield(data{k}, 'inside'),
%           tmpnvox   = numel(data{k}.inside)+numel(data{k}.outside);
%           tmpinside = data{k}.inside(:);
%         end
%       else
%         if strcmp(dimtok{catdim},'pos')
%           tmp       = [tmp;       getsubfield(data{k}, dimtok{catdim})];
%           tmpinside = [tmpinside; data{k}.inside(:)+tmpnvox];
%           tmpnvox   = tmpnvox+numel(data{k}.inside)+numel(data{k}.outside);
%           sortflag  = 0;
%         elseif strcmp(dimtok{catdim}, 'time') || strcmp(dimtok{catdim}, 'freq')
%           tmp       = [tmp(:)' data{k}.(dimtok{catdim})];
%           sortflag  = 1;
%         else
%           tmp       = [tmp(:); data{k}.(dimtok{catdim})];
%           sortflag  = 0;
%         end
%       end
%     end
%     data{1} = setsubfield(data{1}, dimtok{catdim}, tmp);
%     if isfield(data{1}, 'inside'),
%       data{1} = setsubfield(data{1}, 'inside',  tmpinside);
%       data{1} = setsubfield(data{1}, 'outside', setdiff(1:size(data{1}.pos,1)', tmpinside));
%     end
%     
%     %FIXME think about this
%     tryfields = {'dof'};
%   else
%     % no such field as {'label','time','freq','pos'} has to be concatenated
%     sortflag  = 0;
%     tryfields = {'cumsumcnt','cumtapcnt','trialinfo'};
%   end
%   
%   % concatenate the relevant descriptive fields in the data-structure (continued)
%   for k = 1:length(tryfields)
%     try
%       for m = 1:length(data)
%         if m==1,
%           tmpfield = data{m}.(tryfields{k});
%         else
%           tmpfield = [tmpfield; data{m}.(tryfields{k})];
%         end
%       end
%       data{1}.(tryfields{k}) = tmpfield;
%     catch
%     end
%   end
%   % FIXME handle inside in previous loop
%   
%   % FIXME this is ugly: solve it
%   %if issource || isvolume,
%   %  data{1}.dim(catdim) = max(size(tmp));
%   %end
%   
%   % sort concatenated data FIXME this is also ugly and depends on tmp
%   % FIXME if functional data in cell-array no sorting takes place
%   if sortflag && ~iscell(tmp) && ~iscell(data{1}.(param{1})),
%     [srt, ind] = sort(tmp, 2);
%     data{1}.(dimtok{catdim}) = tmp(ind);
%     for k = 1:length(param)
%       try
%         tmp     = data{1}.(param{k});
%       catch
%         continue;
%       end
%       tmp     = permute(tmp, [catdim setdiff(1:length(size(tmp)), catdim)]);
%       tmp     = ipermute(tmp(ind,:,:,:,:), [catdim setdiff(1:length(size(tmp)), catdim)]);
%       data{1}.(param{k}) = tmp;
%     end
%   elseif exist('tmp', 'var') && iscell(tmp)
%     %in this case (ugly!) tmp is probably a cell-array containing functional data
%   end
%   % remove unspecified parameters
%   if ~issource,
%     %rmparam = setdiff(parameterselection('all',data{1}),[param 'pos' 'inside' 'outside' 'freq' 'time']);
%     rmparam = {};
%   else
%     rmparam = setdiff(fieldnames(data{1}), [param(:)' paramdimord(:)' 'pos' 'inside' 'outside' 'dim' 'cfg' 'vol' 'cumtapcnt' 'orilabel' 'time' 'freq']);
%   end
%   for k = 1:length(rmparam)
%     data{1} = rmfield(data{1}, rmparam{k});
%   end
%   
%   % keep the first structure only
%   data        = data{1};
%   dimord      = dimord{1};
%   if ~issource,
%     data.dimord = dimord;
%   else
%     data.([param{1},'dimord']) = dimord;
%   end
%   if isfield(data, 'dim'),
%     data.dim    = dim;
%     %data.dim = size(data.(param{1}));
%   elseif isfield(data, 'dim')
%     data     = rmfield(data, 'dim'); %source data should not contain a dim
%     %FIXME this should be handled by ft_checkdata once the source structure is
%     %unequivocally defined
%   end
%   
% elseif length(data)>1 && israw
%   error('concatenation of several raw data-structures is done by ''ft_appenddata''');
% else
%   % nothing to do
%   data   = data{1};
%   dimord = dimord{1};
% end

% add version information to the configuration
cfg.version.name = mfilename('fullpath');
cfg.version.id = '$Id$';

% add information about the Matlab version used to the configuration
cfg.callinfo.matlab = version();
  
% add information about the function call to the configuration
cfg.callinfo.proctime = toc(ftFuncTimer);
cfg.callinfo.calltime = ftFuncClock;
cfg.callinfo.user = getusername();

% remember the configuration details of the input data
cfg.previous = cell(1,length(varargin));
for i=1:numel(varargin)
  try, cfg.previous{i} = varargin{i}.cfg; end
end

% remember the exact configuration details in the output
freq.cfg = cfg;

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'freq', freq); % use the variable name "data" in the output file
end

% subfunction to check uniqueness of freq bins
function [boolval, faxis] = checkfreq(varargin)

Ndata = numel(varargin);
Nfreq = zeros(1,Ndata);
faxis = zeros(1,0);
for i=1:Ndata
  Nfreq(i) = numel(varargin{i}.freq);
  faxis    = [faxis;varargin{i}.freq(:)];
end
boolval = numel(unique(faxis))==numel(faxis) && ~all(isnan(faxis));
% the second condition is included when the freq is set to dummy nan

% subfunction to check uniqueness of time bins
function [boolval, taxis] = checktime(varargin)

Ndata = numel(varargin);
Ntime = zeros(1,Ndata);
taxis = zeros(1,0);
for i=1:Ndata
  Ntime(i) = numel(varargin{i}.time);
  taxis    = [taxis;varargin{i}.time(:)];
end
boolval = numel(unique(taxis))==numel(taxis) && ~all(isnan(taxis));
% the second condition is included when the time is set to dummy nan

% subfunction to check uniqueness of channel labels
function [boolval, list] = checkchan(varargin)

Ndata = numel(varargin);
Nchan = zeros(1,Ndata);
list  = cell(0,1);
for i=1:Ndata
  Nchan(i) = numel(varargin{i}.label);
  list     = [list;varargin{i}.label(:)];
end
boolval = numel(unique(list))==numel(list);

