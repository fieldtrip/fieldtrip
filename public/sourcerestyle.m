function [output] = sourcerestyle(input, type)

% SOURCERESTYLE converts old style source structures into new style source structures and the
% other way around
%
% Use as:
%   sourcerestyle(input, type)
%    where input is a source structure,
%
% Typically, old style source structures contain
%   avg.XXX or trial.XXX fields
%
% The ne wstyle source structure contains:
%   source.pos
%   source.dim (optional, if the list of positions describes a 3D volume
%   source.XXX the old style subfields in avg/trial
%   source.XXXdimord string how to interpret the respective XXX field:
%     e.g. source.leadfield = cell(1,Npos), source.leadfielddimord = '{pos}_chan_ori'
%          source.mom       = cell(1,Npos), source.momdimord       = '{pos}_ori_rpttap'

if nargin==1,
  type = 'old';
end

fnames = fieldnames(input);
tmp    = cell2mat(strfind(fnames, 'dimord')); %get dimord like fields
if any(tmp>1),
  current = 'new';
elseif any(tmp==1),
  %don't know what to do yet data is JM's own invention
else
  current = 'old';
end

if strcmp(current, type),
  %do nothing
  return
elseif strcmp(current, 'old') && strcmp(type, 'new'),
  %go from old to new

  if isfield(input, 'avg'),
    stuff  = getfield(input, 'avg');
    output = rmfield(input,  'avg');
  elseif isfield(input, 'trial'),
    stuff  = getfield(input, 'trial');
    output = rmfield(input,  'trial');
  else
    %this could occur later in the pipeline, e.g. when doing group statistics using individual subject
    %descriptive statistics
    error('the input does not contain an avg or trial field');
  end
  
  %remove and rename the specified fields if present
  removefields = {'xgrid';'ygrid';'zgrid';'method'};
  renamefields = {'frequency' 'freq'};
  fnames       = fieldnames(output);
  for k = 1:numel(fnames)
    ix = strmatch(fnames{k}, removefields);
    if ~isempty(ix),
      output = rmfield(output, fnames{k});
    end
    ix = strmatch(fnames{k}, renamefields(:,1));
    if ~isempty(ix),
      output = setfield(output, renamefields{ix,2}, ...
                        getfield(output, renamefields{ix,1}));
      output = rmfield(output, fnames{k});
    end
  end

  %put the stuff originally in avg or trial one level up in the structure
  fnames       = fieldnames(stuff);
  npos         = size(input.pos,1);
  for k = 1:numel(fnames)
    tmp = getfield(stuff, fnames{k});
    siz = size(tmp);
    if siz(1) ~= npos && siz(2) ==npos,
      tmp = transpose(tmp);
    end
    output    = setfield(output, fnames{k}, tmp);
    newdimord = createdimord(output, fnames{k}); 
    output    = setfield(output, [fnames{k},'dimord'], newdimord);
  end

elseif strcmp(current, 'new') && strcmp(type, 'old')
  %go from new to old
  error('not implemented yet');
end

function [dimord] = createdimord(output, fname);

tmp = getfield(output, fname);

dimord = '';
dimnum = 1;
if iscell(tmp) && size(output.pos,1)==size(tmp,dimnum)
  dimord = [dimord,'{pos}'];
  dimnum = dimnum + 1;
elseif ~iscell(tmp) && size(output.pos,1)==size(tmp,dimnum)
  dimord = [dimord,'pos'];
  dimnum = dimnum + 1;
end

switch fname
  case 'csd'
    dimord = [dimord,'_ori_ori'];
  case 'filter'
    dimord = [dimord,'_ori_chan']; 
  case 'leadfield'
    dimord = [dimord,'_chan_ori'];
  case 'mom'
  case 'nai'
    if isfield(output, 'freq') && numel(output.freq)==size(tmp,dimnum)
      dimord = [dimord,'_freq'];
    end
  case 'noise'
    if isfield(output, 'freq') && numel(output.freq)==size(tmp,dimnum)
      dimord = [dimord,'_freq'];
    end
  case 'noisecsd'
    dimord = [dimord,'_ori_ori'];
  case 'ori'
    %this field is equivalent to a pos-field
    %FIXME should this be matricized (size is not that big)
  case 'pow'
    if isfield(output, 'freq') && numel(output.freq)==size(tmp,dimnum)
      dimord = [dimord,'_freq'];
    end
  otherwise
    error(sprintf('unknown fieldname %s', fname));
end
