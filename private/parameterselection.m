function [select] = parameterselection(param, data);

% PARAMETERSELECTION selects the parameters that are present as a volume in the data
% add that have a dimension that is compatible with the specified dimensions of the
% volume, i.e. either as a vector or as a 3D volume.
%
% Use as
%   [select] = parameterselection(param, data)
% where
%   param    cell-array, or single string, can be 'all'
%   data     structure with anatomical or functional data
%   select   returns the selected parameters as a cell-array

% Copyright (C) 2005-2008, Robert oostenveld
%
% $Log: parameterselection.m,v $
% Revision 1.14  2008/07/30 07:41:39  roboos
% also detect parameters in source.trial
%
% Revision 1.13  2008/07/25 09:05:02  roboos
% not only check the parameter dimensions against source.dim, but optionally also against source.pos (required if it is an irregular grid)
%
% Revision 1.12  2007/07/31 08:36:27  jansch
% added support for volume data with dimensionality > 3
%
% Revision 1.11  2006/10/02 13:54:14  roboos
% replaced predefined list for 'all' by automatic generated list
%
% Revision 1.10  2006/05/16 10:41:13  roboos
% removed trialA/trialB elements, since typically they cannot be dealt with automatically
%
% Revision 1.9  2006/03/30 07:41:15  roboos
% added fields df and dof in case of selecting 'all'
%
% Revision 1.8  2006/03/29 08:25:57  roboos
% Ensure that it is a cell array, also when empty. Remove empty fields (thanks to Vladimir).
%
% Revision 1.7  2006/02/13 07:45:16  jansch
% added paramters such as lbex, leadfield, trialX.xxx so that this function
% can also be used within source2sparse, and source2full
%
% Revision 1.6  2006/02/09 09:43:20  roboos
% added 'leadfield' as a known volume
%
% Revision 1.5  2006/01/05 13:29:36  roboos
% always return a cell-array (this function used to return a string if only
% one parameter was selected and a cell-array if there were multiple)
%
% Revision 1.4  2005/09/29 00:37:57  roboos
% if 'all' is specified, do not replace complete list, but only the 'all' element
% prevent double occurences in parameter list
%
% Revision 1.3  2005/09/12 09:50:46  jansch
% added ref as a parameter
%
% Revision 1.2  2005/08/19 17:26:27  roboos
% added support for 'pow' ->  will be replaced by 'avg.pow' (similar for other volumes)
% added help documentation
%

if ischar(param)
  param = {param};   % it should be a cell-array
elseif isempty(param)
  param = {};        % even being empty, it should be a cell-array
end

sel = find(strcmp(param, 'all'));
if ~isempty(sel)
  % the old default was a list of all known volume parameters
  % the new default is to try all fields present in the data
  allparam = fieldnames(data);
  % fields can be nested in source.avg
  if isfield(data, 'avg')
    tmp = fieldnames(data.avg);
    for i=1:length(tmp)
      tmp{i} = ['avg.' tmp{i}];
    end
    allparam = cat(1, allparam, tmp);
  end
  % fields can be nested in source.trial
  if isfield(data, 'trial')
    tmp = fieldnames(data.trial);
    for i=1:length(tmp)
      tmp{i} = ['trial.' tmp{i}];
    end
    allparam = cat(1, allparam, tmp);
  end
  param(sel) = [];                          % remove the 'all'
  param      = [param(:)' allparam(:)'];    % add the list of all possible parameters, these will be tested later
else
  % check all specified parameters and give support for some parameters like 'pow' and 'coh'
  % which most often will indicate 'avg.pow' and 'avg.coh'
  for i=1:length(param)
    if ~issubfield(data, param{i}) && issubfield(data, ['avg.' param{i}])
      % replace the parameter xxx by avg.xxx
      param{i} = ['avg.' param{i}];
    end
  end
end

% remove empty fields
param(find(cellfun('isempty', param))) = [];

% ensure that there are no double entries
param = unique(param);

select = {};
for i=1:length(param)
  if issubfield(data, param{i})
    % the field is present, check whether the dimension is correct
    dim = size(getsubfield(data, param{i}));
    if isfield(data, 'dim') && isequal(dim(:), data.dim(:))
      select{end+1} = param{i}; 
    elseif isfield(data, 'dim') && prod(dim)==prod(data.dim)
      select{end+1} = param{i}; 
    elseif isfield(data, 'pos') && prod(dim)==size(data.pos, 1)
      select{end+1} = param{i}; 
    end
  end
end

