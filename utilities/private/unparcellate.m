function varargout = unparcellate(data, parcellation, parameter, parcelparam, varargin)

% UNPARCELLATE performs the reverse of a parcellation, by assigigning each
% parcel's activation to the vertices that contributed to that parcel.
%
% Use as
%
%   fun = unparcellate(data, parcellation, parameter, parcelparam, varargin)
%
% Required inputs:
%
%   data          = structure (or matrix) containing the parcellated functional data
%   parcellation  = structure describing the parcellation, i.e. the parcel
%                   membership for each of the vertices
%   parameter     = string (or cell-array with labels) that specifies the
%                   parameter to be used (if data is a structure) or how to
%                   interpret the rows in the data matrix (if data is a matrix)
%
% Additional inputs are key-value pairs and pertain to bivariate data with
% a 'labelcmb' specified in the input argument 'parameter'.
%
%   avgoverref     = 'yes' (or 'no')
%   directionality = 'both' (or 'inflow'/'outflow')
%
% Outputs:
%   fun = matrix Nvertices x size(data.(parameter),2) (or Nvertices x
%           size(data,2), containing the unparcellated data
%
%   If the input was bivariate data with a labelcmb, an optional second
%   output argument gives a list of the reference parcels.

if isstruct(data)
  tmp = getsubfield(data, parameter);
elseif size(parameter,2)==1
  tmp = data;
  clear data;
  data.label = parameter;
elseif size(parameter,2)==2
  % this contains labelcmb e.g. pairwise granger
  tmp = data;
  clear data;
  data.labelcmb = parameter;
end

if isfield(data, 'label')
  % univariate data or 'chan_chan_xxx' bivariate data
  fun = zeros(size(parcellation.pos, 1), size(tmp, 2))+nan;
  for k = 1:numel(data.label)
    sel = match_str(parcellation.([parcelparam,'label']), data.label{k});
    if ~isempty(sel)
      sel = parcellation.(parcelparam)==sel;
      sel = sel(:); % ensure it is a column vector
      fun(sel,:) = repmat(tmp(k,:), [sum(sel) 1]);
    end
  end
  
  varargout{1} = fun;
  
elseif isfield(data, 'labelcmb')
  % bivariate data
  
  avgoverref     = istrue(ft_getopt(varargin, 'avgoverref', 1));
  directionality = ft_getopt(varargin, 'directionality', 'both');
  
  labelcmb = data.labelcmb;
  
  % remove the square brackets, if present
  for k = 1:numel(labelcmb)
    tok = tokenize(labelcmb{k}, '[');
    labelcmb{k} = tok{1};
  end
  
  % identify which parcel(s) was/were used as seed parcel
  uniq_label = unique(labelcmb(:));
  uniq_count = zeros(size(uniq_label));
  for k = 1:numel(uniq_label)
    uniq_count(k) = sum(strcmp(labelcmb(:), uniq_label{k}));
  end
  
  % assume the labelcmb to be well behaved (i.e. have a full parcellation per seed
  % either or not excluding the seed regions themselves), and to contain quadruplets
  % per combination (2 auto combinations and 2 directional interaction
  % terms)
  nref = sum(uniq_count==max(uniq_count));
  ref  = uniq_label(uniq_count==max(uniq_count));
  
  if strcmp(directionality, 'both')
    % recurse and combine
    fun1 = unparcellate(tmp,parcellation,data.labelcmb,parcelparam,'avgoverref',avgoverref,'directionality','inflow');
    fun2 = unparcellate(tmp,parcellation,data.labelcmb,parcelparam,'avgoverref',avgoverref,'directionality','outflow');
    fun  = [fun1;fun2]; clear fun1 fun2;
    
  else
    
    for k = 1:nref
      if strcmp(directionality, 'outflow')
        sel = strcmp(labelcmb(:,1),ref{k})&~strcmp(labelcmb(:,2),ref{k});
        sel2 = 2;
      elseif strcmp(directionality, 'inflow')
        sel = strcmp(labelcmb(:,2),ref{k})&~strcmp(labelcmb(:,1),ref{k});
        sel2 = 1;
      end
      
      tmpdat   = tmp(sel,:);
      tmplabel = labelcmb(sel,sel2);
      
      % recurse into unparcellate
      tmpfun   = unparcellate(tmpdat, parcellation, tmplabel, parcelparam);
      
      if avgoverref
        if k==1
          fun = tmpfun;
        else
          fun = fun+tmpfun;
        end
      else
        if k==1
          fun = tmpfun;
        else
          fun = cat(1,fun, tmpfun);
        end
      end
      
    end
    
    if avgoverref
      fun = fun./nref;
    end
  end
  
  varargout{1} = fun;
  varargout{2} = ref;
  
end
