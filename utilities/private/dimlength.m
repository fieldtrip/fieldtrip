function [n, fn] = dimlength(data, seldim, fld)

% DIMLENGTH(DATA, SELDIM, FLD) is a helper function to obtain n, the number
% of elements along dimension seldim from the appropriate field from the
% input data containing functional data.
%
% Use als
%   [n, fn] = dimlength(data, seldim, fld)
%
% It can be called with one input argument only, in which case it will
% output two cell arrays containing the size of the functional fields,
% based on the XXXdimord, and the corresponding XXXdimord fields.
%
% When the data contains a single dimord field (everything except source
% data), the cell-arrays in the output only contain one element.
%
% See also FIXSOURCE, CREATEDIMORD

% Copyright (C) 2010, Jan-Mathijs Schoffelen
%
% $Id$

if nargin<3
  fld = 'dimord';
end

if nargin==1
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % call recursively to get the dimensionality of all fields XXX which are
  % accompanied by a XXXdimord field or the general dimord field
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % get all fields of the data structure
  fn = fieldnames(data);
  % get all dimord like fields
  fn    = fn(~cellfun('isempty', strfind(fn, 'dimord')));
  for k = 1:numel(fn)
    fndimord{k} = data.(fn{k});
  end
  
  for k = 1:numel(fndimord)
    dimtok = tokenize(fndimord{k}, '_');
    ndim   = numel(dimtok);
    n{k,1} = zeros(1, ndim);
    for i = 1:ndim
      n{k}(i) = dimlength(data, dimtok{i}, fn{k});
    end
  end
  return
  
else
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % do the actual work
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  switch seldim
    case 'rpt'
      if numel(fld)>6 && isfield(data, fld(1:end-6)),
        % source level data
        dimtok = tokenize(data.(fld), '_');
        tmp    = data.(fld(1:end-6));
        % remove the first cell-dimension
        if iscell(tmp)
          if isfield(data, 'inside'),
            tmp = tmp{data.inside(1)};
          else
            tmp = tmp{1};
          end
          dimtok = dimtok(2:end);
        end
        ix = find(~cellfun('isempty', strfind(dimtok, seldim)));
        n  = size(tmp, ix);
        
      elseif strcmp(data.(fld), 'rpt_pos')
        %FIXME HACK to be fixed
        x = setdiff(fld(data),'inside');
        for k = 1:length(x)
          dims = size(getsubfield(data,x{k}));
          if dims(2)==size(data.pos,1) && numel(dims)==2,
            n = dims(1);
            return
          end
        end
        
      elseif strcmp(data.(fld), 'rpt_pos_freq'),
        %FIXME HACK to be fixed
        x = fld(data);
        for k = 1:length(x)
          dims = size(getsubfield(data,x{k}));
          if dims(2)==size(data.pos,1) && (numel(dims)==2 || dims(3)==length(data.freq)),
            n = dims(1);
            return
          end
        end
        
      elseif strcmp(data.(fld), 'rpt_pos_time'),
        %FIXME HACK to be fixed
        x = fld(data);
        for k = 1:length(x)
          dims = size(getsubfield(data,x{k}));
          if dims(2)==size(data.pos,1) && (numel(dims)==2 || dims(3)==length(data.time)),
            n = dims(1);
            return
          end
        end
        
      elseif strcmp(data.(fld)(1:4), 'rpt_')
        n  = [];
        % generic solution for XXXspctrm
        fnames = fieldnames(data);
        tmp    = find(~cellfun('isempty', strfind(fnames, 'spctrm')));
        for k = 1:numel(tmp)
          n = [n size(data.(fnames{tmp(k)}), 1)];
        end
        % some other possibilities
        if isfield(data, 'cov'),           n = [n size(data.cov,           1)]; end
        %if isfield(data, 'crsspctrm'),     n = [n size(data.crsspctrm,     1)]; end
        %if isfield(data, 'powcovspctrm'),  n = [n size(data.powcovspctrm,  1)]; end
        %if isfield(data, 'powspctrm'),     n = [n size(data.powspctrm,     1)]; end
        if isfield(data, 'trial'),         n = [n size(data.trial,         1)]; end
        %if isfield(data, 'fourierspctrm'), n = [n size(data.fourierspctrm, 1)]; end
        if isfield(data, 'individual'),    n = [n size(data.individual,    1)]; end
        if isfield(data, 'stat'),          n = [n size(data.stat,          1)]; end
        if ~all(n==n(1))
          ft_error('inconsistent number of repetitions for dim "%s"', seldim);
        end
        n = n(1);
        
      else
        %error('cannot determine number of repetitions for dim "%s"', seldim);
        n = nan;
      end
      
    case 'rpttap'
      if numel(fld)>6 && isfield(data, fld(1:end-6)),
        dimtok = tokenize(data.(fld), '_');
        tmp    = data.(fld(1:end-6));
        % remove the first cell-dimension
        if iscell(tmp)
          if isfield(data, 'inside'),
            tmp = tmp{data.inside(1)};
          else
            tmp = tmp{1};
          end
          dimtok = dimtok(2:end);
        end
        ix = find(~cellfun('isempty', strfind(dimtok, seldim)));
        n  = size(tmp, ix);
        
      elseif strcmp(data.(fld)(1:7), 'rpttap_')
        n  = [];
        if isfield(data, 'cov'),           n = [n size(data.cov,           1)]; end
        if isfield(data, 'crsspctrm'),     n = [n size(data.crsspctrm,     1)]; end
        if isfield(data, 'powcovspctrm'),  n = [n size(data.powcovspctrm,  1)]; end
        if isfield(data, 'powspctrm'),     n = [n size(data.powspctrm,     1)]; end
        if isfield(data, 'trial'),         n = [n size(data.trial,         1)]; end
        if isfield(data, 'fourierspctrm'), n = [n size(data.fourierspctrm, 1)]; end
        if ~all(n==n(1))
          ft_error('inconsistent number of repetitions for dim "%s"', seldim);
        end
        n = n(1);
        
      else
        %error('cannot determine number of repetitions for dim "%s"', seldim);
        n = nan;
      end
      
    case 'subj'
      n  = [];
      if isfield(data, 'cov'),           n = [n size(data.cov,           1)]; end
      if isfield(data, 'crsspctrm'),     n = [n size(data.crsspctrm,     1)]; end
      if isfield(data, 'fourierspctrm'), n = [n size(data.fourierspctrm, 1)]; end
      if isfield(data, 'individual'),    n = [n size(data.individual,    1)]; end
      if isfield(data, 'powcovspctrm'),  n = [n size(data.powcovspctrm,  1)]; end
      if isfield(data, 'powspctrm'),     n = [n size(data.powspctrm,     1)]; end
      if isfield(data, 'trial'),         n = [n size(data.trial,         1)]; end
      
      if ~all(n==n(1)), ft_error('inconsistent number of repetitions for dim "%s"', seldim); end
      n = n(1);
      
    case 'chan'
      if ~isfield(data, 'inside'),
        n = length(data.label);
      else
        %error('cannot determine number of repetitions for dim "%s"', seldim);
        n = nan; %FIXME discuss appending label to source-like data
      end
      
    case 'chancmb'
      if ~isfield(data, 'inside'),
        n = size(data.labelcmb, 1);        
      else
        %error('cannot determine number of repetitions for dim "%s"', seldim);
        n = nan; %FIXME discuss appending label to source-like data
      end
      
    case 'freq'
      n = length(data.freq);
      
    case 'time'
      n = length(data.time);
      
    case {'pos' '{pos}'}
      n = size(data.pos,1);
      
    case {'ori'}
      % this is comparable to rpt and to rpttap
      if numel(fld)>6 && isfield(data, fld(1:end-6)),
        dimtok = tokenize(data.(fld), '_');
        tmp    = data.(fld(1:end-6));
        if iscell(tmp)
          if isfield(data, 'inside'),
            tmp = tmp{data.inside(1)};
          else
            tmp = tmp{1};
          end
          dimtok = dimtok(2:end);
        end
        % it can be {pos}_ori_ori, in which case only the first orientation needs to be checked
        ix = find(~cellfun('isempty', strfind(dimtok, seldim)), 1, 'first');
        n  = size(tmp, ix);
      else
        %error('cannot determine number of repetitions for dim "%s"', seldim);
        n = nan;
      end
      
    otherwise
      ft_error('unsupported dim "%s"', seldim);
  end
  
end % if nargin
