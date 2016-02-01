function [params, s_new, avg] = denoise_avg(params, s, state)
% DSS denoising function: Quasiperiodic averaging
%   [params, s_new] = denoise_avg(params, s, state)
%     params                Function specific modifiable parameters
%     params.tr             Trigger indices
%     params.fs             sampling frequency (if not specified
%                           params.begin and params.end assumed to be
%                           indices)  
%     params.begin          how many ms after the trigger the ON state
%                           begins 
%     params.end            how many ms after the trigger the ON state
%                           ends 
%     state                 DSS algorithm state
%     s                     Source signal estimate, matrix of row vector
%                           signals 
%     s_new                 Denoised signal estimate

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id$

if nargin<3 | ~isstruct(state)
    params.name = 'Quasiperiodic averaging with known triggers';
    params.description = '';
    params.param = {'fs','tr','tr_inds','begin','end'};
    params.param_value ={[], [], [], [], []};
    params.param_type = {'scalar','vector','vector','scalar','scalar'};
    params.param_desc = {'sampling frequency','trigger indices','used trigger indices','beginning of the ON mask','end of the ON mask'};
    params.approach = {'pca','defl','symm'};
    params.alpha = {};
    params.beta = {'beta_global'};
    return;
end

if isfield(params, 'fs')
  % sampling frequency available, converting begin and end to ms
  tr_begin = params.tr_begin/1000*params.fs;
  tr_end = params.tr_end/1000*params.fs;
else
  % not available, using as direct indices
  tr_begin = params.tr_begin;
  tr_end = params.tr_end;
end

if length(params.tr)==length(state.X)
  % trigger not as indices but as signal
  error('trigger indices must be extracted from the trigger signal first (automatic extraction not implemented yet)');
end

if ~isfield(params,'tr_inds')
  tr_inds = 1 : length(params.tr);
else
  tr_inds = params.tr_inds;
end

s_new = zeros(size(s));
avg = zeros(size(s,1),tr_end-tr_begin+1);

if (params.tr(end)+tr_end > size(s,2))
  dss_message(state, 3, sprintf('  Last response ends outside the data range. Dropping...\n'));
  tr_inds = tr_inds(1:end-1);
end

% calculating the average
for i = 1 : length(tr_inds)
  avg = avg + s(:,(params.tr(tr_inds(i))+tr_begin):(params.tr(tr_inds(i))+tr_end));
end
avg = avg  / length(params.tr);
% reconstructing the signals
for i = 1 : length(tr_inds)
  inds = (params.tr(tr_inds(i))+tr_begin):(params.tr(tr_inds(i))+tr_end);
  s_new(:,inds)=avg;
end
