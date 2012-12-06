function comprefcorr = ft_componentclassification(cfg, comp, refdata)

% FT_COMPONENTCLASSIFICATION performs a classification of the spatiotemporal
% components
%
% Use as
%   compclass = ft_componentclassification(cfg, comp) 
% where comp is the output of FT_COMPONENTANALYSIS and cfg is a       
% configuration structure that should contain 
%
%   cfg.option1    = value, explain the value here (default = something)
%   cfg.option2    = value, describe the value here and if needed
%                    continue here to allow automatic parsing of the help
%
% The configuration can optionally contain
%   cfg.option3   = value, explain it here (default is automatic)
%
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following options:
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_COMPONENTANALYSIS, FT_TOPOPLOTIC

% Copyright (C) 2011, Jan-Mathijs Schoffelen
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
ftFuncMem   = memtic();

% ensure that the input data is valiud for this function, this will also do 
% backward-compatibility conversions of old data that for example was 
% read from an old *.mat file
comp = ft_checkdata(comp, 'datatype', 'comp', 'feedback', 'yes');
if nargin>2
  refdata = ft_checkdata(refdata, 'datatype', 'raw', 'feedback', 'yes');
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'trackconfig', 'on');

% ensure that the required options are present
cfg = ft_checkconfig(cfg, 'required', 'method');

% ensure that the options are valid
% cfg = ft_checkopt(cfg, 'vartrllen', 'double', {0, 1, 2});
% cfg = ft_checkopt(cfg, 'method', 'char', {'mtm', 'convol'});

% get the options
cfg.inputfile  = ft_getopt(cfg, 'inputfile', '');
cfg.outputfile = ft_getopt(cfg, 'outputfile', '');
method         = ft_getopt(cfg, 'method'); % there is no default

hasdata = (nargin>1);
if ~isempty(cfg.inputfile)
  % the input data should be read from file
  if hasdata
    error('cfg.inputfile should not be used in conjunction with giving input data to this function');
  else
    comp = loadvar(cfg.inputfile, 'data');
  end
end

if strcmp(method, 'template_timeseries') && nargin<=2
  error('for the method ''template_timeseries'' the input to this function should contain the reference time series as a separate input');
end
  
switch method
  case 'template_spectrum'
    cfg                   = ft_checkconfig(cfg,          'required', 'template');
    cfg.template          = ft_checkconfig(cfg.template, 'required', 'spectrum');
    cfg.template.spectrum = ft_checkconfig(cfg.template.spectrum, 'required', {'freq' 'powspctrm'});
    
    % template
    template = cfg.template.spectrum;
    
    % do spectral analysis
    tmpcfg        = [];
    tmpcfg.method = 'mtmfft';
    tmpcfg.taper  = 'hanning';
    tmpcfg.output = 'pow';
    freq          = ft_freqanalysis(tmpcfg, comp);
    
    % do some interpolation here if needed
    % FIXME
    if ~all(template.freq==freq.freq)
    end
    
    % regress
    
    
    x=1;
    
  case 'template_timeseries'
    % check whether the inputs are compatible
    ok = true;
    for k = 1:numel(comp.trial)
      ok = size(comp.trial{k},2) == size(refdata.trial{k},2);
      if ~ok
        error('the input data structures are incompatible because of different trial lengths');
        break;
      end
    end
    Ncomp = numel(comp.label);
    comp  = ft_appenddata([], comp, refdata);
    
    % compute covariance 
    comprefcov = cellcov(comp.trial, [], 2);
    
    % compute correlation
    comprefcorr = comprefcov(1:Ncomp,(Ncomp+1):end)./sqrt(diag(comprefcov(1:Ncomp,1:Ncomp))*diag(comprefcov((Ncomp+1):end,(Ncomp+1):end))');
    
  case 'template_topography'
    error('unknown method of classification');    
  case '1/f'
    error('unknown method of classification');    
  case 'kurtosis'
    error('unknown method of classification');    
  case 'whiteness'
    error('unknown method of classification');    
  case 'something else'
    error('unknown method of classification');    
  otherwise
    error('unknown method of classification');    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% accessing this field here is needed for the configuration tracking
% by accessing it once, it will not be removed from the output cfg
cfg.outputfile;

% get the output cfg
cfg = ft_checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');

% add the version details of this function call to the configuration
cfg.version.name = mfilename('fullpath'); % this is helpful for debugging
cfg.version.id   = '$Id$'; % this will be auto-updated by the revision control system

% add information about the Matlab version used to the configuration
cfg.callinfo.matlab = version();
  
% add information about the function call to the configuration
cfg.callinfo.proctime = toc(ftFuncTimer);
cfg.callinfo.procmem  = memtoc(ftFuncMem);
cfg.callinfo.calltime = ftFuncClock;
cfg.callinfo.user = getusername(); % this is helpful for debugging
fprintf('the call to "%s" took %d seconds and an estimated %d MB\n', mfilename, round(cfg.callinfo.proctime), round(cfg.callinfo.procmem/(1024*1024)));

if hasdata && isfield(comp, 'cfg')
  % remember the configuration details of the input data
  cfg.previous = comp.cfg;
end

% remember the exact configuration details in the output
dataout.cfg = cfg;

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'data', dataout); % use the variable name "data" in the output file
end

%-----cellcov
function [c] = cellcov(x, y, dim, flag)

% [C] = CELLCOV(X, DIM) computes the covariance, across all cells in x along 
% the dimension dim. When there are three inputs, covariance is computed between
% all cells in x and y
% 
% X (and Y) should be linear cell-array(s) of matrices for which the size in at 
% least one of the dimensions should be the same for all cells 

if nargin==2,
  flag = 1;
  dim  = y;
  y    = [];
elseif nargin==3,
  flag = 1;
end

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1),
  error('incorrect input for cellmean');
end

if nargin==1,
  scx1 = cellfun('size', x, 1);
  scx2 = cellfun('size', x, 2);
  if     all(scx2==scx2(1)), dim = 2; %let second dimension prevail
  elseif all(scx1==scx1(1)), dim = 1;
  else   error('no dimension to compute covariance for');
  end
end

if flag,
  mx   = cellmean(x, 2);
  x    = cellvecadd(x, -mx);
  if ~isempty(y),
    my = cellmean(y, 2);
    y  = cellvecadd(y, -my);
  end
end

nx   = max(nx);
nsmp = cellfun('size', x, dim);
if isempty(y), 
  csmp = cellfun(@covc, x, repmat({dim},1,nx), 'UniformOutput', 0);
else
  csmp = cellfun(@covc, x, y, repmat({dim},1,nx), 'UniformOutput', 0);
end
nc   = size(csmp{1});
c    = sum(reshape(cell2mat(csmp), [nc(1) nc(2) nx]), 3)./sum(nsmp); 

function [c] = covc(x, y, dim)

if nargin==2,
  dim = y;
  y   = x;
end

if dim==1,
  c = x'*y;
elseif dim==2,
  c = x*y';
end

%-----cellmean
function [m] = cellmean(x, dim)

% [M] = CELLMEAN(X, DIM) computes the mean, across all cells in x along 
% the dimension dim.
% 
% X should be an linear cell-array of matrices for which the size in at 
% least one of the dimensions should be the same for all cells 

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1),
  error('incorrect input for cellmean');
end

if nargin==1,
  scx1 = cellfun('size', x, 1);
  scx2 = cellfun('size', x, 2);
  if     all(scx2==scx2(1)), dim = 2; %let second dimension prevail
  elseif all(scx1==scx1(1)), dim = 1;
  else   error('no dimension to compute mean for');
  end
end

nx   = max(nx);
nsmp = cellfun('size', x, dim);
ssmp = cellfun(@sum,   x, repmat({dim},1,nx), 'UniformOutput', 0);
m    = sum(cell2mat(ssmp), dim)./sum(nsmp);  

%-----cellstd
function [sd] = cellstd(x, dim, flag)

% [M] = CELLSTD(X, DIM, FLAG) computes the standard deviation, across all cells in x along 
% the dimension dim, normalising by the total number of samples 
% 
% X should be an linear cell-array of matrices for which the size in at 
% least one of the dimensions should be the same for all cells. If flag==1, the mean will
% be subtracted first (default behaviour, but to save time on already demeaned data, it
% can be set to 0).

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1),
  error('incorrect input for cellstd');
end

if nargin<2,
  scx1 = cellfun('size', x, 1);
  scx2 = cellfun('size', x, 2);
  if     all(scx2==scx2(1)), dim = 2; %let second dimension prevail
  elseif all(scx1==scx1(1)), dim = 1;
  else   error('no dimension to compute mean for');
  end
elseif nargin==2,
  flag = 1;
end

if flag,
  m    = cellmean(x, dim);
  x    = cellvecadd(x, -m);
end

nx   = max(nx);
nsmp = cellfun('size', x, dim);
ssmp = cellfun(@sumsq,   x, repmat({dim},1,nx), 'UniformOutput', 0);
sd   = sqrt(sum(cell2mat(ssmp), dim)./sum(nsmp));  

function [s] = sumsq(x, dim)

s = sum(x.^2, dim);

%-----cellvecadd
function [y] = cellvecadd(x, v)

% [Y]= CELLVECADD(X, V) - add vector to all rows or columns of each matrix 
% in cell-array X

% check once and for all to save time
persistent bsxfun_exists;
if isempty(bsxfun_exists); 
    bsxfun_exists=(exist('bsxfun')==5); 
    if ~bsxfun_exists; 
        error('bsxfun not found.');
    end
end

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1),
  error('incorrect input for cellmean');
end

if ~iscell(v),
  v = repmat({v}, nx);
end

sx1 = cellfun('size', x, 1);
sx2 = cellfun('size', x, 2);
sv1 = cellfun('size', v, 1);
sv2 = cellfun('size', v, 2);
if all(sx1==sv1) && all(sv2==1),    
  dim = mat2cell([ones(length(sx2),1) sx2(:)]', repmat(2,nx(1),1), repmat(1,nx(2),1)); 
elseif all(sx2==sv2) && all(sv1==1),
  dim = mat2cell([sx1(:) ones(length(sx1),1)]', repmat(2,nx(1),1), repmat(1,nx(2),1));
elseif all(sv1==1) && all(sv2==1),
  dim = mat2cell([sx1(:) sx2(:)]'', nx(1), nx(2));
else   error('inconsistent input');
end  

y  = cellfun(@bsxfun, repmat({@plus}, nx), x, v, 'UniformOutput', 0);
%y = cellfun(@vplus, x, v, dim, 'UniformOutput', 0);

function y = vplus(x, v, dim)

y = x + repmat(v, dim);

%-----cellvecmult
function [y] = cellvecmult(x, v)

% [Y]= CELLVECMULT(X, V) - multiply vectors in cell-array V
% to all rows or columns of each matrix in cell-array X
% V can be a vector or a cell-array of vectors

% check once and for all to save time
persistent bsxfun_exists;
if isempty(bsxfun_exists); 
    bsxfun_exists=(exist('bsxfun')==5); 
    if ~bsxfun_exists; 
        error('bsxfun not found.');
    end
end

nx = size(x);
if ~iscell(x) || length(nx)>2 || all(nx>1),
  error('incorrect input for cellmean');
end

if ~iscell(v),
  v = repmat({v}, nx);
end

sx1 = cellfun('size', x, 1);
sx2 = cellfun('size', x, 2);
sv1 = cellfun('size', v, 1);
sv2 = cellfun('size', v, 2);
if all(sx1==sv1) && all(sv2==1),    
elseif all(sx2==sv2) && all(sv1==1),
elseif all(sv1==1) && all(sv2==1),
else   error('inconsistent input');
end  

y  = cellfun(@bsxfun, repmat({@times}, nx), x, v, 'UniformOutput', 0);
