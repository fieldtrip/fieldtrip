function [grandavg] = sourcegrandaverage(cfg, varargin);

% SOUREGRANDAVERAGE averages source reconstructions over either multiple
% subjects or conditions. It computes the average and variance for all
% known source parameters. The output can be used in SOURCESTATISTICS
% with the method 'parametric'.
% 
% Alternatively, it can construct an average for multiple input source
% reconstructions in two conditions after randomly reassigning the
% input data over the two conditions. The output then can be used in
% SOURCESTATISTICS with the method 'randomization' or 'randcluster'.
%
% The input source structures should be spatially alligned to each other
% and should have the same positions for the source grid.
% 
% Use as 
%  [grandavg] = SOUREGRANDAVERAGE(cfg, source1, source2, ...) 
%
% where the source structures are obtained from SOURCEANALYSIS or
% from VOLUMENORMALISE, and the configuration can contain the
% following fields:
%   cfg.parameter          = string, describing the functional data to be processed, e.g. 'pow', 'nai' or 'coh'
%   cfg.keepindividual     = 'no' or 'yes'  
%
% See also SOURCEANALYSIS, VOLUMENORMALISE, SOURCESTATISTICS

% Undocumented local options
%  You can also use SOURCEGRANDAVERAGE to compute averages after
% randomizing the assignment of the functional data over two conditions.
% The resulting output can then be used in a statistical test just like
% the randomized single-subject source reconstruction that results from 
% randomization in SOURCEANALYSIS. This involves the following options
%   cfg.randomization      = 'no' or 'yes'
%   cfg.permutation        = 'no' or 'yes'
%   cfg.numrandomization   = number, e.g. 500
%   cfg.numpermutation     = number, e.g. 500 or 'all'
%   cfg.c1                 = list with subjects belonging to condition 1 (or A)
%   cfg.c2                 = list with subjects belonging to condition 2 (or B)

% Copyright (C) 2005, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

fieldtripdefs

cfg = checkconfig(cfg, 'trackconfig', 'on');

% check if the input data is valid for this function
for i=1:length(varargin)
  varargin{i} = checkdata(varargin{i}, 'datatype', {'source', 'volume'}, 'feedback', 'no');
end

% set the defaults
if ~isfield(cfg, 'parameter'),      cfg.parameter = 'pow';     end
if ~isfield(cfg, 'keepindividual'), cfg.keepindividual = 'no'; end
if ~isfield(cfg, 'concatenate'),    cfg.concatenate    = 'no'; end
if ~isfield(cfg, 'randomization'),  cfg.randomization = 'no';  end
if ~isfield(cfg, 'permutation'),    cfg.permutation = 'no';    end
if ~isfield(cfg, 'c1'),             cfg.c1 = [];               end
if ~isfield(cfg, 'c2'),             cfg.c2 = [];               end

if strcmp(cfg.concatenate, 'yes') && strcmp(cfg.keepindividual, 'yes'),
  error('you can either specify cfg.keepindividual or cfg.concatenate to be yes');
end

Nsubject = length(varargin);
Nvoxel   = prod(varargin{1}.dim);
dat      = zeros(Nvoxel, Nsubject);
inside   = zeros(Nvoxel, Nsubject);

if isfield(varargin{1}, 'pos')
  % check that the source locations of each input source reconstruction are the same
  for i=2:Nsubject
    if size(varargin{i}.pos,1)~=size(varargin{1}.pos,1) || any(varargin{i}.pos(:)~=varargin{1}.pos(:))
      error('different grid locations in source reconstructions');
    end
  end
  grandavg.pos = varargin{1}.pos;
end

if isfield(varargin{1}, 'dim')
  % check that the dimensions of each input volume is the same
  for i=2:Nsubject
    if any(varargin{i}.dim(:)~=varargin{1}.dim(:))
      error('different dimensions of the source reconstructions');
    end
  end
  grandavg.dim = varargin{1}.dim;
end

if isfield(varargin{1}, 'xgrid') && isfield(varargin{1}, 'ygrid') && isfield(varargin{1}, 'zgrid')
  % check that the grid locations of each input volume are the same
  for i=2:Nsubject
    if length(varargin{i}.xgrid)~=length(varargin{1}.xgrid) || any(varargin{i}.xgrid~=varargin{1}.xgrid)
      error('different xgrid in source reconstructions');
    elseif length(varargin{i}.ygrid)~=length(varargin{1}.ygrid) || any(varargin{i}.ygrid~=varargin{1}.ygrid)
      error('different ygrid in source reconstructions');
    elseif length(varargin{i}.zgrid)~=length(varargin{1}.zgrid) || any(varargin{i}.zgrid~=varargin{1}.zgrid)
      error('different zgrid in source reconstructions');
    end
  end
  grandavg.xgrid = varargin{1}.xgrid;
  grandavg.ygrid = varargin{1}.ygrid;
  grandavg.zgrid = varargin{1}.zgrid;
end

if isfield(varargin{1}, 'transform')
  % check that the homogenous transformation matrix of each input volume is the same
  for i=2:Nsubject
    if any(varargin{i}.transform(:)~=varargin{1}.transform(:))
      error('different homogenous transformation matrices in source reconstructions');
    end
  end
  grandavg.transform = varargin{1}.transform;
end

% get the source parameter from each input source reconstruction
% get the inside parameter from each input source reconstruction
for i=1:Nsubject
  % TODO this function should use parameterselection
  if issubfield(varargin{i}, ['avg.' cfg.parameter])
    tmp = getsubfield(varargin{i}, ['avg.' cfg.parameter]);
  else
    tmp = getsubfield(varargin{i}, cfg.parameter);
  end
  dat(:,i) = tmp(:);
  tmp = getsubfield(varargin{i}, 'inside');
  inside(tmp,i) = 1;
end
% ensure that voxels that are not in the scanned brain region are excluded from the averaging
dat(~inside) = nan;

if strcmp(cfg.randomization, 'yes') || strcmp(cfg.permutation, 'yes')
  if strcmp(cfg.keepindividual, 'yes')
    error('you cannot keep individual data in combination with randomization or permutation');
  end

  % construct a design vector that contains the condition number 1 or 2
  design = zeros(1,Nsubject);
  design(cfg.c1) = 1;
  design(cfg.c2) = 2;
  if any(design==0)
    error('not all input source structures have been assigned to a condition');
  elseif length(design)~=Nsubject
    error('not enough input source structures given cfg.c1 and cfg.c2');
  end

  % create a matrix with all randomized assignments to the two conditions
  if strcmp(cfg.randomization, 'yes')
    res = zeros(cfg.numrandomization, Nsubject);
    for i=1:cfg.numrandomization
      res(i,:) = design(randperm(Nsubject));
    end
  elseif strcmp(cfg.permutation, 'yes')
    sel1 = find(design==1);
    sel2 = find(design==2);
    if length(sel1)~=length(sel2)
      error('permutation requires that there is an equal number of replications in each conditions')
    end
    res = zeros(cfg.numpermutation, Nsubject);
    for i=1:cfg.numpermutation
      flip = randn(1,length(sel1))>0;
      res(i,sel1) = 1;
      res(i,sel2) = 2;
      res(i,sel1(find(flip))) = 2;
      res(i,sel2(find(flip))) = 1;
    end
  end % randomization or permutation

  % randomize the input source parameter between the two conditions
  clear trialA
  clear trialB
  for i=1:size(res,1)
    selA = find(res(i,:)==1);
    selB = find(res(i,:)==2);
    % create the randomized averaged data
    trialA(i) = setsubfield([], cfg.parameter, nan_mean(dat(:,selA),2));
    trialB(i) = setsubfield([], cfg.parameter, nan_mean(dat(:,selB),2));
  end
  % create the observed average data
  selA = find(design==1);
  selB = find(design==2);
  avgA = setsubfield([], cfg.parameter, nan_mean(dat(:,selA),2));
  avgB = setsubfield([], cfg.parameter, nan_mean(dat(:,selB),2));

  % construct a source structure that can be fed into SOURCESTATISTICS_RANDOMIZATION or SOURCESTATISTICS_RANDCLUSTER
  grandavg.trialA  = trialA;
  grandavg.trialB  = trialB;
  grandavg.avgA    = avgA;
  grandavg.avgB    = avgB;

else
  if strcmp(cfg.concatenate, 'no'),
    % compute a plain average and variance over all input source structures
    grandavg.avg    = setsubfield([], cfg.parameter, nan_mean(dat,2));
    grandavg.var    = setsubfield([], cfg.parameter, nan_std(dat')'.^2);  % nan_std operates over the first dimension
    grandavg.dimord = 'voxel';
  else
    grandavg.avg    = setsubfield([], cfg.parameter, dat);
    grandavg.dimord = 'voxel_freq';
    grandavg.dim    = [grandavg.dim size(dat,2)];
  end 
 
  if strcmp(cfg.keepindividual, 'yes')
    clear trial
    for i=1:Nsubject
      trial(i) = setsubfield([], cfg.parameter, dat(:,i));
    end
    grandavg.trial = trial;
  end
end 

% determine which sources were inside or outside the brain in all subjects
allinside  = find(all( inside,2));
alloutside = find(all(~inside,2));
someinside = find(any( inside,2));
fprintf('%d voxels are inside the brain of all subjects\n',               length(allinside));
fprintf('%d voxels are inside the brain of some, but not all subjects\n', length(someinside));
fprintf('%d voxels are outside the brain of all subjects\n',              length(alloutside));
warning('marking only voxels inside the brain of all subjects as ''inside''');
if strcmp(cfg.concatenate, 'no'),
  grandavg.inside  = find(all( inside,2));
  grandavg.outside = find(any(~inside,2));
  grandavg.df      = sum(inside,2);
else
  grandavg.inside  = find(inside(:));
  grandavg.outside = setdiff([1:prod(size(dat))]', grandavg.inside); 
end

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes'); 

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: sourcegrandaverage.m,v 1.21 2009/01/20 13:01:31 sashae Exp $';
% remember the configuration details of the input data
cfg.previous = [];
for i=1:Nsubject
  try, cfg.previous{i} = varargin{i}.cfg; end
end
% remember the exact configuration details in the output 
grandavg.cfg = cfg;

