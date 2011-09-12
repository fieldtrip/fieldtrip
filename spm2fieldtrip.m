function [data] = spm2fieldtrip(D)

% SPM2FIELDTRIP converts an SPM8 meeg object into a FieldTrip raw data structure
%
% Use as
%   data = spm2fieldtrip(D)
% where D is the SPM8 meeg object which you can load in with SPM_EEG_LOAD
% and where data is a FieldTrip raw data structure as if it were returned
% by FT_PREPROCESSING.
%
% See also FT_PREPROCESSING, SPM_EEG_LOAD

ft_defaults

% record start time and total processing time
ftFuncTimer = tic();
ftFuncClock = clock();
ftFuncMem   = memtic();

if ~ft_hastoolbox('spm8')
  error('this requires a full version of SPM8 on your MATLAB path');
end

if ~isa(D, 'meeg')
  error('this requires an SPM8 "meeg" object as input')
end

% this is how SPM8 represents it
data = D.ftraw(0);

clist      = D.condlist;
conditions = D.conditions;

data.trialinfo = zeros(ntrials,1);
for k = 1:numel(clist)
  fprintf('mapping condition label "%s" to condition code %d\n', clist{k}, k);
  sel=strcmp(clist{k}, conditions);
  data.trialinfo(sel) = k;
end

% FIXME the following is not correct
%
% data.sampleinfo = zeros(ntrials,2);
% for i=1:ntrials
%   data.sampleinfo(i,1) = D.indsample(i);
%   data.sampleinfo(i,2) = D.indsample(i) + D.nsamples;
% end

% add the version details of this function call to the configuration
cfg.version.name = mfilename('fullpath');
cfg.version.id   = '$Id$';

% add information about the Matlab version used to the configuration
cfg.callinfo.matlab = version();

% add information about the function call to the configuration
cfg.callinfo.proctime = toc(ftFuncTimer);
cfg.callinfo.procmem  = memtoc(ftFuncMem);
cfg.callinfo.calltime = ftFuncClock;
cfg.callinfo.user = getusername();
fprintf('the call to "%s" took %d seconds and an estimated %d MB\n', mfilename, round(cfg.callinfo.proctime), round(cfg.callinfo.procmem/(1024*1024)));

% remember the configuration details
data.cfg = cfg;

