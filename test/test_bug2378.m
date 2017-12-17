function test_bug2378

% MEM 1gb
% WALLTIME 00:10:00

% TEST ft_redefinetrial

% This function tests the functionality of ft_redefinetrial, whether it
% keeps the trialinfo in all the different ways of calling it.

%------------------------------------------
%generate data
data = [];
data.fsample = 1000;
data.cfg     = [];

nsmp  = 1000;
nchan = 80;
for k = 1:10
  data.trial{k} = randn(nchan,nsmp);
  data.time{k}  = ((1:nsmp)-1)./data.fsample;
end

% create grad-structure and add to data
grad.pnt  = randn(nchan,3);
grad.ori  = randn(nchan,3);
grad.tra  = eye(nchan);
for k = 1:nchan
  grad.label{k,1} = ['chan',num2str(k,'%03d')];
end
data.grad  = ft_datatype_sens(grad);
data.label = grad.label;
data.trialinfo = (1:10)';
data = ft_checkdata(data, 'hassampleinfo', 'yes');


% call ft_redefinetrial with all different methods
% "WORKS FINE" comments were added on 20 Nov 2013.

% overlapping data
cfg = [];
cfg.length = 0.5;
cfg.overlap = 0;
datanew = ft_redefinetrial(cfg, data);
assert(isfield(datanew, 'trialinfo'));  % WORKS FINE

% select subset of trials:
cfg = [];
cfg.trials = [1,1,1];   
datanew = ft_redefinetrial(cfg, data);
assert(isfield(datanew, 'trialinfo'));  % WORKS FINE

% realign to new reference time point
cfg = [];
cfg.offset = 200; 
% cfg.offset    = single number or Nx1 vector, expressed in samples relative to current t=0
datanew = ft_redefinetrial(cfg, data);
assert(isfield(datanew, 'trialinfo'));  % WORKS FINE

% select specific (shorter) time window with toilim
cfg = [];
cfg.toilim = [0.2 0.4];
datanew = ft_redefinetrial(cfg, data);
assert(isfield(datanew, 'trialinfo'));  % WORKS FINE

% select specific (shorter) time window with beg & end sample
cfg = [];
cfg.begsample = 100;
cfg.endsample = 500;
datanew = ft_redefinetrial(cfg, data);
assert(isfield(datanew, 'trialinfo'));  % WORKS FINE

mytrl = [
  0001 0300 0 1
  0301 0600 0 2
  0601 0900 0 3
  0901 1200 0 4
  0201 1500 0 5
  0501 1800 0 6
  1801 2100 0 7
  ];

% select specific time window with new trl matrix
cfg = [];
cfg.trl = mytrl; % use ft_definetrial.m
datanew = ft_redefinetrial(cfg, data);
assert(isfield(datanew, 'trialinfo'));  % WORKS FINE
