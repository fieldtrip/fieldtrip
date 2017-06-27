function test_ft_freqanalysis(datainfo, writeflag, version)

% MEM 8000mb
% WALLTIME 01:30:00

% TEST test_ft_freqanalysis
% TEST ft_freqanalysis ref_datasets

% writeflag determines whether the output should be saved to disk
% version determines the output directory

if nargin<1
  datainfo = ref_datasets;
end
if nargin<2
  writeflag = 0;
end
if nargin<3
  version = 'latest';
end

for k = 1:numel(datainfo)
  datanew = freqanalysisMtmfft(datainfo(k), writeflag, version, 'fourier', 'yes');
  fname   = fullfile(datainfo(k).origdir,version,'freq',datainfo(k).type,['freq_mtmfft_fourier_trl_',datainfo(k).datatype]);
  load(fname); datanew = rmfield(datanew, 'cfg'); freq = rmfield(freq, 'cfg');
  [ok,msg] = isalmostequal(freq, datanew, 'reltol', 1e-4);
  if ~ok, error('stored and computed data not identical: %s', msg{:}); end
  
  datanew = freqanalysisMtmfft(datainfo(k), writeflag, version, 'powandcsd', 'yes');
  fname   = fullfile(datainfo(k).origdir,version,'freq',datainfo(k).type,['freq_mtmfft_powandcsd_trl_',datainfo(k).datatype]);
  load(fname); datanew = rmfield(datanew, 'cfg'); freq = rmfield(freq, 'cfg');
  [ok,msg] = isalmostequal(freq, datanew, 'reltol', 1e-4);
  if ~ok, error('stored and computed data not identical: %s', msg{:}); end
  
  datanew = freqanalysisMtmfft(datainfo(k), writeflag, version, 'powandcsd', 'no');
  fname   = fullfile(datainfo(k).origdir,version,'freq',datainfo(k).type,['freq_mtmfft_powandcsd_',datainfo(k).datatype]);
  load(fname); datanew = rmfield(datanew, 'cfg'); freq = rmfield(freq, 'cfg');
  [ok,msg] = isalmostequal(freq, datanew, 'reltol', 1e-4);
  if ~ok, error('stored and computed data not identical: %s', msg{:}); end
  
  datanew = freqanalysisMtmfft(datainfo(k), writeflag, version, 'pow', 'no');
  fname = fullfile(datainfo(k).origdir,version,'freq',datainfo(k).type,['freq_mtmfft_',datainfo(k).datatype]);
  load(fname); datanew = rmfield(datanew, 'cfg'); freq = rmfield(freq, 'cfg');
  [ok,msg] = isalmostequal(freq, datanew,'reltol',1e-4);
  if ~ok, error('stored and computed data not identical: %s', msg{:}); end
end

for k = 1:numel(datainfo)
  datanew = freqanalysisMtmconvol(datainfo(k), writeflag, version, 'fourier', 'yes');
  fname   = fullfile(datainfo(k).origdir,version,'freq',datainfo(k).type,['freq_mtmconvol_fourier_trl_',datainfo(k).datatype]);
  load(fname); datanew = rmfield(datanew, 'cfg'); freq = rmfield(freq, 'cfg');
  [ok,msg] = isalmostequal(freq, datanew, 'reltol', 1e-4);
  if ~ok, error('stored and computed data not identical: %s', msg{:}); end
  
  datanew = freqanalysisMtmconvol(datainfo(k), writeflag, version, 'powandcsd', 'yes');
  fname   = fullfile(datainfo(k).origdir,version,'freq',datainfo(k).type,['freq_mtmconvol_powandcsd_trl_',datainfo(k).datatype]);
  load(fname); datanew = rmfield(datanew, 'cfg'); freq = rmfield(freq, 'cfg');
  [ok,msg] = isalmostequal(freq, datanew, 'reltol', 1e-4);
  if ~ok, error('stored and computed data not identical: %s', msg{:}); end
  
  datanew = freqanalysisMtmconvol(datainfo(k), writeflag, version, 'powandcsd', 'no');
  fname   = fullfile(datainfo(k).origdir,version,'freq',datainfo(k).type,['freq_mtmconvol_powandcsd_',datainfo(k).datatype]);
  load(fname); datanew = rmfield(datanew, 'cfg'); freq = rmfield(freq, 'cfg');
  [ok,msg] = isalmostequal(freq, datanew, 'reltol', 1e-4);
  if ~ok, error('stored and computed data not identical: %s', msg{:}); end
  
  datanew = freqanalysisMtmconvol(datainfo(k), writeflag, version, 'pow', 'no');
  fname = fullfile(datainfo(k).origdir,version,'freq',datainfo(k).type,['freq_mtmconvol_',datainfo(k).datatype]);
  load(fname); datanew = rmfield(datanew, 'cfg'); freq = rmfield(freq, 'cfg');
  [ok,msg] = isalmostequal(freq, datanew,'reltol', 1e-4);
  if ~ok
    error('stored and computed data not identical: %s', msg{:});
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [freq] = freqanalysisMtmconvol(dataset, writeflag, version, output, keeptrials)

if isempty(output)
  output = 'pow';
  % output = 'powandcsd';
  % output = 'fourier';
end

if isempty(keeptrials)
  output = 'no';
  % output = 'yes';
end

% the file names should distinguish between the cfg.output and cfg.keeptrials option
postfix = '';
switch output
  case 'pow'
    % don't change
  case 'powandcsd'
    postfix = [postfix 'powandcsd_'];
  case 'fourier'
    postfix = [postfix 'fourier_'];
  otherwise
    error('unexpected output');
end

% the file names should distinguish between the cfg.output and cfg.keeptrials option
switch keeptrials
  case 'no'
    % don't change
  case 'yes'
    postfix = [postfix 'trl_'];
  otherwise
    error('unexpected keeptrials');
end

% --- HISTORICAL --- attempt forward compatibility with function handles
if ~exist('ft_freqanalysis') && exist('freqanalysis')
  eval('ft_freqanalysis = @freqanalysis;');
end

fprintf('testing mtmconvol with datatype=%s, output=%s, keeptrials=%s...\n',...
  dataset.datatype, output, keeptrials);

cfg            = [];
cfg.method     = 'mtmconvol';
cfg.output     = output;
cfg.keeptrials = keeptrials;
cfg.foi        = 2:2:30;
cfg.taper      = 'hanning';
cfg.t_ftimwin  = ones(1,numel(cfg.foi)).*0.5;
cfg.toi        = (250:50:750)./1000;
cfg.polyremoval= 0;
cfg.inputfile = fullfile(dataset.origdir,version,'raw',dataset.type,['preproc_',dataset.datatype]);
if writeflag
  cfg.outputfile = fullfile(dataset.origdir,version,'freq',dataset.type,['freq_mtmconvol_',postfix,dataset.datatype]);
end

if ~strcmp(version, 'latest') && str2num(version)<20100000
  % -- HISTORICAL --- older FieldTrip versions don't support inputfile and outputfile
  load(cfg.inputfile, 'data');
  freq = ft_freqanalysis(cfg, data);
  save(cfg.outputfile, 'freq');
else
  freq = ft_freqanalysis(cfg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [freq] = freqanalysisMtmfft(dataset, writeflag, version, output, keeptrials)

% --- HISTORICAL --- attempt forward compatibility with function handles
if ~exist('ft_freqanalysis') && exist('freqanalysis')
  eval('ft_freqanalysis = @freqanalysis;');
end

if isempty(output)
  output = 'pow';
  % output = 'powandcsd';
  % output = 'fourier';
end

if isempty(keeptrials)
  keeptrials = 'no';
  % keeptrials = 'yes';
end

% the file names should distinguish between the cfg.output and cfg.keeptrials option
postfix = '';
switch output
  case 'pow'
    % don't change
  case 'powandcsd'
    postfix = [postfix 'powandcsd_'];
  case 'fourier'
    postfix = [postfix 'fourier_'];
  otherwise
    error('unexpected output');
end

% the file names should distinguish between the cfg.output and cfg.keeptrials option
switch keeptrials
  case 'no'
    % don't change
  case 'yes'
    postfix = [postfix 'trl_'];
  otherwise
    error('unexpected keeptrials');
end

fprintf('testing mtmfft with datatype=%s, output=%s, keeptrials=%s...\n',...
  dataset.datatype, output, keeptrials);

cfg            = [];
cfg.method     = 'mtmfft';
cfg.output     = output;
cfg.keeptrials = keeptrials;
cfg.foilim     = [0 100];
cfg.taper      = 'hanning';
cfg.polyremoval= 0;
cfg.inputfile  = fullfile(dataset.origdir,version,'raw',dataset.type,['preproc_',dataset.datatype]);
if writeflag
  cfg.outputfile = fullfile(dataset.origdir,version,'freq',dataset.type,['freq_mtmfft_',postfix,dataset.datatype]);
end

if ~strcmp(version, 'latest') && str2num(version)<20100000
  % -- HISTORICAL --- older FieldTrip versions don't support inputfile and outputfile
  load(cfg.inputfile, 'data');
  freq = ft_freqanalysis(cfg, data);
  save(cfg.outputfile, 'freq');
else
  freq = ft_freqanalysis(cfg);
end

