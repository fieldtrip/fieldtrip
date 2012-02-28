function test_ft_freqanalysis(datainfo, writeflag, version)

% TEST test_ft_freqanalysis 
% TEST ft_freqanalysis test_datasets

% writeflag determines whether the output should be saved to disk
% version determines the output directory

if nargin<1
  datainfo = test_datasets;
end
if nargin<2
  writeflag = 0;
end
if nargin<3
  version = 'latest';
end

for k = 1:numel(datainfo)
  datanew = freqanalysisMtmfft(datainfo(k), writeflag, version);

  fname = fullfile(datainfo(k).origdir,version,'freq',datainfo(k).type,['freq_mtmfft_',datainfo(k).datatype]);
  load(fname);
  datanew = rmfield(datanew, 'cfg'); % these are per construction different if writeflag = 0;
  freq    = rmfield(freq,    'cfg');
  assert(isequal(freq, datanew));
end

for k = 1:numel(datainfo)
  datanew = freqanalysisMtmconvol(datainfo(k), writeflag, version);

  fname = fullfile(datainfo(k).origdir,version,'freq',datainfo(k).type,['freq_mtmconvol_',datainfo(k).datatype]);
  load(fname);
  datanew = rmfield(datanew, 'cfg'); % these are per construction different if writeflag = 0;
  freq    = rmfield(freq,    'cfg');
  assert(isequalwithequalnans(freq, datanew));
end

%----------------------
% subfunctions
%----------------------
function [freq] = freqanalysisMtmconvol(dataset, writeflag, version)

% --- HISTORICAL --- attempt forward compatibility with function handles
if ~exist('ft_freqanalysis') && exist('freqanalysis')
  eval('ft_freqanalysis = @freqanalysis;');
end

cfg        = [];
cfg.method = 'mtmconvol';
cfg.output = 'pow';
cfg.foi    = 2:2:30;
cfg.taper  = 'hanning';
cfg.t_ftimwin = ones(1,numel(cfg.foi)).*0.5;
cfg.toi    = (250:50:750)./1000;
cfg.inputfile = fullfile(dataset.origdir,version,'raw',dataset.type,['preproc_',dataset.datatype]);
if writeflag,
  cfg.outputfile = fullfile(dataset.origdir,version,'freq',dataset.type,['freq_mtmconvol_',dataset.datatype]);
end

if ~strcmp(version, 'latest') && str2num(version)<20100000
  % -- HISTORICAL --- older fieldtrip versions don't support inputfile and outputfile
  load(cfg.inputfile, 'data');
  freq = ft_freqanalysis(cfg, data);
  save(cfg.outputfile, 'freq');
else
  freq = ft_freqanalysis(cfg);
end

function [freq] = freqanalysisMtmfft(dataset, writeflag, version)

% --- HISTORICAL --- attempt forward compatibility with function handles
if ~exist('ft_freqanalysis') && exist('freqanalysis')
  eval('ft_freqanalysis = @freqanalysis;');
end

cfg        = [];
cfg.method = 'mtmfft';
cfg.output = 'pow';
cfg.foilim = [0 100];
cfg.taper  = 'hanning';
cfg.inputfile = fullfile(dataset.origdir,version,'raw',dataset.type,['preproc_',dataset.datatype]);
if writeflag,
  cfg.outputfile = fullfile(dataset.origdir,version,'freq',dataset.type,['freq_mtmfft_',dataset.datatype]);
end

if ~strcmp(version, 'latest') && str2num(version)<20100000
  % -- HISTORICAL --- older fieldtrip versions don't support inputfile and outputfile
  load(cfg.inputfile, 'data');
  freq = ft_freqanalysis(cfg, data);
  save(cfg.outputfile, 'freq');
else
  freq = ft_freqanalysis(cfg);
end

