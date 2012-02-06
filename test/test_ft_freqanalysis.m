function test_ft_freqanalysis(datainfo, writeflag)

% TEST test_ft_freqanalysis 
% TEST ft_freqanalysis test_datasets

% the optional writeflag determines whether the output 
% should be saved to disk

if nargin<2
  writeflag = 0;
end

if nargin<1
  datainfo = test_datasets;
end

for k = 1:numel(datainfo)
  datanew = freqanalysisMtmfft(datainfo(k), writeflag);

  fname = [datainfo(k).origdir,'latest/freq/',datainfo(k).type,'freq_mtmfft_',datainfo(k).datatype];
  load(fname);
  datanew = rmfield(datanew, 'cfg'); % these are per construction different if writeflag = 0;
  freq    = rmfield(freq,    'cfg');
  assert(isequal(freq, datanew));
end

for k = 1:numel(datainfo)
  datanew = freqanalysisMtmconvol(datainfo(k), writeflag);

  fname = [datainfo(k).origdir,'latest/freq/',datainfo(k).type,'freq_mtmconvol_',datainfo(k).datatype];
  load(fname);
  datanew = rmfield(datanew, 'cfg'); % these are per construction different if writeflag = 0;
  freq    = rmfield(freq,    'cfg');
  assert(isequalwithequalnans(freq, datanew));
end

%----------------------
% subfunctions
%----------------------
function [freq] = freqanalysisMtmconvol(dataset, writeflag)

cfg        = [];
cfg.method = 'mtmconvol';
cfg.output = 'pow';
cfg.foi    = 2:2:30;
cfg.taper  = 'hanning';
cfg.t_ftimwin = ones(1,numel(cfg.foi)).*0.5;
cfg.toi    = (250:50:750)./1000;
cfg.inputfile  = [dataset.origdir,'latest/raw/',dataset.type,'preproc_',dataset.datatype];
if writeflag,
  cfg.outputfile = [dataset.origdir,'latest/freq/',dataset.type,'freq_mtmconvol_',dataset.datatype];
end
freq = ft_freqanalysis(cfg);

function [freq] = freqanalysisMtmfft(dataset, writeflag)

cfg        = [];
cfg.method = 'mtmfft';
cfg.output = 'pow';
cfg.foilim = [0 100];
cfg.taper  = 'hanning';
cfg.inputfile  = [dataset.origdir,'latest/raw/',dataset.type,'preproc_',dataset.datatype];
if writeflag,
  cfg.outputfile = [dataset.origdir,'latest/freq/',dataset.type,'freq_mtmfft_',dataset.datatype];
end
freq = ft_freqanalysis(cfg);
