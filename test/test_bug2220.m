function test_bug2220(datainfo)

% MEM 2gb
% WALLTIME 00:10:00

% TEST ft_preprocessing ft_preproc_padding preproc


%load('C:\Users\jorhor\Downloads\bugdatafile.mat')
cfg = [];
cfg.baselinewindow = [-.5 0];
cfg.demean = 'yes';

data = [];
data.label  = {'chan'};
data.fsample = 512;

for tr=1:10
  data.time{tr}    = linspace(-1, 2, 3*data.fsample+1);
  data.trial{tr}   = rand(1, length(data.time{tr})) .^ 2;
end

tmpdata = ft_preprocessing(cfg, data);

tbase = tmpdata.time{1} >= cfg.baselinewindow(1) & tmpdata.time{end} <= cfg.baselinewindow(2);

[ok, msg] = isalmostequal(sum(tmpdata.trial{1}(:, tbase), 2), zeros(length(data.label), 1),'abstol',eps*100);

if ~ok
  error('baseline not zero: %s', msg{1});
end

%% test with data on disk
if nargin<1
  datainfo = ref_datasets;
end

for k = 1:numel(datainfo)
  dataset = datainfo(k);
  
  cfg            = [];
  cfg.dataset    = fullfile(dataset.origdir,'original',dataset.type,dataset.datatype,'/',dataset.filename);

  % get header and event information
  if ~isempty(dataset.dataformat)
    hdr   = ft_read_header(cfg.dataset, 'headerformat', dataset.dataformat);
    event = ft_read_event(cfg.dataset, 'eventformat', dataset.dataformat);

    cfg.dataformat   = dataset.dataformat;
    cfg.headerformat = dataset.dataformat;
  else
    hdr   = ft_read_header(cfg.dataset);
    event = ft_read_event(cfg.dataset);
  end

  % create 10 1-second trials to be used as test-case 
  begsample = ((1:10)-1)*round(hdr.Fs) + 1;
  endsample = ((1:10)  )*round(hdr.Fs);
  offset    = zeros(1,10);

  cfg.trl   = [begsample(:) endsample(:) offset(:)];
  sel       = cfg.trl(:,2)<=hdr.nSamples*hdr.nTrials;
  cfg.trl   = cfg.trl(sel,:);

  cfg.continuous     = 'yes';
  cfg.demean         = 'yes';
  cfg.baselinewindow = [0 0.5];
  if strcmp(dataset.type, 'meg') || strcmp(dataset.type, 'eeg')
    cfg.channel      = dataset.type;
  else
    cfg.channel      = 'all';
  end
  data               = ft_preprocessing(cfg);

  tbase = data.time{1} >= cfg.baselinewindow(1) & data.time{end} <= cfg.baselinewindow(2);
 [ok, msg] = isalmostequal(sum(data.trial{1}(1:end, tbase), 2), zeros(length(data.label), 1),'abstol',abs(.0001 * mean(data.trial{1}(:))));

  if ~ok
    error('baseline not zero: %s', msg{1});
  end
  
end
