function cfg = ft_databrowser_multitrial(cfg,data)

% FT_DATABROWSER_MULTITRIAL can be used for visual inspection of data. Artifacts that were
% detected by artifact functions (see FT_ARTIFACT_xxx functions where xxx is the type
% of artifact) are marked. Additionally data pieces can be marked and unmarked as
% artifact by manual selection. The output cfg contains the updated specification of
% the artifacts.

% add one option: how many trials shown on the first screen
cfg.ntrialswin = ft_getopt(cfg, 'ntrialswin', 5);
% retrieve events (ft_databrowser code)
cfg.event = ft_getopt(cfg, 'event'); 
if ~isempty(cfg.event)
  % use the events that the user passed in the configuration
  event = cfg.event;
else
  % fetch the events from the data structure in memory
  event = ft_fetch_event(data);
end

% store original sampleinfo and events
urevent = event;

% add the sample indices as extra channel
for i=1:numel(data.trial)
  smp{i} = data.sampleinfo(i,1):data.sampleinfo(i,2);
end
smp = cat(2,smp{:});

% compute new sampleinfo
begsample = [0 cumsum(cellfun(@length, data.trial(1:end-1)))]' + 1;
endsample = cumsum(cellfun(@length, data.trial))';
% pretend that the trials are pseudo-continuous
data.sampleinfo = [begsample endsample];

% retrieve event samples
event = renamefields(event,'sample','ursample');
for iev = 1:numel(event)
  event(iev).sample = find(smp == event(iev).ursample);
  % events outside of trials aren't found
  if isempty(event(iev).sample)
    event(iev).sample = NaN;
  end
end
% retrieve artifact samples
if isfield(cfg,'artfctdef')
  artnames = fieldnames(cfg.artfctdef);
  for iartfct = 1:numel(artnames)
    art = cfg.artfctdef.(artnames{iartfct}).artifact;
    if isempty(art)
      continue
    end
    cfg.artfctdef.(artnames{iartfct}).urartifact = cfg.artfctdef.(artnames{iartfct}).artifact;
    for iart = 1:numel(art)
      cfg.artfctdef.(artnames{iartfct}).artifact(iart) = find(smp == art(iart));
    end
  end
end

%% add events to highlight trial boundaries
boundarysample = data.sampleinfo(:,1);
boundary = [];
for i = 1:numel(boundarysample)
  boundary(i).type = '';
  boundary(i).value = '';
  boundary(i).sample = boundarysample(i);
  boundary(i).offset = [];
  boundary(i).duration = 0;
  boundary(i).timestamp = boundary(i).sample / data.fsample;
  boundary(i).ursample = NaN;
end
cfg.event = [event boundary];

%%%%%%%%%% forward to databrowser

cfg.continuous = 'yes'; 
cfg.blocksize = (size(cat(2,data.trial{1:cfg.ntrialswin}),2) / data.fsample);

cfgout = ft_databrowser(cfg,data);

%%%%%%%%%% recover previous events
% (these should not have changed in ft_databrowser)
cfgout.event = urevent;

%%%%%%%%%% switch artifacts back to original time frame
if isfield(cfgout,'artfctdef')
  artnames = fieldnames(cfgout.artfctdef);
  for iart = 1:numel(artnames)
    art = cfgout.artfctdef.(artnames{iart}).artifact;
    cfgout.artfctdef.(artnames{iart}).artifact = smp(art);
  end
end

cfg = cfgout;

