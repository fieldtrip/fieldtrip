function test_pull1377b

% MEM 6gb
% WALLTIME 1:00:00
% DEPENDENCY ft_sourceanalysis ft_dipolefitting

load(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/raw/meg/preproc_ctf151.mat'), 'data');
grad = data.grad;
clear data

%%

vol = [];
vol.o = [0 0 4];
vol.r = 12;
vol.unit = 'cm';
vol.type = 'singlesphere';

%%

cfg = [];
cfg.dip.pos = [0 0 9];
cfg.dip.mom = [1 0 0];
cfg.dip.unit = 'cm';
cfg.headmodel = vol;
cfg.grad = grad;
data = ft_dipolesimulation(cfg);

cfg = [];
timelock = ft_timelockanalysis(cfg, data);

cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.output = 'fourier';
freq = ft_freqanalysis(cfg, data);

%%

timelockmethod = {'lcmv', 'sam', 'mne', 'rv', 'music', 'pcc', 'mvl', 'sloreta', 'eloreta'};
timelockmethod = {'lcmv', 'pcc'}; % these are the only one that work

for i=1:numel(timelockmethod)
  tfdata = timelock;
  method = timelockmethod{i};
  test_sourceanalysis
end

freqmethod = {'dics', 'pcc', 'eloreta', 'mne','harmony', 'rv', 'music'};
freqmethod = {'dics', 'pcc', }; % these are the only one that work

for i=1:numel(freqmethod)
  tfdata = freq;
  method = freqmethod{i};
  test_sourceanalysis
end

tfdata = timelock;
test_dipolefitting

% this one does not work
% tfdata = freq;
% test_dipolefitting

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

  function test_sourceanalysis
    
    cfg = [];
    cfg.sourcemodel.pos = [0 0 7];
    cfg.sourcemodel.unit = 'cm';
    cfg.headmodel = vol;
    cfg.method = method;
    cfg.keepleadfield = 'yes';
    cfg.(method).keepfilter = 'yes';
    cfg.(method).lambda = '10%';
    cfg.channel = 'MEG';
    source1 = ft_sourceanalysis(cfg, tfdata);
    
    %%
    % reuse the leadfield
    
    cfg = [];
    cfg.sourcemodel = keepfields(source1, {'pos' 'inside' 'leadfield', 'leadfielddimord', 'label'});
    cfg.headmodel = vol;
    cfg.method = method;
    cfg.(method).lambda = '10%';
    cfg.channel = 'MEG';
    source2a = ft_sourceanalysis(cfg, tfdata);
    
    %%
    % reuse the leadfield, but shuffle channels
    
    cfg = [];
    cfg.sourcemodel = keepfields(source1, {'pos' 'inside' 'leadfield', 'leadfielddimord', 'label'});
    
    % revert the channel order
    cfg.sourcemodel.label  = cfg.sourcemodel.label(end:-1:1,:);
    for i=1:size(cfg.sourcemodel.pos,1)
      cfg.sourcemodel.leadfield{i} = cfg.sourcemodel.leadfield{i}(end:-1:1,:);
    end
    
    cfg.headmodel = vol;
    cfg.method = method;
    cfg.(method).lambda = '10%';
    cfg.channel = 'MEG';
    source2b = ft_sourceanalysis(cfg, tfdata);
    
    %%
    % reuse the filter
    
    cfg = [];
    cfg.sourcemodel = keepfields(source1, {'pos' 'inside'});
    cfg.sourcemodel.label = source1.avg.label;
    cfg.sourcemodel.filter = source1.avg.filter;
    cfg.sourcemodel.filterdimord = source1.avg.filterdimord;
    
    cfg.headmodel = vol;
    cfg.method = method;
    cfg.(method).lambda = '10%';
    cfg.channel = 'MEG';
    source2c = ft_sourceanalysis(cfg, tfdata);
    
    %%
    % reuse the filter, but shuffle channels
    
    cfg = [];
    cfg.sourcemodel = keepfields(source1, {'pos' 'inside' 'filter', 'filterdimord', 'label'});
    cfg.sourcemodel.filter = source1.avg.filter;
    cfg.sourcemodel.filterdimord = source1.avg.filterdimord;
    
    % revert the channel order
    cfg.sourcemodel.label  = cfg.sourcemodel.label(end:-1:1,:);
    for i=1:size(cfg.sourcemodel.pos,1)
      cfg.sourcemodel.filter{i} = cfg.sourcemodel.filter{i}(:,end:-1:1);
    end
    
    cfg.headmodel = vol;
    cfg.method = method;
    cfg.(method).lambda = '10%';
    cfg.channel = 'MEG';
    source2d = ft_sourceanalysis(cfg, tfdata);
    
    %%
    
    if issubfield(source1, 'avg.mom')
      % these should all be the same
      assert(isequal(source1.avg.mom, source2a.avg.mom));
      assert(isequal(source1.avg.mom, source2b.avg.mom));
      assert(isequal(source1.avg.mom, source2c.avg.mom));
      assert(isequal(source1.avg.mom, source2d.avg.mom));
    end
    
    if issubfield(source1, 'avg.pow')
      % these should all be the same
      assert(isequal(source1.avg.pow, source2a.avg.pow));
      assert(isequal(source1.avg.pow, source2b.avg.pow));
      assert(isequal(source1.avg.pow, source2c.avg.pow));
      assert(isequal(source1.avg.pow, source2d.avg.pow));
    end
    
  end % function test_sourceanalysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

  function test_dipolefitting
    
    cfg = [];
    cfg.sourcemodel.pos = [0 0 7];
    cfg.sourcemodel.unit = 'cm';
    cfg.headmodel = vol;
    sourcemodel = ft_prepare_leadfield(cfg, tfdata);
    
    %%
    
    cfg = [];
    cfg.sourcemodel = sourcemodel;
    cfg.headmodel = vol;
    cfg.gridsearch = 'yes';
    cfg.nonlinear = 'no';
    cfg.channel = 'MEG';
    cfg.frequency = [10 11 12]; % only used for frequency data
    source3a = ft_dipolefitting(cfg, tfdata);
    
    %%
    
    cfg = [];
    cfg.sourcemodel = sourcemodel;
    cfg.headmodel = vol;
    cfg.gridsearch = 'no';
    cfg.nonlinear = 'yes';
    cfg.channel = 'MEG';
    cfg.frequency = 10; % only used for frequency data
    source3b = ft_dipolefitting(cfg, tfdata);
    
    %%
    
    try
      cfg = [];
      cfg.sourcemodel = source1;
      cfg.headmodel = vol;
      cfg.channel = 'MEGREF'; % the channels in the precomputed leadfield do not include MEGREF
      cfg.frequency = 10; % only used for frequency data
      source3c = ft_dipolefitting(cfg, tfdata);
      passed = true;
    catch
      passed = false;
    end
    assert(~passed, 'this should have resulted in an error');
  end % function test_dipolefitting


end % function
