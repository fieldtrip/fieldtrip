function test_bug2468

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_sourcedescriptives
% DATA private

% this test function tests the functionality to do projectmom on the output
% of a pcc-beamformer with a different number of components per dipole.

filename = dccnpath('/project/3031000.02/test/bug2468.mat');
load(filename);

% this is a dirty fix which is intended to be fixed upstream in the
% analysis pipeline, i.e. in ft_inverse_pcc
for k = 1:numel(sdata.avg.mom)
  csdlabel{k,1} = repmat({'scandip'}, [1 size(sdata.avg.mom{k},1)]);
end
sdata.avg.csdlabel = csdlabel;

cfg            = [];
cfg.projectmom = 'yes';
sd             = ft_sourcedescriptives(cfg, sdata);

% recompute source level pcc data
cfg                             = [];
cfg.frequency                   = fdata.freq;
cfg.headmodel                   = sourceVol;
cfg.sourcemodel                        = leadfieldModel;
cfg.sourcemodel.filter                 = spatialFilters.avg.filter;
cfg.keeptrials                  = 'no';
cfg.method                      = 'pcc';
cfg.(cfg.method).keepfilter     = 'yes';
cfg.(cfg.method).lambda         = '5%';
cfg.(cfg.method).fixedori       = 'no';
cfg.(cfg.method).keepmom        = 'yes';
sdatanew                        = ft_sourceanalysis(cfg, fdata);

cfg            = [];
cfg.projectmom = 'yes';
sdnew          = ft_sourcedescriptives(cfg, sdatanew);

% now also test the new 'functionality' for the other options of
% ft_sourcedescriptives
cfg            = [];
cfg.projectmom = 'yes';
cfg.keeptrials = 'yes';
sdnew2a        = ft_sourcedescriptives(cfg, sdata);
sdnew2b        = ft_sourcedescriptives(cfg, sdatanew);
cfg.projectmom = 'no';
sdnew3a        = ft_sourcedescriptives(cfg, sdata);
sdnew3b        = ft_sourcedescriptives(cfg, sdatanew);

