function test_ft_sourceanalysis(datainfo, writeflag, version)

% TEST test_ft_sourceanalysis 
% TEST ft_sourceanalysis test_datasets

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

% make vol
vol = [];
vol.o = [0 0 4];
vol.r = 12;
vol.unit = 'cm';
vol.type = 'singlesphere';

% cortical sheet
load('/home/common/matlab/fieldtrip/data/test/corticalsheet.mat');
cfg=[];
cfg.vol = vol;
cfg.grid.pos = corticalsheet.pnt;
% FIXME this should not be needed
cfg.grid.inside = 1:size(corticalsheet.pnt,1);
cfg_corticalsheet = cfg;

% 3D regular grid
cfg=[];
cfg.vol = vol;
cfg.grid.resolution = 1.5;
cfg_3dgrid = cfg;

% small number of sparse regions of interest
cfg=[];
cfg.vol = vol;
cfg.grid.pos = [0 0 5; 1 0 5; -1 0 5; 0 1 5];
cfg_roi = cfg;

for k = 1:numel(datainfo)
  sourcenew = sourceanalysis_pcc(datainfo(k), writeflag, version, cfg_corticalsheet);
  sourcenew = sourceanalysis_pcc(datainfo(k), writeflag, version, cfg_3dgrid);
  sourcenew = sourceanalysis_pcc(datainfo(k), writeflag, version, cfg_roi);

  fname = fullfile(datainfo(k).origdir,version,'source',datainfo(k).type,['source_',datainfo(k).datatype]);
  load(fname);
  sourcenew = rmfield(sourcenew, 'cfg'); % these are per construction different if writeflag = 0;
  source    = rmfield(source, 'cfg');
  assert(isequal(source, sourcenew));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function source = sourceanalysis_pcc(dataset, writeflag, version, initialcfg)
source = [];
% FIXME, here goes the stuff from Arje, Stephen and Jorn

