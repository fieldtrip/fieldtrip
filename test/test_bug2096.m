function test_bug2096

% MEM 4000mb
% WALLTIME 00:20:00

% TEST test_bug2096
% TEST ft_sourcewrite ft_read_cifti ft_write_cifti

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2096'));

%%
source = [];
source.dim    = [11 12 13];
[X Y Z]       = ndgrid(1:source.dim(1), 1:source.dim(2), 1:source.dim(3));
source.transform = eye(4);
source.transform(1,4) = -(source.dim(1)+1)/2;
source.transform(2,4) = -(source.dim(2)+1)/2;
source.transform(3,4) = -(source.dim(3)+1)/2;
source.pos    = ft_warp_apply(source.transform, [X(:) Y(:) Z(:)]);
source.pow    = (1:prod(source.dim))';
source.dimord = 'pos';

cfg = [];
cfg.filetype  = 'cifti';
cfg.parameter = 'pow';
cfg.filename  = 'test_bug2096';
ft_sourcewrite(cfg, source);

source1 = ft_read_cifti('test_bug2096.pow.dscalar.nii')

cfg = [];
cfg.filetype  = 'cifti';
cfg.parameter = 'pow';
cfg.filename  = 'test_bug2096b';
ft_sourcewrite(cfg, source1);

source2 = ft_read_cifti('test_bug2096b.pow.dscalar.nii')

% assert(isequal(source , source1)); % source1 has more details
% assert(isequal(source1, source2)); % numerical differences

%%
[pnt, tri] = icosahedron;

source = [];
source.pos    = pnt;
source.tri    = tri;
source.pow    = (1:size(pnt,1))';
source.dimord = 'pos';
% source.BrainStructure = ones(1, size(pnt,1));
% source.BrainStructurelabel = {'CORTEX'};
source.BrainStructure = [1 1 1 1 1 1 2 2 2 2 2 2];
source.BrainStructurelabel = {'CORTEX_LEFT', 'CORTEX_RIGHT'};

cfg = [];
cfg.filetype  = 'cifti';
cfg.parameter = 'pow';
cfg.filename  = 'test_bug2096';
ft_sourcewrite(cfg, source);

source1 = ft_read_cifti('test_bug2096.pow.dscalar.nii')

cfg = [];
cfg.filetype  = 'cifti';
cfg.parameter = 'pow';
cfg.filename  = 'test_bug2096b';
ft_sourcewrite(cfg, source1);

source2 = ft_read_cifti('test_bug2096b.pow.dscalar.nii')

% assert(isequal(source , source1));
% assert(isequal(source1, source2));


%%
source = [];
source.dim = [5 6 7];
[X Y Z] = ndgrid(1:source.dim(1), 1:source.dim(2), 1:source.dim(3));
source.transform = eye(4);
source.pos     = ft_warp_apply(source.transform, [X(:) Y(:) Z(:)]);
source.imagcoh = randn(prod(source.dim));
source.dimord  = 'pos_pos';

cfg = [];
cfg.filetype  = 'cifti';
cfg.parameter = 'imagcoh';
cfg.filename  = 'test_bug2096';
ft_sourcewrite(cfg, source);

source1 = ft_read_cifti('test_bug2096.imagcoh.dconn.nii');

cfg = [];
cfg.filetype  = 'cifti';
cfg.parameter = 'imagcoh';
cfg.filename  = 'test_bug2096b';
ft_sourcewrite(cfg, source1);

source2 = ft_read_cifti('test_bug2096b.imagcoh.dconn.nii')

% assert(isequal(source , source1));
% assert(isequal(source2, source1));


%%
source = [];
source.dim = [5 6 7];
[X Y Z] = ndgrid(1:source.dim(1), 1:source.dim(2), 1:source.dim(3));
source.transform = eye(4);
source.pos = ft_warp_apply(source.transform, [X(:) Y(:) Z(:)]);
source.time = 1:10;
source.timeseries = randn(prod(source.dim), length(source.time));
source.dimord = 'pos_time';

cfg = [];
cfg.filetype = 'cifti';
cfg.parameter = 'timeseries';
cfg.filename = 'test_bug2096';
cfg.precision = 'single';
ft_sourcewrite(cfg, source);

source1 = ft_read_cifti('test_bug2096.timeseries.dtseries.nii');

%% test the dscalar output
[pnt, tri] = icosahedron;
pntL = pnt; pntL(:,1) = pntL(:,1) - 1; % shift along X
pntR = pnt; pntR(:,1) = pntR(:,1) + 1; % shift along X

source = [];
source.pos    = [pntL; pntR];
source.tri    = [tri; tri+12];
source.BrainStructure = [1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2];
source.BrainStructurelabel = {'CORTEX_LEFT', 'CORTEX_RIGHT'};
source.activity = zeros(size(source.pos,1), 1);
source.activity(:,1) = 1:24;
source.dimord = 'pos';

cfg = [];
cfg.filetype  = 'cifti';
cfg.parameter = 'activity';
cfg.filename  = 'test_bug2096';
ft_sourcewrite(cfg, source);

source1 = ft_read_cifti('test_bug2096.activity.dscalar.nii')
ft_plot_mesh(source1, 'vertexcolor', source1.activity(:,1), 'edgecolor', 'none')

%% test the dtsetries output
[pnt, tri] = icosahedron;
pntL = pnt; pntL(:,1) = pntL(:,1) - 1; % shift along X
pntR = pnt; pntR(:,1) = pntR(:,1) + 1; % shift along X

source = [];
source.pos    = [pntL; pntR];
source.tri    = [tri; tri+12];
source.BrainStructure = [1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2];
source.BrainStructurelabel = {'CORTEX_LEFT', 'CORTEX_RIGHT'};
source.time = 1:10;
source.timeseries = zeros(size(source.pos,1), length(source.time));
for i=1:size(source.timeseries,2)
  source.timeseries(:,i) = 1:24;
end
source.dimord = 'pos_time';

cfg = [];
cfg.filetype  = 'cifti';
cfg.parameter = 'timeseries';
cfg.filename  = 'test_bug2096';
ft_sourcewrite(cfg, source);

source1 = ft_read_cifti('test_bug2096.timeseries.dtseries.nii')
ft_plot_mesh(source1, 'vertexcolor', source1.timeseries(:,1), 'edgecolor', 'none')

%% test the dconn output
[pnt, tri] = icosahedron;
pntL = pnt; pntL(:,1) = pntL(:,1) - 1; % shift along X
pntR = pnt; pntR(:,1) = pntR(:,1) + 1; % shift along X

source = [];
source.pos    = [pntL; pntR];
source.tri    = [tri; tri+12];
source.imagcoh = rand(size(source.pos,1));
source.dimord = 'pos_pos';
source.BrainStructure = [1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2];
source.BrainStructurelabel = {'CORTEX_LEFT', 'CORTEX_RIGHT'};

cfg = [];
cfg.filetype  = 'cifti';
cfg.parameter = 'imagcoh';
cfg.filename  = 'test_bug2096';
ft_sourcewrite(cfg, source);

source1 = ft_read_cifti('test_bug2096.imagcoh.dconn.nii')

%% version 1
cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2096/cifti1'));


% cii1 = ft_read_cifti('DenseConnectome.dconn.nii');            % this one is disabled because it is 10GB large
% cii2 = ft_read_cifti('DenseTimeSeries.dtseries.nii');         % this one is disabled because the file contains an internal error (number of greynodes is not consistent with size of data)
% cii3 = ft_read_cifti('ParcellatedTimeSeries.ptseries.nii');   % this one is disabled because the code cannot deal with cifti-1 parcels
cii4 = ft_read_cifti('BOLD_REST2_LR.dtseries.nii');
cii5 = ft_read_cifti('BOLD_REST2_LR_Atlas.dtseries.nii');

%% version 2
cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2096/cifti2'));


cii1 = ft_read_cifti('ones.dscalar.nii');
cii2 = ft_read_cifti('Conte69.MyelinAndCorrThickness.32k_fs_LR.dscalar.nii');
cii3 = ft_read_cifti('Conte69.MyelinAndCorrThickness.32k_fs_LR.dtseries.nii');
% cii4 = ft_read_cifti('Conte69.MyelinAndCorrThickness.32k_fs_LR.ptseries.nii');
cii5 = ft_read_cifti('Conte69.parcellations_VGD11b.32k_fs_LR.dlabel.nii');

%% release data
cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2096/hcp_release/fsaverage_LR32k'));

filename = {
  '177746.ArealDistortion.32k_fs_LR.dscalar.nii'
  '177746.aparc.32k_fs_LR.dlabel.nii'
  '177746.BA.32k_fs_LR.dlabel.nii'
  '177746.aparc.a2009s.32k_fs_LR.dlabel.nii'
  '177746.MyelinMap.32k_fs_LR.dscalar.nii'
  '177746.corrThickness.32k_fs_LR.dscalar.nii'
  '177746.MyelinMap_BC.32k_fs_LR.dscalar.nii'
  '177746.curvature.32k_fs_LR.dscalar.nii'
  '177746.SmoothedMyelinMap.32k_fs_LR.dscalar.nii'
  '177746.sulc.32k_fs_LR.dscalar.nii'
  '177746.SmoothedMyelinMap_BC.32k_fs_LR.dscalar.nii'
  '177746.thickness.32k_fs_LR.dscalar.nii'
  };

datafield = {
  'arealdistortion'
  'aparc'
  'ba'
  'aparc_a2009s'
  'myelinmap'
  'corrthickness'
  'myelinmap_bc'
  'curvature'
  'smoothedmyelinmap'
  'sulc'
  'smoothedmyelinmap_bc'
  'thickness'
  };

for i=1:length(filename)
  disp(filename{i});
  source = ft_read_cifti(filename{i}, 'representation', 'source');
  figure
  ft_plot_mesh(source, 'vertexcolor', source.(datafield{i}), 'edgecolor', 'none');
  title(datafield{i});
end

%% MEG specific development data from DVE
cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2096/hcp_devel/TestParcelsForMEG'));

filename = {
  '3T_Q1-Q6related468_MSMsulc_d100_ts2_Znet2.pconn.nii'
  '3T_Q1-Q6related468_MSMsulc_d100_ts3_Znet2.pconn.nii'
  '3T_Q1-Q6related468_MSMsulc_d25_ts2_Znet2.pconn.nii'
  '3T_Q1-Q6related468_MSMsulc_d25_ts3_Znet2.pconn.nii'
  'HCP_Q1-Q6_R210mgtr_MSMRSN_R468_MSMsulc_d25_melodic_IC_ftb.32k_fs_LR.pconn.nii'
  %   'Q1-Q6_R440.L.midthickness.32k_fs_LR.surf.gii'
  %   'Q1-Q6_R440.L.very_inflated.32k_fs_LR.surf.gii'
  %   'Q1-Q6_R440.R.midthickness.32k_fs_LR.surf.gii'
  %   'Q1-Q6_R440.R.very_inflated.32k_fs_LR.surf.gii'
  'Q1-Q6_R440.sulc.32k_fs_LR.dscalar.nii'
  'Q1-Q6_RelatedParcellation210_mgtr_MSMRSNOrig3_d26_DR_DeDrift_Q1-Q6related468_MSMsulc_d25_melodic_IC_ftb_thr200sqmm_distance9.0.32k_fs_LR.ptseries.nii'
  'groupICA_3T_Q1-Q6related468_MSMsulc_d25_melodic_IC_ftb_thr200sqmm_distance9.0.dlabel.nii'
  'melodic_IC_ftb.dlabel.nii'
  %   'Q1-Q6_R440_ForMEG-Parcels.32k_fs_LR.wb.spec'
  %   'R468_pconn_31aug14.scene'
  };

for i=1:length(filename)
  disp(filename{i});
  source = ft_read_cifti(filename{i}, 'representation', 'source');
end
