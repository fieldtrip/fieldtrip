function test_bug2096

% MEM 4000mb
% WALLTIME 00:20:00

% TEST ft_sourcewrite ft_read_cifti ft_write_cifti

% needed for the dccnpath function, since we will change directory later on
addpath(fileparts(mfilename('fullpath')));

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2096'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% general purpose tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all

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
cfg.filename  = 'test_bug2096.pow';
ft_sourcewrite(cfg, source);

source1 = ft_read_cifti('test_bug2096.pow.dscalar.nii');

cfg = [];
cfg.filetype  = 'cifti';
cfg.parameter = 'pow';
cfg.filename  = 'test_bug2096b.pow';
ft_sourcewrite(cfg, source1);

source2 = ft_read_cifti('test_bug2096b.pow.dscalar.nii');

% assert(isequal(source , source1)); % source1 has more details
% assert(isequal(source1, source2)); % numerical differences

%%
[pnt, tri] = icosahedron;

source = [];
source.pos    = pnt;
source.tri    = tri;
source.pow    = (1:size(pnt,1))';
source.dimord = 'pos';
% source.brainstructure = ones(1, size(pnt,1));
% source.brainstructurelabel = {'CORTEX'};
source.brainstructure = [1 1 1 1 1 1 2 2 2 2 2 2];
source.brainstructurelabel = {'CORTEX_LEFT', 'CORTEX_RIGHT'};

cfg = [];
cfg.filetype  = 'cifti';
cfg.parameter = 'pow';
cfg.filename  = 'test_bug2096';
ft_sourcewrite(cfg, source);

source1 = ft_read_cifti('test_bug2096.pow.dscalar.nii');

cfg = [];
cfg.filetype  = 'cifti';
cfg.parameter = 'pow';
cfg.filename  = 'test_bug2096b';
ft_sourcewrite(cfg, source1);

source2 = ft_read_cifti('test_bug2096b.pow.dscalar.nii');

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
cfg.filename  = 'test_bug2096.imagcoh';
ft_sourcewrite(cfg, source);

source1 = ft_read_cifti('test_bug2096.imagcoh.dconn.nii');

cfg = [];
cfg.filetype  = 'cifti';
cfg.parameter = 'imagcoh';
cfg.filename  = 'test_bug2096b.imagcoh';
ft_sourcewrite(cfg, source1);

source2 = ft_read_cifti('test_bug2096b.imagcoh.dconn.nii');

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
cfg.filename = 'test_bug2096.timeseries';
cfg.precision = 'single';
ft_sourcewrite(cfg, source);

source1 = ft_read_cifti('test_bug2096.timeseries.dtseries.nii');

%%
try
  parcellation = ft_read_atlas(dccnpath('/home/common/matlab/fieldtrip/template/atlas/aal/ROI_MNI_V4.nii'));
catch
  parcellation = ft_read_atlas(fullfile(getenv('HOME'), '/matlab/fieldtrip/template/atlas/aal/ROI_MNI_V4.nii'));
end
source = ft_checkdata(parcellation, 'datatype', 'source');
source = removefields(source, {'tissue', 'tissuelabel'});
source.pow = randn(prod(parcellation.dim), 1);
source.pow(parcellation.tissue==0) = nan;
source.dimord = 'pos';

cfg = [];
cfg.funparameter = 'pow';
ft_sourceplot(cfg, source);

cfg = [];
cfg.filetype  = 'cifti';
cfg.parameter = 'pow';
cfg.filename  = 'test_bug2096.pow';
ft_sourcewrite(cfg, source);

source1 = ft_read_cifti('test_bug2096.pow.dscalar.nii');

cfg = [];
cfg.filetype  = 'cifti';
cfg.parameter = 'pow';
cfg.filename  = 'test_bug2096b.pow';
ft_sourcewrite(cfg, source1);

source2 = ft_read_cifti('test_bug2096b.pow.dscalar.nii');

cfg = [];
cfg.parameter = 'pow';
sourcep = ft_sourceparcellate(cfg, source, parcellation);

cfg = [];
cfg.funparameter = 'pow';
ft_sourceplot(cfg, sourcep);

assert(isfield(sourcep, 'pow'));
assert(isfield(sourcep, 'brainordinate'));

cfg = [];
cfg.filetype  = 'cifti';
cfg.parameter = 'pow';
cfg.parcellation = 'tissue';
cfg.filename  = 'test_bug2096.pow';
ft_sourcewrite(cfg, sourcep);

sourcep1 = ft_read_cifti('test_bug2096.pow.pscalar.nii');

cfg = [];
cfg.filetype  = 'cifti';
cfg.parameter = 'pow';
cfg.filename  = 'test_bug2096b.pow';
ft_sourcewrite(cfg, sourcep1);

sourcep2 = ft_read_cifti('test_bug2096b.pow.pscalar.nii');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% specific tests for dscalar, dtseries, dconn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all

%% test the dscalar output
[pnt, tri] = icosahedron;
pntL = pnt; pntL(:,1) = pntL(:,1) - 1; % shift along X
pntR = pnt; pntR(:,1) = pntR(:,1) + 1; % shift along X

source = [];
source.pos    = [pntL; pntR];
source.tri    = [tri; tri+12];
source.brainstructure = [1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2];
source.brainstructurelabel = {'CORTEX_LEFT', 'CORTEX_RIGHT'};
source.activity = zeros(size(source.pos,1), 1);
source.activity(:,1) = 1:24;
source.dimord = 'pos';

cfg = [];
cfg.filetype  = 'cifti';
cfg.parameter = 'activity';
cfg.filename  = 'test_bug2096.activity';
ft_sourcewrite(cfg, source);

source1 = ft_read_cifti('test_bug2096.activity.dscalar.nii');
ft_plot_mesh(source1, 'vertexcolor', source1.activity(:,1), 'edgecolor', 'none')

%% test the dtsetries output
[pnt, tri] = icosahedron;
pntL = pnt; pntL(:,1) = pntL(:,1) - 1; % shift along X
pntR = pnt; pntR(:,1) = pntR(:,1) + 1; % shift along X

source = [];
source.pos    = [pntL; pntR];
source.tri    = [tri; tri+12];
source.brainstructure = [1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2];
source.brainstructurelabel = {'CORTEX_LEFT', 'CORTEX_RIGHT'};
source.time = 1:10;
source.timeseries = zeros(size(source.pos,1), length(source.time));
for i=1:size(source.timeseries,2)
  source.timeseries(:,i) = 1:24;
end
source.dimord = 'pos_time';

cfg = [];
cfg.filetype  = 'cifti';
cfg.parameter = 'timeseries';
cfg.filename  = 'test_bug2096.timeseries';
ft_sourcewrite(cfg, source);

source1 = ft_read_cifti('test_bug2096.timeseries.dtseries.nii');
ft_plot_mesh(source1, 'vertexcolor', source1.timeseries(:,1), 'edgecolor', 'none')

%% test the dconn output
[pnt, tri] = icosahedron;
pntL = pnt; pntL(:,1) = pntL(:,1) - 1; % shift along X
pntR = pnt; pntR(:,1) = pntR(:,1) + 1; % shift along X

source = [];
source.pos    = [pntL; pntR];
source.tri    = [tri; tri+12];
source.brainstructure = [1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2];
source.brainstructurelabel = {'CORTEX_LEFT', 'CORTEX_RIGHT'};
source.imagcoh = rand(24,24);
source.dimord = 'pos_pos';

cfg = [];
cfg.filetype  = 'cifti';
cfg.parameter = 'imagcoh';
cfg.filename  = 'test_bug2096.imagcoh';
ft_sourcewrite(cfg, source);

source1 = ft_read_cifti('test_bug2096.imagcoh.dconn.nii');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% specific tests for pscalar, ptseries, pconn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all

%% test the pscalar output
[pnt, tri] = icosahedron;
pntL = pnt; pntL(:,1) = pntL(:,1) - 1; % shift along X
pntR = pnt; pntR(:,1) = pntR(:,1) + 1; % shift along X

source = [];
source.brainordinate.pos    = [pntL; pntR];
source.brainordinate.tri    = [tri; tri+12];
source.brainordinate.brainstructure = [1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2];
source.brainordinate.brainstructurelabel = {'CORTEX_LEFT', 'CORTEX_RIGHT'};
source.brainordinate.parcellation =   [1 1 1 2 2 2 3 3 3 4 4 4 4 4 4 5 5 5 6 6 6 7 7 7];
source.brainordinate.parcellationlabel = {'parcel1' 'parcel2' 'parcel3' 'parcel4' 'parcel5' 'parcel6' 'parcel7'};
source.activity = zeros(length(source.brainordinate.parcellationlabel) , 1);
source.activity(:,1) = 1:7;
source.dimord = 'chan';
source.label  = source.brainordinate.parcellationlabel;

figure
ft_plot_mesh(source.brainordinate, 'vertexcolor', source.brainordinate.parcellation(:), 'edgecolor', 'none')

cfg = [];
cfg.filetype  = 'cifti';
cfg.parameter = 'activity';
cfg.filename  = 'test_bug2096.activity';
ft_sourcewrite(cfg, source);

source1 = ft_read_cifti('test_bug2096.activity.pscalar.nii');

figure
ft_plot_mesh(source1.brainordinate, 'vertexcolor', source1.brainordinate.parcellation(:), 'edgecolor', 'none')

%% test the ptsetries output
[pnt, tri] = icosahedron;
pntL = pnt; pntL(:,1) = pntL(:,1) - 1; % shift along X
pntR = pnt; pntR(:,1) = pntR(:,1) + 1; % shift along X

source = [];
source.brainordinate.pos    = [pntL; pntR];
source.brainordinate.tri    = [tri; tri+12];
source.brainordinate.brainstructure = [1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2];
source.brainordinate.brainstructurelabel = {'CORTEX_LEFT', 'CORTEX_RIGHT'};
source.brainordinate.parcellation =   [1 1 1 2 2 2 3 3 3 4 4 4 4 4 4 5 5 5 6 6 6 7 7 7];
source.brainordinate.parcellationlabel = {'parcel1' 'parcel2' 'parcel3' 'parcel4' 'parcel5' 'parcel6' 'parcel7'};
source.time = 1:10;
source.timeseries = zeros(7, length(source.time));
for i=1:size(source.timeseries,2)
  source.timeseries(:,i) = 1:7;
end
source.dimord = 'chan_time';
source.label  = source.brainordinate.parcellationlabel;

figure
ft_plot_mesh(source.brainordinate, 'vertexcolor', source.brainordinate.parcellation(:), 'edgecolor', 'none')

cfg = [];
cfg.filetype  = 'cifti';
cfg.parameter = 'timeseries';
cfg.filename  = 'test_bug2096.timeseries';
ft_sourcewrite(cfg, source);

source1 = ft_read_cifti('test_bug2096.timeseries.ptseries.nii');

figure
ft_plot_mesh(source1.brainordinate, 'vertexcolor', source1.brainordinate.parcellation, 'edgecolor', 'none')

%% test the pconn output
[pnt, tri] = icosahedron;
pntL = pnt; pntL(:,1) = pntL(:,1) - 1; % shift along X
pntR = pnt; pntR(:,1) = pntR(:,1) + 1; % shift along X

source = [];
source.brainordinate.pos    = [pntL; pntR];
source.brainordinate.tri    = [tri; tri+12];
source.brainordinate.brainstructure = [1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2];
source.brainordinate.brainstructurelabel = {'CORTEX_LEFT', 'CORTEX_RIGHT'};
source.brainordinate.parcellation =   [1 1 1 2 2 2 3 3 3 4 4 4 4 4 4 5 5 5 6 6 6 7 7 7];
source.brainordinate.parcellationlabel = {'parcel1' 'parcel2' 'parcel3' 'parcel4' 'parcel5' 'parcel6' 'parcel7'};
source.imagcoh = rand(7,7);
source.dimord = 'chan_chan';
source.label  = source.brainordinate.parcellationlabel;

figure
ft_plot_mesh(source.brainordinate, 'vertexcolor', source.brainordinate.parcellation(:), 'edgecolor', 'none')

cfg = [];
cfg.filetype  = 'cifti';
cfg.parameter = 'imagcoh';
cfg.filename  = 'test_bug2096.imagcoh';
ft_sourcewrite(cfg, source);

source1 = ft_read_cifti('test_bug2096.imagcoh.pconn.nii');

figure
ft_plot_mesh(source1.brainordinate, 'vertexcolor', source1.brainordinate.parcellation, 'edgecolor', 'none')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% specific tests reading from external files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all

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
cii4 = ft_read_cifti('Conte69.MyelinAndCorrThickness.32k_fs_LR.ptseries.nii');
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

% the ones starting with x177746 have been confirmed by inspecting the XML section of the files
% the files contain <MapName>177746_aparc</MapName>
% the x in front is needed to make it a valid fieldname and results from "fixname"
datafield = {
  'arealdistortion'
  'x177746_aparc'
  'x177746_ba'
  'x177746_aparc_a2009s'
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
  source = ft_read_cifti(filename{i});
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
  source = ft_read_cifti(filename{i});
end


