function test_bug1295

% MEM 1gb
% WALLTIME 00:03:39

% first show the issue: read in dicom test data with the old and new code

filename = '/home/common/matlab/fieldtrip/data/test/latest/mri/dicom/ERIVDBER_030731_R.OOSTERVELD.MR.PAUGAA_ANATOMICAL-3D.2.99.2003.7.31.11.19.16.15000.53831976.IMA';

mri1 = ft_read_mri(filename, 'format', 'dicom_old');
mri2 = ft_read_mri(filename, 'format', 'dicom');

% realign
cfg = [];
cfg.fiducial.nas = [107  88  33];
cfg.fiducial.lpa = [118 146 115];
cfg.fiducial.rpa = [124  27 124];
tmp1 = ft_volumerealign(cfg, mri1);

% realign
cfg = [];
cfg.fiducial.nas = [ 88 107  33];
cfg.fiducial.lpa = [149 119 116];
cfg.fiducial.rpa = [ 24 122 119];
tmp2 = ft_volumerealign(cfg, mri2);

% check the handedness of the transformation matrices
T1 = tmp1.transform;
T2 = tmp2.transform;

dir1 = cross(T1(1,1:3),T1(2,1:3));
dir2 = cross(T2(1,1:3),T2(2,1:3));

isrighthanded1 = sign(dir1*T1(3,1:3)')==1
isrighthanded2 = sign(dir2*T2(3,1:3)')==1

% now do the realignment with an additional zpoint
% this is new functionality as of feb 2012
cfg = [];
cfg.fiducial.nas = [107  88  33];
cfg.fiducial.lpa = [118 146 115];
cfg.fiducial.rpa = [124  27 124];
cfg.fiducial.zpoint = [43 88 107];

tmp1b = ft_volumerealign(cfg, mri1);
T1b   = tmp1b.transform;
dir1b = cross(T1b(1,1:3),T1b(2,1:3));

isrighthanded1b = sign(dir1b*T1b(3,1:3)')==1




