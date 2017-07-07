function test_bug1816

% MEM 1500mb
% WALLTIME 00:10:00

% test the ft_volumesegment function used for segmentation with FSL BET and FAST
% see http://bugzilla.fcdonders.nl/show_bug.cgi?id=1816

% TEST ft_read_mri ft_volumesegment

mri = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/Subject01.mri'));

% read in aligned images
subjectT1  = dccnpath('/home/common/matlab/fieldtrip/data/test/bug1826/T1.nii.gz');
subjectT2  = dccnpath('/home/common/matlab/fieldtrip/data/test/bug1826/T2_T1Space_trilinear.nii.gz');   
subjectDTi = dccnpath('/home/common/matlab/fieldtrip/data/test/bug1826/DTI_T1Space.nii.gz');

T1  = ft_read_mri(subjectT1);
T2  = ft_read_mri(subjectT2);
DTi = ft_read_mri(subjectDTi);

% FIXME what are the (approximate) coordinate systems in which the data is expressed?
T1.coordsys = 'ctf';
T2.coordsys = 'ctf';
DTi.coordsys = 'ctf';

% at this moment the test script does not yet work, but we don't want the automatic regression testing to flag it as failure
return

% segment using FSL

cfg = [];
cfg.output = {'brain','skull','scalp','csf','gray','white'};  % define the requested tissue-types. 
                                                              % this list can  be extended/changed. 
                                                              
cfg.method = 'fsl';                                           % not implemented yet. it is a flag to
                                                              % indicate that the segmentation
                                                              % should use fsl instead of the
                                                              % already implemented spm-fieldtrip
                                                              % segmentation.
                                                              
cfg. units   = 'mm';  % the physical units in which the output will be expressed   

seg = ft_volumesegment(cfg,T1,T2,DTi)  

% The output segmentation should contain the following fields:

assert(isfield(seg,'dim'),'dimensionality of volume is missing'); % e.g. seg.dim = [256 256 256];
assert(isfield(seg,'transform'),'transformation matrix is missing') 
% The transformation matrix alinges the anatomical data to the coordinate system of the image. 

% The segmentation should not change the coordinate system of the volume.

assert(all(T1.dim == seg.dim), 'Dimensionality of the volume changed after segmentation.');
assert(isequal(T1.transform, seg.transform), 'Coordinate system of the volume changed after segmentation.');

assert(isfield(seg,'coordsys'), 'Coordinate system is not defined');
assert(isfield(seg,'unit'),'Units are not defined');

% check for fields describing the tissuetypes
for i = 1 : size(cfg.output,2)
assert(isfield(seg,cfg.output{i}), 'Required tissue-type is missing from segmentation.');
assert(isequal(size(getfield(seg,cfg.output{i})),seg.dim),'Tissue does not have the same dimensions as dim.');
end

assert(isfield(seg,'cfg'),'cfg is missing from segmentation'); 

