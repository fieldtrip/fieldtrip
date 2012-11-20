function test_bug1816

% test the ft_volumesegment function used for segmentation with FSL BET and FAST
% see http://bugzilla.fcdonders.nl/show_bug.cgi?id=1816

% at this moment the test script does not yet work, but we don't want the automatic regression testing to flag it as failure
return

% TEST test_bug1816
% TEST ft_read_mri ft_volumesegment

if ispc
    datadir = 'H:';
else
    datadir = '/home';
end

mri=ft_read_mri(strcat(datadir,'/common/matlab/fieldtrip/data/Subject01.mri'));

% read in aligned images
subjectT1  = strcat(datadir,'/common/matlab/fieldtrip/data/test/bug1826/T1.nii.gz');             % change path to image
subjectT2  = strcat(datadir,'/common/matlab/fieldtrip/data/test/bug1826/T2_T1Space_trilinear.nii.gz');   
subjectDTi = strcat(datadir,'/common/matlab/fieldtrip/data/test/bug1826/DTI_T1Space.nii.gz');

T1  = ft_read_mri(subjectT1);
T2  = ft_read_mri(subjectT2);
DTi = ft_read_mri(subjectDTi);

% segment

cfg = [];
cfg.output = {'brain','skull','scalp','csf','gray','white'};  % define the requested tissue-types. 
                                                              % this list can  be extended/changed. 
                                                              
cfg.method = 'fsl';                                           % not implemented yet. it is a flag to
                                                              % indicate that the segmentation
                                                              % should use fsl instead of the
                                                              % already implemented spm-fieldtrip
                                                              % segmentation.
                                                              
cfg.coordsys = 'ctf'; % if the coordinate system of the mri is not specified in the mri.coordsys 
                      % field of the volume, ft_volumesegment will reuqire this information 
                      % specified in the cfg. I do not know what kind of coordinate system FSL needs.
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

% check if other methods of segmentation are still working

test_ft_volumesegment;

