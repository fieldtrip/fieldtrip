function start_mars(T1,T2,TPM,TCM,numGaussian,runningMode,beta,convergence,saveBiasField,saveNative,saveWarp,saveDeformField)
% start_mars(T1,T2,TPM,TCM,numGaussian,runningMode,beta,convergence,saveBiasField,saveNative,saveWarp,saveDeformField)
%
% Function interface to start the toolbox
% 'MARS (Morphologically and Anatomically accuRate Segmentation)'
%
% It is recommended to use the graphic user interface (GUI) of the MARS. If
% you have a lot of images needed to be segmented, it will be easier to use
% this function interface. However, please note that this function cannot
% control all the options of MARS. Most code in this function is generated
% from the MARS GUI. For complete control of MARS options, one is recommended
% to use the GUI to config all the options first and then generate the code
% from there and change this function accordingly.
% 
% To start the MARS GUI, start SPM8 GUI, click the "Batch" button, and from
% there, go to SPM-->Tools-->MARS. If you do not see "MARS" under Tools,
% please put the MARS folder under SPM8_DIR/toolbox/
%
% If you have a bunch of images to segment, just put the file names
% vertically into T1 and/or T2 using strvcat(), and make sure each
% T2 corresponds to its corresponding T1, e.g.:
% T1 = strvcat('/home/user/segProject/img1-t1.nii','/home/user/segProject/img2-t1.nii');
% T2 = strvcat('/home/user/segProject/img1-t2.nii','/home/user/segProject/img2-t2.nii');
% The program will process them one by one. However, all the images in the
% list will be segmented by the same TPM and/or TCM. If you want them to be
% processed by different TPM/TCM, call this function for each image separately,
% each time using a different TPM/TCM.
%
% Examples:
%
% start_mars
% run MARS with all options default, the program will ask you to choose
% input T1 image(s)
%
% start_mars('/home/user/segProject/img1-t1.nii')
% segment T1 image 'img1-t1.nii', using all default options
%
% start_mars('/home/user/segProject/img1-t1.nii','/home/user/segProject/img1-t2.nii')
% segment T1 image 'img1-t1.nii' together with T2 'img1-t2.nii', using all default options
%
% start_mars('/home/user/segProject/img1-t1.nii',[],[],[],[],'New Segment with MRF cleanup')
% segment T1 image 'img1-t1.nii', using all default options, but only run
% New Segment with SPM8's MRF-cleanup (no MARS iterations)
%
% start_mars('/home/user/segProject/img1-t1.nii',[],[],[],[],'New Segment')
% segment T1 image 'img1-t1.nii', using all default options, but only run
% New Segment (no MARS iterations, no SPM8's MRF-cleanup)
%
% See the help texts in the GUI for all the info related to the input arguments.
%
% Adapted from Chris Rorden, 2011.11
% Yu (Andy) Huang, 2015.06
% Neural Engineering Lab, Dept. of Biomedical Engineering, City College of New York
% yhuang16@citymail.cuny.edu

if nargin < 1 || isempty(T1) % no input images
    T1 = spm_select(inf,'image','Select images for MARS');
end

if nargin < 2 || isempty(T2) % no T2 specified
    T2 = [];
end

if nargin < 3 || isempty(TPM) % no TPM specified
    TPM = fullfile(spm('Dir'),'toolbox','mars','eTPM.nii'); % default TPM for MARS
end

if nargin < 4 || isempty(TCM) % no TCM specified
    TCM = fullfile(spm('Dir'),'toolbox','mars','rTCM_BW20_S1.mat'); % default TCM for MARS
end

if nargin < 5 || isempty(numGaussian) % not specified number of clusters in the Gaussian mixure
    numGaussian = [2 2 2 3 4 2];
end

if nargin < 6 || isempty(runningMode) % running mode not specified
    mode = 2;
else
    switch runningMode
        case 'MARS'
            mode = 2;
        case 'New Segment with MRF cleanup'
            mode = 1;
        case 'New Segment'
            mode = 0;
        otherwise
            error('unrecognized option for MARS running mode! Please specify one of the followings: MARS; New Segment with MRF cleanup; New Segment');
    end
end

if nargin < 7 || isempty(beta) % beta not specified
    beta = 0.1;
end

if nargin < 8 || isempty(convergence) % convergence threshold for MARS not specified
    convergence = 0.1;
end

if nargin < 9 || isempty(saveBiasField) % save bias field or not
    bf = [0 0];
else
    switch saveBiasField
        case 'Save Nothing'
            bf = [0 0];
        case 'Save Bias Corrected'
            bf = [0 1];
        case 'Save Bias Field'
            bf = [1 0];
        case 'Save Field and Corrected'
            bf = [1 1];
        otherwise
            error('unrecognized option for bias field! Please specify one of the followings: Save Nothing; Save Bias Corrected; Save Bias Field; Save Field and Corrected');
    end
end

if nargin < 10 || isempty(saveNative) % save segmentation in native space or not
    n = [1 0];
else
    switch saveNative
        case 'None'
            n = [0 0];
        case 'Native Space'
            n = [1 0];
        case 'DARTEL Imported'
            n = [0 1];
        case 'Native + DARTEL Imported'
            n = [1 1];
        otherwise
            error('unrecognized option for saving segmentation results! Please specify one of the followings: None; Native Space; DARTEL Imported; Native + DARTEL Imported');
    end
end

if nargin < 11 || isempty(saveWarp) % save segmentation in warped/normalized space or not
    w = [0 0];
else
    switch saveWarp
        case 'None'
            w = [0 0];
        case 'Modulated'
            w = [0 1];
        case 'Unmodulated'
            w = [1 0];
        case 'Modulated + Unmodulated'
            w = [1 1];
        otherwise
            error('unrecognized option for saving warped segmentation results! Please specify one of the followings: None; Modulated; Unmodulated; Modulated + Unmodulated');
    end
end

if nargin < 12 || isempty(saveDeformField) % save deformation field or not
    df = [0 0];
else
    switch saveDeformField
        case 'None'
            df = [0 0];
        case 'Inverse'
            df = [1 0];
        case 'Forward'
            df = [0 1];
        case 'Inverse + Forward'
            df = [1 1];
        otherwise
            error('unrecognized option for deformation field! Please specify one of the followings: None; Inverse; Forward; Inverse + Forward');
    end
end

if isempty(T1)
    error('No image selected. Please select an image to segment!');
end

[ptht,namt,extt] = spm_fileparts(deblank(TPM(1,:)));
Tem = [ptht,filesep,namt,extt];
spm_jobman('initcfg');

for i=1:size(T1,1)
    
    ref = deblank(T1(i,:));
    [pth,nam,ext] = spm_fileparts(ref);
    
    matlabbatch{1}.spm.tools.mars.channel(1).vols = {ref};  % image to be segmented
    matlabbatch{1}.spm.tools.mars.channel(1).biasreg = 0.0001; % P(beta)
    matlabbatch{1}.spm.tools.mars.channel(1).biasfwhm = 60;
    matlabbatch{1}.spm.tools.mars.channel(1).write = bf;
    
    if ~isempty(T2) %T2 specified
        ref2 = deblank(T2(i,:));
        [pth2,nam2,ext2] = spm_fileparts(ref2);
        if ~isempty(TCM), fprintf('Using %s and %s to segment T1 and T2 images: %s %s\n', TPM, TCM, ref, ref2);
        else fprintf('Using %s to segment T1 and T2 images: %s %s\n', TPM, ref, ref2);
        end
        
        matlabbatch{1}.spm.tools.mars.channel(2).vols = {ref2}; % the 2nd image aiding segmentation
        matlabbatch{1}.spm.tools.mars.channel(2).biasreg = 0.0001;
        matlabbatch{1}.spm.tools.mars.channel(2).biasfwhm = 60;
        matlabbatch{1}.spm.tools.mars.channel(2).write = bf;
    else
        if ~isempty(TCM), fprintf('Using %s and %s to segment T1 %s\n', TPM, TCM, ref);
        else fprintf('Using %s to segment T1 %s\n', TPM, ref);
        end
    end
    
    
    matlabbatch{1}.spm.tools.mars.tissue(1).tpm = {[Tem,',1']}; % TPM tissue #1
    matlabbatch{1}.spm.tools.mars.tissue(1).ngaus = numGaussian(1); % number of Gaussians to model this tissue
    matlabbatch{1}.spm.tools.mars.tissue(1).native = n;
    matlabbatch{1}.spm.tools.mars.tissue(1).warped = w;
    matlabbatch{1}.spm.tools.mars.tissue(2).tpm = {[Tem,',2']}; % TPM tissue #2
    matlabbatch{1}.spm.tools.mars.tissue(2).ngaus = numGaussian(2);
    matlabbatch{1}.spm.tools.mars.tissue(2).native = n;
    matlabbatch{1}.spm.tools.mars.tissue(2).warped = w;
    matlabbatch{1}.spm.tools.mars.tissue(3).tpm = {[Tem,',3']}; % TPM tissue #3
    matlabbatch{1}.spm.tools.mars.tissue(3).ngaus = numGaussian(3);
    matlabbatch{1}.spm.tools.mars.tissue(3).native = n;
    matlabbatch{1}.spm.tools.mars.tissue(3).warped = w;
    matlabbatch{1}.spm.tools.mars.tissue(4).tpm = {[Tem,',4']}; % TPM tissue #4
    matlabbatch{1}.spm.tools.mars.tissue(4).ngaus = numGaussian(4);
    matlabbatch{1}.spm.tools.mars.tissue(4).native = n;
    matlabbatch{1}.spm.tools.mars.tissue(4).warped = w;
    matlabbatch{1}.spm.tools.mars.tissue(5).tpm = {[Tem,',5']}; % TPM tissue #5
    matlabbatch{1}.spm.tools.mars.tissue(5).ngaus = numGaussian(5);
    matlabbatch{1}.spm.tools.mars.tissue(5).native = n;
    matlabbatch{1}.spm.tools.mars.tissue(5).warped = w;
    matlabbatch{1}.spm.tools.mars.tissue(6).tpm = {[Tem,',6']}; % TPM tissue #6
    matlabbatch{1}.spm.tools.mars.tissue(6).ngaus = numGaussian(6);
    matlabbatch{1}.spm.tools.mars.tissue(6).native = n;
    matlabbatch{1}.spm.tools.mars.tissue(6).warped = w;
    
    matlabbatch{1}.spm.tools.mars.warp.reg = 4; % P(alpha)
    matlabbatch{1}.spm.tools.mars.warp.affreg = 'mni'; % template type for initial affine registration. Important!
    matlabbatch{1}.spm.tools.mars.warp.samp = 3; % downsample distance for Unified Segmentation Algorithm
    matlabbatch{1}.spm.tools.mars.warp.write = df;
    
    matlabbatch{1}.spm.tools.mars.marsOptions.marsMode = mode; % run which algorithm?
    matlabbatch{1}.spm.tools.mars.marsOptions.tcm = {TCM}; % TCM (global, local or regional)
    matlabbatch{1}.spm.tools.mars.marsOptions.beta = beta; % Beta for adjusting the strength of the local TCM
    matlabbatch{1}.spm.tools.mars.marsOptions.marsConvergence = convergence; % User-defined convergence threshold for MARS iterations
    
    
    spm_jobman('run',matlabbatch);
end % for each image
