To install MARS, unzip the "mars" folder to SPM8_DIR/toolbox/

To get started, you can do either of the followings:

(1) start SPM8 GUI, click the "Batch" button, and from there, go to SPM-->Tools-->MARS. Once you're in MARS GUI, you can config all input arguments to start using MARS.

(2) use the script start_mars.m, below are some examples:

start_mars
run MARS with all options default, the program will ask you to choose input T1 image(s)

start_mars('/home/user/segProject/img1-t1.nii')
segment T1 image 'img1-t1.nii', using all default options

start_mars('/home/user/segProject/img1-t1.nii','/home/user/segProject/img1-t2.nii')
segment T1 image 'img1-t1.nii' together with T2 'img1-t2.nii', using all default options

start_mars('/home/user/segProject/img1-t1.nii',[],[],[],[],'New Segment with MRF cleanup')
segment T1 image 'img1-t1.nii', using all default options, but only run New Segment with SPM8's MRF-cleanup (no MARS iterations)

start_mars('/home/user/segProject/img1-t1.nii',[],[],[],[],'New Segment')
segment T1 image 'img1-t1.nii', using all default options, but only run New Segment (no MARS iterations, no SPM8's MRF-cleanup)

To get help, you can type "help start_mars" in the command window, or go to MARS GUI for complete help text under each option.

Info about the prior maps (TPM/TCM):

TPM_BW20.nii: TPM built from the BrainWeb 20-subject dataset.

lTCM_BW20.mat: local TCM built from the BrainWeb 20-subject dataset.

gTCM_BW20_S1.mat: global TCM trained from the first subject in BrainWeb 20-subject dataset.

rTCM_BW20_S1.mat: regional TCM from gTCM_BW20_S1.mat.

eTPM.nii: TPM built from 26 real human heads by Dr. Chris Rorden (see http://neuralengr.com/segment/).

gTCM_H1.mat: global TCM trained from one real human head.

rTCM_H1.mat: regional TCM from gTCM_H1.mat.

The default TPM and TCM of MARS is eTPM.nii and rTCM_BW20_S1.mat, respectively.

Note that there is no local TCM built from real human heads, simply because we did not have large enough prior dataset that includes at least 20 human head MRIs with each pixel segmented out.

Note:

(1) MARS was not designed to segment heads with lesions.

(2) Segmentation results from MARS are continuous probability maps. If you want to build volume conductor model out of these results, you need to use any meshing software that is able to generate finite element mesh from non-binary segmentation data; otherwise you have to binarize the segmentation to proceed. If you have to binarize MARS results, I recommend you to use my previous toolbox (mysegment.m) for that purpose, which is available at http://neuralengr.com/segment/

(3) MARS has been tested on the following OS: Ubuntu 14.04 64-bit; Windows 7 64-bit; OS X 10.10 Yosemite 64-bit.

(4) MARS has not been made compatible to SPM12, now can only run under SPM8.

Yu (Andy) Huang, Jan 2016
