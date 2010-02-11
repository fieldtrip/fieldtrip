Sparse Logistic Regression ToolBox ver1.0beta 
(Updated on 2009/07/20 by Okito Yamashita)

*** Please see README.pdf for more detailed explanation ****

The Sparse Logistic Regression toolbox (SLR toolbox hereafter) is a suite of MATLAB functions for solving classification problems. It provides one of solutions for binary or multi-class classification problem. The unique feature is parameters of the classifier are learned in a sparse way, resulting in automatic feature selection while learning weight parameters in the classifier. This feature may be appropriate to classification problems of neuro-imaging data where only a limited number of training data (from several tens to several hundreds) can be used to classify a rather high dimensional feature vector (over one thousand). Furthermore SLR releases users from the time-consuming feature selection task preceding classification although results may be suboptimal.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SLR toolbox is a suite of MATLAB functions and scripts. MATLAB, a commercial engineering mathematics package, is required to use SLR toolbox. A couple of functions require the optimization toolbox (see below). Codes in the toolbox were written for MATLAB ver7.0.1 or later under UNIX. This toolbox has originally been developed by Okito Yamashita in ATR Computational Neuroscience laboratories for personal use.

To get installed the toolbox, you just download and unzip the file (SLR1.0beta.zip) wherever you like. You may also download two test-data acquired from two real experiments in order to see how binary and multi-class classification problems are solved (this is optional). Please start from demo functions ('demo_*.m') to learn how the functions in SLR toolbox works.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Referencing the toolbox
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
When using this tool for a paper please refer to the following paper: 

Yamashita O, Sato MA, Yoshioka T, Tong F, Kamitani Y (2008).
Sparse estimation automatically selects voxels relevant for the decoding of fMRI activity patterns. Neuroimage. Oct 1;42(4):1414-29.

The above manuscript contains basics of SLR (SLR-LAP) and applications to fMRI decoding.
 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Feedback & Bug report  
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Any feedback and bug report are welcome. Please keep contact with me (oyamashi@atr.jp). I would like to respond as quickly as possible.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Licencing & Copy Right
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SLR toolbox is free but copyright software, distributed under the terms of the GNU General Public Licence as published by the Free Software Foundation. Further details on "copyleft" can be found at http://www.gnu.org/copyleft/. No formal support or maintenance is provided or implied.

-------------------------------
 Acknowledgements
-------------------------------
This toolbox is brought to you by ATR Computational Neuroscience laboratories in Kyoto. This research was supported in part by the NICT, Honda Research Institute, the SCOPE, SOUMU, the Nissan Science Foundation, and grants from the National Eye Institute to FT (R01 EY017082 and R01 EY14202). 

