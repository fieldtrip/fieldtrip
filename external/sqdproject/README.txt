Date: 14 February 2010
Version:  3.0 beta
Creator: Shantanu Ray, ISR, University of Maryland
Maintained by: the University of Maryland lab of Jonathan Z. Simon, jzsimon@umd.edu
Most Recent Version location: http://www.isr.umd.edu/Labs/CSSL/simonlab/resources.html

Contents:
Information Files:
 README.txt: this file
Basic Usage Files:
 sqdread.m: Matlab script to read sqd files, both MEG data and other information
 sqdwrite.m: Matlab script to write MEG data to sqd files
 sqdgettriggers.m: Detects triggers needed to epoch data
 sqdmakeepochs.m: Breaks data into epochs (and/or averages) & saves to Matlab data files.
 writechanloc.m: Matlab script to write sensor location for other uses
Demo Files:
 trigger_and_epoching_examples.m: how to epoch your MEG data within Matlab instead of MEG160.
 sqddemo.m: Demonstration script for sqdread,sqdwrite
 testsqd.sqd: Demo sqd file; Used in examples in the help
 sqdtestdata.mat: Demo matlab data file;  Used in examples in the help
Other Files and Folders not typically called by the user:
 @acqparam: Object which saves the acquisition paramters
 @chanhandle: Object which contains data for each channel
 @sqdhandle: Object which serves as a information container for sqd info

Installation instructions:
1. Remove files of the previous toolbox from MATLAB path. This is far more important than you might think.
2. Please copy all of the above files and folders to a place you can add to your MATLAB path

For instructions on how to use this toolbox:
1. In MATLAB,sqddemo for demonstration on how to use sqdread and sqdwrite.
2. or	help sqdread
	help sqdwrite
	help writechanloc
	
	
Changes incorporated in 3.0 beta (2/14/2010)
* New file trigger_and_epoching_examples.m to demonstrate new trigger and epoching features.
* Additions and changes to sqdmakeepochs.m

Changes incorporated in 2.9.9 beta (2/10/2010)
* Really fixed bug in gain multiplier. The data used to have the wrong amplitude if the  
  OutputGain was not 100 (the factor off was (OutputGain/100)^2). [The incorrect factor
  in sqdread was reversed in sqdwrite, so using the two serially, e.g. for denoising,
  had no impact, even before the bug was fixed.]
* Fixed bug in reading epoched files ("acquistion type 3"). Previously only the first
  epoch was read in. Now all epochs are read in (or fewer, if fewer samples are requested).
  If all epochs (trials) are requested (the default), then the data is returned in a 3 
  dimensional matrix, in the form sample x channel x epoch. In this form, the mean 
  across epochs can be calculated using, e.g. mean(data,3). If fewer samples are requested,
  the epoched data is concatenated, giving the usual 2 dimensional data matrix of the form
  samples x channel.
* Added new function file sqdgettriggers.m. This function automatically detects triggers so that
  the data can be epoched.
* Added new function file sqdmakeepochs.m. This function uses the detected triggers to epoch (and
  optinally, average) the data and save the epochs in Matlab data files.

Changes incorporated in 2.2 beta (12/11/2008)
* Bugs in sqdwrite when overwriting data fixed.

Changes incorporated in 2.1 beta and earlier
* Faster sqdread and getdata
* Fixed Bug in getdata - erroneous rawdata to double conversion for trigger files
* Corrected the help for sqdread and putdata
* Faster sqdwrite and putdata
* Bug in gain multiplier. There was an error in the gain-factor in the 
   conversion formula that was being used to convert the data from raw format
   to double.  The data might have wrong amplitudes if the values of 
   Amplifier-InputGain and Amplifier-OutputGain were anything
   other than 2.00 and 100.00 respectively. This bug has been duly corrected.
   *Note, this bug was not fixed correctly until version 3.0 beta*
* SQDREAD  now also returns:
   a. Patient Information
   b. Sensor Location
* New function WRITECHANLOC has been added which can write the sensor location
   to a file to be used in conjunction with EEGLAB toolbox in MATLAB
* Working examples have been included in the MATLAB help
* The object structure has been changed with the addition of a new CHANHANDLE class

Bugs still open:
Not all information contained in the "info" object is available to the user.
