% FieldTrip is a Matlab toolbox for MEG and EEG analysis that is being
% developed at the Centre for Cognitive Neuroimaging of the Donders
% Institute for Brain, Cognition and Behaviour together with collaborating
% institutes. The development of FieldTrip is supported by funding from the
% BrainGain consortium. The FieldTrip software is released as open source
% under the GNU general public license.
%
% The toolbox includes algorithms for simple and advanced analysis of MEG
% and EEG data, such as time-frequency analysis, source reconstruction
% using dipoles, distributed sources and beamformers and non-parametric
% statistical testing. It supports the data formats of all major MEG
% systems (CTF, Neuromag, BTi) and of the most popular EEG systems, and new
% formats can be added easily. FieldTrip contains high-level functions that
% you can use to construct your own analysis protocol in Matlab.
% Furthermore, it easily allows developers to incorporate low-level
% algorithms for new EEG/MEG analysis methods.
%
% The FieldTrip software is free but copyrighted software, distributed
% under the terms of the GNU General Public Licence as published by
% the Free Software Foundation (either version 2, or at your option
% any later version). See the file COPYING for more details.
% 
% The functions in this toolbox are copyrighted by their respective authors:
%   Robert Oostenveld, DCCN, FCDC, SMI, MBFYS
%   Jan-Matthijs Schoffelen, CCNi, FCDC
%   Pascal Fries, FCDC
%   Markus Bauer, FCDC
%   Ole Jensen, FCDC
%   Markus Siegel, FCDC, UKE
%   Jens Schwarzbach, FCDC
%   Eric Maris, DCC, FCDC
%   Ingrid Nieuwenhuis, DCCN, FCDC
%   Saskia Haegens, DCCN, FCDC
%
% Copyrights (C) 2008-2009, Donders Institute for Brain, Cognition and Behaviour, The Netherlands (DCCN, DCC, DCN)
% Copyrights (C) 2008-2009, Centre for Cognitive Neuroimaging in Glasgow, United Kingdom (CCNi)
% Copyrights (C) 2003-2008, F.C. Donders Centre, University Nijmegen, The Netherlands (FCDC)
% Copyrights (C) 2004-2007, Nijmegen Institute for Cognition and Information, The Netherlands (NICI)
% Copyrights (C) 2004-2005, Universitatsklinikum Hamburg-Eppendorf, Germany (UKE)
% Copyrights (C) 2003-2004, Center for Sensory Motor Interaction, University Aalborg, Denmark (SMI)
% Copyrights (C) 1999-2003, Department of Medical Physics, Radboud University Nijmegen, The Netherlands (MBFYS)
%
% The FieldTrip toolbox depend on functions from other toolboxes to do a
% large part of the actual work, such as reading data from binary files and
% forward and inverse modelling of the EEG/MEG. These low-level functions
% are contained in the private subdirectory. These other toolboxes on which
% the framework depends are copyrighted by their respective authors, see
% each individual matlab file for the details.
%
% Unauthorised copying and distribution of functions that are not
% explicitely covered by the GPL is not allowed!
%
% Below is an overview of the most important FieldTrip functions, sorted by
% category. You can get more details on a function by typing "help functionname"
% in Matlab.
%
% Preprocessing and reading data
%   definetrial
%   rejectartifact
%   rejectvisual
%   preprocessing
%   appenddata
%   resampledata
%   channelrepair
%   recodeevent
%   redefinetrial
%   read_fcdc_header
%   read_fcdc_data
%   read_fcdc_event
%   read_fcdc_mri
%
% Event-Related Fields or Potentials
%   timelockanalysis
%   timelockgrandaverage
%   timelockstatistics
%   singleplotER
%   topoplotER
%   multiplotER
%
% Frequency and Time-Frequency analysis
%   freqanalysis
%   freqanalysis_mtmfft
%   freqanalysis_mtmwelch
%   freqanalysis_mtmconvol
%   freqanalysis_wltconvol
%   freqanalysis_tfr
%   freqgrandaverage
%   freqdescriptives
%   freqstatistics
%   singleplotTFR
%   topoplotTFR
%   multiplotTFR
%
% Source analysis
%   dipolefitting
%   dipolesimulation
%   sourceanalysis
%   sourcegrandaverage
%   sourcedescriptives
%   sourcestatistics
%   sourceplot
%   sourceinterpolate
%   prepare_localspheres
%   prepare_singleshell
%   prepare_bemmodel
%   prepare_leadfield
%   prepare_atlas
%   volumelookup
%
% Statistical analysis
%   timelockstatistics
%   freqstatistics
%   sourcestatistics
%
% Plotting and display of data
%   prepare_layout
%   layoutplot
%   topoplot
%   topoplotER
%   topoplotTFR
%   multiplotER
%   multiplotTFR
%   singleplotER
%   singleplotTFR
%   sourceplot
%   clusterplot

% $Log: Contents.m,v $
% Revision 1.11  2009/02/04 10:14:30  roboos
% updated
%
% Revision 1.10  2009/02/04 10:10:17  roboos
% updated copyrights, added list with main function names
%
% Revision 1.9  2007/05/30 06:53:29  roboos
% added Eric and Ingrid
%
% Revision 1.8  2006/07/04 16:14:24  roboos
% removed obsolete function list, clarified and cleaned up copyrights
%
% Revision 1.7  2005/06/02 12:23:20  roboos
% updated some of the descriptions, it is still not complete
%
% Revision 1.6  2004/02/04 14:00:13  roberto
% added Jens to authors and included his functions
%
% Revision 1.5  2004/01/27 10:20:08  roberto
% added combineplanar
%
% Revision 1.4  2004/01/21 16:09:03  roberto
% updated documentation
%
% Revision 1.3  2004/01/14 10:34:21  roberto
% added Markus Siegels functions and listed him as author
%
% Revision 1.2  2003/11/28 21:37:26  roberto
% updated help and function list
%
% Revision 1.1  2003/11/07 16:20:14  roberto
% updated contents with new functions
%

