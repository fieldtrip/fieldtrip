% FieldTrip
% Version unknown www.fieldtriptoolbox.org DD-MM-YYYY
%
% FieldTrip is the MATLAB software toolbox for MEG and EEG analysis that is being
% developed by a team of researchers at the Donders Institute for Brain, Cognition and
% Behaviour in Nijmegen, the Netherlands in close collaboration with collaborating
% institutes.
%
% The toolbox offers advanced analysis methods of MEG, EEG, and invasive
% electrophysiological data, such as time-frequency analysis, source reconstruction
% using dipoles, distributed sources and beamformers and non-parametric statistical
% testing. It supports the [data formats](/faq/dataformat) of all major MEG systems
% (CTF, Neuromag/Elekta/Megin, BTi/4D, Yokogawa/Ricoh) and of most popular EEG systems,
% and new formats can be added easily. FieldTrip contains high-level functions that you
% can use to construct your own analysis protocols in MATLAB. Furthermore, it easily
% allows methods researchers to incorporate new methods for EEG/MEG analysis.
%
% The FieldTrip software is free but copyrighted software, distributed under the terms
% of the GNU General Public Licence as published by the Free Software Foundation (either
% version 2, or at your option any later version). See the file COPYING for more
% details.
%
% The functions in this toolbox are copyrighted by their respective authors:
%   Robert Oostenveld, DCCN, FCDC, SMI, MBFYS
%   Jan-Mathijs Schoffelen, CCNi, FCDC
%   Pascal Fries, FCDC
%   Markus Bauer, FCDC
%   Ole Jensen, FCDC
%   Markus Siegel, FCDC, UKE
%   Jens Schwarzbach, FCDC
%   Eric Maris, DCC, FCDC
%   Ingrid Nieuwenhuis, DCCN, FCDC
%   Saskia Haegens, DCCN, FCDC
%   Vladimir Litvak, UCL
%   Eelke Spaak, DCCN
%   Jorn Horschich, DCCN
%   Roemer van der Meij, DCC
%   Lilla Magyari, MPI, DCCN
%   and many others ...
%
% Copyright (C) 2008-2018, Donders Institute, Radboud University, The Netherlands (DCCN, DCC, DCN, DCMN)
% Copyright (C) 2014-2018, Karolinska Institute, Stockholm, Sweden (NatMEG)
% Copyright (C) 2012-2016, Max Planck Institute for Psycholinguistics, The Netherlands (MPI)
% Copyright (C) 2010-2013, Swammerdam Institute for Life Sciences, University of Amsterdam (SILS)
% Copyright (C) 2008-2009, Centre for Cognitive Neuroimaging in Glasgow, United Kingdom (CCNi)
% Copyright (C) 2009-2009, Netherlands Institute for Neuroscience (NIN)
% Copyright (C) 2003-2008, F.C. Donders Centre, Radboud University Nijmegen, The Netherlands (FCDC)
% Copyright (C) 2004-2007, Nijmegen Institute for Cognition and Information, The Netherlands (NICI)
% Copyright (C) 2004-2005, Universitatsklinikum Hamburg-Eppendorf, Germany (UKE)
% Copyright (C) 2003-2004, Center for Sensory Motor Interaction, University Aalborg, Denmark (SMI)
% Copyright (C) 1999-2003, Department of Medical Physics, Katholieke Universiteit Nijmegen, The Netherlands (MBFYS)
%
% The FieldTrip toolbox depend on functions from other toolboxes to do a large part of
% the actual work, such as reading data from binary files and forward and inverse
% modelling of the EEG/MEG. These low-level functions are contained in the private
% subdirectory. These other toolboxes on which the framework depends are copyrighted by
% their respective authors, see each individual MATLAB file for the details.
%
% Unauthorised copying and distribution of functions that are not explicitely covered by
% the GPL is not allowed!
%
% Below is a non-exhaustive overview of some of the important FieldTrip functions, sorted by category.
% You can get more details on a function by typing 'help functionname' in MATLAB.
%
% Preprocessing, reading and converting data
%   ft_definetrial
%   ft_artifact_eog
%   ft_artifact_jump
%   ft_artifact_muscle
%   ft_rejectartifact
%   ft_rejectvisual
%   ft_preprocessing
%   ft_appenddata
%   ft_resampledata
%   ft_channelrepair
%   ft_recodeevent
%   ft_redefinetrial
%
% Event-Related Fields or Potentials
%   ft_timelockanalysis
%   ft_timelockgrandaverage
%   ft_timelockstatistics
%   ft_singleplotER
%   ft_topoplotER
%   ft_multiplotER
%
% Frequency and Time-Frequency analysis
%   ft_freqanalysis
%   ft_freqbaseline
%   ft_freqgrandaverage
%   ft_freqdescriptives
%   ft_freqstatistics
%   ft_singleplotTFR
%   ft_topoplotTFR
%   ft_multiplotTFR
%
% Source analysis
%   ft_dipolefitting
%   ft_dipolesimulation
%   ft_sourceanalysis
%   ft_sourcegrandaverage
%   ft_sourcedescriptives
%   ft_sourcestatistics
%   ft_sourceparcellate
%   ft_sourceplot
%   ft_sourceinterpolate
%   ft_prepare_leadfield
%   ft_volumelookup
%   ft_volumenormalise
%   ft_volumesegment
%
% Statistical analysis
%   ft_timelockstatistics
%   ft_freqstatistics
%   ft_sourcestatistics
%   ft_statfun_actvsblT
%   ft_statfun_depsamplesFunivariate
%   ft_statfun_depsamplesFmultivariate
%   ft_statfun_depsamplesT
%   ft_statfun_depsamplesregrT
%   ft_statfun_diff
%   ft_statfun_diff_itc
%   ft_statfun_indepsamplesF
%   ft_statfun_indepsamplesT
%   ft_statfun_indepsamplesZcoh
%   ft_statfun_indepsamplesregrT
%   ft_statfun_mean
%   ft_statfun_pooledT
%   ft_statfun_roc
%
% Plotting and display of data
%   ft_clusterplot
%   ft_layoutplot
%   ft_movieplotER
%   ft_movieplotTFR
%   ft_multiplotCC
%   ft_multiplotER
%   ft_multiplotTFR
%   ft_mvaranalysis
%   ft_neighbourplot
%   ft_prepare_layout
%   ft_singleplotER
%   ft_singleplotTFR
%   ft_sourceplot
%   ft_topoplotER
%   ft_topoplotIC
%   ft_topoplotTFR
