% FieldTrip
% Version unknown www.fieldtriptoolbox.org DD-MM-YYYY
%
% FieldTrip is the MATLAB toolbox for MEG and EEG analysis that is being developed at
% the Centre for Cognitive Neuroimaging of the Donders Institute for Brain, Cognition
% and Behaviour together with collaborating institutes. The FieldTrip software is
% released as open source under the GNU general public license.
%
% The FieldTrip toolbox includes algorithms for simple and advanced analysis of MEG and
% EEG data, such as time-frequency analysis, source reconstruction using dipoles,
% distributed sources and beamformers and non-parametric statistical testing. It
% supports the data formats of all major MEG systems (CTF, Neuromag, BTi) and of the
% most popular EEG systems, and new formats can be added easily. FieldTrip contains
% high-level functions that you can use to construct your own analysis protocol in
% MATLAB. Furthermore, it easily allows developers to incorporate low-level algorithms
% for new EEG/MEG analysis methods.
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
%   and others ...
%
% Copyrights (C) 2008-2017, Donders Institute for Brain, Cognition and Behaviour, Radboud University Nijmegen, The Netherlands (DCCN, DCC, DCN)
% Copyrights (C) 2011-2017, Max Planck Institute for Psycholynguistics, Nijmegen, The Netherlands (MPI)
% Copyrights (C) 2007-2017, The Wellcome Trust Centre for Neuroimaging, University College London, UK (UCL)
% Copyrights (C) 2008-2010, Centre for Cognitive Neuroimaging, Glasgow, United Kingdom (CCNi)
% Copyrights (C) 2003-2008, F.C. Donders Centre, Radboud University Nijmegen, The Netherlands (FCDC)
% Copyrights (C) 2004-2007, Nijmegen Institute for Cognition and Information, Radboud University Nijmegen, The Netherlands (NICI)
% Copyrights (C) 2004-2005, Universitatsklinikum Hamburg-Eppendorf, Germany (UKE)
% Copyrights (C) 2003-2004, Center for Sensory Motor Interaction, University Aalborg, Denmark (SMI)
% Copyrights (C) 1999-2003, Department of Medical Physics, Radboud University Nijmegen, The Netherlands (MBFYS)
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
% Below is an overview of the most important FieldTrip functions, sorted by category.
% You can get more details on a function by typing "help functionname" in MATLAB.
%
% Preprocessing and reading data
%   ft_definetrial
%   ft_rejectartifact
%   ft_rejectvisual
%   ft_preprocessing
%   ft_appenddata
%   ft_resampledata
%   ft_channelrepair
%   ft_recodeevent
%   ft_redefinetrial
%   ft_databrowser
%   ft_read_header
%   ft_read_data
%   ft_read_event
%   ft_read_mri
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
%   ft_sourceplot
%   ft_sourceinterpolate
%   ft_prepare_localspheres
%   ft_prepare_singleshell
%   ft_prepare_bemmodel
%   ft_prepare_leadfield
%   ft_prepare_atlas
%   ft_volumelookup
%
% Statistical analysis
%   ft_timelockstatistics
%   ft_freqstatistics
%   ft_sourcestatistics
%
% Plotting and display of data
%   ft_prepare_layout
%   ft_layoutplot
%   ft_singleplotER
%   ft_multiplotER
%   ft_topoplotER
%   ft_singleplotTFR
%   ft_multiplotTFR
%   ft_topoplotTFR
%   ft_topoplotCC
%   ft_multiplotCC
%   ft_sourceplot
%   ft_clusterplot

% $Id$
% $Rev$

