function [status] = ft_hastoolbox(toolbox, autoadd, silent)

% FT_HASTOOLBOX tests whether an external toolbox is installed. Optionally it will
% try to determine the path to the toolbox and install it automatically.
%
% Use as
%   [status] = ft_hastoolbox(toolbox, autoadd, silent)
%
% autoadd = -1 means that it will check and give an error when not yet installed
% autoadd =  0 means that it will check and give a warning when not yet installed
% autoadd =  1 means that it will check and give an error if it cannot be added
% autoadd =  2 means that it will check and give a warning if it cannot be added
% autoadd =  3 means that it will check but remain silent if it cannot be added
%
% silent = 0 means that it will give some feedback about adding the toolbox
% silent = 1 means that it will not give feedback

% Copyright (C) 2005-2022, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% this function is called many times in FieldTrip and associated toolboxes
% use efficient handling if the same toolbox has been investigated before
% persistent previous previouspath
%
% if ~isequal(previouspath, path)
%   previous = [];
% end
%
% if isempty(previous)
%   previous = struct;
% elseif isfield(previous, fixname(toolbox))
%   status = previous.(fixname(toolbox));
%   return
% end

% Note to developers: please do NOT use ft_warning and ft_error
% inside this function, but rather the normal warning and error.

% these are for speeding up subsequent calls with the same input arguments
persistent previous_argin previous_argout

if isdeployed
  % it is not possible to check the presence of functions or change the path in a compiled application
  status = true;
  return
end

if nargin<2
  % default is not to add the path automatically
  autoadd = 0;
end

if nargin<3
  % default is not to be silent
  silent = 0;
end

current_argin = {toolbox, autoadd, silent};
if isequal(current_argin, previous_argin)
  % return the previous output from cache
  status = previous_argout{1};
  return
end


% this points the user to the website where he/she can download the toolbox
url = {
  '14888-MUTUAL-INFORMATION-COMPUTATION'  'see http://www.mathworks.com/matlabcentral/fileexchange/14888-mutual-information-computation'
  '35625-INFORMATION-THEORY-TOOLBOX'      'see http://www.mathworks.com/matlabcentral/fileexchange/35625-information-theory-toolbox'
  '4D-VERSION'                            'contact Christian Wienbruch'
  'AFNI'                                  'see http://afni.nimh.nih.gov'
  'BAYESFACTOR'                           'see https://klabhub.github.io/bayesFactor'
  'BCI2000'                               'see http://bci2000.org'
  'BCT'                                   'see http://www.brain-connectivity-toolbox.net/'
  'BEMCP'                                 'contact Christophe Phillips'
  'BEOWULF'                               'see http://robertoostenveld.nl, or contact Robert Oostenveld'
  'BESA'                                  'see http://www.besa.de/downloads/matlab/ and get the "BESA MATLAB Readers"'
  'BIOSIG'                                'see http://biosig.sourceforge.net'
  'BRAINSTORM'                            'see http://neuroimage.ucs.edu/brainstorm'
  'BRAINSUITE'                            'see http://brainsuite.bmap.ucla.edu/processing/additional-tools/'
  'BRAINVISA'                             'see http://brainvisa.info'
  'BREWERMAP'                             'see https://nl.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer--attractive-and-distinctive-colormaps'
  'BSMART'                                'see http://www.brain-smart.org'
  'CCA'                                   'see http://www.imt.liu.se/~magnus/cca or contact Magnus Borga'
  'CELLFUNCTION'                          'see https://github.com/schoffelen/cellfunction'
  'CMOCEAN'                               'see https://nl.mathworks.com/matlabcentral/fileexchange/57773-matplotlib-perceptually-uniform-colormaps'
  'COLORCET'                              'see https://www.peterkovesi.com/matlabfns/index.html#colour'
  'COMM'                                  'see http://www.mathworks.com/products/communications'
  'COMPILER'                              'see http://www.mathworks.com/products/compiler'
  'CONNECTIVITY'                          'see http://www.fieldtriptoolbox.org'
  'CPD'                                   'see https://sites.google.com/site/myronenko/research/cpd'
  'DATAHASH'                              'see http://www.mathworks.com/matlabcentral/fileexchange/31272'
  'DENOISE'                               'see http://lumiere.ens.fr/Audition/adc/meg, or contact Alain de Cheveigne'
  'DIPOLI'                                'see ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip/external'
  'DISTCOMP'                              'see http://www.mathworks.nl/products/parallel-computing/'
  'DSS'                                   'see http://www.cis.hut.fi/projects/dss'
  'DUNEURO'                               'see http://duneuro.org/ and https://www.fieldtriptoolbox.org/workshop/ohbm2018/'
  'EEG'                                   'see http://eeg.sourceforge.net'
  'EEGLAB'                                'see http://www.sccn.ucsd.edu/eeglab'
  'EEGSF'                                 'see http://eeg.sourceforge.net'  % alternative name
  'EEPROBE'                               'see http://www.ant-neuro.com, or contact Maarten van der Velde'
  'EGI_MFF_V2'                            'see http://www.egi.com/ or contact either Phan Luu or Colin Davey at EGI'
  'EZC3D'                                 'see https://github.com/pyomeca/ezc3d'
  'FASTICA'                               'see http://www.cis.hut.fi/projects/ica/fastica'
  'FILEEXCHANGE'                          'see http://www.mathworks.com/matlabcentral/fileexchange/'
  'FILEIO'                                'see http://www.fieldtriptoolbox.org'
  'FNS'                                   'see http://hhvn.nmsu.edu/wiki/index.php/FNS'
  'FORWARD'                               'see http://www.fieldtriptoolbox.org'
  'FREESURFER'                            'see http://surfer.nmr.mgh.harvard.edu/fswiki'
  'GCMI'                                  'see https://github.com/robince/gcmi'
  'GIFTI'                                 'see http://www.artefact.tk/software/matlab/gifti'
  'GTEC'                                  'see http://www.gtec.at'
  'HOMER3'                                'see https://github.com/BUNPC/Homer3 and https://github.com/fNIRS/snirf_homer3'
  'IBTB'                                  'see Magri et al. BMC Neurosci 2009, 10:81'
  'ICASSO'                                'see http://www.cis.hut.fi/projects/ica/icasso'
  'IMAGES'                                'see http://www.mathworks.com/products/image'  % Mathworks refers to this as IMAGES
  'INVERSE'                               'see http://www.fieldtriptoolbox.org'
  'ISO2MESH'                              'see http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?Home or contact Qianqian Fang'
  'ITAB'                                  'contact Stefania Della Penna'
  'JSONIO'                                'see https://github.com/gllmflndn/JSONio'
  'JSONLAB'                               'see https://se.mathworks.com/matlabcentral/fileexchange/33381-jsonlab--a-toolbox-to-encode-decode-json-files'
  'JNIFTI'                                'see https://github.com/NeuroJSON/jnifti'
  'LAGEXTRACTION'                         'see https://github.com/agramfort/eeglab-plugin-ieee-tbme-2010'
  'MARS'                                  'see http://www.parralab.org/mars'
  'MATLAB2BESA'                           'see http://www.besa.de/downloads/matlab/ and get the "MATLAB to BESA Export functions"'
  'MATNWB'                                'see https://neurodatawithoutborders.github.io/matnwb/'
  'MATPLOTLIB'                            'see https://nl.mathworks.com/matlabcentral/fileexchange/62729-matplotlib-perceptually-uniform-colormaps'
  'MAYO_MEF'                              'see https://github.com/MultimodalNeuroimagingLab/mef_reader_fieldtrip and https://msel.mayo.edu/codes.html'
  'MEG-CALC'                              'this is a commercial toolbox from Neuromag, see http://www.neuromag.com'
  'MEG-PD'                                'see http://www.kolumbus.fi/kuutela/programs/meg-pd'
  'MENTAT'                                'see http://robertoostenveld.nl, or contact Robert Oostenveld'
  'MFFMATLABIO'                           'see https://github.com/arnodelorme/mffmatlabio'
  'MISC'                                  'various functions that were downloaded from http://www.mathworks.com/matlabcentral/fileexchange and elsewhere'
  'MNE'                                   'see https://mne.tools/stable/overview/matlab.html'
  'MRI'                                   'see http://eeg.sourceforge.net'  % alternative name
  'MRTRIX'                                'see https://mrtrix.org'
  'MVPA-LIGHT'                            'see https://github.com/treder/MVPA-Light'
  'MYSQL'                                 'see http://www.mathworks.com/matlabcentral/fileexchange/8663-mysql-database-connector'
  'NETCDF'                                'see http://www.mathworks.com/matlabcentral/fileexchange/15177'
  'NEURALYNX_V3'                          'see http://www.urut.ch/new/serendipity/index.php?/pages/nlxtomatlab.html'
  'NEURALYNX_V6'                          'see https://neuralynx.com/software/category/matlab-netcom-utilities/ and take the version from Neuralynx'
  'NEURONE'                               'see http://www.megaemg.com/support/unrestricted-downloads'
  'NEUROSHARE'                            'see http://www.neuroshare.org'
  'NLXNETCOM'                             'see http://www.neuralynx.com'
  'NPMK'                                  'see https://github.com/BlackrockMicrosystems/NPMK'
  'NWAY'                                  'see http://www.models.kvl.dk/source/nwaytoolbox'
  'OPENMEEG'                              'see http://openmeeg.github.io and http://www.fieldtriptoolbox.org/faq/how_do_i_install_the_openmeeg_binaries'
  'OPTIM'                                 'see http://www.mathworks.com/products/optim'
  'PEER'                                  'see http://www.fieldtriptoolbox.org'
  'PEER'                                  'see http://www.fieldtriptoolbox.org/development/peer'
  'PLEXON'                                'available from http://www.plexon.com/assets/downloads/sdk/ReadingPLXandDDTfilesinMatlab-mexw.zip'
  'PLOT2SVG'                              'see http://www.mathworks.com/matlabcentral/fileexchange/7401-scalable-vector-graphics-svg-export-of-figures'
  'PLOTTING'                              'see http://www.fieldtriptoolbox.org'
  'PLOTTING'                              'see http://www.fieldtriptoolbox.org'
  'PREPROC'                               'see http://www.fieldtriptoolbox.org'
  'PRTOOLS'                               'see http://www.prtools.org'
  'REALTIME'                              'see http://www.fieldtriptoolbox.org'
  'RICOH_MEG_READER'                      'contact Ricoh engineers'
  'SIGNAL'                                'see http://www.mathworks.com/products/signal'
  'SIMBIO'                                'see https://www.mrt.uni-jena.de/simbio/index.php/Main_Page'
  'SON2'                                  'see http://www.kcl.ac.uk/depsta/biomedical/cfnr/lidierth.html, or contact Malcolm Lidierth'
  'SPECEST'                               'see http://www.fieldtriptoolbox.org'
  'SPIKE'                                 'see http://www.fieldtriptoolbox.org'
  'SPLINES'                               'see http://www.mathworks.com/products/splines'
  'SPM'                                   'see http://www.fil.ion.ucl.ac.uk/spm'
  'SPM12'                                 'see http://www.fil.ion.ucl.ac.uk/spm'
  'SPM2'                                  'see http://www.fil.ion.ucl.ac.uk/spm'
  'SPM5'                                  'see http://www.fil.ion.ucl.ac.uk/spm'
  'SPM8'                                  'see http://www.fil.ion.ucl.ac.uk/spm'
  'SPM99'                                 'see http://www.fil.ion.ucl.ac.uk/spm'
  'SQDPROJECT'                            'see http://www.isr.umd.edu/Labs/CSSL/simonlab'
  'TCP_UDP_IP'                            'see http://www.mathworks.com/matlabcentral/fileexchange/345, or contact Peter Rydesaeter'
  'TOOLBOX_GRAPH'                         'see http://www.mathworks.com/matlabcentral/fileexchange/5355-toolbox-graph or contact Gabriel Peyre'
  'VISION'                                'see https://nl.mathworks.com/products/computer-vision/'
  'VGRID'                                 'see http://www.rheinahrcampus.de/~medsim/vgrid/manual.html'
  'VIDEOMEG'                              'see https://github.com/andreyzhd/VideoMEG'
  'WAVEFRONT'                             'see http://mathworks.com/matlabcentral/fileexchange/27982-wavefront-obj-toolbox'
  'XDF'                                   'see https://github.com/xdf-modules/xdf-Matlab'
  'XML4MAT'                               'see http://www.mathworks.com/matlabcentral/fileexchange/6268-xml4mat-v2-0'
  'XSENS'                                 'see https://www.xsens.com/motion-capture and http://www.fieldtriptoolbox.org/getting_started/xsens/'
  'XUNIT'                                 'see http://www.mathworks.com/matlabcentral/fileexchange/22846-matlab-xunit-test-framework'
  'YOKOGAWA'                              'this is deprecated, please use YOKOGAWA_MEG_READER instead'
  'YOKOGAWA_MEG_READER'                   'contact Ricoh engineers'
  };

% determine whether the toolbox is installed
toolbox = upper(toolbox);

% In case SPMxUP not available, allow to use fallback toolbox
fallback_toolbox='';

switch toolbox
  case 'AFNI'
    dependency= {'BrikLoad', 'BrikInfo'};
  case 'DSS'
    dependency= {'denss', 'dss_create_state'};
  case 'EEGLAB'
    dependency = 'runica';
  case 'NWAY'
    dependency = 'parafac';
  case 'SPM'
    dependency = 'spm'; % any version of SPM is fine
  case 'SPM99'
    dependency = {'spm', get_spm_version()==99};
  case 'SPM2'
    dependency = {'spm', get_spm_version()==2};
  case 'SPM2UP' % version 2 or later, but not SPM 9X
    dependency = {'spm', get_spm_version()>=2, get_spm_version()<95};
    % this is to avoid crashes when trying to add SPM to the path
    fallback_toolbox = 'SPM2';
  case 'SPM5'
    dependency = {'spm', get_spm_version()==5};
  case 'SPM5UP' % version 5 or later, but not SPM 9X
    dependency = {'spm', get_spm_version()>=5, get_spm_version()<95};
    % this is to avoid crashes when trying to add SPM to the path
    fallback_toolbox = 'SPM5';
  case 'SPM8'
    dependency = {'spm', get_spm_version()==8};
  case 'SPM8UP' % version 8 or later, but not SPM 9X
    dependency = {'spm', get_spm_version()>=8, get_spm_version()<95};
    % this is to avoid crashes when trying to add SPM to the path
    fallback_toolbox = 'SPM8';
  case 'SPM12'
    dependency = {'spm', get_spm_version()==12};
  case 'SPM12UP' % version 12 or later, but not SPM 9X
    dependency = {'spm', get_spm_version()>=12, get_spm_version()<95};
    % this is to avoid crashes when trying to add SPM to the path
    fallback_toolbox = 'SPM12';
  case 'MEG-PD'
    dependency = {'rawdata', 'channames'};
  case 'MEG-CALC'
    dependency = {'megmodel', 'megfield', 'megtrans'};
  case 'MVPA-LIGHT'
    dependency = {'mv_classify','train_lda'};
  case 'BIOSIG'
    dependency = {'sopen', 'sread'};
  case 'EEG'
    dependency = {'ctf_read_res4', 'ctf_read_meg4'};
  case 'EEGSF'  % alternative name
    dependency = {'ctf_read_res4', 'ctf_read_meg4'};
  case 'MRI'    % other functions in the mri section
    dependency = {'avw_hdr_read', 'avw_img_read'};
  case 'NEUROSHARE'
    dependency = {'ns_OpenFile', 'ns_SetLibrary', 'ns_GetAnalogData'};
  case 'ARTINIS'
    dependency = {'read_artinis_oxy3'};
  case 'BESA'
    dependency = {'readBESAavr', 'readBESAelp', 'readBESAswf'};
  case 'MATLAB2BESA'
    dependency = {'besa_save2Avr', 'besa_save2Elp', 'besa_save2Swf'};
  case 'EEPROBE'
    dependency = {'read_eep_avr', 'read_eep_cnt'};
  case 'YOKOGAWA'
    dependency = @()hasyokogawa('16bitBeta6');
  case 'YOKOGAWA12BITBETA3'
    dependency = @()hasyokogawa('12bitBeta3');
  case 'YOKOGAWA16BITBETA3'
    dependency = @()hasyokogawa('16bitBeta3');
  case 'YOKOGAWA16BITBETA6'
    dependency = @()hasyokogawa('16bitBeta6');
  case 'YOKOGAWA_MEG_READER'
    dependency = @()hasyokogawa('1.5');
  case 'RICOH_MEG_READER'
    dependency = @()hasricoh('1.0');
  case 'BEOWULF'
    dependency = {'evalwulf', 'evalwulf', 'evalwulf'};
  case 'MENTAT'
    dependency = {'pcompile', 'pfor', 'peval'};
  case 'SON2'
    dependency = {'SONFileHeader', 'SONChanList', 'SONGetChannel'};
  case '4D-VERSION'
    dependency  = {'read4d', 'read4dhdr'};
  case {'STATS', 'STATISTICS'}
    dependency = {has_license('statistics_toolbox'), 'betacdf', 'raylcdf', 'unidcdf'};      % also check the availability of a toolbox license
  case {'OPTIM', 'OPTIMIZATION'}
    dependency = {has_license('optimization_toolbox'), 'fminunc', 'optimset'};              % also check the availability of a toolbox license
  case {'SPLINES', 'CURVE_FITTING'}
    dependency = {has_license('curve_fitting_toolbox'), 'smooth', 'fit'};                   % also check the availability of a toolbox license
  case 'COMM'
    dependency = {has_license('communication_toolbox'), 'de2bi', 'fskmod', 'pskmod'};       % also check the availability of a toolbox license
  case 'SIGNAL'
    dependency = {has_license('signal_toolbox'), 'window', 'hanning'};                      % also check the availability of a toolbox license
  case 'IMAGES'
    dependency = {has_license('image_toolbox'), 'imerode', 'imdilate'};                     % also check the availability of a toolbox license
  case 'VISION'
    dependency = {has_license('video_and_image_blockset'), 'pointCloud', 'pcnormals'};      % also check the availability of a toolbox license
  case {'DCT', 'DISTCOMP'}
    dependency = {has_license('distrib_computing_toolbox'), 'parpool', 'batch'};            % also check the availability of a toolbox license
  case 'COMPILER'
    dependency = {has_license('compiler'), 'mcc', 'mcr'};                                   % also check the availability of a toolbox license
  case 'FASTICA'
    dependency = 'fpica';
  case 'BRAINSTORM'
    dependency = 'process_fooof';
  case 'DENOISE'
    dependency = {'tsr', 'sns'};
  case 'CTF'
    dependency = {'getCTFBalanceCoefs', 'getCTFdata'};
  case 'BCI2000'
    dependency  = {'load_bcidat'};
  case 'NLXNETCOM'
    dependency = {'MatlabNetComClient', 'NlxConnectToServer', 'NlxGetNewCSCData'};
  case 'DIPOLI'
    dependency = {'dipoli.maci', 'file'};
  case 'MNE'
    dependency = {'fiff_read_meas_info', 'fiff_setup_read_raw'};
  case 'TCP_UDP_IP'
    dependency = {'pnet', 'pnet_getvar', 'pnet_putvar'};
  case 'BEMCP'
    dependency = {'bem_Cij_cog', 'bem_Cij_lin', 'bem_Cij_cst'};
  case 'OPENMEEG'
    dependency = {'om_save_tri'};
  case 'PRTOOLS'
    dependency = {'prversion', 'dataset', 'svc'};
  case 'ITAB'
    dependency = {'lcReadHeader', 'lcReadData'};
  case 'BSMART'
    dependency = 'bsmart';
  case 'FREESURFER'
    dependency = {'MRIread', 'vox2ras_0to1'};
  case 'FNS'
    dependency = 'elecsfwd';
  case 'SIMBIO'
    dependency = {'calc_stiff_matrix_val', 'sb_transfer'};
  case 'VGRID'
    dependency = 'vgrid';
  case 'GIFTI'
    dependency = 'gifti';
  case 'XML4MAT'
    dependency = {'xml2struct', 'xml2whos'};
  case 'SQDPROJECT'
    dependency = {'sqdread', 'sqdwrite'};
  case 'BCT'
    dependency = {'macaque71.mat', 'motif4funct_wei'};
  case 'CCA'
    dependency = {'ccabss'};
  case 'MFFMATLABIO'
    dependency = {'eegplugin_mffmatlabio', 'mff_getobj'};
  case 'EGI_MFF'
    dependency = {'mff_getObject', 'mff_getSummaryInfo'};
  case 'TOOLBOX_GRAPH'
    dependency = 'toolbox_graph';
  case 'NETCDF'
    dependency = {'netcdf'};
  case 'MYSQL'
    % this only consists of a single mex file
    dependency = has_mex('mysql');
  case 'ISO2MESH'
    dependency = {'vol2surf', 'qmeshcut'};
  case 'DATAHASH'
    dependency = {'DataHash'};
  case 'IBTB'
    dependency = {'binr','information', 'eqpop', 'eqspace', 'ceqspace', 'gseqspace'};
  case 'ICASSO'
    dependency = {'icassoEst'};
  case 'XUNIT'
    dependency = {'initTestSuite', 'runtests'};
  case 'PLEXON'
    dependency = {'plx_adchan_gains', 'mexPlex'};
  case '35625-INFORMATION-THEORY-TOOLBOX'
    dependency = {'conditionalEntropy', 'entropy', 'jointEntropy',...
      'mutualInformation' 'nmi' 'nvi' 'relativeEntropy'};
  case '14888-MUTUAL-INFORMATION-COMPUTATION'
    dependency = {'condentropy', 'demo_mi', 'estcondentropy.cpp', 'estjointentropy.cpp', 'estpa.cpp', 'findjointstateab.cpp', 'makeosmex.m', 'mutualinfo.m', 'condmutualinfo.m', 'entropy.m', 'estentropy.cpp', 'estmutualinfo.cpp', 'estpab.cpp', 'jointentropy.m' 'mergemultivariables.m' };
  case 'PLOT2SVG'
    dependency = {'plot2svg.m', 'simulink2svg.m'};
  case 'BRAINSUITE'
    dependency = {'readdfs.m', 'writedfc.m'};
  case 'BRAINVISA'
    dependency = {'loadmesh.m', 'plotmesh.m', 'savemesh.m'};
  case 'NEURALYNX_V6'
    dependency = has_mex('Nlx2MatCSC');
  case 'NEURALYNX_V3'
    dependency = has_mex('Nlx2MatCSC_v3');
  case 'NPMK'
    dependency = {'OpenNSx' 'OpenNEV'};
  case 'VIDEOMEG'
    dependency = {'comp_tstamps' 'load_audio0123', 'load_video123'};
  case 'WAVEFRONT'
    dependency = {'write_wobj' 'read_wobj'};
  case 'NEURONE'
    dependency = {'readneurone' 'readneuronedata' 'readneuroneevents'};
  case 'BREWERMAP'
    dependency = {'brewermap' 'brewermap_view'};
  case 'GTEC'
    dependency = {'ghdf5read' 'ghdf5fileimport'};
  case 'MARS'
    dependency = {'spm_mars_mrf'};
  case 'LAGEXTRACTION'
    dependency = {'extractlag' 'perform_realign'};
  case 'JSONLAB'
    dependency = {'loadjson' 'savejson'};
  case 'JNIFTI'
    dependency = {'loadjnifti' 'savejnifti'};
  case 'PLOTLY'
    dependency = {'fig2plotly' 'savejson'};
  case 'JSONIO'
    dependency = {'jsonread', 'jsonwrite', 'jsonread.mexa64'};
  case 'CPD'
    dependency = {'cpd', 'cpd_affine', 'cpd_P'};
  case 'XDF'
    dependency = {'load_xdf', 'load_xdf_innerloop'};
  case 'MRTRIX'
    dependency = {'read_mrtrix'};
  case 'BAYESFACTOR'
    dependency = {'bf.ttest', 'bf.ttest2'};
  case 'EZC3D'
    dependency = {'ezc3dRead', 'ezc3dWrite'};
  case 'GCMI'
    dependency = {'copnorm' 'mi_gg'};
  case 'XSENS'
    dependency = {'load_mvnx'};
  case 'MAYO_MEF' % MEF 2.1 and MEF 3.0
    dependency = {'MEFFieldTrip_2p1', 'MEFFieldTrip_3p0'};
  case 'MATNWB'
    dependency = {'nwbRead', 'generateCore'};
  case 'MATPLOTLIB'
    dependency = {'cividis', 'inferno', 'magma', 'plasma', 'tab10', 'tab20', 'tab20b', 'tab20c', 'twilight', 'viridis'};
  case 'CMOCEAN'
    dependency = {'cmocean'};
  case 'COLORCET'
    dependency = {'colorcet'};
  case 'FILEEXCHANGE'
    dependency = is_subdir_in_fieldtrip_path('/external/fileexchange');
  case 'HOMER3'
    dependency = {'SnirfClass' 'DataClass' 'AuxClass' 'MeasListClass' 'MetaDataTagsClass' 'ProbeClass' 'StimClass'};
  case 'DUNEURO'
    dependency = {'duneuro_meeg', 'duneuro_function', 'compute_B_primary'};

    % the following are FieldTrip modules or toolboxes
  case 'FILEIO'
    dependency = {'ft_read_header', 'ft_read_data', 'ft_read_event', 'ft_read_sens'};
  case 'FORWARD'
    dependency = {'ft_compute_leadfield', 'ft_prepare_vol_sens'};
  case 'INVERSE'
    dependency = {'ft_inverse_dics', 'ft_inverse_dipolefit', 'ft_inverse_lcmv', 'ft_inverse_mne', 'ft_inverse_pcc'};
  case 'PLOTTING'
    dependency = {'ft_plot_topo', 'ft_plot_mesh', 'ft_plot_matrix'};
  case 'QSUB'
    dependency = {'qsubcellfun', 'qsubfeval', 'qsubget'};
  case 'PEER'
    dependency = {'peercellfun', 'peerfeval', 'peerget'};
  case 'ENGINE'
    dependency = {'enginecellfun', 'enginefeval', 'engineget'};
  case 'CONNECTIVITY'
    dependency = {'ft_connectivity_corr', 'ft_connectivity_granger'};
  case 'SPIKE'
    dependency = {'ft_spiketriggeredaverage', 'ft_spiketriggeredspectrum'};
  case 'CELLFUNCTION'
    dependency = {'cellmean', 'cellvecadd', 'cellcat'};
  case 'SPECEST'
    dependency = {'ft_specest_mtmconvol', 'ft_specest_mtmfft', 'ft_specest_wavelet'};
  case 'PREPROC'
    dependency = {'ft_preproc_detrend', 'ft_preproc_baselinecorrect', 'ft_preproc_bandpassfilter', 'ft_preproc_bandstopfilter'};
  case {'REALTIME', 'STATFUN', 'TRIALFUN', 'TEMPLATE/LAYOUT', 'TEMPLATE/ANATOMY', 'TEMPLATE/ATLAS', 'TEMPLATE/DEWAR', 'TEMPLATE/HEADMODEL', 'TEMPLATE/ELECTRODE', 'TEMPLATE/NEIGHBOURS', 'TEMPLATE/SOURCEMODEL'}
    dependency = is_subdir_in_fieldtrip_path(toolbox);
  otherwise
    if ~silent, ft_warning('cannot determine whether the %s toolbox is present', toolbox); end
    dependency = false;
end

status = is_present(dependency);
if ~status && ~isempty(fallback_toolbox)
  % in case of SPMxUP
  toolbox = fallback_toolbox;
end

% try to determine the path of the requested toolbox and add it
if ~status && autoadd>0

  % for core FieldTrip modules
  prefix = fileparts(which('ft_defaults'));
  if ~status
    status = myaddpath(fullfile(prefix, lower(toolbox)), silent);
  end

  % for external FieldTrip modules
  prefix = fullfile(fileparts(which('ft_defaults')), 'external');
  if ~status
    status = myaddpath(fullfile(prefix, lower(toolbox)), silent);
    licensefile = [lower(toolbox) '_license'];
    if status && exist(licensefile, 'file')
      % this will execute openmeeg_license, mne_license and duneuro_license
      % which display the license on screen for three seconds
      feval(licensefile);
    end
  end

  % for contributed FieldTrip extensions
  prefix = fullfile(fileparts(which('ft_defaults')), 'contrib');
  if ~status
    status = myaddpath(fullfile(prefix, lower(toolbox)), silent);
    licensefile = [lower(toolbox) '_license'];
    if status && exist(licensefile, 'file')
      % this will execute openmeeg_license, mne_license and artinis_license
      % which display the license on screen for a few seconds
      feval(licensefile);
    end
  end

  % for linux computers in the Donders Centre for Cognitive Neuroimaging
  prefix = '/home/common/matlab';
  if ~status && is_folder(prefix)
    status = myaddpath(fullfile(prefix, lower(toolbox)), silent);
  end

  % for windows computers in the Donders Centre for Cognitive Neuroimaging
  prefix = 'h:\common\matlab';
  if ~status && is_folder(prefix)
    status = myaddpath(fullfile(prefix, lower(toolbox)), silent);
  end

  % use the MATLAB subdirectory in your homedirectory, this works on linux and mac
  prefix = fullfile(getenv('HOME'), 'matlab');
  if ~status && is_folder(prefix)
    status = myaddpath(fullfile(prefix, lower(toolbox)), silent);
  end

  if ~status
    % the toolbox is not on the path and cannot be added
    sel = find(strcmp(url(:,1), toolbox));
    if ~isempty(sel)
      msg = sprintf('the %s toolbox is not installed, %s', toolbox, url{sel, 2});
    else
      msg = sprintf('the %s toolbox is not installed', toolbox);
    end
    if autoadd==1
      error(msg);
    elseif autoadd==2
      warning(msg);
    else
      % fail silently
    end
  end

elseif ~status && autoadd<0
  % the toolbox is not on the path and should not be added
  sel = find(strcmp(url(:,1), toolbox));
  if ~isempty(sel)
    msg = sprintf('the %s toolbox is not installed, %s', toolbox, url{sel, 2});
  else
    msg = sprintf('the %s toolbox is not installed', toolbox);
  end
  error(msg);
end

% this function is called many times in FieldTrip and associated toolboxes
% use efficient handling if the same toolbox has been investigated before
% if status
%  previous.(fixname(toolbox)) = status;
% end

% remember the previous path, allows us to determine on the next call
% whether the path has been modified outise of this function
% previouspath = path;


% remember the current input and output arguments, so that they can be
% reused on a subsequent call in case the same input argument is given
current_argout = {status};
previous_argin  = current_argin;
previous_argout = current_argout;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = myaddpath(toolbox, silent)
global ft_default

if ~is_folder(toolbox)
  % search for a case-insensitive match, this is needed for MVPA-Light
  [p, f] = fileparts(toolbox);
  dirlist = dir(p);
  sel = strcmpi({dirlist.name}, f);
  if sum(sel)==1
    toolbox = fullfile(p, dirlist(sel).name);
  end
end

if isdeployed
  ft_warning('cannot change path settings for %s in a compiled application', toolbox);
  status = true;
elseif is_folder(toolbox)
  if ~silent
    ft_warning('off','backtrace');
    ft_warning('adding %s toolbox to your MATLAB path', toolbox);
    ft_warning('on','backtrace');
  end
  if any(~cellfun(@isempty, regexp(lower(toolbox), {'spm2$', 'spm5$', 'spm8$', 'spm12$'})))
    % SPM needs to be added with all its subdirectories
    addpath(genpath(toolbox));
    % check whether the mex files are compatible
    check_spm_mex;
  elseif ~isempty(regexp(lower(toolbox), 'mvpa-light$', 'once'))
    % this comes with its own startup script
    addpath(fullfile(toolbox, 'startup'))
    startup_MVPA_Light;
  elseif ~isempty(regexp(lower(toolbox), 'ibtb', 'once'))
    % this needs to be added with all its subdirectories
    addpath(genpath(toolbox));
  else
    addpath(toolbox);
  end
  % remember the toolbox that was just added to the path, it will be cleaned up by FT_POSTAMBLE_HASTOOLBOX
  if ~isfield(ft_default, 'toolbox') || ~isfield(ft_default.toolbox, 'cleanup')
    ft_default.toolbox.cleanup = {};
  end
  ft_default.toolbox.cleanup{end+1} = toolbox;
  status = true;
elseif (~isempty(regexp(toolbox, 'spm2$', 'once')) || ~isempty(regexp(toolbox, 'spm5$', 'once')) || ~isempty(regexp(toolbox, 'spm8$', 'once')) || ~isempty(regexp(toolbox, 'spm12$', 'once'))) && exist([toolbox 'b'], 'dir')
  % the final release version of SPM is not available, add the beta version instead
  status = myaddpath([toolbox 'b'], silent);
else
  status = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function path = unixpath(path)
%path(path=='\') = '/'; % replace backward slashes with forward slashes
path = strrep(path,'\','/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = hasfunction(funname, toolbox)
try
  % call the function without any input arguments, which probably is inapropriate
  feval(funname);
  % it might be that the function without any input already works fine
  status = true;
catch
  % either the function returned an error, or the function is not available
  % availability is influenced by the function being present and by having a
  % license for the function, i.e. in a concurrent licensing setting it might
  % be that all toolbox licenses are in use
  m = lasterror;
  if strcmp(m.identifier, 'MATLAB:license:checkouterror')
    if nargin>1
      ft_warning('the %s toolbox is available, but you don''t have a license for it', toolbox);
    else
      ft_warning('the function ''%s'' is available, but you don''t have a license for it', funname);
    end
    status = false;
  elseif strcmp(m.identifier, 'MATLAB:UndefinedFunction')
    status = false;
  else
    % the function seems to be available and it gave an unknown error,
    % which is to be expected with inappropriate input arguments
    status = true;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = is_subdir_in_fieldtrip_path(toolbox_name)
fttrunkpath = unixpath(fileparts(which('ft_defaults')));
fttoolboxpath = fullfile(fttrunkpath, lower(toolbox_name));
needle   = [pathsep fttoolboxpath pathsep];
haystack = [pathsep path() pathsep];
status   = contains(haystack, needle);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = has_mex(name)
full_name=[name '.' mexext];
status = (exist(full_name, 'file')==3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = get_spm_version()
if ~is_present('spm')
  v=NaN;
  return
end

version_str = spm('ver');
token = regexp(version_str,'(\d*)','tokens');
v = str2num([token{:}{:}]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = check_spm_mex()
status = true;
try
  % this will always result in an error
  spm_conv_vol
catch
  me = lasterror;
  % any error is ok, except when due to an invalid MEX file
  status = ~isequal(me.identifier, 'MATLAB:mex:ErrInvalidMEXFile');
end
if ~status
  % SPM8 mex file issues are common on macOS with recent MATLAB versions
  ft_warning('the SPM mex files are incompatible with your platform, see http://bit.ly/2OGF6US');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = has_license(toolbox_name)
% NOTE: this explicitly checks out a license, which may be suboptimal in
% terms of license use. Consider using the option 'test', but this needs to
% be checked with respect to backward compatibility first
status = license('checkout', toolbox_name)==1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = is_present(dependency)
if iscell(dependency)
  % use recursion
  status = all(cellfun(@is_present,dependency));
elseif islogical(dependency)
  % boolean
  status = all(dependency);
elseif ischar(dependency)
  % name of a function
  status = is_function_present_in_search_path(dependency);
elseif isa(dependency, 'function_handle')
  status = dependency();
else
  assert(false,'this should not happen');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = is_function_present_in_search_path(function_name)
w = which(function_name);

% must be in path and not a variable
status = ~isempty(w) && ~isequal(w, 'variable');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ISFOLDER is needed for versions prior to 2017b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = is_folder(dirpath)
tf = exist(dirpath,'dir') == 7;
