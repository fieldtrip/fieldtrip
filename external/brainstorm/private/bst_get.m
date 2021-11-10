function [argout1, argout2, argout3, argout4, argout5] = bst_get( varargin )
% BST_GET: Get a Brainstorm structure.
% This function is used to abstract the way that these structures are stored.
%
% USAGE :
% ====== DIRECTORIES ==================================================================
%    - bst_get('UserDir')               : User home directory
%    - bst_get('BrainstormHomeDir')     : Application directory of brainstorm
%    - bst_get('BrainstormUserDir')     : User home directory for brainstorm (<home>/.brainstorm/)
%    - bst_get('BrainstormTmpDir')      : User brainstorm temporary directory (Default: <home>/.brainstorm/tmp/)
%    - bst_get('BrainstormTmpDir', isForcedDefault)   : User DEFAULT brainstorm temporary directory (<home>/.brainstorm/tmp/)
%    - bst_get('BrainstormDocDir')      : Doc folder folder of the Brainstorm distribution (may vary in compiled versions)
%    - bst_get('BrainstormDefaultsDir') : Defaults folder folder of the Brainstorm distribution (may vary in compiled versions)
%    - bst_get('BrainstormPluginDir')   : User home directory for brainstorm (<home>/.brainstorm/)
%    - bst_get('UserReportsDir')        : User reports directory (<home>/.brainstorm/reports/)
%    - bst_get('UserMexDir')            : User temporary directory (<home>/.brainstorm/mex/)
%    - bst_get('UserProcessDir')        : User custom processes directory (<home>/.brainstorm/process/)
%    - bst_get('UserDefaultsDir')       : User defaults directory (<home>/.brainstorm/defaults/)
%    - bst_get('UserPluginsDir')        : User plugins directory (<home>/.brainstorm/plugins/)
%    - bst_get('BrainstormDbFile')      : User brainstorm.mat file (<home>/.brainstorm/brainstorm.mat)
%    - bst_get('BrainstormDbDir')       : User database directory (contains all the brainstorm protocols)
%    - bst_get('DirDefaultSubject')     : Directory name of the default subject
%    - bst_get('DirDefaultStudy')       : Directory name of the default study for each subject
%    - bst_get('DirAnalysisInter')      : Directory name of the inter-subject analysis study 
%    - bst_get('DirAnalysisIntra')      : Directory name of the intra-subject analysis study (for each subject)
%    - bst_get('AnatomyDefaults')       : Get the contents of directory bstDir/defaults/anatomy
%    - bst_get('MniAtlasDefaults')      : Get the contents of directory bstDir/defaults/mniatlas
%    - bst_get('EegDefaults')           : Get the contents of directory bstDir/defaults/eeg
%    - bst_get('LastUsedDirs')          : Structure with all the last used directories (last used)
%    - bst_get('OsType', isMatlab=1)    : Get a string that describes the operating system (if isMatlab=1 return the Matlab/JVM platform, else return the real host system)
%    - bst_get('FileFilters', DataType) : Get the list of import filters for a specific data type
%    - bst_get('PluginCustomPath')      : Full custom path to all plugins
%    - bst_get('BrainSuiteDir')         : Full path to a local installation of BrainSuite
%    - bst_get('SpmTpmAtlas')           : Full path to the SPM atlas TPM.nii
%    - bst_get('PythonExe')             : Path to the python executable
%
% ====== PROTOCOLS ====================================================================
%    - bst_get('iProtocol')             : Indice of current protocol 
%    - bst_get('Protocol', ProtocolName): Return the indice of the protocol ProtocolName, or [] if it doesn't exist
%    - bst_get('ProtocolInfo')          : Definition structure for current protocol
%    - bst_get('ProtocolSubjects')      : Subjects list for current protocol
%    - bst_get('ProtocolStudies')       : Studies list for current protocol
%    - bst_get('isProtocolLoaded')      : 1 if the protocol has been loaded, 0 else
%    - bst_get('isProtocolModified')    : 1 if the protocol has been modified, 0 else
%
% ====== STUDIES ======================================================================
%    - bst_get('Study', StudyFileName)  : Get one study in current protocol with its file name
%    - bst_get('Study', iStudies)       : Get one or more studies
%    - bst_get('Study')                 : Get current study in current protocol
%    - bst_get('StudyCount')            : Get number of studies in the current protocol
%    - bst_get('StudyWithSubject',   SubjectFile)          : Find studies associated with a given subject file (WITHOUT the system studies ('intra_subject', 'default_study'))
%    - bst_get('StudyWithSubject',   ..., 'intra_subject') : Find studies ... INCLUDING 'intra_subject' study
%    - bst_get('StudyWithSubject',   ..., 'default_study') : Find studies ... INCLUDING 'default_study' study
%    - bst_get('StudyWithCondition', ConditionPath)        : Find studies for a given condition path
%    - bst_get('ChannelStudiesWithSubject', iSubjects)     : Get all the studies where there should be a channel file for a list of subjects
%    - bst_get('AnalysisIntraStudy', iSubject)    : Get the default analysis study for target subject
%    - bst_get('AnalysisInterStudy')              : Get the default analysis study for inter-subject analysis
%    - bst_get('DefaultStudy', iSubject)          : Get the default study for target subject (by subject indice)
%    - bst_get('DefaultStudy')                    : Get the global default study (common to all subjects)
%    - bst_get('DefaultStudy', SubjectFile)       : Get the default study for target subject (by filename)
%    - bst_get('ChannelFile', ChannelFile)        : Find a channel file in current protocol
%    - bst_get('ChannelFileForStudy', StudyFile/DataFile)  : Find a channel file in current protocol
%    - bst_get('ChannelForStudy',   iStudies)     : Return current Channel struct for target study
%    - bst_get('ChannelModalities', ChannelFile)  : Return displayable modalities for a Channel file
%    - bst_get('ChannelModalities', DataFile)     : Return displayable modalities for a Data/Results/Timefreq... file
%    - bst_get('ChannelDevice', ChannelFile)      : Return acquisistion device for a channel file
%    - bst_get('TimefreqDisplayModalities', TfFile)     : Get displayable modalities for a TF file based on recordings
%    - bst_get('HeadModelForStudy',   iStudy)           : Return current HeadModel struct for target study
%    - bst_get('HeadModelFile', HeadModelFile)          : Find a HeadModel file in current protocol
%    - bst_get('HeadModelFile', HeadModelFile, iStudies): Find a HeadModel file in current protocol
%    - bst_get('NoiseCovFile', NoiseCovFile)            : Find a noise covariance file file in current protocol
%    - bst_get('NoiseCovFile', NoiseCovFile, iStudies)  : Find a noise covariance file in current protocol
%    - bst_get('NoiseCovFile', DataCovFile)             : Find a data covariance file file in current protocol
%    - bst_get('NoiseCovFile', DataCovFile, iStudies)   : Find a data covariance file file in current protocol
%    - bst_get('DataFile',    DataFile)                 : Find a DataFile in current protocol
%    - bst_get('DataFile',    DataFile, iStudies)       : Find a DataFile in current protocol
%    - bst_get('DataForDataList', iStudy, DataListName) : Find all the DataFiles grouped by a data list
%    - bst_get('DataForStudy', iStudy)                  : Find all the Data files that are dependent on the channel/headmodel of a given study
%    - bst_get('DataForStudies', iStudies)
%    - bst_get('DataForChannelFile', ChannelFile)       : Find all the DataFiles that use the given ChannelFile
%    - bst_get('ResultsFile', ResultsFile)              : Find a ResultsFile in current protocol
%    - bst_get('ResultsFile', ResultsFile, iStudies)    : Find a ResultsFile in input studies
%    - bst_get('ResultsForDataFile', DataFile, iStudies): Find all results computed based on DataFile
%    - bst_get('StatFile', StatFile)                    : Find a StatFile in current protocol
%    - bst_get('StatFile', StatFile, iStudies)          : Find a StatFile in input studies
%    - bst_get('StatForDataFile',   DataFile, iStudies)
%    - bst_get('StatForDataFile',   DataFile)
%    - bst_get('TimefreqFile',      TimefreqFile)
%    - bst_get('TimefreqFile',      TimefreqFile, iStudies)
%    - bst_get('TimefreqForFile',   FileName, iStudies)   : Find all timefreq files computed based on target file
%    - bst_get('TimefreqForKernel', KernelFile) 
%    - bst_get('DipolesFile',       DipolesFile)
%    - bst_get('DipolesFile',       DipolesFile, iStudies)
%    - bst_get('DipolesForFile',    FileName, iStudies)   : Find all dipoles files computed based on target file
%    - bst_get('DipolesForKernel',  KernelFile) 
%    - bst_get('MatrixFile',        MatrixFile)
%    - bst_get('MatrixFile',        MatrixFile, iStudies)
%    - bst_get('MatrixForMatrixList', iStudy, MatrixListName)
%    - bst_get('AnyFile',           AnyFile)
%    - bst_get('AnyFile',           AnyFile, iStudies)
%    - bst_get('RelatedDataFile',   FileName)
%    - bst_get('RelatedDataFile',   FileName, iStudies)
%    - bst_get('GetFilenames')
%
% ====== SUBJECTS ======================================================================
%    - bst_get('Subject', SubjectFileName, isRaw) : Find a subject in current protocol with its file name
%    - bst_get('Subject', SubjectName, isRaw)     : Find a subject in current protocol with its name
%    - bst_get('Subject', iSubject)               : Get a subject (normal or default if iSubject==0)
%    - bst_get('Subject')                         : Get current subject in current protocol
%    - bst_get('SubjectCount')                    : Get number of studies in the current protocol
%    - bst_get('NormalizedSubjectName')           : Name of the subject with a normalized anatomy
%    - bst_get('NormalizedSubject')               : Get groupd analysis subject for the current protocol
%    - bst_get('ConditionsForSubject', SubjectFile)           : Find all conditions for a given subject
%    - bst_get('SurfaceFile',          SurfaceFile)           : Find a surface in current protocol
%    - bst_get('SurfaceFileByType',    iSubject,    SurfaceType) : Find surfaces with given type for subject #i (default only)
%    - bst_get('SurfaceFileByType',    SubjectFile, SurfaceType) : Find surfaces with given type for subject SubjectFile (default only)
%    - bst_get('SurfaceFileByType',    SurfaceName, SurfaceType) : Find surfaces with given type for subject that also has surface SurfaceName (default only)
%    - bst_get('SurfaceFileByType',    MriName,     SurfaceType) : Find surfaces with given type for subject that also has MRI MriName  (default only)
%    - bst_get('SurfaceFileByType',    ...,  ..., isDefaultOnly) : If 0, return all the surfaces of the given type, instead of only the default surface
%    - bst_get('MriFile',              MriFile)               : Find a MRI in current protocol
% 
% ====== GUI =================================================================
%    - bst_get('BstControls')    : Return main Brainstorm GUI structure
%    - bst_get('BstFrame')       : Return main Brainstorm JFrame
%    - bst_get('isGUI')          : Return 1 if the Brainstorm interface is displayed
%    - bst_get('GuiLevel')       : Return GUI level:  -1=server, 0=nogui, 1=normal, 2=autopilot
%    - bst_get('ScreenDef')      : Get screens configuration
%    - bst_get('DecorationSize') : Get dimensions of the windows decorations
%    - bst_get('Layout')         : Configuration of the main Brainstorm window
%    - bst_get('Layout', prop)   : Get one property in the layout properties
%    - bst_get('PanelContainer')                : Display list of registered panel containers
%    - bst_get('PanelContainer', ContainerName) : Get a panel container handle
%    - bst_get('Panel')                         : Display list of registered panels
%    - bst_get('Panel',         PanelName)      : Find a panel with its name
%    - bst_get('PanelControls', PanelName)      : Get the controls of a panel
%    - bst_get('Clipboard')       : Nodes that were copied from the interface
%    - bst_get('FigFont')         : Standard font size displayed in the figures
%    - bst_get('Font')            : Create a Java font, scaled for the operating system
%
% ====== CONFIGURATION =================================================================
%    - bst_get('Version')               : Brainstorm version
%    - bst_get('MatlabVersion')         : Matlab version (version number * 100, eg. 801)
%    - bst_get('MatlabReleaseName')     : Matlab version (release name, eg. "R2014a")
%    - bst_get('JavaVersion')           : Java version
%    - bst_get('isJavacomponent')       : Returns 1 if javacomponent is available (Matlab < 2019b), 0 otherwise
%    - bst_get('SystemMemory')          : Amount of memory available, in Mb
%    - bst_get('ByteOrder')             : {'l','b'} - Byte order used to read and save binary files 
%    - bst_get('TSDisplayMode')         : {'butterfly','column'}
%    - bst_get('ElectrodeConfig', Modality)  : Structure describing the current display values for SEEG/ECOG/EEG contacts
%    - bst_get('AutoUpdates')                : {0,1} - If 1, check automatically for updates at startup
%    - bst_get('ForceMatCompression')   : {0,1} - If 1, always save mat-files using the v7 format instead of v6
%    - bst_get('IgnoreMemoryWarnings')  : {0,1} - If 1, do not display memory warnings at the Brainstorm startup
%    - bst_get('ExpertMode')            : {0,1} - If 1, show advanced options that regular user do not see
%    - bst_get('DisplayGFP')            : {0,1} - If 1, the GFP is displayed on all the time series figures
%    - bst_get('DownsampleTimeSeries')  : {0,1,...} - If > 0, downsample dense time series for faster display
%    - bst_get('DisableOpenGL')         : {0,1,2} - If 1, do not use OpenGL renderer; if 2, use software OpenGL
%    - bst_get('InterfaceScaling')      : {100,125,150,...} - Scales the Brainstorm GUI by a fixed factor
%    - bst_get('GraphicsSmoothing')     : {0,1} - If 1, uses the graphics smoothing (Matlab >= 2014b)
%    - bst_get('SystemCopy')            : {0,1} - If 1, uses the system calls mv/cp instead of movefile/copyfile (Linux only)
%    - bst_get('JOGLVersion')           : {0,1,2}, Detect the current version of JOGL available in Matlab
%    - bst_get('DefaultFormats')        : Default formats for importing/exporting data, channels, ... (last used)
%    - bst_get('BFSProperties')         : Conductivities and thicknesses for 3-shell spherical forward model
%    - bst_get('ImportDataOptions')     : Import options for recordings
%    - bst_get('ImportEegRawOptions')   : Importation options for RAW EEG format
%    - bst_get('BugReportOptions')      : Bug reporter options
%    - bst_get('DefaultSurfaceDisplay') : Default display options for surfaces (smooth, data threshold, sulci map)
%    - bst_get('MagneticExtrapOptions') : Structure with the options for magnetic field extrapolation
%    - bst_get('DefaultFreqBands')
%    - bst_get('TimefreqOptions_morlet')
%    - bst_get('TimefreqOptions_fft')
%    - bst_get('TimefreqOptions_psd')
%    - bst_get('TimefreqOptions_hilbert')
%    - bst_get('OpenMEEGOptions')
%    - bst_get('DuneuroOptions')
%    - bst_get('GridOptions_headmodel')
%    - bst_get('GridOptions_dipfit')
%    - bst_get('UniformizeTimeSeriesScales') : {0,1} - If 1, the Y-axis of all the time series figures have the scale
%    - bst_get('FlipYAxis')                  : {0,1} - If 1, the recordings are displayed with the negative values UP
%    - bst_get('AutoScaleY')                 : {0,1} - If 1, the axis limits are updated when the figure is updated
%    - bst_get('ShowXGrid')                  : {0,1} - If 1, show the XGrid in the time series figures
%    - bst_get('ShowYGrid')                  : {0,1} - If 1, show the YGrid in the time series figures
%    - bst_get('ShowZeroLines')              : {0,1} - If 1, show the Y=0 lines in the columns view
%    - bst_get('ShowEventsMode')             : {'dot','line','none'}
%    - bst_get('Resolution')                 : [resX,resY] fixed resolutions for X and Y axes
%    - bst_get('FixedScaleY', Modality)      : Struct with the scales to impose on the recordings for the selected modality
%    - bst_get('XScale', XScale)             : {'log', 'linear'}
%    - bst_get('YScale', YScale)             : {'log', 'linear'}
%    - bst_get('UseSigProcToolbox')       : Use Matlab's Signal Processing Toolbox when available
%    - bst_get('RawViewerOptions', sFile) : Display options for RAW recordings, adapt for specific file
%    - bst_get('RawViewerOptions')        : Display options for RAW recordings
%    - bst_get('TopoLayoutOptions')       : Display options for 2DLayout display
%    - bst_get('StatThreshOptions')       : Options for online statistical thresholding
%    - bst_get('ContactSheetOptions')     : Display options for contact sheets
%    - bst_get('ProcessOptions')          : Options related with the data processing
%    - bst_get('CustomColormaps')         : Gets the list of user defined colormaps
%    - bst_get('MriOptions')              : Configuration for MRI display
%    - bst_get('DigitizeOptions')         : Digitizer options
%    - bst_get('ReadOnly')                : Read only interface
%    - bst_get('NodelistOptions')         : Structure with the options for file selection in the Process1 and Process2 panels
%    - bst_get('ResizeFunction')          : Get the appropriate resize function
%    - bst_get('groot')                   : Get the root graphic object
%    - bst_get('JFrame', hFig)            : Get the underlying java frame for a Matlab figure
%    - bst_get('LastPsdDisplayFunction')  : Display option of measure for spectrum (log, power, magnitude, etc.)
%    - bst_get('PlotlyCredentials')       : Get the credentials and URL to connect to plot.ly server
%    - bst_get('ExportBidsOptions')       : Additional metadata for BIDS export
%
% SEE ALSO bst_set

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2020 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Francois Tadel, 2008-2021
%          Martin Cousineau, 2017

%% ==== PARSE INPUTS ====
global GlobalData;
if ((nargin >= 1) && ischar(varargin{1}))
    contextName = varargin{1};
else
    return
end
% Initialize returned variable
argout1 = [];
argout2 = [];
argout3 = [];
argout4 = [];
argout5 = [];

% Get required context structure
switch contextName
%% ==== SUBFUNCTIONS =====
    case 'findFileInStudies'
        [argout1, argout2, argout3] = findFileInStudies(varargin{2:end});
        
%% ==== BRAINSTORM CONFIGURATION ====
    case 'Version'
        argout1 = GlobalData.Program.Version;
        
    case 'MatlabVersion'
        if ~exist('version','builtin')
            Version = 601;
        else
            % Get version number
            str_vers = version();
            vers = sscanf(str_vers, '%d.%d%*s');
            if isempty(vers) || any(~isnumeric(vers)) || any(isnan(vers))
                vers = 0;
            end
            Version = 100 * vers(1) + vers(2);
            % Get release name
            ipar = [find(str_vers == '('), find(str_vers == ')')];
        end
        argout1 = Version;
        
    case 'MatlabReleaseName'
        if ~exist('version','builtin')
            Release = 'Matlab6';
        else
            % Get version number
            str_vers = version();
            % Get release name
            ipar = [find(str_vers == '('), find(str_vers == ')')];
            Release = str_vers(ipar(1)+1:ipar(2)-1);
        end
        argout1 = Release;

    case 'JavaVersion'
        strver = char(java.lang.System.getProperty('java.version'));
        iDot = find(strver == '.');
        if (length(iDot) < 2)
            argout1 = str2num(strver);
        else 
            argout1 = str2num(strver(1:iDot(2)-1));
        end
        
    case 'isJavacomponent'
        % After Matlab 2019b, javacomponent() and JavaFrame property have been deprecated
        argout1 = (bst_get('MatlabVersion') <= 906);
        
    case 'SystemMemory'
        maxvar = [];
        totalmem = [];
        if ispc && (bst_get('MatlabVersion') >= 706)
            try
                % Get memory info
                usermem  = memory();
                maxvar   = round(usermem.MaxPossibleArrayBytes / 1024 / 1024);
                totalmem = round(usermem.MemAvailableAllArrays / 1024 / 1024);
            catch
                % Whatever...
            end
        end
        argout1 = maxvar;
        argout2 = totalmem;
            
    case 'BrainstormHomeDir'
        argout1 = GlobalData.Program.BrainstormHomeDir;
    case 'BrainstormDbDir'
        argout1 = GlobalData.DataBase.BrainstormDbDir;
    case 'UserDir'
        try
            if ispc
                userDir = getenv('USERPROFILE');
            else
                userDir = char(java.lang.System.getProperty('user.home'));
            end
        catch
            userDir = '';
        end
        if isempty(userDir)
            userDir = bst_get('BrainstormHomeDir');
        end
        argout1 = userDir;
        
    case 'BrainstormUserDir'
        bstUserDir = bst_fullfile(bst_get('UserDir'), '.brainstorm');
        if ~isdir(bstUserDir)
            res = mkdir(bstUserDir);
            if ~res
                error(['Cannot create Brainstorm user directory: "' bstUserDir '".']); 
            end
        end
        argout1 = bstUserDir;
        
    case 'BrainstormTmpDir'
        tmpDir = '';
        isForcedDefault = ((nargin >= 2) && varargin{2});
        % Default folder: userdir/tmp
        defDir = bst_fullfile(bst_get('BrainstormUserDir'), 'tmp');
        % If temporary directory is set in the preferences
        if ~isForcedDefault && isfield(GlobalData, 'Preferences') && isfield(GlobalData.Preferences, 'BrainstormTmpDir')
            tmpDir = GlobalData.Preferences.BrainstormTmpDir;
        end 
        % Else: use directory userdir/tmp
        if isempty(tmpDir)
            tmpDir = defDir;
        end
        % Create directory if it does not exist yet
        if ~isdir(tmpDir)
            res = mkdir(tmpDir);
            if ~res && ~strcmp(tmpDir, defDir)
                disp(['BST> Cannot create Brainstorm temporary directory: "' tmpDir '".']); 
                disp(['BST> Using default directory instead: "' defDir '".']);
                % Create temporary folder
                tmpDir = defDir;
                if ~isdir(tmpDir)
                    res = mkdir(tmpDir);
                else
                    res = 1;
                end
            end
            % Error: cannot create any temporary folder
            if ~res
                error(['Cannot create Brainstorm temporary directory: "' tmpDir '".']); 
            end
        end
        argout1 = tmpDir;
        
    case 'BrainstormDocDir'
        docDir = bst_fullfile(GlobalData.Program.BrainstormHomeDir, 'doc');
        if ~exist(docDir, 'file')
            % Matlab compiler >= 2018b stores 'doc' under 'bst_javabuil'
            docDir = bst_fullfile(GlobalData.Program.BrainstormHomeDir, 'bst_javabuil', 'doc');
            if ~exist(docDir, 'file')
                docDir = '';
                disp('BST> Could not find "doc" folder.');
                disp(['BST> BrainstormHomeDir = ' GlobalData.Program.BrainstormHomeDir]);
            end
        end
        argout1 = docDir;

    case 'BrainstormDefaultsDir'
        defaultsDir = bst_fullfile(GlobalData.Program.BrainstormHomeDir, 'defaults');
        if ~exist(defaultsDir, 'file')
            % Matlab compiler >= 2018b stores 'defaults' under 'bst_javabuil'
            defaultsDir = bst_fullfile(GlobalData.Program.BrainstormHomeDir, 'bst_javabuil', 'defaults');
            if ~exist(defaultsDir, 'file')
                defaultsDir = '';
                disp('BST> Could not find "defaults" folder.');
                disp(['BST> BrainstormHomeDir = ' GlobalData.Program.BrainstormHomeDir]);
            end
        end
        argout1 = defaultsDir;
        
    case 'UserReportsDir'
        reportDir = bst_fullfile(bst_get('BrainstormUserDir'), 'reports');
        if ~isdir(reportDir)
            res = mkdir(reportDir);
            if ~res
                error(['Cannot create user reports directory: "' reportDir '".']); 
            end
        end
        argout1 = reportDir;
        
    case 'UserMexDir'
        mexDir = bst_fullfile(bst_get('BrainstormUserDir'), 'mex');
        if ~isdir(mexDir)
            res = mkdir(mexDir);
            if ~res
                error(['Cannot create Brainstorm mex-files directory: "' mexDir '".']); 
            end
        end
        argout1 = mexDir;
        
    case 'UserProcessDir'
        processDir = bst_fullfile(bst_get('BrainstormUserDir'), 'process');
        if ~isdir(processDir)
            res = mkdir(processDir);
            if ~res
                error(['Cannot create Brainstorm custom processes directory: "' processDir '".']); 
            end
        end
        argout1 = processDir;
        
    case 'UserDefaultsDir'
        defDir = bst_fullfile(bst_get('BrainstormUserDir'), 'defaults');
        defDirAnat = bst_fullfile(defDir, 'anatomy');
        defDirEeg  = bst_fullfile(defDir, 'eeg');
        if ~isdir(defDir)
            res = mkdir(defDir);
            if ~res
                error(['Cannot create user templates directory: "' defDir '".']); 
            end
        end
        if ~isdir(defDirAnat)
            res = mkdir(defDirAnat);
            if ~res
                error(['Cannot create user templates directory: "' defDirAnat '".']); 
            end
        end
        if ~isdir(defDirEeg)
            res = mkdir(defDirEeg);
            if ~res
                error(['Cannot create user templates directory: "' defDirEeg '".']); 
            end
        end
        argout1 = defDir;
        
    case 'UserPluginsDir'
        pluginsDir = bst_fullfile(bst_get('BrainstormUserDir'), 'plugins');
        if ~isdir(pluginsDir)
            res = mkdir(pluginsDir);
            if ~res
                error(['Cannot create plugins directory: "' pluginsDir '".']); 
            end
        end
        argout1 = pluginsDir;
        
    case 'BrainstormDbFile'
        argout1 = bst_fullfile(bst_get('BrainstormUserDir'), 'brainstorm.mat');

%% ==== PROTOCOL ====
    case 'iProtocol'
        if isempty(GlobalData.DataBase.iProtocol)
            argout1 = 0;
        else
            argout1 = GlobalData.DataBase.iProtocol;
        end
        
    case 'Protocol'
        if ~isempty(GlobalData.DataBase.ProtocolInfo)
            argout1 = find(strcmpi({GlobalData.DataBase.ProtocolInfo.Comment}, varargin{2}));
        else
            argout1 = [];
        end
        
    case {'ProtocolInfo', 'ProtocolSubjects', 'ProtocolStudies', 'isProtocolLoaded', 'isProtocolModified'}
        argout2 = GlobalData.DataBase.iProtocol;
        % No protocol: return empty matrix
        if isempty(argout2) || ~isnumeric(argout2) || argout2 == 0
            return;
        end
        % Check index integrity
        if ((argout2 <= 0) || (argout2 > length(GlobalData.DataBase.ProtocolInfo))), warning('Brainstorm:InvalidIndex', 'Invalid index'), return, end
        % Get requested protocol structure
        argout1 = GlobalData.DataBase.(contextName)(argout2);

        
%% ==== STUDY ====
    % Usage: [sStudy, iStudy] = bst_get('Study', StudyFileName)
    %        [sStudy, iStudy] = bst_get('Study')
    %        [sStudy, iStudy] = bst_get('Study', iStudies)
    case 'Study'
        if isempty(GlobalData.DataBase.iProtocol) || (GlobalData.DataBase.iProtocol == 0)
            return;
        end
        % Get list of current protocols
        ProtocolInfo    = GlobalData.DataBase.ProtocolInfo(GlobalData.DataBase.iProtocol);
        ProtocolStudies = GlobalData.DataBase.ProtocolStudies(GlobalData.DataBase.iProtocol);
        % Get list of current protocol studies
        if isempty(ProtocolStudies) || isempty(ProtocolInfo)
            return;
        end
        % ===== PARSE INPUTS =====
        if (nargin < 2)
            % Call: bst_get('Study');
            iStudies = ProtocolInfo.iStudy;
            StudyFileName = [];
        elseif (isnumeric(varargin{2})) 
            iStudies = varargin{2};
            StudyFileName = [];
        elseif (ischar(varargin{2}))
            iStudies = [];
            StudyFileName = strrep(varargin{2}, ProtocolInfo.STUDIES, '');
        end
        % Indices
        iAnalysisStudy = -2;    % CANNOT USE -1 => DISABLES SEARCH FUNCTIONS
        iDefaultStudy  = -3;
        % Indices > 0: normal studies indiced in ProtocolStudies.Study array
            
        % ===== GET STUDY BY INDEX =====
        % Call: bst_get('Study', iStudies);
        if ~isempty(iStudies)
            argout1 = repmat(ProtocolStudies.DefaultStudy, 0);
            % Get analysis study
            iTargetAnalysis = find(iStudies == iAnalysisStudy);
            if ~isempty(iTargetAnalysis)
                argout1(iTargetAnalysis) = repmat(ProtocolStudies.AnalysisStudy, size(iTargetAnalysis));
                argout2(iTargetAnalysis) = repmat(iAnalysisStudy, size(iTargetAnalysis));
            end
            % Get default study
            iTargetDefault = find(iStudies == iDefaultStudy);
            if ~isempty(iTargetDefault) 
                argout1(iTargetDefault) = repmat(ProtocolStudies.DefaultStudy, size(iTargetDefault));
                argout2(iTargetDefault) = repmat(iDefaultStudy, size(iTargetDefault));
            end
            % Get normal studies
            iTargetNormal = find((iStudies >= 1) & (iStudies <= length(ProtocolStudies.Study)));
            if ~isempty(iTargetNormal)
                argout1(iTargetNormal) = ProtocolStudies.Study(iStudies(iTargetNormal));
                argout2(iTargetNormal) = iStudies(iTargetNormal);
            end
            % Error
            if isempty(argout1)
                GlobalData.DataBase.ProtocolInfo(GlobalData.DataBase.iProtocol).iStudy = [];
            end
            
        % ===== GET STUDY BY FILENAME =====
        % Call: bst_get('Study', StudyFileName);
        elseif ~isempty(StudyFileName)
            % NORMAL STUDY
            iStudy = find(file_compare({ProtocolStudies.Study.FileName}, StudyFileName), 1);
            % If a study is found : return it
            if ~isempty(iStudy)
                argout1 = ProtocolStudies.Study(iStudy);
                argout2   = iStudy;
            % DEFAULT STUDY
            elseif ~isempty(ProtocolStudies.DefaultStudy) && file_compare({ProtocolStudies.DefaultStudy.FileName}, StudyFileName)
                argout1 = ProtocolStudies.DefaultStudy;
                argout2   = iDefaultStudy;
            % ANALYSIS STUDY
            elseif ~isempty(ProtocolStudies.AnalysisStudy) && file_compare({ProtocolStudies.AnalysisStudy.FileName}, StudyFileName)
                argout1 = ProtocolStudies.AnalysisStudy;
                argout2   = iAnalysisStudy;
            end
        else
            return
        end
 
        
%% ==== STUDY WITH SUBJECT FILE ====
    % Usage : [sStudies, iStudies] = bst_get('StudyWithSubject', SubjectFile) : WITHOUT the system studies ('intra_subject', 'default_study')
    %         [sStudies, iStudies] = bst_get(..., 'intra_subject', 'default_study') : WITH the system studies: 'intra_subject' | 'default_study'
    case 'StudyWithSubject'
        % Parse inputs
        if (nargin < 2) || ~ischar(varargin{2})
            error('Invalid call to bst_get()');
        end
        if (nargin > 2)
            IntraStudies   = any(strcmpi(varargin(3:end), 'intra_subject'));
            DefaultStudies = any(strcmpi(varargin(3:end), 'default_study'));
        else
            IntraStudies   = 0;
            DefaultStudies = 0;
        end
        SubjectFile = {varargin{2}};
        if isempty(GlobalData.DataBase.iProtocol) || (GlobalData.DataBase.iProtocol == 0)
            return;
        end
        % Get list of current protocol description
        ProtocolInfo    = GlobalData.DataBase.ProtocolInfo(GlobalData.DataBase.iProtocol);
        ProtocolStudies = GlobalData.DataBase.ProtocolStudies(GlobalData.DataBase.iProtocol);
        if isempty(ProtocolStudies) || isempty(ProtocolInfo)
            return;
        end
        
        % Get default subject
        sDefaultSubject = bst_get('Subject', 0);
        % If SubjectFile is the default subject filename
        if ~isempty(sDefaultSubject) && ~isempty(sDefaultSubject.FileName) && file_compare( SubjectFile{1}, sDefaultSubject.FileName)
            % Get all the subjects files that use default anatomy
            ProtocolSubjects = GlobalData.DataBase.ProtocolSubjects(GlobalData.DataBase.iProtocol);
            iSubjectUseDefaultAnat = find([ProtocolSubjects.Subject.UseDefaultAnat]);
            if isempty(iSubjectUseDefaultAnat)
                return
            end
            SubjectFile = {ProtocolSubjects.Subject(iSubjectUseDefaultAnat).FileName};
            % Also updates inter-subject node
            isInterSubject = 1;
        else
            isInterSubject = 0;
        end
        % Search all the current protocol's studies
        iStudies = [];
        for i=1:length(SubjectFile)
            iStudies = [iStudies, find(file_compare({ProtocolStudies.Study.BrainStormSubject}, SubjectFile{i}))];
        end
        % Return results
        if ~isempty(iStudies)
            % Remove "analysis_intra" and "default_study" studies from list
            if ~IntraStudies
                iStudies(strcmpi({ProtocolStudies.Study(iStudies).Name}, bst_get('DirAnalysisIntra'))) = [];
            end
            if ~DefaultStudies
                iStudies(strcmpi({ProtocolStudies.Study(iStudies).Name}, bst_get('DirDefaultStudy'))) = [];
            end
            % Return studies
            argout1 = ProtocolStudies.Study(iStudies);
            argout2   = iStudies;
        else
            argout1 = repmat(db_template('Study'), 0);
            argout2   = [];
        end
        % Add inter-subject node, if needed
        if isInterSubject
            [sInterStudy, iInterStudy] = bst_get('AnalysisInterStudy');
            argout1 = [argout1, sInterStudy];
            argout2   = [argout2,   iInterStudy];
        end
              
        
%% ==== STUDY WITH CONDITION PATH ====
    % USAGE:  [sStudies, iStudies] = bst_get('StudyWithCondition', ConditionPath)
    %
    % INPUT: ConditionPath
    %    - 'SubjectName/conditionName'  : Target condition for the specified subject
    %    - 'SubjectName/@intra'         : Intra-subject condition for the subject
    %    - 'SubjectName/@default_study' : Default condition for the subject (where the subject's shared files are stored)
    %    - '*/conditionName'            : Target condition for all the subjects
    %    - '@inter'                     : Inter-subject condition
    %    - '@default_study'             : Protocol's default condition (where the protocol's shared files are stored)
    
    case 'StudyWithCondition'
        % Parse inputs
        if (nargin ~= 2) || ~ischar(varargin{2})
            error('Invalid call to bst_get()');
        end
        ConditionPath = varargin{2};
        % Get list of current protocol description
        ProtocolInfo    = GlobalData.DataBase.ProtocolInfo(GlobalData.DataBase.iProtocol);
        ProtocolStudies = GlobalData.DataBase.ProtocolStudies(GlobalData.DataBase.iProtocol);
        if isempty(ProtocolStudies) || isempty(ProtocolInfo)
            return;
        end

        % ConditionPath = @inter
        if strcmpi(ConditionPath, '@inter')
            iStudy     = -2;
            argout2   = iStudy;
            argout1 = bst_get('Study', iStudy);
        % ConditionPath = @default_study
        elseif strcmpi(ConditionPath, '@default_study')
            iStudy     = -3;
            argout2   = iStudy;
            argout1 = bst_get('Study', iStudy);
        % ConditionPath = SubjectName/ConditionName
        else 
            % Get subject and condition names
            condSplit = str_split(ConditionPath);
            if (length(condSplit) ~= 2)
                error('Invalid condition path.');
            end
            SubjectName = condSplit{1};
            ConditionName = condSplit{2};

            % If first element is '*', search for condition in all the studies
            if (SubjectName(1) == '*')
                iStudies = 1:length(ProtocolStudies.Study);
            % Else : search for condition only in studies that are linked to the subject specified in the ConditionPath
            else
                % Get subject
                sSubject = bst_get('Subject', SubjectName, 1);
                if isempty(sSubject)
                    return;
                end
                iStudies = find(file_compare({ProtocolStudies.Study.BrainStormSubject}, sSubject.FileName));
            end
            % Nothing to search
            if isempty(iStudies)
                return
            end

            % Search all the current protocol's studies
            iStudies = iStudies(cellfun(@(c)isequal(ConditionName, c), [ProtocolStudies.Study(iStudies).Condition]));
            % Return results
            if ~isempty(iStudies)
                % Remove "analysis_intra" and "default_study" studies from list
                if ~strcmpi(ConditionName, '@intra')
                    iStudies(strcmpi({ProtocolStudies.Study(iStudies).Condition}, bst_get('DirAnalysisIntra'))) = [];
                end
                if ~strcmpi(ConditionName, '@default_study')
                    iStudies(strcmpi({ProtocolStudies.Study(iStudies).Condition}, bst_get('DirDefaultStudy'))) = [];
                end
                % Sort by subject
                if (length(iStudies) > 1)
                    SubjNameList = cell(1,length(iStudies));
                    % For each study, get subject name
                    for i = 1:length(iStudies)
                        sSubject = bst_get('Subject', ProtocolStudies.Study(iStudies(i)).BrainStormSubject);
                        SubjNameList{i} = sSubject.Name;
                    end
                    % Sort subjects names
                    [sortSubjList, iSort] = sort(SubjNameList);
                    % Apply same sorting to studies
                    iStudies = iStudies(iSort);
                end
                % Return studies
                argout1 = ProtocolStudies.Study(iStudies);
                argout2   = iStudies;
            else
                argout1 = repmat(db_template('Study'), 0);
                argout2   = [];
            end
        end

%% ==== CHANNEL STUDIES WITH SUBJECT ====
    % Usage: iStudies = bst_get('ChannelStudiesWithSubject', iSubjects, 'NoIntra')
    case 'ChannelStudiesWithSubject'
        % Parse inputs
        if (nargin >= 2) && isnumeric(varargin{2})
            iSubjects = varargin{2};
        else
            error('Invalid call to bst_get()');
        end
        if (nargin == 3) && strcmpi(varargin{3}, 'NoIntra')
            NoIntra = 1;
        else
            NoIntra = 0;
        end
        % Process all subjects
        iStudies = [];
        for i=1:length(iSubjects)
            iSubject = iSubjects(i);
            sSubject = bst_get('Subject', iSubject, 1);
            % No subject: error
            if isempty(sSubject) 
                continue
            % If subject uses default channel file    
            elseif (sSubject.UseDefaultChannel ~= 0)
                % Get default study for this subject
                [tmp___, iStudiesNew] = bst_get('DefaultStudy', iSubject);
                iStudies = [iStudies, iStudiesNew];
            % Else: get all the studies belonging to this subject
            else
                if NoIntra
                    [tmp___, iStudiesNew] = bst_get('StudyWithSubject', sSubject.FileName);
                else
                    [tmp___, iStudiesNew] = bst_get('StudyWithSubject', sSubject.FileName, 'intra_subject');
                end
                iStudies = [iStudies, iStudiesNew];
            end
        end
        argout1 = iStudies;
    
        
%% ==== STUDIES COUNT ====
    % Usage: [nbStudies] = bst_get('StudyCount')
    case 'StudyCount'
        % Nothing
        if isempty(GlobalData.DataBase.iProtocol) || (GlobalData.DataBase.iProtocol == 0)
            argout1 = 0;
            return;
        end
        % Get list of current protocol studies
        ProtocolStudies = GlobalData.DataBase.ProtocolStudies(GlobalData.DataBase.iProtocol);
        argout1 = length(ProtocolStudies.Study);

%% ==== SUBJECTS COUNT ====
    % Usage: [nbSubjects] = bst_get('SubjectCount')
    case 'SubjectCount'
        if isempty(GlobalData.DataBase.iProtocol) || (GlobalData.DataBase.iProtocol == 0)
            argout1 = 0;
            return;
        end
        % Get list of current protocol studies
        ProtocolSubjects = GlobalData.DataBase.ProtocolSubjects(GlobalData.DataBase.iProtocol);
        argout1 = length(ProtocolSubjects.Subject);
        
%% ==== NORMALIZED SUBJECT ====
    case 'NormalizedSubject'
        % Get normalized subject name
        normSubjName = 'Group_analysis';
        % Try to get normalized subject
        [sNormSubj, iNormSubj] = bst_get('Subject', normSubjName, 0);
        % If normalized subject does not exist: create it
        if isempty(sNormSubj)
            % Always use default anatomy
            UseDefaultAnat = 1;
            % If all the subjects use a global channel file: Use global default as well
            if ~isempty(GlobalData.DataBase.ProtocolSubjects(GlobalData.DataBase.iProtocol).Subject) && ...
                all([GlobalData.DataBase.ProtocolSubjects(GlobalData.DataBase.iProtocol).Subject.UseDefaultChannel] == 2)
                UseDefaultChannel = 2;
            else
                UseDefaultChannel = 1;
            end    
            % Create subject
            [sNormSubj, iNormSubj] = db_add_subject(normSubjName, [], UseDefaultAnat, UseDefaultChannel);
            % Get full subject structure
            [sNormSubj, iNormSubj] = bst_get('Subject', normSubjName, 0);
        end
        argout1 = sNormSubj;
        argout2 = iNormSubj;
        
        
%% ==== ANALYSIS STUDY (INTRA) ====
    % Usage: [sAnalStudy, iAnalStudy] = bst_get('AnalysisIntraStudy', iSubject) 
    case 'AnalysisIntraStudy'
        % Parse inputs
        if (nargin == 2) && isnumeric(varargin{2})
            iSubject = varargin{2};
        else
            error('Invalid call to bst_get()');
        end
        % Get subject
        sSubject = bst_get('Subject', iSubject, 1);
        % Get studies related to subject
        [sSubjStudies, iSubjStudies] = bst_get('StudyWithSubject', sSubject.FileName, 'intra_subject');
        % Look for the 'AnalysisIntra' study
        iFound = find(cellfun(@(c)ismember(bst_get('DirAnalysisIntra'), c), {sSubjStudies.Condition}));
        iAnalStudy = iSubjStudies(iFound);
        sAnalStudy = sSubjStudies(iFound);
        % Return found structure
        argout1 = sAnalStudy;
        argout2 = iAnalStudy;        
        
        
%% ==== ANALYSIS STUDY (INTER) ====
    % Usage: [sAnalStudyInter, iAnalStudyInter] = bst_get('AnalysisInterStudy') 
    case 'AnalysisInterStudy'
        iAnalStudyInter = -2;
        [argout1, argout2] = bst_get('Study', iAnalStudyInter);
        
       
%% ==== DEFAULT STUDY ====
    % Usage: [sDefaulStudy, iDefaultStudy] = bst_get('DefaultStudy', iSubject)
    %        [sDefaulStudy, iDefaultStudy] = bst_get('DefaultStudy')           : iSubject=0
    %        [sDefaulStudy, iDefaultStudy] = bst_get('DefaultStudy', SubjectFile)
    case 'DefaultStudy'
        if isempty(GlobalData.DataBase.iProtocol) || (GlobalData.DataBase.iProtocol == 0)
            return;
        end
        % Parse inputs
        if (nargin == 1)
            iSubject = 0;
        elseif (nargin == 2) && isnumeric(varargin{2})
            iSubject = varargin{2};
        elseif (nargin == 2) && ischar(varargin{2})
            SubjectFile = varargin{2};
            % Get subject attached to study
            [sSubject, iSubject] = bst_get('Subject', SubjectFile, 1);
            if isempty(sSubject) || ~sSubject.UseDefaultChannel
                return;
            end
        else
            error('Invalid call to bst_get()');
        end
        % === DEFAULT SUBJECT ===
        % => Return global default study
        if (iSubject == 0)
            % Get protocol's studies
            ProtocolStudies = GlobalData.DataBase.ProtocolStudies(GlobalData.DataBase.iProtocol);
            % Return Global default study
            argout1 = ProtocolStudies.DefaultStudy;
            argout2   = -3;
        % === NORMAL SUBJECT ===
        else
            % Get subject
            sSubject = bst_get('Subject', iSubject, 1);
            % === GLOBAL DEFAULT STUDY ===
            if sSubject.UseDefaultChannel == 2
                % Get protocol's studies
                ProtocolStudies = GlobalData.DataBase.ProtocolStudies(GlobalData.DataBase.iProtocol);
                % Return Global default study
                argout1 = ProtocolStudies.DefaultStudy;
                argout2   = -3;
            % === SUBJECT'S DEFAULT STUDY ===
            elseif sSubject.UseDefaultChannel == 1
                % Get studies related to subject
                [sSubjStudies, iSubjStudies] = bst_get('StudyWithSubject', sSubject.FileName, 'default_study');
                % Look for the 'DefaultStudy' study
                iFound = find(cellfun(@(c)ismember(bst_get('DirDefaultStudy'), c), {sSubjStudies.Condition}));
                iDefaultStudy = iSubjStudies(iFound);
                sDefaultStudy = sSubjStudies(iFound);
                % Return found structure
                argout1 = sDefaultStudy;
                argout2   = iDefaultStudy;        
            end
        end
        
        
        
%% ==== SUBJECT ====
    % Usage : [sSubject, iSubject] = bst_get('Subject', iSubject,        isRaw)
    %         [sSubject, iSubject] = bst_get('Subject', SubjectFileName, isRaw);
    %         [sSubject, iSubject] = bst_get('Subject', SubjectName,     isRaw);
    %         [sSubject, iSubject] = bst_get('Subject');
    % If isRaw is set: force to return the real brainstormsubject description
    % (ignoring wether it uses protocol's default anatomy or not)
    case 'Subject' 
        if isempty(GlobalData.DataBase.iProtocol) || (GlobalData.DataBase.iProtocol == 0)
            return;
        end
         % Get list of current protocol subjects
        ProtocolSubjects = GlobalData.DataBase.ProtocolSubjects(GlobalData.DataBase.iProtocol);
        sSubject = [];
        SubjectName = [];
        SubjectFileName = [];
        % ISRAW parameter
        if (nargin < 3)
            isRaw = 0;
        else
            isRaw = varargin{3};
        end
        % Call: bst_get('subject', iSubject, isRaw);
        if (nargin >= 2) && isnumeric(varargin{2})
            iSubject = varargin{2};
            if (iSubject > length(ProtocolSubjects.Subject))
                error('Invalid subject indice.');
            end
            % If required subject is default subject (iSubject = 0)
            if (iSubject == 0)
                % Default subject available
                if ~isempty(ProtocolSubjects.DefaultSubject)
                    sSubject = ProtocolSubjects.DefaultSubject;
                % Default subject not available
                else
                    return
                end
            % Normal subject 
            else
                sSubject = ProtocolSubjects.Subject(iSubject);
            end
            
        % Call: bst_get('subject', SubjectFileName, isRaw);
        % Call: bst_get('subject', SubjectName,     isRaw);
        elseif (nargin >= 2) && isempty(varargin{2})
            % If study name is empty: use DefaultSubject
            SubjectFileName = ProtocolSubjects.DefaultSubject.FileName;
        elseif (nargin >= 2) && (ischar(varargin{2}))
            [fName, fBase, fExt] = bst_fileparts(varargin{2});
            % Argument is a Matlab .mat filename
            if strcmpi(fExt, '.mat')
                SubjectFileName = varargin{2};
            % Else : assume argument is a directory
            else
                SubjectName = file_standardize(varargin{2});
            end
                            
        % Call: bst_get('subject');   => looking for current subject 
        elseif (nargin < 2)
            % Get current subject filename in current study
            sStudy = bst_get('Study');
            if isempty(sStudy)
                return
            end
            SubjectFileName = sStudy.BrainStormSubject;
            % If study's subject is not defined, get DefaultSubject
            if isempty(SubjectFileName) && ~isempty(ProtocolSubjects.DefaultSubject)
                SubjectFileName = ProtocolSubjects.DefaultSubject.FileName;
            end
        else
            error('Invalid call to bst_get()');
        end

        % If Subject is defined by its filename/name
        if isempty(sSubject)
            % Look in Default Subject
            if ~isempty(ProtocolSubjects.DefaultSubject) && (file_compare(ProtocolSubjects.DefaultSubject.FileName, SubjectFileName) ...
                                                          || strcmpi(ProtocolSubjects.DefaultSubject.Name, SubjectName))
                sSubject = ProtocolSubjects.DefaultSubject;
                iSubject = 0;
            % If not found : find target subject file name in normal subjects
            elseif ~isempty(SubjectFileName)
                iSubject = find(file_compare({ProtocolSubjects.Subject.FileName}, SubjectFileName), 1);
                sSubject = ProtocolSubjects.Subject(iSubject);
            elseif ~isempty(SubjectName)
                iSubject = find(file_compare({ProtocolSubjects.Subject.Name}, SubjectName), 1);
                sSubject = ProtocolSubjects.Subject(iSubject);
            else
                error('Subject name not specified.');
            end
        end
        
        % Return found subject
        if ~isempty(iSubject) && ~isempty(sSubject)
            % If subject uses default subject
            if sSubject.UseDefaultAnat && ~isRaw && ~isempty(ProtocolSubjects.DefaultSubject) && ~isempty(ProtocolSubjects.DefaultSubject.FileName)
                % Return default subject (WITH REAL SUBJECT'S NAME)
                argout1                   = ProtocolSubjects.DefaultSubject;
                argout1.Name              = sSubject.Name;
                argout1.UseDefaultAnat    = sSubject.UseDefaultAnat;
                argout1.UseDefaultChannel = sSubject.UseDefaultChannel;
                argout2                     = iSubject;
            % Else, return found subject
            else
                argout1  = sSubject;
                argout2    = iSubject;
            end
        end

        
%% ==== SURFACE FILE ====
    % Usage : [sSubject, iSubject, iSurface] = bst_get('SurfaceFile', SurfaceFile)
    case 'SurfaceFile'
        % No protocol in database
        if isempty(GlobalData) || isempty(GlobalData.DataBase) || isempty(GlobalData.DataBase.iProtocol) || (GlobalData.DataBase.iProtocol == 0)
            return;
        end
        % Get list of current protocol subjects
        ProtocolSubjects = GlobalData.DataBase.ProtocolSubjects(GlobalData.DataBase.iProtocol);
        if isempty(ProtocolSubjects)
            return
        end;
        
        % Parse inputs
        if (nargin == 2)
            SurfaceFile = varargin{2};
        else
            error('Invalid call to bst_get().');
        end

        % Remove SUBJECTS path from SurfaceFile
        SurfaceFile = file_short(SurfaceFile);
        % Look for surface file in DefaultSubject
        if ~isempty(ProtocolSubjects.DefaultSubject)
            % Find the first surface that matches the SurfaceFile
            iSurface = find(file_compare(SurfaceFile, {ProtocolSubjects.DefaultSubject.Surface.FileName}), 1);
            % If a surface was found in default subject : return it
            if ~isempty(iSurface)
                argout1  = ProtocolSubjects.DefaultSubject;
                argout2    = 0;
                argout3 = iSurface;
                return
            end
        end
        % Look for surface file in all the surfaces of all subjects
        for iSubj = 1:length(ProtocolSubjects.Subject)
            % Find the first surface that matches the SurfaceFile
            iSurface = find(file_compare(SurfaceFile, {ProtocolSubjects.Subject(iSubj).Surface.FileName}), 1);
            % If a surface was found in current subject : return it
            if ~isempty(iSurface)
                argout1  = ProtocolSubjects.Subject(iSubj);
                argout2    = iSubj;
                argout3 = iSurface;
                return
            end
        end
            
        
%% ==== SURFACE FILE BY TYPE ====
    % Usage : [sSurface, iSurface] = bst_get('SurfaceFileByType', iSubject,    SurfaceType)
    %         [sSurface, iSurface] = bst_get('SurfaceFileByType', SubjectFile, SurfaceType)
    %         [sSurface, iSurface] = bst_get('SurfaceFileByType', SurfaceFile, SurfaceType)
    %         [sSurface, iSurface] = bst_get('SurfaceFileByType', MriFile,     SurfaceType)
    %         [sSurface, iSurface] = bst_get('SurfaceFileByType', ...,         SurfaceType, isDefaultOnly)
    case 'SurfaceFileByType'
        % By default: return only the default surfaces of the category
        if (nargin < 4)
            isDefaultOnly = 1;
        else
            isDefaultOnly = varargin{4};
        end
        % Get subject
        if isempty(varargin{2})
            % Get default subject
            sSubject = bst_get('Subject', 0);
        elseif ischar(varargin{2})
            FileName = varargin{2};
            sSubject = bst_get('AnyFile', FileName);
        else
            iSubject = varargin{2};
            sSubject = bst_get('Subject', iSubject);
        end
        % Error handling
        if isempty(sSubject)
            disp('BST> Warning: Subject not found.');
            return;
        elseif isempty(sSubject.Surface)
            return;
        end
        SurfaceType = varargin{3};
        
        % === RETURN ONLY DEFAULTS ===
        if isDefaultOnly
            % Look for required surface type
            field = ['i' SurfaceType];
            if ~isfield(sSubject, field) || isempty(sSubject.(field))
                return
            end
            argout1 = sSubject.Surface(sSubject.(field));
            argout2 = sSubject.(field); 
        % === RETURN ALL THE SURFACES ===
        else
            % Build the list of tagged surfaces
            fileTag = ['_' lower(SurfaceType)];
            iTargetList = find(cellfun(@(c)~isempty(strfind(c, fileTag)), {sSubject.Surface.FileName}));
            % Put the default cortex on top of the list
            iDefaults = intersect([sSubject.iCortex, sSubject.iScalp, sSubject.iInnerSkull, sSubject.iOuterSkull, sSubject.iFibers, sSubject.iFEM], iTargetList);
            if ~isempty(iDefaults)
                iTargetList = [iDefaults, setdiff(iTargetList, iDefaults)];
            end
            % Return all cortex surfaces
            argout1 = sSubject.Surface(iTargetList);
            argout2 = iTargetList; 
        end
        
        
%% ==== MRI FILE ====
    % Usage : [sSubject, iSubject, iMri] = bst_get('MriFile', MriFile)
    case 'MriFile'
        % No protocol in database
        if isempty(GlobalData.DataBase.iProtocol) || (GlobalData.DataBase.iProtocol == 0)
            return;
        end
        % Get list of current protocol subjects
        ProtocolSubjects = GlobalData.DataBase.ProtocolSubjects(GlobalData.DataBase.iProtocol);
        if isempty(ProtocolSubjects)
            return
        end
        
        % Parse inputs
        if (nargin == 2)
            MriFile = varargin{2};
        else
            error('Invalid call to bst_get().');
        end

        % Remove SUBJECTS path from MriFile
        MriFile = file_short(MriFile);
        % Look for MRI file in DefaultSubject
        if ~isempty(ProtocolSubjects.DefaultSubject)
            % Find the first MRI that matches the MriFile
            iMri = find(file_compare(MriFile, {ProtocolSubjects.DefaultSubject.Anatomy.FileName}), 1);
            % If a MRI was found in default subject : return it
            if ~isempty(iMri)
                argout1  = ProtocolSubjects.DefaultSubject;
                argout2    = 0;
                argout3 = iMri;
                return
            end
        end
        % Look for MRI file in all the MRIs of all subjects
        for iSubj = 1:length(ProtocolSubjects.Subject)
            % Find the first MRI that matches the MriFile
            iMri = find(file_compare(MriFile, {ProtocolSubjects.Subject(iSubj).Anatomy.FileName}), 1);
            % If a MRI was found in current subject : return it
            if ~isempty(iMri)
                argout1  = ProtocolSubjects.Subject(iSubj);
                argout2    = iSubj;
                argout3 = iMri;
                return
            end
        end
        
        
%% ==== CHANNEL FILE ====
    % Usage: [sStudy, iStudy, iChannel] = bst_get('ChannelFile', ChannelFile)
    case 'ChannelFile'
        % No protocol in database
        if isempty(GlobalData.DataBase.iProtocol) || (GlobalData.DataBase.iProtocol == 0)
            return;
        end
        ProtocolInfo = GlobalData.DataBase.ProtocolInfo(GlobalData.DataBase.iProtocol);
        % Parse inputs
        if (nargin == 2)
            ChannelFile = varargin{2};
            ChannelFile = strrep(ChannelFile, ProtocolInfo.STUDIES, '');
        else
            error('Invalid call to bst_get().');
        end
        % Look for Channel file in all the surfaces of all subjects
        [argout1, argout2, argout3] = findFileInStudies('Channel', 'FileName', ChannelFile);
        
        
%% ==== CHANNEL FILE FOR STUDY ====
    % Usage: [ChannelFile, sStudy, iStudy] = bst_get('ChannelFileForStudy', StudyFile/DataFile)
    case 'ChannelFileForStudy'
        % Parse inputs
        if (nargin == 2)
            StudyFile = varargin{2};
        else
            error('Invalid call to bst_get().');
        end
        % Get study in database
        [sStudy, iStudy] = bst_get('Study', StudyFile);
        % If data file instead on Study file
        if isempty(sStudy)
            [sStudy, iStudy] = bst_get('AnyFile', StudyFile);
        end
        sChannel = bst_get('ChannelForStudy', iStudy);
        if ~isempty(sChannel)
            argout1 = sChannel.FileName;
            argout2 = sStudy;
            argout3 = iStudy;
        else
            argout1 = [];
        end
        
        
%% ==== CHANNEL STRUCT FOR STUDY ====
    % Usage: [sChannel, iChanStudy] = bst_get('ChannelForStudy', iStudies)
    case 'ChannelForStudy'
        % Parse inputs
        if (nargin == 2)
            iStudies = varargin{2};
        else
            error('Invalid call to bst_get().');
        end
        iChanStudies = [];
        sListChannel = [];
        for i = 1:length(iStudies)           
            % Get study 
            iStudy = iStudies(i);
            sStudy = bst_get('Study', iStudy);
            if isempty(sStudy)
                continue;
            end
            iChanStudy = iStudy;
            % === Analysis-Inter node ===
            iAnalysisInter      = -2;
            iGlobalDefaultStudy = -3;
            if (iStudy == iAnalysisInter)
                % If no channel file is defined in 'Analysis-intra' node: look in 
                if isempty(sStudy.Channel)
                    % Get global default study
                    sStudy = bst_get('Study', iGlobalDefaultStudy);
                    iChanStudy = iGlobalDefaultStudy;
                end
            % === All other nodes ===
            else
                % Get subject attached to study
                [sSubject, iSubject] = bst_get('Subject', sStudy.BrainStormSubject, 1);
                if isempty(sSubject)
                    return;
                end
                % Subject uses default channel/headmodel
                if (sSubject.UseDefaultChannel ~= 0)
                    [sStudy, iChanStudy] = bst_get('DefaultStudy', iSubject);
                    if isempty(sStudy)
                        return
                    end
                end
            end
            iChanStudies = [iChanStudies, iChanStudy];
            sListChannel = [sListChannel, sStudy.Channel];
        end
        % Return Channel structure
        argout1 = sListChannel;
        argout2 = iChanStudies;

%% ==== CHANNEL MODALITIES =====
    % Usage: [Modalities, DispMod, DefMod] = bst_get('ChannelModalities', ChannelFile)
    %        [Modalities, DispMod, DefMod] = bst_get('ChannelModalities', DataFile/ResultsFile/TimefreqFile...)
    case 'ChannelModalities'
        % Get input file
        [sStudy, iStudy, iItem, DataType, sItem] = bst_get('AnyFile', varargin{2});
        % If channel file
        if strcmpi(DataType, 'channel')
            sChannel = sItem;
        else
            sChannel = bst_get('ChannelForStudy', iStudy);
        end
        % Return modalities
        if ~isempty(sChannel)
            % Get all modalities
            if ~isempty(sChannel.DisplayableSensorTypes)
                % Return the displayable sensors on top of the list
                argout1 = cat(2, sChannel.DisplayableSensorTypes, setdiff(sChannel.Modalities, sChannel.DisplayableSensorTypes));
            else
                argout1 = sChannel.Modalities;
            end
            % Get only sensors that have spatial representation
            argout2 = sChannel.DisplayableSensorTypes;
            % Default candidates
            if ~isempty(argout2)
                defList = argout2;
            else
                defList = argout1;
            end
            if isempty(defList)
                return;
            end
            % Remove EDF and BDF from the default list
            defList = setdiff(defList, {'EDF','BDF','KDF'});
            % Get default modality
            if ismember('SEEG', defList)
                argout3 = 'SEEG';
            elseif ismember('ECOG', defList)
                argout3 = 'ECOG';
            elseif any(ismember({'MEG','MEG GRAD','MEG MAG'}, defList))
                argout3 = 'MEG';
            elseif ismember('EEG', defList)
                argout3 = 'EEG';
            elseif ismember('NIRS', defList)
                argout3 = 'NIRS';
            else
                argout3 = defList{1};
            end
            % Place the default on top of the lists
            if ismember(argout3, argout1)
                argout1(strcmpi(argout1, argout3)) = [];
                argout1 = cat(2, argout3, argout1);
            end
            if ismember(argout3, argout2)
                argout2(strcmpi(argout2, argout3)) = [];
                argout2 = cat(2, argout3, argout2);
            end
        else
            argout1 = [];
        end
        

%% ==== TIMEFREQ DISPLAY MODALITIES ====
    % Usage: DisplayMod = bst_get('TimefreqDisplayModalities', TimefreqFile)
    case 'TimefreqDisplayModalities'
        TimefreqFile = varargin{2};
        % Load sensor names from file
        TimefreqMat = in_bst_timefreq(TimefreqFile, 0, 'RowNames');
        % Get channel file
        ChannelFile = bst_get('ChannelFileForStudy', TimefreqFile);
        if isempty(ChannelFile)
            return;
        end
        % Load channel file
        ChannelMat = load(file_fullpath(ChannelFile), 'Channel');
        % Get channels that are present in the file
        [tmp__,I,J] = intersect({ChannelMat.Channel.Name}, TimefreqMat.RowNames);
        FileMod = unique({ChannelMat.Channel(I).Type});
        % Check if only one type of gradiometer is selected
        if isequal(FileMod, {'MEG GRAD'}) && all(cellfun(@(c)(c(end)=='2'), {ChannelMat.Channel(I).Name}))
            argout1 = {'MEG GRAD2'};
        elseif isequal(FileMod, {'MEG GRAD'}) && all(cellfun(@(c)(c(end)=='3'), {ChannelMat.Channel(I).Name}))
            argout1 = {'MEG GRAD3'};
        % Keep only the modalities that can be displayed (as topography)
        else
            argout1 = intersect(FileMod, {'MEG','MEG GRAD','MEG MAG','EEG','ECOG','SEEG','NIRS'});
        end
        
        
%% ==== CHANNEL DEVICE ====
    % Usage: Device = bst_get('ChannelDevice', ChannelFile)
    case 'ChannelDevice'
        ChannelFile = varargin{2};
        if ~isempty(strfind(ChannelFile, 'vectorview306'))
            Device = 'Vectorview306';
        elseif ~isempty(strfind(ChannelFile, 'ctf_acc1'))
            Device = 'CTF';
        elseif ~isempty(strfind(ChannelFile, '4d_acc1'))
            Device = '4D';
        elseif ~isempty(strfind(ChannelFile, 'babysquid'))
            Device = 'BabySQUID';
        elseif ~isempty(strfind(ChannelFile, 'babymeg'))
            Device = 'BabyMEG';
        elseif ~isempty(strfind(ChannelFile, 'kit'))
            Device = 'KIT';
        elseif ~isempty(strfind(ChannelFile, 'ricoh'))
            Device = 'RICOH';
        elseif ~isempty(strfind(ChannelFile, 'kriss'))
            Device = 'KRISS';
        elseif ~isempty(strfind(ChannelFile, 'nirsbrs'))
            Device = 'NIRS-BRS';
        else
            Device = '';
        end
        argout1 = Device;
        

%% ==== HEADMODEL STRUCT FOR STUDY ====
    % Usage: [sHeadModel] = bst_get('HeadModelForStudy', iStudy)
    case 'HeadModelForStudy'
        % Parse inputs
        if (nargin == 2)
            iStudy = varargin{2};
        else
            error('Invalid call to bst_get().');
        end
        % Get study 
        sStudy = bst_get('Study', iStudy);
        % === Analysis-Inter node ===
        iAnalysisInter      = -2;
        iGlobalDefaultStudy = -3;
        if (iStudy == iAnalysisInter)
            % If no channel file is defined in 'Analysis-intra' node: look in 
            if isempty(sStudy.iHeadModel)
                % Get global default study
                sStudy = bst_get('Study', iGlobalDefaultStudy);
            end
        % === All other nodes ===
        else
            % Get subject attached to study
            [sSubject, iSubject] = bst_get('Subject', sStudy.BrainStormSubject, 1);
            if isempty(sSubject)
                return;
            end
            % Subject uses default channel/headmodel
            if (sSubject.UseDefaultChannel ~= 0)
                sStudy = bst_get('DefaultStudy', iSubject);
                if isempty(sStudy)
                    return
                end
            end
        end
        % Return HeadModel structure
        if ~isempty(sStudy.iHeadModel)
            argout1 = sStudy.HeadModel(sStudy.iHeadModel(1));
        else
            argout1 = [];
        end

        
%% ==== HEADMODEL FILE ====
    % Usage: [sStudy, iStudy, iHeadModel] = bst_get('HeadModelFile', HeadModelFile, iStudies)
    %        [sStudy, iStudy, iHeadModel] = bst_get('HeadModelFile', HeadModelFile)
    case 'HeadModelFile'
        % No protocol in database
        if isempty(GlobalData.DataBase.iProtocol) || (GlobalData.DataBase.iProtocol == 0)
            return;
        end
        % Input #2: HeadModelFile
        ProtocolInfo = GlobalData.DataBase.ProtocolInfo(GlobalData.DataBase.iProtocol);
        HeadModelFile = varargin{2};
        HeadModelFile = strrep(HeadModelFile, ProtocolInfo.STUDIES, '');
        % Input #3: iStudies
        if (nargin < 3)
            iStudies = [];
        else
            iStudies = varargin{3};
        end
        % Look for surface file in all the surfaces of all subjects
        [argout1, argout2, argout3] = findFileInStudies('HeadModel', 'FileName', HeadModelFile, iStudies);
  
%% ==== NOISECOV FILE ====
    % Usage: [sStudy, iStudy, iNoiseCov] = bst_get('NoiseCovFile', NoiseCovFile, iStudies)
    %        [sStudy, iStudy, iNoiseCov] = bst_get('NoiseCovFile', NoiseCovFile)
    % Usage: [sStudy, iStudy, iNoiseCov] = bst_get('DataCovFile', NoiseCovFile, iStudies)
    %        [sStudy, iStudy, iNoiseCov] = bst_get('DataCovFile', NoiseCovFile)
    case {'NoiseCovFile', 'DataCovFile'}
        % No protocol in database
        if isempty(GlobalData.DataBase.iProtocol) || (GlobalData.DataBase.iProtocol == 0)
            return;
        end
        % Input #2: NoiseCovFile
        ProtocolInfo = GlobalData.DataBase.ProtocolInfo(GlobalData.DataBase.iProtocol);
        NoiseCovFile = varargin{2};
        NoiseCovFile = strrep(NoiseCovFile, ProtocolInfo.STUDIES, '');
        % Input #3: iStudies
        if (nargin < 3)
            iStudies = [];
        else
            iStudies = varargin{3};
        end
        % Look for surface file in all the surfaces of all subjects
        [argout1, argout2, argout3] = findFileInStudies('NoiseCov', 'FileName', NoiseCovFile, iStudies);

        
%% ==== DATA FILE ====
    % Usage: [sStudy, iStudy, iData] = bst_get('DataFile', DataFile, iStudies)
    %        [sStudy, iStudy, iData] = bst_get('DataFile', DataFile)
    case 'DataFile'
        % No protocol in database
        if isempty(GlobalData.DataBase.iProtocol) || (GlobalData.DataBase.iProtocol == 0)
            return;
        end
        % Input #2: DataFile
        ProtocolInfo = GlobalData.DataBase.ProtocolInfo(GlobalData.DataBase.iProtocol);
        DataFile = varargin{2};
        DataFile = strrep(DataFile, ProtocolInfo.STUDIES, '');
        % Input #3: iStudies
        if (nargin < 3)
            iStudies = [];
        else
            iStudies = varargin{3};
        end
        % Look for file in all the studies
        [argout1, argout2, argout3] = findFileInStudies('Data', 'FileName', DataFile, iStudies);
        

%% ==== DATA FOR DATA LIST ====
    % Usage: [iFoundData] = bst_get('DataForDataList', iStudy, DataListName)
    case 'DataForDataList'
        iStudy = varargin{2};
        DataListName = varargin{3};
        % Get study structure
        sStudy = bst_get('Study', iStudy);
        % Get all the data files held by this datalist
        listComments = cellfun(@(c)deblank(str_remove_parenth(c)), {sStudy.Data.Comment}, 'UniformOutput', 0);
        iFoundData = find(strcmp(listComments, DataListName));
        % Return found data files
        argout1 = iFoundData;

%% ==== MATRIX FOR MATRIX LIST ====
    % Usage: [iFoundMatrix] = bst_get('MatrixForMatrixList', iStudy, MatrixListName)
    case 'MatrixForMatrixList'
        iStudy = varargin{2};
        MatrixListName = varargin{3};
        % Get study structure
        sStudy = bst_get('Study', iStudy);
        % Get all the matrix files held by this datalist
        listComments = cellfun(@(c)deblank(str_remove_parenth(c)), {sStudy.Matrix.Comment}, 'UniformOutput', 0);
        iFoundMatrix = find(strcmp(listComments, MatrixListName));
        % Return found matrix files
        argout1 = iFoundMatrix;
        
        
%% ==== DATA FOR STUDY (INCLUDING SHARED STUDIES) ====
    % Usage: [iStudies, iDatas] = bst_get('DataForStudy', iStudy)
    case 'DataForStudy'
        % Get target study
        iStudy = varargin{2};
        sStudy = bst_get('Study', iStudy);
        isDefaultStudy  = strcmpi(sStudy.Name, bst_get('DirDefaultStudy'));
        isGlobalDefault = (iStudy == -3);
        
        % If study is the global default study
        sStudies = [];
        iStudies = [];
        if isGlobalDefault
            % Get all the subjects of the protocol
            nbSubjects = bst_get('SubjectCount');
            for iSubject = 1:nbSubjects
                sSubject = bst_get('Subject', iSubject, 1);
                if sSubject.UseDefaultChannel
                    [tmp_sStudies, tmp_iStudies] = bst_get('StudyWithSubject', sSubject.FileName);
                    sStudies = [sStudies, tmp_sStudies];
                    iStudies = [iStudies, tmp_iStudies];
                end
            end
        % Else, if study is a subject's default study (ie. channel file is shared by all studies of one subject)
        elseif isDefaultStudy
            % Get all the subject's studies
            [sStudies, iStudies] = bst_get('StudyWithSubject', sStudy.BrainStormSubject, 'intra_subject', 'default_study');
        else
            % Normal: one channel per condition
            sStudies = sStudy;
            iStudies = iStudy;
        end
        % Get all the DataFiles for all these studies
        for i = 1:length(sStudies)
            nData = length(sStudies(i).Data);
            argout1 = [argout1, repmat(iStudies(i), [1,nData])];
            argout2   = [argout2, 1:nData];
        end
       
        
%% ==== DATA FOR STUDIES (INCLUDING SHARED STUDIES) ====
    % Usage: [iStudies, iDatas] = bst_get('DataForStudies', iStudies)
    case 'DataForStudies'
        iStudies = varargin{2};
        for i = 1:length(iStudies)
            [tmp_iStudies, tmp_iDatas] = bst_get('DataForStudy', iStudies(i));
            argout1 = [argout1, tmp_iStudies];
            argout2 = [argout2, tmp_iDatas];
        end
        
%% ==== DATA FILE FOR CHANNEL FILE ====
    % Usage: DataFiles = bst_get('DataForChannelFile', ChannelFile)
    case 'DataForChannelFile'
        ChannelFile = varargin{2};
        DataFiles = {};
        % Get study for the given channel file
        [sStudy, iStudy] = bst_get('ChannelFile', ChannelFile);
        if isempty(sStudy)
            return;
        end
        % Get dependent data files
        [iStudies, iDatas] = bst_get('DataForStudy', iStudy);
        % Get all the Data filenames
        for i = 1:length(iStudies)
            sStudy = bst_get('Study', iStudies(i));
            DataFiles = cat(2, DataFiles, {sStudy.Data(iDatas(i)).FileName});
        end
        argout1 = DataFiles;
                
        
%% ==== RESULTS FILE ====
    % Usage: [sStudy, iStudy, iResult] = bst_get('ResultsFile', ResultsFile, iStudies)
    %        [sStudy, iStudy, iResult] = bst_get('ResultsFile', ResultsFile)
    case 'ResultsFile'
        % No protocol in database
        if isempty(GlobalData.DataBase.iProtocol) || (GlobalData.DataBase.iProtocol == 0)
            return;
        end
        % Input #2: ResultsFile
        ProtocolInfo = GlobalData.DataBase.ProtocolInfo(GlobalData.DataBase.iProtocol);
        ResultsFile = varargin{2};
        ResultsFile = strrep(ResultsFile, ProtocolInfo.STUDIES, '');
        % Input #3: iStudies
        if (nargin < 3)
            iStudies = [];
        else
            iStudies = varargin{3};
        end
        % Look for surface file in all the surfaces of all subjects
        [argout1, argout2, argout3] = findFileInStudies('Result', 'FileName', ResultsFile, iStudies);
       
        
%% ==== RESULTS FOR DATA FILE ====
    % Usage: [sStudy, iStudy, iResults] = bst_get('ResultsForDataFile', DataFile)           : search the whole protocol
    % Usage: [sStudy, iStudy, iResults] = bst_get('ResultsForDataFile', DataFile, iStudies) : search only the specified studies
    case 'ResultsForDataFile'
        % No protocol in database
        if isempty(GlobalData.DataBase.iProtocol) || (GlobalData.DataBase.iProtocol == 0)
            return;
        end
        ProtocolInfo = GlobalData.DataBase.ProtocolInfo(GlobalData.DataBase.iProtocol);
        % Input #2: DataFile
        DataFile = varargin{2};
        DataFile = strrep(DataFile, ProtocolInfo.STUDIES, '');
        % Determine in which studies to search for ResultsFile
        if (nargin >= 3)
            % Studies specified in argument
            iStudy = varargin{3};
        else
            % Get study in which DataFile is located
            [sStudy, iStudy] = bst_get('DataFile', DataFile);
            if isempty(iStudy)
                return;
            end
        end
        % Search selected studies
        [argout1, argout2, argout3] = findFileInStudies('Result', 'DataFile', DataFile, iStudy);


%% ==== STAT FILE ====
    % Usage: [sStudy, iStudy, iData] = bst_get('StatFile', StatFile, iStudies)
    %        [sStudy, iStudy, iData] = bst_get('StatFile', StatFile)
    case 'StatFile'
        % No protocol in database
        if isempty(GlobalData.DataBase.iProtocol) || (GlobalData.DataBase.iProtocol == 0)
            return;
        end
        % Input #2: SataFile
        ProtocolInfo = GlobalData.DataBase.ProtocolInfo(GlobalData.DataBase.iProtocol);
        StatFile = varargin{2};
        StatFile = strrep(StatFile, ProtocolInfo.STUDIES, '');
        % Input #3: iStudies
        if (nargin < 3)
            iStudies = [];
        else
            iStudies = varargin{3};
        end
        % Look for surface file in all the surfaces of all subjects
        [argout1, argout2, argout3] = findFileInStudies('Stat', 'FileName', StatFile, iStudies);
        
        
%% ==== STAT FOR DATA FILE ====
    % Usage: [sStudy, iStudy, iResults] = bst_get('StatForDataFile', DataFile)           : search the whole protocol
    % Usage: [sStudy, iStudy, iResults] = bst_get('StatForDataFile', DataFile, iStudies) : search only the specified studies
    case 'StatForDataFile'
        % No protocol in database
        if isempty(GlobalData.DataBase.iProtocol) || (GlobalData.DataBase.iProtocol == 0)
            return;
        end
        ProtocolInfo = GlobalData.DataBase.ProtocolInfo(GlobalData.DataBase.iProtocol);
        % Parse inputs
        if (nargin >= 2)
            DataFile = varargin{2};
            DataFile = strrep(DataFile, ProtocolInfo.STUDIES, '');
        else
            error('Invalid call to bst_get().');
        end
        % Determine in which studies to search for ResultsFile
        if (nargin >= 3)
            % Studies specified in argument
            iStudies = varargin{3};
        else
            % Get study in which DataFile is located
            [sStudies, iStudies] = bst_get('DataFile', DataFile);
            if isempty(iStudies)
                return;
            end
        end
        % Search selected studies
        [argout1, argout2, argout3] = findFileInStudies('Stat', 'DataFile', DataFile, iStudies);
     
%% ==== TIMEFREQ FILE ====
    % Usage: [sStudy, iStudy, iTimefreq] = bst_get('TimefreqFile', TimefreqFile, iStudies)
    %        [sStudy, iStudy, iTimefreq] = bst_get('TimefreqFile', TimefreqFile)
    case 'TimefreqFile'
        % No protocol in database
        if isempty(GlobalData.DataBase.iProtocol) || (GlobalData.DataBase.iProtocol == 0)
            return;
        end
        % Input #2: TimefreqFile
        ProtocolInfo = GlobalData.DataBase.ProtocolInfo(GlobalData.DataBase.iProtocol);
        TimefreqFile = varargin{2};
        TimefreqFile = strrep(TimefreqFile, ProtocolInfo.STUDIES, '');
        % Remove optional RefRowName
        iPipe = find(TimefreqFile == '|', 1);
        if ~isempty(iPipe)
            TimefreqFile = TimefreqFile(1:iPipe-1);
        end
        % Input #3: iStudies
        if (nargin < 3)
            iStudies = [];
        else
            iStudies = varargin{3};
        end
        % Look for surface file in all the surfaces of all subjects
        [argout1, argout2, argout3] = findFileInStudies('Timefreq', 'FileName', TimefreqFile, iStudies);
        
%% ==== TIMEFREQ FOR FILE ====
    % Usage: [sStudy, iStudy, iTimefreq] = bst_get('TimefreqForFile', FileName, iStudies) : search only the specified studies
    %        [sStudy, iStudy, iTimefreq] = bst_get('TimefreqForFile', FileName)           : search the whole protocol
    case 'TimefreqForFile'
        % No protocol in database
        if isempty(GlobalData.DataBase.iProtocol) || (GlobalData.DataBase.iProtocol == 0)
            return;
        end
        ProtocolInfo = GlobalData.DataBase.ProtocolInfo(GlobalData.DataBase.iProtocol);
        % Parse inputs
        if (nargin >= 2)
            FileName = varargin{2};
            FileName = strrep(FileName, ProtocolInfo.STUDIES, '');
        else
            error('Invalid call to bst_get().');
        end
        % Get study in which file is located
        if (nargin >= 3)
            iStudies = varargin{3};
            [sStudy, iStudy, iFile, DataType] = bst_get('AnyFile', FileName, iStudies);
        else
            [sStudy, iStudy, iFile, DataType] = bst_get('AnyFile', FileName);
        end
        % If file was not found
        if isempty(iStudy)
            return;
        end
        % Search direct dependent files
        [tmp, tmp, iTf] = findFileInStudies('Timefreq', 'DataFile', FileName, iStudy);
        % Data files: get all the depending results, and then all the timefreq for those results
        if strcmpi(DataType, 'data')
            [tmp, tmp, iResults] = bst_get('ResultsForDataFile', FileName, iStudy);
            for i = 1:length(iResults);
                % Search selected studies
                [tmp, tmp, iTfRes] = findFileInStudies('Timefreq', 'DataFile', sStudy.Result(iResults(i)).FileName, iStudy);
                if ~isempty(iTfRes)
                    iTf = [iTf iTfRes];
                end
            end
        end
        % Return results
        if ~isempty(iTf)
            argout1 = sStudy;
            argout2 = iStudy;
            argout3 = iTf;
        end
        
        
%% ==== DIPOLES FOR FILE ====
    % Usage: [sStudy, iStudy, iDipoles] = bst_get('DipolesForFile', FileName, iStudies) : search only the specified studies
    %        [sStudy, iStudy, iDipoles] = bst_get('DipolesForFile', FileName)           : search the whole protocol
    case 'DipolesForFile'
        % No protocol in database
        if isempty(GlobalData.DataBase.iProtocol) || (GlobalData.DataBase.iProtocol == 0)
            return;
        end
        ProtocolInfo = GlobalData.DataBase.ProtocolInfo(GlobalData.DataBase.iProtocol);
        % Parse inputs
        if (nargin >= 2)
            FileName = varargin{2};
            FileName = strrep(FileName, ProtocolInfo.STUDIES, '');
        else
            error('Invalid call to bst_get().');
        end
        % Get study in which file is located
        if (nargin >= 3)
            iStudies = varargin{3};
            [sStudy, iStudy, iFile, DataType] = bst_get('AnyFile', FileName, iStudies);
        else
            [sStudy, iStudy, iFile, DataType] = bst_get('AnyFile', FileName);
        end
        % If file was not found
        if isempty(iStudy)
            return;
        end
        % Search direct dependent files
        [tmp, tmp, iDip] = findFileInStudies('Dipoles', 'DataFile', FileName, iStudy);
        % Data files: get all the depending results, and then all the timefreq for those results
        if strcmpi(DataType, 'data')
            [tmp, tmp, iResults] = bst_get('ResultsForDataFile', FileName, iStudy);
            for i = 1:length(iResults);
                % Search selected studies
                [tmp, tmp, iDipRes] = findFileInStudies('Dipoles', 'DataFile', sStudy.Result(iResults(i)).FileName, iStudy);
                if ~isempty(iDipRes)
                    iDip = [iDip, iDipRes];
                end
            end
        end
        % Return results
        if ~isempty(iDip)
            argout1 = sStudy;
            argout2 = iStudy;
            argout3 = iDip;
        end
        
        
%% ==== TIMEFREQ FOR KERNEL ====
    % Find all the timefreq files dependent from links due to a given kernel
    % Usage: [sStudy, iStudy, iTimefreq] = bst_get('TimefreqForKernel', KernelFile) 
    case 'TimefreqForKernel'
        sFoundStudy = [];
        iFoundStudy = [];
        iFoundTimefreq = [];
        % No protocol in database
        if isempty(GlobalData.DataBase.iProtocol) || (GlobalData.DataBase.iProtocol == 0)
            return;
        end
        % Get study in which file is located
        KernelFile = varargin{2};
        [sStudy, iStudy, iFile, DataType] = bst_get('ResultsFile', KernelFile);
        if isempty(iStudy)
            return;
        end
        % Get all the data files relative with this kernel
        [iDepStudies, iDepDatas] = bst_get('DataForStudy', iStudy);
        % Keep only once each study
        iDepStudies = unique(iDepStudies);
        % Process all the dependent studies
        for iSt = 1:length(iDepStudies)
            % Get the study structure
            sDepStudy = bst_get('Study', iDepStudies(iSt));
            % Process each timefreq file separately
            for iTf = 1:length(sDepStudy.Timefreq)
                DataFile = sDepStudy.Timefreq(iTf).DataFile;
                % Keep only the files that are linked to LINKS
                if isempty(DataFile) || (length(DataFile) < 5) || ~isequal(DataFile(1:5), 'link|')
                    continue;
                end
                % Split link
                splitFile = str_split(DataFile, '|');
                % If the kernel is found: add it to the found list
                if file_compare(splitFile{2}, KernelFile)
                    sFoundStudy = [sFoundStudy sDepStudy];
                    iFoundStudy = [iFoundStudy, iDepStudies(iSt)];
                    iFoundTimefreq = [iFoundTimefreq, iTf];
                end
            end
        end
        % Return findings
        argout1 = sFoundStudy;
        argout2 = iFoundStudy;
        argout3 = iFoundTimefreq;
        
        
%% ==== DIPOLES FOR KERNEL ====
    % Find all the dipoles files dependent from links due to a given kernel
    % Usage: [sStudy, iStudy, iDipoles] = bst_get('TimefreqForKernel', KernelFile) 
    case 'DipolesForKernel'
        sFoundStudy = [];
        iFoundStudy = [];
        iFoundDipoles = [];
        % No protocol in database
        if isempty(GlobalData.DataBase.iProtocol) || (GlobalData.DataBase.iProtocol == 0)
            return;
        end
        % Get study in which file is located
        KernelFile = varargin{2};
        [sStudy, iStudy, iFile, DataType] = bst_get('ResultsFile', KernelFile);
        if isempty(iStudy)
            return;
        end
        % Get all the data files relative with this kernel
        [iDepStudies, iDepDatas] = bst_get('DataForStudy', iStudy);
        % Keep only once each study
        iDepStudies = unique(iDepStudies);
        % Process all the dependent studies
        for iSt = 1:length(iDepStudies)
            % Get the study structure
            sDepStudy = bst_get('Study', iDepStudies(iSt));
            % Process each timefreq file separately
            for iDip = 1:length(sDepStudy.Dipoles)
                DataFile = sDepStudy.Dipoles(iDip).DataFile;
                % Keep only the files that are linked to LINKS
                if isempty(DataFile) || (length(DataFile) < 5) || ~isequal(DataFile(1:5), 'link|')
                    continue;
                end
                % Split link
                splitFile = str_split(DataFile, '|');
                % If the kernel is found: add it to the found list
                if file_compare(splitFile{2}, KernelFile)
                    sFoundStudy = [sFoundStudy sDepStudy];
                    iFoundStudy = [iFoundStudy, iDepStudies(iSt)];
                    iFoundDipoles = [iFoundDipoles, iDip];
                end
            end
        end
        % Return findings
        argout1 = sFoundStudy;
        argout2 = iFoundStudy;
        argout3 = iFoundDipoles;
        
        
%% ==== DIPOLES FILE ====
    % Usage: [sStudy, iStudy, iDipole] = bst_get('DipolesFile', DipolesFile, iStudies)
    %        [sStudy, iStudy, iDipole] = bst_get('DipolesFile', DipolesFile)
    case 'DipolesFile'
        % No protocol in database
        if isempty(GlobalData.DataBase.iProtocol) || (GlobalData.DataBase.iProtocol == 0)
            return;
        end
        % Input #2: DipolesFile
        ProtocolInfo = GlobalData.DataBase.ProtocolInfo(GlobalData.DataBase.iProtocol);
        DipolesFile = varargin{2};
        DipolesFile = strrep(DipolesFile, ProtocolInfo.STUDIES, '');
        % Input #3: iStudies
        if (nargin < 3)
            iStudies = [];
        else
            iStudies = varargin{3};
        end
        % Look for surface file in all the surfaces of all subjects
        [argout1, argout2, argout3] = findFileInStudies('Dipoles', 'FileName', DipolesFile, iStudies);
        
        
%% ==== MATRIX FILE ====
    % Usage: [sStudy, iStudy, iDipole] = bst_get('MatrixFile', MatrixFile, iStudies)
    %        [sStudy, iStudy, iDipole] = bst_get('MatrixFile', MatrixFile)
    case 'MatrixFile'
        % No protocol in database
        if isempty(GlobalData.DataBase.iProtocol) || (GlobalData.DataBase.iProtocol == 0)
            return;
        end
        % Input #2: MatrixFile
        ProtocolInfo = GlobalData.DataBase.ProtocolInfo(GlobalData.DataBase.iProtocol);
        MatrixFile = varargin{2};
        MatrixFile = strrep(MatrixFile, ProtocolInfo.STUDIES, '');
        % Input #3: iStudies
        if (nargin < 3)
            iStudies = [];
        else
            iStudies = varargin{3};
        end
        % Look for surface file in all the surfaces of all subjects
        [argout1, argout2, argout3] = findFileInStudies('Matrix', 'FileName', MatrixFile, iStudies);
        
%% ==== IMAGE FILE ====
    % Usage: [sStudy, iStudy, iDipole] = bst_get('ImageFile', ImageFile, iStudies)
    %        [sStudy, iStudy, iDipole] = bst_get('ImageFile', ImageFile)
    case 'ImageFile'
        % No protocol in database
        if isempty(GlobalData.DataBase.iProtocol) || (GlobalData.DataBase.iProtocol == 0)
            return;
        end
        % Input #2: ImageFile
        ProtocolInfo = GlobalData.DataBase.ProtocolInfo(GlobalData.DataBase.iProtocol);
        ImageFile = varargin{2};
        ImageFile = strrep(ImageFile, ProtocolInfo.STUDIES, '');
        % Input #3: iStudies
        if (nargin < 3)
            iStudies = [];
        else
            iStudies = varargin{3};
        end
        % Look for surface file in all the surfaces of all subjects
        [argout1, argout2, argout3] = findFileInStudies('Image', 'FileName', ImageFile, iStudies);
        
        
%% ==== ANY FILE ====
    % Usage: [sStudy, iStudy, iFile, DataType, sItem] = bst_get('AnyFile', FileName, iStudies)
    %        [sStudy, iStudy, iFile, DataType, sItem] = bst_get('AnyFile', FileName)
    case 'AnyFile'
        % Input #2: FileName
        FileName = varargin{2};
        if isempty(FileName)
            return
        end
        % Input #3: iStudies
        if (nargin < 3)
            iStudies = [];
        else
            iStudies = varargin{3};
        end
        % Get data format
        fileType = file_gettype(FileName);
        if isempty(fileType)
            error('File type is not recognized.');
        end
        sItem = [];
        % Get information related with this file
        switch (fileType)
            % ===== FUNCTIONAL =====
            case 'channel'
                [sStudy, iStudy] = bst_get('ChannelFile', FileName);
                iItem = 1;
                if (nargout >= 5) && ~isempty(sStudy)
                    sItem = sStudy.Channel;
                end
            case 'headmodel'
                [sStudy, iStudy, iItem] = bst_get('HeadModelFile', FileName);
                if (nargout >= 5) && ~isempty(sStudy)
                    sItem = sStudy.HeadModel(iItem);
                end
            case 'noisecov'
                [sStudy, iStudy, iItem] = bst_get('NoiseCovFile', FileName);
                if (nargout >= 5) && ~isempty(sStudy)
                    sItem = sStudy.NoiseCov(iItem);
                end
            case 'ndatacov'
                [sStudy, iStudy, iItem] = bst_get('DataCovFile', FileName);
                if (nargout >= 5) && ~isempty(sStudy)
                    sItem = sStudy.NoiseCov(iItem);
                end
            case 'data'
                [sStudy, iStudy, iItem] = bst_get('DataFile', FileName, iStudies);
                if (nargout >= 5) && ~isempty(sStudy)
                    sItem = sStudy.Data(iItem);
                end
            case {'results', 'link'}
                [sStudy, iStudy, iItem] = bst_get('ResultsFile', FileName, iStudies);
                if (nargout >= 5) && ~isempty(sStudy)
                    sItem = sStudy.Result(iItem);
                end
            case {'presults', 'pdata','ptimefreq','pmatrix'}
                [sStudy, iStudy, iItem] = bst_get('StatFile', FileName, iStudies);
                if (nargout >= 5) && ~isempty(sStudy)
                    sItem = sStudy.Stat(iItem);
                end
            case 'dipoles'
                [sStudy, iStudy, iItem] = bst_get('DipolesFile', FileName, iStudies);
                if (nargout >= 5) && ~isempty(sStudy)
                    sItem = sStudy.Dipoles(iItem);
                end
            case 'timefreq'
                % Remove optional RefRowName
                iPipe = find(FileName == '|', 1);
                if ~isempty(iPipe)
                    FileName = FileName(1:iPipe-1);
                end
                [sStudy, iStudy, iItem] = bst_get('TimefreqFile', FileName, iStudies);
                if (nargout >= 5) && ~isempty(sStudy)
                    sItem = sStudy.Timefreq(iItem);
                end
            case 'matrix'
                [sStudy, iStudy, iItem] = bst_get('MatrixFile', FileName, iStudies);
                if (nargout >= 5) && ~isempty(sStudy)
                    sItem = sStudy.Matrix(iItem);
                end
            case 'brainstormstudy'
                [sStudy, iStudy] = bst_get('Study', FileName);
                iItem = 0;
                if (nargout >= 5) && ~isempty(sStudy)
                    sItem = sStudy;
                end
            case {'image', 'video', 'videolink'}
                [sStudy, iStudy, iItem] = bst_get('ImageFile', FileName, iStudies);
                if (nargout >= 5) && ~isempty(sStudy)
                    sItem = sStudy.Image(iItem);
                end
            % ===== ANATOMY =====
            case {'cortex','scalp','innerskull','outerskull','tess','fibers','fem'}
                [sStudy, iStudy, iItem] = bst_get('SurfaceFile', FileName);
                if (nargout >= 5) && ~isempty(sStudy)
                    sItem = sStudy.Surface(iItem);
                end
            case 'subjectimage'
                [sStudy, iStudy, iItem] = bst_get('MriFile', FileName);
                if (nargout >= 5) && ~isempty(sStudy)
                    sItem = sStudy.Anatomy(iItem);
                end
            case 'brainstormsubject'
                [sStudy, iStudy] = bst_get('Subject', FileName);
                iItem = 0;
                if (nargout >= 5) && ~isempty(sStudy)
                    sItem = sStudy;
                end
            otherwise
                error('File type is not recognized.');
        end 
        argout1 = sStudy;
        argout2 = iStudy;
        argout3 = iItem;
        argout4 = fileType;
        if (nargout >= 5)
            argout5 = sItem;
        end
       
        
%% ==== GET RELATED DATA FILE ====
    % Usage: DataFile = bst_get('RelatedDataFile', FileName, iStudies)
    %        DataFile = bst_get('RelatedDataFile', FileName)
    case 'RelatedDataFile'
        % Input #2: FileName
        FileName = varargin{2};
        % Input #3: iStudies
        if (nargin < 3)
            iStudies = [];
        else
            iStudies = varargin{3};
        end
        % Get file in database
        [sStudy, iStudy, iFile, fileType] = bst_get('AnyFile', FileName, iStudies);
        % If this data file does not belong to any study
        if isempty(sStudy)
            return;
        end
        % Get associated data file
        switch (fileType)
            case 'data'
                RelatedDataFile = sStudy.Data(iFile).FileName;
            case {'pdata','presults','ptimefreq','pmatrix'}
                RelatedDataFile = sStudy.Stat(iFile).DataFile;
            case {'results', 'link'}
                RelatedDataFile = sStudy.Result(iFile).DataFile;
            case 'dipoles'
                RelatedDataFile = sStudy.Dipoles(iFile).DataFile;
            case 'timefreq'
                RelatedDataFile = sStudy.Timefreq(iFile).DataFile;
            otherwise
                RelatedDataFile = '';
        end
        % If related file is results: get related data file
        if ~isempty(RelatedDataFile)
            relFileType = file_gettype(RelatedDataFile);
            if ismember(relFileType, {'link','results'})
                RelatedDataFile = bst_get('RelatedDataFile', RelatedDataFile, iStudy);
            end
        end
        % Return file
        argout1 = RelatedDataFile;
        
%% ==== ALL CONDITIONS FOR ONE SUBJECT ====
    % Usage: [Conditions] =  bst_get('ConditionsForSubject', SubjectFile)
    case 'ConditionsForSubject'
        % Parse inputs
        if (nargin == 2)
            SubjectFile = varargin{2};
        else
            error('Invalid call to bst_get().');
        end
        % Get list of studies associated with subject
        sStudies = bst_get('StudyWithSubject', SubjectFile);
        % Get Conditions for each study
        Conditions = {};
        for i = 1:length(sStudies)
            % Test if the condition of this study was not added previously
            isNewCondition = 1;
            for iCond = 1:length(Conditions)
                % If new condition is found 
                % (and excludes DirAnalysisIntra and DirDefaultSubject from list)
                if isempty(sStudies(i).Condition) || ...
                   isequal(sStudies(i).Condition, Conditions(iCond)) || ...
                   strcmpi(sStudies(i).Condition{1}, bst_get('DirAnalysisIntra')) || ...
                   strcmpi(sStudies(i).Condition{1}, bst_get('DirDefaultSubject'))
                    isNewCondition = 0;
                    break;
                end
            end
            % If Condition is not added yet : add it to the list
            if isNewCondition && ~isempty(sStudies(i).Condition)
                Conditions{end+1} = sStudies(i).Condition{1};
            end
        end
        % Return conditions list
        argout1 = Conditions;
        
        
%% ==== ANATOMY DEFAULTS ====
    % Returns the list of all the anatomy defaults (distributed with the software + user defined)
    case 'AnatomyDefaults'
        % Parse inputs
        if (nargin == 2)
            AnatName = varargin{2};
        else
            AnatName = [];
        end
        % Get templates from the brainstorm3 folder
        progDir   = bst_fullfile(bst_get('BrainstormDefaultsDir'), 'anatomy');
        progFiles = dir(progDir);
        % Get templates from the user folder
        userDir   = bst_fullfile(bst_get('UserDefaultsDir'), 'anatomy');
        userFiles = dir(userDir);
        % Combine the two lists
        AllFiles = cat(2, cellfun(@(c)bst_fullfile(progDir,c), {progFiles.name}, 'UniformOutput', 0), ...
                          cellfun(@(c)bst_fullfile(userDir,c), setdiff({userFiles.name}, {progFiles.name}), 'UniformOutput', 0));
        % Initialize list of defaults
        sTemplates = repmat(struct('FilePath',[],'Name',[]), 0);
        % Find all the valid defaults (.zip files or subdirectory with a brainstormsubject.mat in it)
        for i = 1:length(AllFiles)
            % Decompose file name
            [fPath, fBase, fExt] = bst_fileparts(AllFiles{i});
            % Entry is a directory W/ a name that does not start with a '.' 
            if isempty(fBase) || strcmpi(fBase(1),'.') || (~isempty(fExt) && ~strcmpi(fExt, '.zip'))
                continue;
            end
            % If it's a folder: check for a brainstormsubject file
            if isdir(AllFiles{i})
                bstFiles = dir(bst_fullfile(AllFiles{i}, 'brainstormsubject*.mat'));
                if (length(bstFiles) == 1)
                    sTemplates(end+1).FilePath = AllFiles{i};
                    sTemplates(end).Name = fBase;
                end
            % If it's a zip file
            elseif isequal(fExt, '.zip')
                sTemplates(end+1).FilePath = AllFiles{i};
                sTemplates(end).Name = fBase;
            end
        end
        % Get defaults from internet 
        if ~ismember('icbm152', lower({sTemplates.Name}))
            sTemplates(end+1).FilePath = 'http://neuroimage.usc.edu/bst/getupdate.php?t=ICBM152_2019';
            sTemplates(end).Name = 'ICBM152';
        end
        if ~ismember('icbm152_2019', lower({sTemplates.Name}))
            sTemplates(end+1).FilePath = 'http://neuroimage.usc.edu/bst/getupdate.php?t=ICBM152_2019';
            sTemplates(end).Name = 'ICBM152_2019';
        end
        if ~ismember('icbm152_brainsuite_2016', lower({sTemplates.Name}))
            sTemplates(end+1).FilePath = 'http://neuroimage.usc.edu/bst/getupdate.php?t=ICBM152_BrainSuite_2016';
            sTemplates(end).Name = 'ICBM152_BrainSuite_2016';
        end
        if ~ismember('colin27_2016', lower({sTemplates.Name}))
            sTemplates(end+1).FilePath = 'http://neuroimage.usc.edu/bst/getupdate.php?t=Colin27_2016';
            sTemplates(end).Name = 'Colin27_2016';
        end
        if ~ismember('colin27_brainsuite_2016', lower({sTemplates.Name}))
            sTemplates(end+1).FilePath = 'http://neuroimage.usc.edu/bst/getupdate.php?t=Colin27_BrainSuite_2016';
            sTemplates(end).Name = 'Colin27_BrainSuite_2016';
        end
        if ~ismember('bci-dni_brainsuite_2020', lower({sTemplates.Name}))
            sTemplates(end+1).FilePath = 'http://neuroimage.usc.edu/bst/getupdate.php?t=BCI-DNI_BrainSuite_2020';
            sTemplates(end).Name = 'BCI-DNI_BrainSuite_2020';
        end
        if ~ismember('uscbrain_brainsuite_2020', lower({sTemplates.Name}))
            sTemplates(end+1).FilePath = 'http://neuroimage.usc.edu/bst/getupdate.php?t=USCBrain_BrainSuite_2020';
            sTemplates(end).Name = 'USCBrain_BrainSuite_2020';
        end
        if ~ismember('fsaverage_2020', lower({sTemplates.Name}))
            sTemplates(end+1).FilePath = 'http://neuroimage.usc.edu/bst/getupdate.php?t=FSAverage_2020';
            sTemplates(end).Name = 'FsAverage_2020';
        end
        if ~ismember('kabdebon_7w', lower({sTemplates.Name}))
            sTemplates(end+1).FilePath = 'http://neuroimage.usc.edu/bst/getupdate.php?t=Kabdebon_7w';
            sTemplates(end).Name = 'Kabdebon_7w';
        end
        if ~ismember('oreilly_0.5m_2021', lower({sTemplates.Name}))
            sTemplates(end+1).FilePath = 'http://neuroimage.usc.edu/bst/getupdate.php?t=Oreilly_0.5m_2021';
            sTemplates(end).Name = 'Oreilly_0.5m_2021';
        end
        if ~ismember('oreilly_1m_2021', lower({sTemplates.Name}))
            sTemplates(end+1).FilePath = 'http://neuroimage.usc.edu/bst/getupdate.php?t=Oreilly_1m_2021';
            sTemplates(end).Name = 'Oreilly_1m_2021';
        end
        if ~ismember('oreilly_2m_2021', lower({sTemplates.Name}))
            sTemplates(end+1).FilePath = 'http://neuroimage.usc.edu/bst/getupdate.php?t=Oreilly_2m_2021';
            sTemplates(end).Name = 'Oreilly_2m_2021';
        end
        if ~ismember(lower({sTemplates.Name}), 'oreilly_3m_2021')
            sTemplates(end+1).FilePath = 'http://neuroimage.usc.edu/bst/getupdate.php?t=Oreilly_3m_2021';
            sTemplates(end).Name = 'Oreilly_3m_2021';
        end
        if ~ismember('oreilly_4.5m_2021', lower({sTemplates.Name}))
            sTemplates(end+1).FilePath = 'http://neuroimage.usc.edu/bst/getupdate.php?t=Oreilly_4.5m_2021';
            sTemplates(end).Name = 'Oreilly_4.5m_2021';
        end
        if ~ismember('oreilly_6m_2021', lower({sTemplates.Name}))
            sTemplates(end+1).FilePath = 'http://neuroimage.usc.edu/bst/getupdate.php?t=Oreilly_6m_2021';
            sTemplates(end).Name = 'Oreilly_6m_2021';
        end
        if ~ismember('oreilly_7.5m_2021', lower({sTemplates.Name}))
            sTemplates(end+1).FilePath = 'http://neuroimage.usc.edu/bst/getupdate.php?t=Oreilly_7.5m_2021';
            sTemplates(end).Name = 'Oreilly_7.5m_2021';
        end
        if ~ismember('oreilly_9m_2021', lower({sTemplates.Name}))
            sTemplates(end+1).FilePath = 'http://neuroimage.usc.edu/bst/getupdate.php?t=Oreilly_9m_2021';
            sTemplates(end).Name = 'Oreilly_9m_2021';
        end
        if ~ismember('oreilly_10.5m_2021', lower({sTemplates.Name}))
            sTemplates(end+1).FilePath = 'http://neuroimage.usc.edu/bst/getupdate.php?t=Oreilly_10.5m_2021';
            sTemplates(end).Name = 'Oreilly_10.5m_2021';
        end
        if ~ismember('oreilly_12m_2021', lower({sTemplates.Name}))
            sTemplates(end+1).FilePath = 'http://neuroimage.usc.edu/bst/getupdate.php?t=Oreilly_12m_2021';
            sTemplates(end).Name = 'Oreilly_12m_2021';
        end
        if ~ismember('oreilly_15m_2021', lower({sTemplates.Name}))
            sTemplates(end+1).FilePath = 'http://neuroimage.usc.edu/bst/getupdate.php?t=Oreilly_15m_2021';
            sTemplates(end).Name = 'Oreilly_15m_2021';
        end
        if ~ismember('oreilly_18m_2021', lower({sTemplates.Name}))
            sTemplates(end+1).FilePath = 'http://neuroimage.usc.edu/bst/getupdate.php?t=Oreilly_18m_2021';
            sTemplates(end).Name = 'Oreilly_18m_2021';
        end
        if ~ismember('oreilly_24m_2021', lower({sTemplates.Name}))
            sTemplates(end+1).FilePath = 'http://neuroimage.usc.edu/bst/getupdate.php?t=Oreilly_24m_2021';
            sTemplates(end).Name = 'Oreilly_24m_2021';
        end
        % If a specific template was requested
        if ~isempty(AnatName)
            iAnat = find(strcmpi({sTemplates.Name}, AnatName));
            sTemplates = sTemplates(iAnat);
        end
        % Sort in alphabetical order
        if ~isempty(sTemplates)
            [tmp__, I] = sort_nat({sTemplates(2:end).Name});
            sTemplates = sTemplates([1, I+1]);
        end
        % Return defaults list
        argout1 = sTemplates;
        
        
%% ==== MNI ATLASES ====
    % Returns the list of all the available MNI atlases
    case 'MniAtlasDefaults'
        % Get templates from the brainstorm3 folder
        mniDir   = bst_fullfile(bst_get('UserDefaultsDir'), 'mniatlas');
        mniFiles = dir(bst_fullfile(mniDir, '*.nii.gz'));
        mniFiles = cellfun(@(c)bst_fullfile(mniDir,c), {mniFiles.name}, 'UniformOutput', 0);
        % Initialize list of defaults
        sTemplates = repmat(struct('FilePath',[],'Name',[],'Info',[]), 0);
        % Find all the valid defaults (.zip files or subdirectory with a brainstormsubject.mat in it)
        for i = 1:length(mniFiles)
            % Decompose file name
            [fPath, fBase, fExt] = bst_fileparts(mniFiles{i});
            % Keep only files with .nii and .nii.gz extensions
            if ~isempty(fBase) && (fBase(1) ~= '.') && ~isempty(fExt) && strcmpi(fExt, '.gz')
                sTemplates(end+1).FilePath = mniFiles{i};
                sTemplates(end).Name = strrep(fBase, '.nii', '');
                sTemplates(end).Info = '';
            end
        end
        % Sort in alphabetical order
        if ~isempty(sTemplates)
            [tmp__, I] = sort_nat(lower({sTemplates.Name}));
            sTemplates = sTemplates(I);
        end
        
        % Get defaults from internet
        if ~ismember('aal2', lower({sTemplates.Name}))
            sTemplates(end+1).FilePath = 'http://neuroimage.usc.edu/bst/getupdate.php?t=mni_AAL2';
            sTemplates(end).Name = 'AAL2';
            sTemplates(end).Info = 'https://www.gin.cnrs.fr/en/tools/aal/';
        end
        if ~ismember('aal3', lower({sTemplates.Name}))
            sTemplates(end+1).FilePath = 'http://neuroimage.usc.edu/bst/getupdate.php?t=mni_AAL3';
            sTemplates(end).Name = 'AAL3';
            sTemplates(end).Info = 'https://www.gin.cnrs.fr/en/tools/aal/';
        end
        if ~ismember('aicha', lower({sTemplates.Name}))
            sTemplates(end+1).FilePath = 'http://neuroimage.usc.edu/bst/getupdate.php?t=mni_AICHA';
            sTemplates(end).Name = 'AICHA';
            sTemplates(end).Info = 'https://www.gin.cnrs.fr/en/tools/aicha';
        end
        if ~ismember('brainnetome', lower({sTemplates.Name}))
            sTemplates(end+1).FilePath = 'http://neuroimage.usc.edu/bst/getupdate.php?t=mni_Brainnetome';
            sTemplates(end).Name = 'Brainnetome';
            sTemplates(end).Info = 'http://atlas.brainnetome.org/';
        end
        if ~ismember('brainnetome_leaddbs', lower({sTemplates.Name}))
            sTemplates(end+1).FilePath = 'http://neuroimage.usc.edu/bst/getupdate.php?t=mni_Brainnetome_leaddbs';
            sTemplates(end).Name = 'Brainnetome_leaddbs';
            sTemplates(end).Info = 'http://atlas.brainnetome.org/';
        end
        if ~ismember('brodmann', lower({sTemplates.Name}))
            sTemplates(end+1).FilePath = 'http://neuroimage.usc.edu/bst/getupdate.php?t=mni_Brodmann';
            sTemplates(end).Name = 'Brodmann';
            sTemplates(end).Info = 'https://people.cas.sc.edu/rorden/mricro/lesion.html#brod';
        end
        if ~ismember('hammers83', lower({sTemplates.Name}))
            sTemplates(end+1).FilePath = 'http://neuroimage.usc.edu/bst/getupdate.php?t=mni_Hammers';
            sTemplates(end).Name = 'Hammers';
            sTemplates(end).Info = 'http://brain-development.org/brain-atlases/adult-brain-atlases/';
        end
        if ~ismember('neuromorphometrics', lower({sTemplates.Name}))
            sTemplates(end+1).FilePath = 'http://neuroimage.usc.edu/bst/getupdate.php?t=mni_Neuromorphometrics';
            sTemplates(end).Name = 'Neuromorphometrics';
            sTemplates(end).Info = 'https://search.kg.ebrains.eu/instances/Dataset/ef48c5e9-6b3c-4d5a-a9a9-e678fe10bdf6';
        end
        if ~ismember('julich-brain-v25', lower({sTemplates.Name}))
            sTemplates(end+1).FilePath = 'http://neuroimage.usc.edu/bst/getupdate.php?t=mni_Julich-Brain-v25';
            sTemplates(end).Name = 'Julich-Brain-v25';
            sTemplates(end).Info = 'https://search.kg.ebrains.eu/instances/Dataset/ef48c5e9-6b3c-4d5a-a9a9-e678fe10bdf6';
        end
        if ~ismember('schaefer2018_100_7net', lower({sTemplates.Name}))
            sTemplates(end+1).FilePath = 'http://neuroimage.usc.edu/bst/getupdate.php?t=mni_Schaefer2018';
            sTemplates(end).Name = 'Schaefer2018';
            sTemplates(end).Info = 'https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal';
        end
        % Return defaults list
        argout1 = sTemplates;
        
        
%% ==== EEG DEFAULTS ====
    % Returns an array of struct(fullpath, name) of all the Brainstorm eeg nets defaults
    % Usage: EegDefaults = bst_get('EegDefaults')
    %        EegDefaults = bst_get('EegDefaults', TemplateName=[], SetupName=[])
    case 'EegDefaults'
        % Parse inputs
        if (nargin >= 3)
            SetupName = varargin{3};
        else
            SetupName = [];
        end
        if (nargin >= 2)
            TemplateName = varargin{2};
        else
            TemplateName = [];
        end
        % Get templates from the brainstorm3 folder
        progDir   = bst_fullfile(bst_get('BrainstormDefaultsDir'), 'eeg');
        progFiles = dir(bst_fullfile(progDir, '*'));
        % Get templates from the user folder
        userDir   = bst_fullfile(bst_get('UserDefaultsDir'), 'eeg');
        userFiles = dir(bst_fullfile(userDir, '*'));
        % Combine the two lists
        dirList = cat(2, cellfun(@(c)bst_fullfile(progDir,c), {progFiles.name}, 'UniformOutput', 0), ...
                          cellfun(@(c)bst_fullfile(userDir,c), setdiff({userFiles.name}, {progFiles.name}), 'UniformOutput', 0));
        % Initialize list of folders
        fullDefaultsList = repmat(struct('contents','', 'name',''), 0);
        % For each template directory
        for iDir = 1:length(dirList)
            % Decompose file name
            [fPath, fBase, fExt] = bst_fileparts(dirList{iDir});
            % Entry is a not a folder, or starts with a "."
            if ~isdir(dirList{iDir}) || isempty(fBase) || strcmpi(fBase(1),'.')
                continue;
            end
            % Skip if it is not the requested template
            if ~isempty(TemplateName) && ~strcmpi(fBase, TemplateName)
                continue;
            end
            % Get files list
            fileList = dir(bst_fullfile(dirList{iDir}, 'channel*.mat'));
            defaultsList = repmat(struct('fullpath','', 'name',''), 0);
            % Find all the valid defaults (channel files)
            for iFile = 1:length(fileList)
                [tmp__, baseName] = bst_fileparts(fileList(iFile).name);
                baseName = strrep(baseName, 'channel_', '');
                baseName = strrep(baseName, '_channel', '');
                baseName = strrep(baseName, '_', ' ');
                % Skip if it is not the requested template
                if ~isempty(SetupName) && ~strcmpi(baseName, SetupName)
                    continue;
                end
                % Add to list of templates
                iNewDefault = length(defaultsList) + 1;
                defaultsList(iNewDefault).fullpath = bst_fullfile(dirList{iDir}, fileList(iFile).name);
                defaultsList(iNewDefault).name = baseName;
            end
            % Add files list to defaults list
            if ~isempty(defaultsList)
                 fullDefaultsList(end + 1) = struct('contents', defaultsList, ...
                                                    'name',     fBase);
            end
        end
        % Return defaults list
        argout1 = fullDefaultsList;
        
        
%% ==== GET FILENAMES ====
    case 'GetFilenames'
        iStudies = varargin{2};
        iItems = varargin{3};
        DataType = varargin{4};
        FileNames = cell(1, length(iStudies));
        argout1 = {};
        for i = 1:length(iStudies)
            % Get study definition
            sStudy = bst_get('Study', iStudies(i));
            if isempty(sStudy)
                continue;
            end
            % Recordings or sources
            switch (DataType)
                case 'data'
                    if (iItems(i) > length(sStudy.Data))
                        return;
                    end
                    FileNames{i} = sStudy.Data(iItems(i)).FileName;
                case 'results'
                    if (iItems(i) > length(sStudy.Result))
                        return;
                    end
                    FileNames{i} = sStudy.Result(iItems(i)).FileName;
                case 'timefreq'
                    if (iItems(i) > length(sStudy.Timefreq))
                        return;
                    end
                    FileNames{i} = sStudy.Timefreq(iItems(i)).FileName;
                case 'matrix'
                    if (iItems(i) > length(sStudy.Matrix))
                        return;
                    end
                    FileNames{i} = sStudy.Matrix(iItems(i)).FileName;
                case {'pdata','presults','ptimfreq'}
                    if (iItems(i) > length(sStudy.Stat))
                        return;
                    end
                    FileNames{i} = sStudy.Stat(iItems(i)).FileName;
            end
        end
        argout1 = FileNames;

        
%% ==== GUI ====
    case 'BstFrame'
        if isempty(GlobalData) || isempty(GlobalData.Program.GUI) || isempty(GlobalData.Program.GUI.mainWindow)
            argout1 = [];
        else
            argout1 = GlobalData.Program.GUI.mainWindow.jBstFrame;
        end
    case 'BstControls'
        if isempty(GlobalData) || isempty(GlobalData.Program) || isempty(GlobalData.Program.GUI) || isempty(GlobalData.Program.GUI.mainWindow)
            argout1 = [];
        else
            argout1 = GlobalData.Program.GUI.mainWindow;
        end
    case 'isGUI'
        if isempty(GlobalData) || isempty(GlobalData.Program) || ~isfield(GlobalData.Program, 'GuiLevel')
            argout1 = [];
        else
            argout1 = (GlobalData.Program.GuiLevel >= 1);
        end
    case 'GuiLevel'
        if isempty(GlobalData) || isempty(GlobalData.Program) || ~isfield(GlobalData.Program, 'GuiLevel')
            argout1 = [];
        else
            argout1 = GlobalData.Program.GuiLevel;
        end
    case 'ScreenDef'
        if isempty(GlobalData) || isempty(GlobalData.Program) || ~isfield(GlobalData.Program, 'ScreenDef')
            argout1 = [];
        else
            argout1 = GlobalData.Program.ScreenDef;
        end
    case 'DecorationSize'
        if isempty(GlobalData) || isempty(GlobalData.Program) || ~isfield(GlobalData.Program, 'DecorationSize')
            argout1 = [];
        else
            argout1 = GlobalData.Program.DecorationSize;
        end
    case 'Layout'
        % Default or current layout structure
        if ~isfield(GlobalData, 'Preferences') || ~isfield(GlobalData.Preferences, 'Layout') || ~((nargin == 1) || isfield(GlobalData.Preferences.Layout, varargin{2})) || ~isfield(GlobalData.Preferences.Layout, 'MainWindowPos')
            GlobalData.Preferences.Layout = db_template('Layout');
        end
        % Structure or property call
        if (nargin == 2) && ischar(varargin{2}) && isfield(GlobalData.Preferences.Layout, varargin{2})
            argout1 = GlobalData.Preferences.Layout.(varargin{2});
        elseif (nargin == 1)
            argout1 = GlobalData.Preferences.Layout;
        else
            error('Invalid call to bst_get.');
        end

    case 'ByteOrder'
        if isfield(GlobalData, 'Preferences') && isfield(GlobalData.Preferences, 'ByteOrder')
            argout1 = GlobalData.Preferences.ByteOrder;
        else
            argout1 = 'l';
        end

    case 'UniformizeTimeSeriesScales'
        if isfield(GlobalData, 'Preferences') && isfield(GlobalData.Preferences, 'UniformizeTimeSeriesScales')
            argout1 = GlobalData.Preferences.UniformizeTimeSeriesScales;
        else
            argout1 = 1;
        end 
        
    case 'FlipYAxis'
        if isfield(GlobalData, 'Preferences') && isfield(GlobalData.Preferences, 'FlipYAxis')
            argout1 = GlobalData.Preferences.FlipYAxis;
        else
            argout1 = 0;
        end 
        
    case 'AutoScaleY'
        if isfield(GlobalData, 'Preferences') && isfield(GlobalData.Preferences, 'AutoScaleY')
            argout1 = GlobalData.Preferences.AutoScaleY;
        else
            argout1 = 1;
        end
        
    case 'ShowXGrid'
        if isfield(GlobalData, 'Preferences') && isfield(GlobalData.Preferences, 'ShowXGrid')
            argout1 = GlobalData.Preferences.ShowXGrid;
        else
            argout1 = 0;
        end 
        
    case 'ShowYGrid'
        if isfield(GlobalData, 'Preferences') && isfield(GlobalData.Preferences, 'ShowYGrid')
            argout1 = GlobalData.Preferences.ShowYGrid;
        else
            argout1 = 0;
        end
        
    case 'ShowZeroLines'
        if isfield(GlobalData, 'Preferences') && isfield(GlobalData.Preferences, 'ShowZeroLines')
            argout1 = GlobalData.Preferences.ShowZeroLines;
        else
            argout1 = 1;
        end
        
    case 'Resolution'
        if isfield(GlobalData, 'Preferences') && isfield(GlobalData.Preferences, 'Resolution')
            argout1 = GlobalData.Preferences.Resolution;
        else
            argout1 = [0 0];
        end 
        
    case 'FixedScaleY'
        if isfield(GlobalData, 'Preferences') && isfield(GlobalData.Preferences, 'FixedScaleY') && isfield(GlobalData.Preferences.FixedScaleY, varargin{2}) && ~isempty(GlobalData.Preferences.FixedScaleY.(varargin{2}))
            argout1 = GlobalData.Preferences.FixedScaleY.(varargin{2});
        else
            argout1 = [];
        end
        
    case 'XScale'
        if isfield(GlobalData, 'Preferences') && isfield(GlobalData.Preferences, 'XScale')
            argout1 = GlobalData.Preferences.XScale;
        else
            argout1 = 'linear';
        end
        
    case 'YScale'
        if isfield(GlobalData, 'Preferences') && isfield(GlobalData.Preferences, 'YScale')
            argout1 = GlobalData.Preferences.YScale;
        else
            argout1 = 'linear';
        end        
        
    case 'ShowEventsMode'
        if isfield(GlobalData, 'Preferences') && isfield(GlobalData.Preferences, 'ShowEventsMode')
            argout1 = GlobalData.Preferences.ShowEventsMode;
        else
            argout1 = 'dot';
        end

    case 'AutoUpdates'
        if isfield(GlobalData, 'Preferences') && isfield(GlobalData.Preferences, 'AutoUpdates')
            argout1 = GlobalData.Preferences.AutoUpdates;
        else
            argout1 = 1;
        end 
        
    case 'ForceMatCompression'
        if isfield(GlobalData, 'Preferences') && isfield(GlobalData.Preferences, 'ForceMatCompression')
            argout1 = GlobalData.Preferences.ForceMatCompression;
        else
            argout1 = 0;
        end 
        
    case 'IgnoreMemoryWarnings'
        if isfield(GlobalData, 'Preferences') && isfield(GlobalData.Preferences, 'IgnoreMemoryWarnings')
            argout1 = GlobalData.Preferences.IgnoreMemoryWarnings;
        else
            argout1 = 0;
        end 
        
    case 'SystemCopy'
        if isfield(GlobalData, 'Preferences') && isfield(GlobalData.Preferences, 'SystemCopy')
            argout1 = GlobalData.Preferences.SystemCopy;
        else
            argout1 = 0;
        end
        
    case 'ExpertMode'
        if isfield(GlobalData, 'Preferences') && isfield(GlobalData.Preferences, 'ExpertMode')
            argout1 = GlobalData.Preferences.ExpertMode;
        else
            argout1 = 0;
        end 
        
    case 'DisplayGFP'
        if isfield(GlobalData, 'Preferences') && isfield(GlobalData.Preferences, 'DisplayGFP')
            argout1 = GlobalData.Preferences.DisplayGFP;
        else
            argout1 = 1;
        end
        
    case 'DownsampleTimeSeries'
        if isfield(GlobalData, 'Preferences') && isfield(GlobalData.Preferences, 'DownsampleTimeSeries')
            if (GlobalData.Preferences.DownsampleTimeSeries == 1)
                GlobalData.Preferences.DownsampleTimeSeries = 5;
            end
            argout1 = GlobalData.Preferences.DownsampleTimeSeries;
        else
            argout1 = 5;
        end
        
    case 'GraphicsSmoothing'
        if isfield(GlobalData, 'Preferences') && isfield(GlobalData.Preferences, 'GraphicsSmoothing')
            argout1 = GlobalData.Preferences.GraphicsSmoothing;
        else
            argout1 = 5;
        end
        
    case 'DisableOpenGL'
        if isfield(GlobalData, 'Preferences') && isfield(GlobalData.Preferences, 'DisableOpenGL')
            argout1 = GlobalData.Preferences.DisableOpenGL;
        else
            argout1 = 0;
        end

    case 'InterfaceScaling'
        if isfield(GlobalData, 'Preferences') && isfield(GlobalData.Preferences, 'InterfaceScaling')
            argout1 = GlobalData.Preferences.InterfaceScaling;
        else
            % Get screen resolution
            if isfield(GlobalData, 'Program') && isfield(GlobalData.Program, 'ScreenDef') && isfield(GlobalData.Program.ScreenDef, 'javaPos') && ~isempty(GlobalData.Program.ScreenDef(1).javaPos)
                AvailableRes = [100 125 150 200 250 300 400];
                iRes = bst_closest(GlobalData.Program.ScreenDef(1).javaPos.width * 100 / 1920, AvailableRes);
                argout1 = AvailableRes(iRes);
            else
                argout1 = 100;
            end
        end
        
    case 'JOGLVersion'
        % If JOGL1 is available
        if exist('javax.media.opengl.GLCanvas', 'class')
            argout1 = 1;
        % If JOGL2 is available
        elseif exist('javax.media.opengl.awt.GLCanvas', 'class')
            argout1 = 2;
        % No JOGL available
        else
            argout1 = 0;
        end

    case 'TSDisplayMode'
        if isfield(GlobalData, 'Preferences') && isfield(GlobalData.Preferences, 'TSDisplayMode')
            argout1 = GlobalData.Preferences.TSDisplayMode;
        else
            argout1 = 'butterfly';
        end 

    case 'PluginCustomPath'
        if isfield(GlobalData, 'Preferences') && isfield(GlobalData.Preferences, 'PluginCustomPath') && ~isempty(GlobalData.Preferences.PluginCustomPath)
            argout1 = GlobalData.Preferences.PluginCustomPath;
        else
            argout1 = [];
        end
        
    case 'BrainSuiteDir'
        if isfield(GlobalData, 'Preferences') && isfield(GlobalData.Preferences, 'BrainSuiteDir') && ~isempty(GlobalData.Preferences.BrainSuiteDir)
            if isdir(GlobalData.Preferences.BrainSuiteDir) && file_exist(bst_fullfile(GlobalData.Preferences.BrainSuiteDir, 'bdp'))
                argout1 = GlobalData.Preferences.BrainSuiteDir;
            else
                argout1 = [];
            end
        else
            argout1 = [];
        end
        
    case 'SpmTpmAtlas'
        % Get template file
        tpmUser = bst_fullfile(bst_get('BrainstormUserDir'), 'defaults', 'spm', 'TPM.nii');
        if file_exist(tpmUser)
            argout1 = tpmUser;
            disp(['SPM12 template found: ' tpmUser]);
            return;
        end
        % If it does not exist: check in brainstorm3 folder
        tpmDistrib = bst_fullfile(bst_get('BrainstormHomeDir'), 'defaults', 'spm', 'TPM.nii');
        if file_exist(tpmDistrib)
            argout1 = tpmDistrib;
            disp(['SPM12 template found: ' tpmDistrib]);
            return;
        end
        % If it does not exist: check in spm12 folder
        PlugSpm = bst_plugin('GetInstalled', 'spm12');
        if ~isempty(PlugSpm)
            tpmSpm = bst_fullfile(PlugSpm.Path, PlugSpm.SubFolder, 'tpm', 'TPM.nii');
            if file_exist(tpmSpm)
                argout1 = tpmSpm;
                disp(['SPM12 template found: ' tpmSpm]);
                return;
            end
        else
            tpmSpm = '';
        end
        % Not found...
        disp('SPM12 template not found in any of the following folders:');
        disp([' - ' tpmUser]);
        disp([' - ' tpmDistrib]);
        if ~isempty(tpmSpm)
            disp([' - ' tpmSpm]);
        end
        % Return the preferred location: .brainstorm/defaults/spm/TPM.nii
        argout1 = tpmUser;
        
    case 'PythonExe'
        % Get saved value
        if isfield(GlobalData, 'Preferences') && isfield(GlobalData.Preferences, 'PythonExe') && ~isempty(GlobalData.Preferences.PythonExe)
            if file_exist(GlobalData.Preferences.PythonExe)
                argout1 = GlobalData.Preferences.PythonExe;
            else
                disp(['BST> Error: Python executable not found: ' GlobalData.Preferences.PythonExe]);
                argout1 = [];
            end
        else
            argout1 = [];
        end
        % If not defined in Brainstorm, but set in Matlab
        if isempty(argout1)
            [pyVer, PythonExe] = bst_python_ver();
            if ~isempty(PythonExe) && file_exist(PythonExe)
                disp(['BST> Found Python executable: ' PythonExe]);
                argout1 = PythonExe;
                bst_set('PythonExe', PythonExe);
            end
        end
        
    case 'ElectrodeConfig'
        % Get modality
        Modality = varargin{2};
        if isempty(Modality) || ~ismember(Modality, {'EEG','ECOG','SEEG','ECOG+SEEG'})
            disp(['GET> Invalid modality: ' Modality]);
            Modality = 'EEG';
        end
        % Value was saved previously
        if isfield(GlobalData, 'Preferences') && isfield(GlobalData.Preferences, 'ElectrodeConfig') && isfield(GlobalData.Preferences.ElectrodeConfig, Modality) && isfield(GlobalData.Preferences.ElectrodeConfig.(Modality), 'ContactDiameter')
            argout1 = GlobalData.Preferences.ElectrodeConfig.(Modality);
        % Get default value
        else
            switch (Modality)
                case 'EEG'
                    ElectrodeConfig.Type            = 'eeg';
                    ElectrodeConfig.ContactDiameter = 0.010;
                    ElectrodeConfig.ContactLength   = 0.002;
                    ElectrodeConfig.ElecDiameter    = [];
                    ElectrodeConfig.ElecLength      = [];
                case 'ECOG'
                    ElectrodeConfig.Type            = 'ecog';
                    ElectrodeConfig.ContactDiameter = 0.004;
                    ElectrodeConfig.ContactLength   = 0.001;
                    ElectrodeConfig.ElecDiameter    = 0.0005;
                    ElectrodeConfig.ElecLength      = [];
                case {'SEEG','ECOG+SEEG'}
                    ElectrodeConfig.Type            = 'seeg';
                    ElectrodeConfig.ContactDiameter = 0.0008;
                    ElectrodeConfig.ContactLength   = 0.002;
                    ElectrodeConfig.ElecDiameter    = 0.0007;
                    ElectrodeConfig.ElecLength      = 0.070;
            end
            argout1 = ElectrodeConfig;
        end
        
    case 'UseSigProcToolbox'
        % In a parfor loop: GlobalData is empty => Check only if the toolbox is installed (ignore user preferences) 
        if isempty(GlobalData) || ~isfield(GlobalData, 'Program') || ~isfield(GlobalData.Program, 'HasSigProcToolbox')
            argout1 = exist('fir2', 'file');
        else
            % Save the result of the check for the SigProc tb
            if isempty(GlobalData.Program.HasSigProcToolbox)
                % Check if Signal Processing Toolbox is installed
                GlobalData.Program.HasSigProcToolbox = (exist('fir2', 'file') == 2);
            end
            % Return user preferences
            if ~GlobalData.Program.HasSigProcToolbox
                argout1 = 0;
            elseif isfield(GlobalData, 'Preferences') && isfield(GlobalData.Preferences, 'UseSigProcToolbox')
                argout1 = GlobalData.Preferences.UseSigProcToolbox;
            else
                argout1 = 1;
            end
        end

    case 'CustomColormaps'
        if isfield(GlobalData, 'Preferences') && isfield(GlobalData.Preferences, 'CustomColormaps') && ~isempty(GlobalData.Preferences.CustomColormaps)
            argout1 = GlobalData.Preferences.CustomColormaps;
        else
            argout1 = repmat(struct('Name', '', 'CMap', []), 0);
        end
        
    case 'BFSProperties'
        if isfield(GlobalData, 'Preferences') && isfield(GlobalData.Preferences, 'BFSProperties') && ~isempty(GlobalData.Preferences.BFSProperties)
            argout1 = GlobalData.Preferences.BFSProperties;
        else
            argout1 = [.33 .0042 .33 .88 .93];
        end
        
    case 'LastUsedDirs'
        defPref = struct(...
            'ImportData',      '', ...
            'ImportChannel',   '', ...
            'ImportAnat',      '', ...
            'ImportMontage',   '', ...
            'ExportChannel',   '', ...
            'ExportData',      '', ...
            'ExportAnat',      '', ...
            'ExportProtocol',  '', ...
            'ExportImage',     '', ...
            'ExportScript',    '', ...
            'ExportMontage',   '');
        argout1 = FillMissingFields(contextName, defPref);
        % Check that all folders are valid
        fields = fieldnames(argout1);
        for i = 1:length(fields)
            if ~ischar(argout1.(fields{i})) || ~file_exist(argout1.(fields{i}))
                argout1.(fields{i}) = '';
            end
        end
        
    case 'DefaultFormats'
        defPref = struct(...
            'AnatIn',      'FreeSurfer', ...
            'ChannelIn',   '', ...
            'ChannelOut',  '', ...
            'DataIn',      'CTF', ...
            'DataOut',     '', ...
            'DipolesIn',   '', ...
            'DipolesOut',  '', ...
            'ImageOut',    '', ...
            'EventsIn',    '', ...
            'EventsOut',   '', ...
            'MriIn',       '', ...
            'MriOut',      'Nifti1', ...
            'NoiseCovIn',  '', ...
            'NoiseCovOut', '', ...
            'ResultsIn',   '', ...
            'ResultsOut',  '', ...
            'SpmOut',      '', ...
            'SspIn',       '', ...
            'SspOut',      '', ...
            'SurfaceIn',   '', ...
            'SurfaceOut',  '', ...
            'LabelIn',     '', ...
            'LabelOut',    '', ...
            'TimefreqIn',  '', ...
            'TimefreqOut', '', ...
            'MatrixIn',    '', ...
            'MatrixOut',   '', ...
            'MontageIn',   '', ...
            'MontageOut',  '', ...
            'FibersIn',    '');
        argout1 = FillMissingFields(contextName, defPref);

    case 'OsType'
        switch (mexext)
            case 'mexglx',    argout1 = 'linux32';
            case 'mexa64',    argout1 = 'linux64';
            case 'mexmaci',   argout1 = 'mac32';
            case 'mexmaci64', argout1 = 'mac64';
            case 'mexs64',    argout1 = 'sol64';
            case 'mexw32',    argout1 = 'win32';
            case 'mexw64',    argout1 = 'win64';
            otherwise,        error('Unsupported extension.');
        end
        % CALL: bst_get('OsType', isMatlab=0)
        if (nargin >= 2) && isequal(varargin{2}, 0)
            if strcmpi(argout1, 'win32') && (~isempty(strfind(java.lang.System.getProperty('java.home'), '(x86)')) || ~isempty(strfind(java.lang.System.getenv('ProgramFiles(x86)'), '(x86)')))
                argout1 = 'win64';
            end
        end
        
    case 'ImportDataOptions'
        defPref = db_template('ImportOptions');
        argout1 = FillMissingFields(contextName, defPref);
        
    case 'RawViewerOptions'
        defPref =  struct(...
            'PageDuration',   3, ...
            'RemoveBaseline', 'all', ...
            'UseCtfComp',     1, ...
            'Shortcuts',      []);
        defPref.Shortcuts = {...
            '1', 'event1', 'simple', []; ...   % Key, event name, event type (simple,extended,page), epoch time
            '2', 'event2', 'simple', []; ... 
            '3', 'event3', 'simple', []; ...
            '4', 'event4', 'simple', []; ...
            '5', 'event5', 'simple', []; ...
            '6', 'event6', 'simple', []; ...
            '7', 'event7', 'simple', []; ...
            '8', 'event8', 'simple', []; ...
            '9', 'event9', 'simple', []};
        argout1 = FillMissingFields(contextName, defPref);
        % If invalid PageDuration: reset to default
        if (argout1.PageDuration <= 0.1)
            argout1.PageDuration = defPref.PageDuration;
        end
        % If old shortcuts: reset to defaults
        if any(size(argout1.Shortcuts) ~= size(defPref.Shortcuts))
            disp('BST> Warning: Reset keyboard shortcuts to include new options.');
            argout1.Shortcuts = defPref.Shortcuts;
            bst_set('RawViewerOptions', argout1);
        end

    case 'MontageOptions'
        defPref = struct('Shortcuts', []);
        defPref.Shortcuts = {
           %'a', []; ... Note: A is reserved for All channels
            'b', []; ...
            'c', []; ...
            'd', []; ...
            'e', []; ...
            'f', []; ...
            'g', []; ...
            'h', []; ...
            'i', []; ...
            'j', []; ...
            'k', []; ...
            'l', []; ...
            'm', []; ...
            'n', []; ...
            'o', []; ...
            'p', []; ...
            'q', []; ...
            'r', []; ...
            's', []; ...
            't', []; ...
            'u', []; ...
            'v', []; ...
            'w', []; ...
            'x', []; ...
            'y', []; ...
            'z', []; ...
            };
        argout1 = FillMissingFields(contextName, defPref);
        
    case 'TopoLayoutOptions'
        defPref = struct(...
            'TimeWindow',      [], ...
            'WhiteBackground', 0, ...
            'ShowRefLines',    1, ...
            'ShowLegend',      1, ...
            'ContourLines',    10);
        argout1 = FillMissingFields(contextName, defPref);
        
    case 'StatThreshOptions'
        defPref = struct(...
            'pThreshold',    .05, ...
            'durThreshold',  0, ...
            'Correction',    'fdr', ...
            'Control',       [1 2 3]);
        argout1 = FillMissingFields(contextName, defPref);
        % Make sure that Control is not a string (previous brainstorm version)
        if ischar(argout1.Control)
            argout1.Control = defPref.Control;
        end
        % Make sure that 'no' is used instead of 'none' (previous brainstorm version)
        if strcmpi(argout1.Correction, 'none')
            argout1.Correction = 'no';
        end
        
    case 'ContactSheetOptions'
        defPref = struct(...
            'nImages',    20, ...
            'TimeRange',  [], ...
            'SkipVolume', 0.2);
        argout1 = FillMissingFields(contextName, defPref);
        
    case 'ProcessOptions'
        defPref = struct(...
            'SavedParam',    struct(), ...
            'MaxBlockSize',  100 / 8 * 1024 * 1024, ...   % 100Mb
            'LastMaxBlockSize',  100 / 8 * 1024 * 1024);  % 100Mb
        argout1 = FillMissingFields(contextName, defPref);
        
    case 'ImportEegRawOptions'
        defPref = struct(...
            'isCanceled',        0, ...
            'BaselineDuration',  0, ...
            'SamplingRate',      1000, ...
            'MatrixOrientation', 'channelXtime', ... % {'channelXtime', 'timeXchannel'}
            'VoltageUnits',      'V', ...            % {'\muV', 'mV', 'V'}
            'SkipLines',         0, ...
            'nAvg',              1, ...
            'isChannelName',     0);                 % 1 if the first entry contains the channel name
        argout1 = FillMissingFields(contextName, defPref);
        
    case 'BugReportOptions'
        defPref = struct(...
            'isEnabled',  0, ...
            'SmtpServer', 'mailhost.chups.jussieu.fr', ...
            'UserEmail',  '');
        argout1 = FillMissingFields(contextName, defPref);
        
    case 'DefaultSurfaceDisplay'
        defPref = struct(...
            'SurfShowSulci',    1, ...
            'SurfSmoothValue',  0, ...
            'DataThreshold',    0.5, ...
            'SizeThreshold',    1, ...
            'DataAlpha',        0);
        argout1 = FillMissingFields(contextName, defPref);
        
    case 'MagneticExtrapOptions'
        defPref = struct(...
            'ForceWhitening', 0, ...
            'EpsilonValue',   0.0001);
        argout1 = FillMissingFields(contextName, defPref);
        
    case 'DefaultFreqBands'
        argout1 = {...
            'delta',  '2, 4',   'mean'; ...
            'theta',  '5, 7',   'mean'; ...
            'alpha',  '8, 12',  'mean'; ...
            'beta',   '15, 29', 'mean'; ...
            'gamma1', '30, 59', 'mean'; ...
            'gamma2', '60, 90', 'mean'};
        
    case 'TimefreqOptions_morlet'
        defPref.isTimeBands     = 0;
        defPref.isFreqBands     = 0;
        defPref.isFreqLog       = 0;
        defPref.TimeBands       = {};
        defPref.Freqs           = '1:1:60';
        defPref.FreqsLog        = '1:40:150';
        defPref.FreqBands       = bst_get('DefaultFreqBands');
        defPref.Measure         = 'power';
        defPref.SaveKernel      = 0;
        defPref.Output          = 'all';
        defPref.RemoveEvoked    = 0;
        defPref.ClusterFuncTime = 'after';
        defPref.MorletFc        = 1;
        defPref.MorletFwhmTc    = 3;
        argout1 = FillMissingFields(contextName, defPref);
        if isempty(argout1.Freqs)
            argout1.Freqs = defPref.Freqs;
        end
        if ~isempty(argout1.FreqBands) && ((size(argout1.FreqBands,2) ~= 3) || ~all(cellfun(@ischar, argout1.FreqBands(:))) || any(cellfun(@(c)isempty(strtrim(c)), argout1.FreqBands(:))))
            argout1.FreqBands = defPref.FreqBands;
        end

    case 'TimefreqOptions_hilbert'
        defPref.isTimeBands     = 0;
        defPref.isFreqBands     = 1;
        defPref.isFreqLog       = 0;
        defPref.TimeBands       = {};
        defPref.Freqs           = [];
        defPref.FreqsLog        = [];
        defPref.FreqBands       = bst_get('DefaultFreqBands');
        defPref.Measure         = 'power';
        defPref.SaveKernel      = 0;
        defPref.Output          = 'all';
        defPref.RemoveEvoked    = 0;
        defPref.ClusterFuncTime = 'after';
        argout1 = FillMissingFields(contextName, defPref);
        if isempty(argout1.Freqs)
            argout1.Freqs = defPref.Freqs;
        end
        if ~isempty(argout1.FreqBands) && (size(argout1.FreqBands,2) == 3) && ~ischar(argout1.FreqBands{1,2})
            argout1.FreqBands = defPref.FreqBands;
        end
        
    case 'TimefreqOptions_plv'
        defPref.isTimeBands     = 0;
        defPref.isFreqBands     = 1;
        defPref.isFreqLog       = 0;
        defPref.TimeBands       = {};
        defPref.Freqs           = [];
        defPref.FreqsLog        = [];
        defPref.FreqBands       = bst_get('DefaultFreqBands');
        defPref.Measure         = 'other';
        defPref.SaveKernel      = 0;
        defPref.Output          = 'all';
        defPref.ClusterFuncTime = 'after';
        argout1 = FillMissingFields(contextName, defPref);
        if isempty(argout1.Freqs)
            argout1.Freqs = defPref.Freqs;
        end
        if ~isempty(argout1.FreqBands) && ~ischar(argout1.FreqBands{1,2})
            argout1.FreqBands = defPref.FreqBands;
        end
        
    case 'TimefreqOptions_fft'
        defPref.isTimeBands     = 0;
        defPref.isFreqBands     = 0;
        defPref.isFreqLog       = 0;
        defPref.TimeBands       = {};
        defPref.Freqs           = [];
        defPref.FreqsLog        = [];
        defPref.FreqBands       = bst_get('DefaultFreqBands');
        defPref.Measure         = 'power';
        defPref.Output          = 'all';
        defPref.ClusterFuncTime = 'after';
        argout1 = FillMissingFields(contextName, defPref);
        if isempty(argout1.Freqs)
            argout1.Freqs = defPref.Freqs;
        end
        if ~isempty(argout1.FreqBands) && ~ischar(argout1.FreqBands{1,2})
            argout1.FreqBands = defPref.FreqBands;
        end
        
    case 'TimefreqOptions_psd'
        defPref.isTimeBands     = 0;
        defPref.isFreqBands     = 0;
        defPref.isFreqLog       = 0;
        defPref.TimeBands       = {};
        defPref.Freqs           = [];
        defPref.FreqsLog        = [];
        defPref.FreqBands       = bst_get('DefaultFreqBands');
        defPref.Measure         = 'power';
        defPref.Output          = 'all';
        defPref.ClusterFuncTime = 'after';
        argout1 = FillMissingFields(contextName, defPref);
        if isempty(argout1.Freqs)
            argout1.Freqs = defPref.Freqs;
        end
        if ~isempty(argout1.FreqBands) && ~ischar(argout1.FreqBands{1,2})
            argout1.FreqBands = defPref.FreqBands;
        end
    
    case 'ExportBidsOptions'
        defPref.ProjName    = [];
        defPref.ProjID      = [];
        defPref.ProjDesc    = [];
        defPref.Categories  = [];
        defPref.JsonDataset = ['{' 10 '  "License": "PD"' 10 '}'];
        defPref.JsonMeg     = ['{' 10 '  "TaskDescription": "My task"' 10 '}'];
        argout1 = FillMissingFields(contextName, defPref);        
        
    case 'OpenMEEGOptions'
        defPref.BemFiles     = {};
        defPref.BemNames     = {'Scalp', 'Skull', 'Brain'};
        defPref.BemCond      = [1, 0.0125, 1];
        defPref.BemSelect    = [1 1 1];
        defPref.isAdjoint    = 0;
        defPref.isAdaptative = 1;
        defPref.isSplit      = 0;
        defPref.SplitLength  = 4000;
        argout1 = FillMissingFields(contextName, defPref);
        
    case 'DuneuroOptions'
        defPref = duneuro_defaults();        
        argout1 = FillMissingFields(contextName, defPref);
        
    case 'GridOptions_dipfit'
        defPref = struct(...
            'Method',         'isotropic', ...
            'nLayers',        17, ...
            'Reduction',      3, ...
            'nVerticesInit',  4000, ...
            'Resolution',     0.020, ...
            'FileName',       '');        
        argout1 = FillMissingFields(contextName, defPref);
        
    case 'GridOptions_headmodel'
        defPref = struct(...
            'Method',         'isotropic', ...
            'nLayers',        17, ...
            'Reduction',      3, ...
            'nVerticesInit',  4000, ...
            'Resolution',     0.005, ...
            'FileName',       '');
        argout1 = FillMissingFields(contextName, defPref);
        
    case 'MriOptions'
        defPref = struct(...
            'isRadioOrient',    0, ...
            'isMipAnatomy',     0, ...
            'isMipFunctional',  0, ...
            'OverlaySmooth',    0, ...
            'InterpDownsample', 3, ...
            'DistanceThresh',   6, ...
            'UpsampleImage',    0, ...
            'DefaultAtlas',     []);
        argout1 = FillMissingFields(contextName, defPref);
        
    case 'DigitizeOptions'
        defPref = struct(...
            'ComPort',      'COM1', ...
            'ComRate',      9600, ...
            'ComByteCount', 94, ...  % 47 bytes * 2 receivers
            'UnitType',     'fastrak', ...
            'PatientId',    'S001', ...
            'nFidSets',     2, ...
            'isBeep',       1, ...
            'isMEG',        1, ...
            'isSimulate',   0, ...
            'Montages',     [...
                struct('Name',   'No EEG', ...
                       'Labels', []), ...
                struct('Name',   'Default', ...
                       'Labels', [])], ...
            'iMontage',     1);
        argout1 = FillMissingFields(contextName, defPref);
    
    case 'ConnectGraphOptions'
        defPref = struct(...
            'LobeFullLabel', 1, ... 
            'TextDisplayMode', [1 2], ...
            'LabelSize', 7, ...  
            'NodeSize', 5, ...         
            'LinkSize', 1.5, ...           
            'BgColor', [0 0 0], ...        
            'HierarchyNodeIsVisible', 1);
        % If we have an additional argument, get the default values
        if nargin > 1
            argout1 = defPref;
        % Otherwise, get the saved values
        else
            savedValues = FillMissingFields(contextName, defPref);
            
            % if any of the fields are [], replace by default value
            % do it here to avoid touching the common FillMissingFields
            % function, as other tools may actually want to set [] as desired property           
            fields = fieldnames(savedValues);
            for i=1:numel(fields)
                if(isempty(savedValues.(fields{i})))
                    savedValues.(fields{i}) = defPref.(fields{i});
                end
            end
            argout1 = savedValues;
        end
    
    case 'NodelistOptions'
        defPref = struct(...
            'String', '', ...         % What to search for
            'Target', 'Comment', ...  % What field to search for: {'FileName', 'Comment'}
            'Action', 'Select');      % What to do with the filtered files: {'Select', 'Exclude'}
        argout1 = FillMissingFields(contextName, defPref);
        
    case 'ReadOnly'
        if isfield(GlobalData.DataBase, 'isReadOnly')
            argout1 = GlobalData.DataBase.isReadOnly;
        else
            argout1 = 0;
        end
        
    case 'LastPsdDisplayFunction'
        if isfield(GlobalData, 'Preferences') && isfield(GlobalData.Preferences, 'LastPsdDisplayFunction')
            argout1 = GlobalData.Preferences.LastPsdDisplayFunction;
        else
            argout1 = [];
        end

    case 'PlotlyCredentials'
        try
            creds = loadplotlycredentials();
            argout1 = creds.username;
            argout2 = creds.api_key;
        catch
            argout1 = '';
            argout2 = '';
        end
        
        try
            config = loadplotlyconfig();
            argout3 = config.plotly_domain;
        catch
            argout3 = '';
        end

    case 'KlustersExecutable'
        if isfield(GlobalData, 'Preferences') && isfield(GlobalData.Preferences, 'KlustersExecutable')
            argout1 = GlobalData.Preferences.KlustersExecutable;
        else
            argout1 = [];
        end
        
        
%% ===== FILE FILTERS =====
    case 'FileFilters'
        switch lower(varargin{2})
            case 'mri'
                argout1 = {...
                    {'.img'},          'MRI: Analyze (*.img/*.hdr)',           'Analyze'; ...
                    {'.ima'},          'MRI: BrainVISA GIS (*.ima/*.dim)',     'GIS'; ...
                    {'.ima'},          'MRI: BrainVISA GIS (*.ima/*.dim)',     'GIS'; ...
                    {'.mri'},          'MRI: CTF (*.mri)',                     'CTF'; ...
                    {'.mat'},          'MRI: FieldTrip (*.mat)',               'FT-MRI'; ...
                    {'.mgh','.mgz'},   'MRI: MGH (*.mgh,*.mgz)',               'MGH'; ...
                    {'.mnc', '.mni'},  'MRI: MNI (*.mnc,*.mni)',               'MINC'; ...
                    {'.nii','.gz'},    'MRI: NIfTI-1 (*.nii;*.nii.gz)',        'Nifti1'; ...
                    {'_subjectimage'}, 'MRI: Brainstorm (*subjectimage*.mat)', 'BST'; ...
                    {'*'},             'MRI: DICOM (SPM converter)',           'DICOM-SPM'; ...
                    {'.mri', '.fif', '.img', '.ima', '.nii', '.mgh', '.mgz', '.mnc', '.mni', '.gz', '_subjectimage'}, 'All MRI files (subject space)', 'ALL'; ...
                    {'.mri', '.fif', '.img', '.ima', '.nii', '.mgh', '.mgz', '.mnc', '.mni', '.gz', '_subjectimage'}, 'All MRI files (MNI space)',     'ALL-MNI'; ...
                    {'.mri', '.fif', '.img', '.ima', '.nii', '.mgh', '.mgz', '.mnc', '.mni', '.gz', '_subjectimage'}, 'Volume atlas (subject space)',  'ALL-ATLAS'; ...
                    {'.mri', '.fif', '.img', '.ima', '.nii', '.mgh', '.mgz', '.mnc', '.mni', '.gz', '_subjectimage'}, 'Volume atlas (MNI space)',      'ALL-MNI-ATLAS'; ...
                   };
            case 'mriout'
                argout1 = {...
                    {'.img'}, 'MRI: Analyze (*.img/*.hdr)',         'Analyze'; ...
                    {'.ima'}, 'MRI: BrainVISA GIS (*.ima/*.dim)',   'GIS'; ...
                    {'.mri'}, 'MRI: CTF (*.mri)',                   'CTF'; ...
                    {'.mat'}, 'MRI: FieldTrip (*.mat)',             'FT-MRI'; ...
                    {'.nii'}, 'MRI: NIfTI-1 (*.nii)',               'Nifti1'...
                   };
            case 'anatin'
                argout1 = {...
                    {'.folder'}, 'FreeSurfer', 'FreeSurfer-fast'; ...
                    {'.folder'}, 'FreeSurfer + Volume atlases', 'FreeSurfer'; ...
                    {'.folder'}, 'FreeSurfer + Volume atlases + Thickness', 'FreeSurfer+Thick'; ...
                    {'.folder'}, 'BrainSuite', 'BrainSuite-fast'; ...
                    {'.folder'}, 'BrainSuite + Volume atlases', 'BrainSuite'; ...
                    {'.folder'}, 'BrainVISA', 'BrainVISA'; ...
                    {'.folder'}, 'CAT12', 'CAT12'; ...
                    {'.folder'}, 'CAT12 + Thickness', 'CAT12+Thick'; ...
                    {'.folder'}, 'CIVET', 'CIVET'; ...
                    {'.folder'}, 'CIVET + Thickness', 'CIVET+Thick'; ...
                    {'.folder'}, 'HCP MEG/anatomy (pipeline v3)', 'HCPv3'; ...
                    {'.folder'}, 'SimNIBS', 'SimNIBS'; ...
                   };
            case 'source4d'
                argout1 = {...
                    {'.folder'}, 'NIfTI-1 (*.nii)',                  'Nifti1';...
                    {'.folder'}, 'Analyze (*.img/*.hdr)',            'Analyze'; ...
                    {'.folder'}, 'Matlab 4D matrix (*voltime*.mat)', 'BST'; ...
                   };
                    
            case 'surface'
                argout1 = {...
                    {'.mesh'},  'BrainVISA (*.mesh)',      'MESH'; ...
                    {'_tess', '_head', '_scalp', '_brain', '_cortex', '_innerskull', '_outerskull'}, 'Brainstorm (*.mat)', 'BST'; ...
                    {'.dfs'},   'BrainSuite (*.dfs)',      'DFS'; ...
                    {'.dsgl'},  'BrainSuite old (*.dsgl)', 'DSGL'; ...
                    {'.bd0','.bd1','.bd2','.bd3','.bd4','.bd5','.bd6','.bd7','.bd8','.bd9', ...
                     '.s00','.s01','.s02','.s03','.s04','.s05','.s06','.s07','.s08','.s09'}, ...
                                'Curry BEM (*.db*;*.s0*)', 'CURRY-BEM';
                    {'.vtk'},   'FSL: VTK (*.vtk)',        'VTK'; ...
                    {'*'},      'FreeSurfer (*.*)',        'FS';
                    {'.off'},   'Geomview OFF (*.off)',    'OFF'; ...
                    {'.gii'},   'GIfTI / MRI coordinates (*.gii)', 'GII'; ...
                    {'.gii'},   'GIfTI / MNI coordinates (*.gii)', 'GII-MNI'; ...
                    {'.gii'},   'GIfTI / World coordinates (*.gii)', 'GII-WORLD'; ...
                    {'.fif'},   'MNE (*.fif)',             'FIF'; ...
                    {'.obj'},   'MNI OBJ (*.obj)',         'MNIOBJ'; ...
                    {'.msh'},   'SimNIBS Gmsh4 (*.msh)',   'SIMNIBS'; ...
                    {'.tri'},   'TRI (*.tri)',             'TRI'; ...
                    {'.mri', '.fif', '.img', '.ima', '.nii', '.mgh', '.mgz', '.mnc', '.mni', '.gz', '_subjectimage'}, 'Volume mask or atlas (subject space)', 'MRI-MASK'; ...
                    {'.mri', '.fif', '.img', '.ima', '.nii', '.mgh', '.mgz', '.mnc', '.mni', '.gz'},                  'Volume mask or atlas (MNI space)',     'MRI-MASK-MNI'; ...
                    {'.nwbaux'},   'Neurodata Without Borders (*.nwbaux)', 'NWB'; ...
                    {'*'},      'All surface files (*.*)', 'ALL'; ...
                   };
               
            case 'surfaceout'
                argout1 = {...
                    {'.mesh'}, 'BrainVISA (*.mesh)',   'MESH'; ...
                    {'.dfs'},  'BrainSuite (*.dfs)',   'DFS'; ...
                    {'.fs'},   'FreeSurfer (*.fs)',    'FS'
                    {'.off'},  'Geomview OFF (*.off)', 'OFF'; ...
                    {'.gii'},  'GIfTI (*.gii)',        'GII'; ...
                    {'.tri'},  'TRI (*.tri)',          'TRI'; ...
                   };
                
            case 'data'
                argout1 = {...
                     {'.*'},                 'MEG/EEG: 4D-Neuroimaging/BTi (*.*)',   '4D'; ...
                     {'.meg4','.res4'},      'MEG/EEG: CTF (*.ds;*.meg4;*.res4)',    'CTF'; ...
                     {'.fif'},               'MEG/EEG: Elekta-Neuromag (*.fif)',     'FIF'; ...
                     {'.mat'},               'MEG/EEG: FieldTrip (*.mat)',           'FT-TIMELOCK'; ...
                     {'.raw'},               'MEG/EEG: ITAB (*.raw)',                'ITAB'; ...
                     {'.kdf'},               'MEG/EEG: KRISS MEG (*.kdf)',           'KDF'; ...
                     {'.mrk','.sqd','.con','.raw','.ave'},  'MEG/EEG: Ricoh (*.sqd;*.con;*.raw;*.ave;*.mrk)', 'RICOH'; ...
                     {'.mat'},               'MEG/EEG: SPM (*.mat/.dat)',            'SPM-DAT'; ...
                     {'.mrk','.sqd','.con','.raw','.ave'},  'MEG/EEG: Yokogawa/KIT (*.sqd;*.con;*.raw;*.ave;*.mrk)', 'KIT'; ...
                     {'.hdf5'},              'MEG/EEG: York Instruments MEGSCAN (.hdf5)', 'MEGSCAN-HDF5'; ...
                     {'.bst'},               'MEG/EEG: Brainstorm binary (*.bst)',   'BST-BIN'; ...
                     {'.adicht'},            'EEG: ADInstruments LabChart (*.adicht)', 'EEG-ADICHT'; ...
                     {'.msr'},               'EEG: ANT ASA (*.msr)',                 'EEG-ANT-MSR'; ...
                     {'.cnt','.avr'},        'EEG: ANT EEProbe (*.cnt;*.avr)',       'EEG-ANT-CNT'; ...
                     {'*'},                  'EEG: ASCII text (*.*)',                'EEG-ASCII'; ...
                     {'.raw'},               'EEG: Axion AxIS (*.raw)',              'EEG-AXION'; ...
                     {'.bdf'},               'EEG: BDF (*.bdf)',                     'EEG-BDF'; ...
                     {'.avr','.mux','.mul'}, 'EEG: BESA exports (*.avr;*.mul;*.mux)', 'EEG-BESA'; ...
                     {'.ns1','.ns2','.ns3','.ns4','.ns5','.ns6'}, 'EEG: Blackrock NeuroPort (*.nsX/*.nev)', 'EEG-BLACKROCK';
                     {'.eeg','.dat'},        'EEG: BrainAmp (*.eeg;*.dat)',          'EEG-BRAINAMP'; ...
                     {'.txt'},               'EEG: BrainVision Analyzer (*.txt)',    'EEG-BRAINVISION'; ...
                     {'.sef','.ep','.eph'},  'EEG: Cartool (*.sef;*.ep;*.eph)',      'EEG-CARTOOL'; ...
                     {'.dat','.cdt'},        'EEG: Curry (*.dat;*.cdt)',             'EEG-CURRY'; ...
                     {'.smr','.son'},        'EEG: CED Spike2 old 32bit (*.smr;*.son)',  'EEG-SMR'; ...
                     {'.smr','.smrx'},       'EEG: CED Spike2 new 64bit (*.smr;*.smrx)', 'EEG-SMRX'; ...
                     {'.rda'},               'EEG: Compumedics ProFusion Sleep (*.rda)',  'EEG-COMPUMEDICS-PFS'; ...
                     {'.bin'},               'EEG: Deltamed Coherence-Neurofile (*.bin)', 'EEG-DELTAMED'; ...
                     {'.edf','.rec'},        'EEG: EDF / EDF+ (*.rec;*.edf)',        'EEG-EDF'; ...
                     {'.set'},               'EEG: EEGLAB (*.set)',                  'EEG-EEGLAB'; ...
                     {'.raw'},               'EEG: EGI Netstation RAW (*.raw)',      'EEG-EGI-RAW'; ...
                     {'.mff','.bin'},        'EEG: EGI-Philips (*.mff)',             'EEG-EGI-MFF'; ...
                     {'.edf'},               'EEG: EmotivPRO (*.edf)',               'EEG-EMOTIV'; ...
                     {'.erp','.hdr'},        'EEG: ERPCenter (*.hdr;*.erp)',         'EEG-ERPCENTER'; ...
                     {'.erp'},               'EEG: ERPLab (*.erp)',                  'EEG-ERPLAB'; ...
                     {'.mat','.hdf5'},       'EEG: g.tec Matlab (*.mat,*.hdf5)',     'EEG-GTEC'; ...
                     {'.rhd','.rhs'},        'EEG: Intan (*.rhd,*.rhs)',             'EEG-INTAN'; ...
                     {'.mb2'},               'EEG: MANSCAN (*.mb2)',                 'EEG-MANSCAN'; ...
                     {'.trc'},               'EEG: Micromed (*.trc)',                'EEG-MICROMED'; ...
                     {'.mat'},               'EEG: Matlab matrix (*.mat)',           'EEG-MAT'; ...
                     {'.csv'},               'EEG: Muse (*.csv)',                    'EEG-MUSE-CSV'; ...
                     {'.ncs'},               'EEG: Neuralynx (*.ncs)',               'EEG-NEURALYNX'; ...
                     {'.bin'},               'EEG: NeurOne session folder',          'EEG-NEURONE'; ...
                     {'.cnt','.avg','.eeg','.dat'}, 'EEG: Neuroscan (*.cnt;*.eeg;*.avg;*.dat)', 'EEG-NEUROSCAN'; ...
                     {'.eeg','.dat'},        'EEG: NeuroScope (*.eeg;*.dat)',        'EEG-NEUROSCOPE'; ...
                     {'.e'},                 'EEG: Nicolet (*.e)',                   'EEG-NICOLET'; ...
                     {'.eeg'},               'EEG: Nihon Kohden (*.eeg)',            'EEG-NK'; ...
                     {'.dat'},               'EEG: Open Ephys flat binary (*.dat)',  'EEG-OEBIN'; ...
                     {'.plx','.pl2'},        'EEG: Plexon (*.plx;*.pl2)',            'EEG-PLEXON'; ...
                     {'.ns1','.ns2','.ns3','.ns4','.ns5','.ns6'}, 'EEG: Ripple Trellis (*.nsX/*.nev)', 'EEG-RIPPLE'; ...
                     {'.h5'},                'EEG: The Virtual Brain (*_TimeSeriesEEG.h5)', 'EEG-TVB'; ...
                     {'.csv'},               'EEG: Wearable Sensing (*.csv)',        'EEG-WS-CSV'; ...
                     {'.nirs'},              'NIRS: Brainsight (*.nirs)',            'NIRS-BRS'; ...
                     {'.bnirs','.jnirs','.snirf'}, 'NIRS: SNIRF (*.snirf)',          'NIRS-SNIRF'; ...
                     {'.edf'},               'EyeLink eye tracker (*.edf)',          'EYELINK'; ...
                    };
            case 'raw'
                argout1 = {...
                     {'.*'},                 'MEG/EEG: 4D-Neuroimaging/BTi (*.*)',   '4D'; ...
                     {'.meg4','.res4'},      'MEG/EEG: CTF (*.ds;*.meg4;*.res4)',    'CTF'; ...
                     {'.fif'},               'MEG/EEG: Elekta-Neuromag (*.fif)',     'FIF'; ...
                     {'.mat'},               'MEG/EEG: FieldTrip (*.mat)',           'FT-TIMELOCK'; ...
                     {'.raw'},               'MEG/EEG: ITAB (*.raw)',                'ITAB'; ...
                     {'.kdf'},               'MEG/EEG: KRISS MEG (*.kdf)',           'KDF'; ...
                     {'.mrk','.sqd','.con','.raw','.ave'},  'MEG/EEG: Ricoh (*.sqd;*.con;*.raw;*.ave;*.mrk)', 'RICOH'; ...
                     {'.mat'},               'MEG/EEG: SPM (*.mat/.dat)',            'SPM-DAT'; ...
                     {'.mrk','.sqd','.con','.raw','.ave'},  'MEG/EEG: Yokogawa/KIT (*.sqd;*.con;*.raw;*.ave;*.mrk)', 'KIT'; ...
                     {'.hdf5'},              'MEG/EEG: York Instruments MEGSCAN (.hdf5)', 'MEGSCAN-HDF5'; ...
                     {'.bst'},               'MEG/EEG: Brainstorm binary (*.bst)',   'BST-BIN'; ...
                     {'.adicht'},            'EEG: ADInstruments LabChart (*.adicht)', 'EEG-ADICHT'; ...
                     {'.msr'},               'EEG: ANT ASA (*.msr)',                 'EEG-ANT-MSR'; ...
                     {'.cnt','.avr'},        'EEG: ANT EEProbe (*.cnt;*.avr)',       'EEG-ANT-CNT'; ...
                     {'*'},                  'EEG: ASCII text (*.*)',                'EEG-ASCII'; ...
                     {'.raw'},               'EEG: Axion AxIS (*.raw)',              'EEG-AXION'; ...
                     {'.bdf'},               'EEG: BDF (*.bdf)',                     'EEG-BDF'; ...
                     {'.avr','.mux','.mul'}, 'EEG: BESA exports (*.avr;*.mul;*.mux)', 'EEG-BESA'; ...
                     {'.ns1','.ns2','.ns3','.ns4','.ns5','.ns6'}, 'EEG: Blackrock NeuroPort (*.nsX/*.nev)', 'EEG-BLACKROCK';
                     {'.eeg','.dat'},        'EEG: BrainAmp (*.eeg;*.dat)',          'EEG-BRAINAMP'; ...
                     {'.txt'},               'EEG: BrainVision Analyzer (*.txt)',    'EEG-BRAINVISION'; ...
                     {'.sef','.ep','.eph'},  'EEG: Cartool (*.sef;*.ep;*.eph)',      'EEG-CARTOOL'; ...
                     {'.smr','.son'},        'EEG: CED Spike2 old 32bit (*.smr;*.son)',  'EEG-SMR'; ...
                     {'.smr','.smrx'},       'EEG: CED Spike2 new 64bit (*.smr;*.smrx)', 'EEG-SMRX'; ...
                     {'.rda'},               'EEG: Compumedics ProFusion Sleep (*.rda)',  'EEG-COMPUMEDICS-PFS'; ...
                     {'.dat','.cdt'},        'EEG: Curry (*.dat;*.cdt)',             'EEG-CURRY'; ...
                     {'.bin'},               'EEG: Deltamed Coherence-Neurofile (*.bin)', 'EEG-DELTAMED'; ...
                     {'.edf','.rec'},        'EEG: EDF / EDF+ (*.rec;*.edf)',        'EEG-EDF'; ...
                     {'.set'},               'EEG: EEGLAB (*.set)',                  'EEG-EEGLAB'; ...
                     {'.raw'},               'EEG: EGI Netstation RAW (*.raw)',      'EEG-EGI-RAW'; ...
                     {'.mff','.bin'},        'EEG: EGI-Philips (*.mff)',             'EEG-EGI-MFF'; ...
                     {'.edf'},               'EEG: EmotivPRO (*.edf)',               'EEG-EMOTIV'; ...
                     {'.mat','.hdf5'},       'EEG: g.tec Matlab (*.mat,*.hdf5)',     'EEG-GTEC'; ...
                     {'.rhd','.rhs'},        'EEG: Intan (*.rhd,*.rhs)',             'EEG-INTAN'; ...
                     {'.mb2'},               'EEG: MANSCAN (*.mb2)',                 'EEG-MANSCAN'; ...
                     {'.mat'},               'EEG: Matlab matrix (*.mat)',           'EEG-MAT'; ...
                     {'.csv'},               'EEG: Muse (*.csv)',                    'EEG-MUSE-CSV'; ...
                     {'.trc'},               'EEG: Micromed (*.trc)',                'EEG-MICROMED'; ...
                     {'.ncs'},               'EEG: Neuralynx (*.ncs)',               'EEG-NEURALYNX'; ...
                     {'.nwb'},               'EEG: Neurodata Without Borders (*.nwb)','NWB'; ...
                     {'.bin'},               'EEG: NeurOne session folder',          'EEG-NEURONE'; ...
                     {'.cnt','.avg','.eeg','.dat'}, 'EEG: Neuroscan (*.cnt;*.eeg;*.avg;*.dat)', 'EEG-NEUROSCAN'; ...
                     {'.eeg','.dat'},        'EEG: NeuroScope (*.eeg;*.dat)',        'EEG-NEUROSCOPE'; ...
                     {'.e'},                 'EEG: Nicolet (*.e)',                   'EEG-NICOLET'; ...
                     {'.eeg'},               'EEG: Nihon Kohden (*.eeg)',            'EEG-NK'; ...
                     {'.dat'},               'EEG: Open Ephys flat binary (*.dat)',  'EEG-OEBIN'; ...
                     {'.plx','.pl2'},        'EEG: Plexon (*.plx;.pl2)'              'EEG-PLEXON'; ...
                     {'.ns1','.ns2','.ns3','.ns4','.ns5','.ns6'}, 'EEG: Ripple Trellis (*.nsX/*.nev)', 'EEG-RIPPLE'; ...
                     {'.h5'},                'EEG: The Virtual Brain (*_TimeSeriesEEG.h5)', 'EEG-TVB'; ...
                     {'.tbk'},               'EEG: Tucker Davis Technologies (*.tbk)',    'EEG-TDT'; ...
                     {'.csv'},               'EEG: Wearable Sensing (*.csv)',        'EEG-WS-CSV'; ...
                     {'.trc','.eeg','.e','.bin','.rda','.edf','.bdf'}, 'SEEG: Deltamed/Micromed/NK/Nicolet/BrainAmp/EDF', 'SEEG-ALL'; ...
                     {'.trc','.eeg','.e','.bin','.rda','.edf','.bdf'}, 'ECOG: Deltamed/Micromed/NK/Nicolet/BrainAmp/EDF', 'ECOG-ALL'; ...
                     {'.nirs'},              'NIRS: Brainsight (*.nirs)',            'NIRS-BRS'; ...
                     {'.bnirs','.jnirs','.snirf'}, 'NIRS: SNIRF (*.snirf)',          'NIRS-SNIRF'; ...
                     {'.edf'},               'EyeLink eye tracker (*.edf)',          'EYELINK'; ...
                    };

            case 'dataout'
                argout1 = {...
                    {'.bst'},   'MEG/EEG: Brainstorm binary (*.bst)',          'BST-BIN'; ...
                    {'.mat'},   'MEG/EEG: FieldTrip timelock (*.mat)',         'FT-TIMELOCK'; ...
                    {'.mat'},   'MEG/EEG: SPM (*.mat/.dat)',                   'SPM-DAT'; ...
                    {'.eeg'},   'EEG: BrainVision BrainAmp (*.eeg)',           'EEG-BRAINAMP'; ...
                    {'.eph'},   'EEG: Cartool EPH (*.eph)',                    'EEG-CARTOOL-EPH'; ...
                    {'.raw'},   'EEG: EGI NetStation RAW (*.raw)',             'EEG-EGI-RAW'; ...
                    {'.edf'},   'EEG: European Data Format (*.edf)',           'EEG-EDF'; ...
                    {'.snirf'}, 'NIRS: SNIRF (*.snirf)',                       'NIRS-SNIRF'; ...
                    {'.txt'},   'ASCII: Space-separated, fixed column size (*.txt)', 'ASCII-SPC'; ...
                    {'.txt'},   'ASCII: Space-separated with header, fixed column size (*.txt)', 'ASCII-SPC-HDR'; ...
                    {'.tsv'},   'ASCII: Tab-separated (*.tsv)',                'ASCII-TSV'; ...
                    {'.tsv'},   'ASCII: Tab-separated with header (*.tsv)',    'ASCII-TSV-HDR'; ...
                    {'.tsv'},   'ASCII: Tab-separated with header transposed (*.tsv)',   'ASCII-TSV-HDR-TR'; ...
                    {'.csv'},   'ASCII: Comma-separated (*.csv)',              'ASCII-CSV'; ...
                    {'.csv'},   'ASCII: Comma-separated with header (*.csv)',  'ASCII-CSV-HDR'; ...
                    {'.csv'},   'ASCII: Comma-separated with header transposed (*.csv)', 'ASCII-CSV-HDR-TR'; ...
                    {'.xlsx'},  'Microsoft Excel (*.xlsx)',                    'EXCEL'; ...
                    {'.xlsx'},  'Microsoft Excel transposed (*.xlsx)',         'EXCEL-TR'; ...
                    {'_timeseries'}, 'Brainstorm matrix (*timeseries*.mat)',   'BST'; ...
                   };
            case 'rawout'
                argout1 = {...
                    {'.bst'},   'MEG/EEG: Brainstorm binary (*.bst)',  'BST-BIN'; ...
                    {'.mat'},   'MEG/EEG: SPM (*.mat/.dat)',           'SPM-DAT'; ...
                    {'.eeg'},   'EEG: BrainVision BrainAmp (*.eeg)',   'EEG-BRAINAMP'; ...
                    {'.raw'},   'EEG: EGI NetStation RAW (*.raw)',     'EEG-EGI-RAW'; ...
                    {'.edf'},   'EEG: European Data Format (*.edf)',   'EEG-EDF'; ...
                    {'.snirf'}, 'NIRS: SNIRF (*.snirf)',               'NIRS-SNIRF'; ...
                   };
            case 'events'
                argout1 = {...
                    {'.trg'},          'ANT EEProbe (*.trg)',           'ANT'; ...
                    {'.mrk'},          'AnyWave (*.mrk)',               'ANYWAVE'; ...
                    {'.evt'},          'BESA (*.evt)',                  'BESA'; ...
                    {'.tsv'},          'BIDS events: onset, duration, trial_type, channel (*.tsv)', 'BIDS'; ...
                    {'.vmrk'},         'BrainVision BrainAmp (*.vmrk)', 'BRAINAMP'; ...
                    {'_events'},       'Brainstorm (events*.mat)',      'BST'; ...
                    {'.mrk'},          'Cartool (*.mrk)',               'CARTOOL'; ...
                    {'.mrk'},          'CTF MarkerFile (*.mrk)',        'CTF'; ...
                    {'.cef'},          'Curry (*.cef)',                 'CURRY'; ...
                    {'.eve','.fif'},   'Elekta-Neuromag MNE (*.eve;*.fif)',    'FIF'; ...
                    {'.evl','.txt'},   'Elekta-Neuromag Graph (*.evl;*.txt)',  'GRAPH'; ...
                    {'.txt','.mat'},   'FieldTrip trial definition (*.txt;*.mat)', 'TRL'; ...
                    {'.trg'},          'KRISS MEG (*.trg)',             'KDF'; ...
                    {'.evt'},          'Micromed (*.evt)',              'MICROMED'; ...
                    {'.ev2'},          'Neuroscan (*.ev2)',             'NEUROSCAN'; ...
                    {'.txt'},          'Nicolet export (*.txt)',        'NICOLET'; ...
                    {'timestamps.npy'},'Open Ephys (timestamps.npy)', 'OEBIN'; ...
                    {'.log'},          'Presentation (*.log)',          'PRESENTATION'; ...
                    {'.mrk','.sqd','.con','.raw','.ave'},   'Ricoh (*.mrk;*.sqd;*.con;*.raw;*.ave)', 'RICOH'; ...
                    {'.txt'},          'XLTEK export (*.txt)',          'XLTEK'; ...
                    {'.mrk','.sqd','.con','.raw','.ave'},   'Yokogawa/KIT (*.mrk;*.sqd;*.con;*.raw;*.ave)', 'KIT'; ...
                    {'.*'},            'Array of times (*.mat;*.*)',    'ARRAY-TIMES'; ...
                    {'.*'},            'Array of samples (*.mat;*.*)',  'ARRAY-SAMPLES'; ...
                    {'.txt','.csv'},   'CSV text file: label, time, duration (*.txt;*.csv)', 'CSV-TIME'; ...
                    {'.*'},            'CTF Video Times (.txt)',        'CTFVIDEO'; ...
                    };
            case 'channel'
                argout1 = {...
                    {'.*'},                        'MEG/EEG: 4D-Neuroimaging/BTi (*.*)',  '4D'; ...
                    {'.meg4','.res4'},             'MEG/EEG: CTF (*.ds;*.meg4;*.res4)',   'CTF' ; ...
                    {'.fif'},                      'MEG/EEG: Elekta-Neuromag (*.fif)',    'FIF'; ...
                    {'.kdf'},                      'MEG/EEG: KRISS MEG (*.kdf)',          'KDF'; ...
                    {'.raw'},                      'MEG/EEG: ITAB (*.raw)',               'ITAB'; ...
                    {'.mrk','.sqd','.con','.raw','.ave'},  'MEG/EEG: Ricoh (*.sqd;*.con;*.raw;*.ave;*.mrk)', 'RICOH'; ...
                    {'.mrk','.sqd','.con','.raw','.ave'},  'MEG/EEG: Yokogawa/KIT (*.sqd;*.con;*.raw;*.ave;*.mrk)', 'KIT'; ...
                    {'.hdf5'},                     'MEG/EEG: York Instruments MEGSCAN (.hdf5)', 'MEGSCAN-HDF5'; ...
                    {'.bst'},                      'MEG/EEG: Brainstorm binary (*.bst)',   'BST-BIN'; ...
                    {'_channel'},                  'MEG/EEG: Brainstorm (channel*.mat)',  'BST'; ...
                    {'.elc'},                      'EEG: ANT ASA/Xensor (*.elc)',         'XENSOR'; ...
                    {'.sfp','.elp','.ela','.eps'}, 'EEG: BESA (*.sfp;*.elp;*.eps/*.ela)', 'BESA'; ...
                    {'.bvef','.bvct','.txt'},      'EEG: BrainVision electrode file (*.bvef,*.bvct,*.txt)', 'BRAINVISION'; ...
                    {'.tsv'},                      'EEG: BIDS electrodes.tsv, subject space mm (*.tsv)',    'BIDS-OTHER-MM'; ...
                    {'.tsv'},                      'EEG: BIDS electrodes.tsv, MNI space mm (*.tsv)',        'BIDS-MNI-MM'; ...
                    {'.els','.xyz'},               'EEG: Cartool (*.els;*.xyz)',          'CARTOOL'; ...
                    {'.eeg'},                      'EEG: MegDraw (*.eeg)',                'MEGDRAW'; ...
                    {'.res','.rs3','.pom'},        'EEG: Curry (*.res;*.rs3;*.pom)',      'CURRY'; ...
                    {'.ced','.xyz','.set'},        'EEG: EEGLAB (*.ced;*.xyz;*.set)',     'EEGLAB'; ...
                    {'.elc'},                      'EEG: EETrak (*.elc)',                 'EETRAK'; ...
                    {'.sfp'},                      'EEG: EGI (*.sfp)',                    'EGI'; ...
                    {'coordinates.xml'},           'EEG: EGI-Philips (coordinates.xml)',  'MFF'; ...
                    {'.elp'},                      'EEG: EMSE (*.elp)',                   'EMSE'; ...
                    {'.pts','.csv'},               'EEG: IntrAnat, subject space (*.pts;*.csv)', 'INTRANAT'; ...
                    {'.pts','.csv'},               'EEG: IntrAnat, MNI space (*.pts;*.csv)',     'INTRANAT_MNI'; ...
                    {'.csv'},                      'EEG: Localite (*.csv)',                      'LOCALITE'; ...
                    {'.dat','.tri','.txt','.asc'}, 'EEG: Neuroscan (*.dat;*.tri;*.txt;*.asc)',   'NEUROSCAN'; ...
                    {'.pos','.pol','.elp','.txt'}, 'EEG: Polhemus (*.pos;*.pol;*.elp;*.txt)',    'POLHEMUS'; ...
                    {'.csv'},                      'EEG: SimNIBS (*.csv)',             'SIMNIBS'; ...
                    {'.h5'},                       'EEG: The Virtual Brain (*_SensorsEEG.h5)',    'TVB'; ...
                    {'*'},                         'EEG: ASCII: Name,XYZ (*.*)',       'ASCII_NXYZ'; ...
                    {'*'},                         'EEG: ASCII: Name,XYZ_MNI (*.*)',   'ASCII_NXYZ_MNI'; ...
                    {'*'},                         'EEG: ASCII: Name,XYZ_World (*.*)', 'ASCII_NXYZ_WORLD'; ...
                    {'*'},                         'EEG: ASCII: Name,XY (*.*)',        'ASCII_NXY'; ...
                    {'*'},                         'EEG: ASCII: XYZ (*.*)',            'ASCII_XYZ'; ...
                    {'*'},                         'EEG: ASCII: XYZ_MNI (*.*)',        'ASCII_XYZ_MNI'; ...
                    {'*'},                         'EEG: ASCII: XYZ_World (*.*)',      'ASCII_XYZ_WORLD'; ...
                    {'*'},                         'EEG: ASCII: XY (*.*)',             'ASCII_XY'; ...
                    {'*'},                         'EEG: ASCII: XYZ,Name (*.*)',       'ASCII_XYZN'; ...
                    {'*'},                         'EEG: ASCII: XYZ_MNI,Name (*.*)',   'ASCII_XYZN_MNI'; ...
                    {'*'},                         'EEG: ASCII: XYZ_World,Name (*.*)', 'ASCII_XYZN_WORLD'; ...
                    {'*'},                         'EEG: ASCII: Name,Theta,Phi (*.*)', 'ASCII_NTP'; ...
                    {'*'},                         'EEG: ASCII: Theta,Phi (*.*)',      'ASCII_TP'; ...
                    };
            case 'channelout'
                argout1 = {...
                    {'.pos'}, 'EEG+Headshape: Polhemus (*.pos)',     'POLHEMUS'; ...
                    {'.eeg'}, 'Headshape: MegDraw (*.eeg)',          'MEGDRAW'; ...
                    {'.pos'}, 'Headshape: Polhemus (*.pos)',         'POLHEMUS-HS'; ...
                    {'.txt'}, 'Headshape: ASCII: XYZ (*.txt)',       'ASCII_XYZ-HS'; ...
                    {'.txt'}, 'Headshape: ASCII: XYZ_World (*.txt)', 'ASCII_XYZ_WORLD-HS'; ...
                    {'.txt'}, 'Headshape: ASCII: Name,XYZ (*.txt)',  'ASCII_NXYZ-HS'; ...
                    {'.txt'}, 'Headshape: ASCII: Name,XYZ_World (*.txt)',  'ASCII_NXYZ_WORLD-HS'; ...
                    {'.txt'}, 'Headshape: ASCII: XYZ,Name (*.txt)',  'ASCII_XYZN-HS'; ...
                    {'.txt'}, 'Headshape: ASCII: XYZ_World,Name (*.txt)',  'ASCII_XYZN_WORLD-HS'; ...
                    {'.sfp'}, 'EEG: BESA (*.sfp)',                   'BESA-SFP'; ...
                    {'.elp'}, 'EEG: BESA (*.elp)',                   'BESA-ELP'; ...
                    {'.xyz'}, 'EEG: Cartool (*.xyz)',                'CARTOOL-XYZ'; ...
                    {'.res'}, 'EEG: Curry (*.res)',                  'CURRY-RES'; ...
                    {'.xyz'}, 'EEG: EEGLAB (*.xyz)',                 'EEGLAB-XYZ'; ...
                    {'.sfp'}, 'EEG: EGI (*.sfp)',                    'EGI'; ...
                    {'.txt'}, 'EEG: ASCII: XYZ (*.txt)',             'ASCII_XYZ-EEG'; ...
                    {'.txt'}, 'EEG: ASCII: XYZ_MNI (*.*)',           'ASCII_XYZ_MNI-EEG'; ...
                    {'.txt'}, 'EEG: ASCII: XYZ_World (*.*)',         'ASCII_XYZ_WORLD-EEG'; ...
                    {'.txt'}, 'EEG: ASCII: Name,XYZ (*.txt)',        'ASCII_NXYZ-EEG'; ...
                    {'.txt'}, 'EEG: ASCII: Name,XYZ_MNI (*.*)',      'ASCII_NXYZ_MNI-EEG'; ...
                    {'.txt'}, 'EEG: ASCII: Name,XYZ_World (*.*)',    'ASCII_NXYZ_WORLD-EEG'; ...
                    {'.txt'}, 'EEG: ASCII: XYZ,Name (*.txt)',        'ASCII_XYZN-EEG'; ...
                    {'.txt'}, 'EEG: ASCII: XYZ_MNI,Name (*.*)',      'ASCII_XYZN_MNI-EEG'; ...
                    {'.txt'}, 'EEG: ASCII: XYZ_World,Name (*.*)',    'ASCII_XYZN_WORLD-EEG'; ...
                    {'.txt'}, 'NIRS: Brainsight (*.txt)',            'BRAINSIGHT-TXT'; ...
                    };
            case 'labelin'
                argout1 = {...
                    {'.dfs'},    'BrainSuite atlas (*.dfs)',      'DFS'; ...
                    {'.annot'},  'FreeSurfer atlas (*.annot)',    'FS-ANNOT'; ...
                    {'.label'},  'FreeSurfer ROI, single scout (*.label)',    'FS-LABEL-SINGLE'; ...
                    {'.label'},  'FreeSurfer ROI, probability map (*.label)', 'FS-LABEL'; ...
                    {'.gii'},    'GIfTI texture (*.gii)',         'GII-TEX'; ...
                    {'.dset'},   'SUMA atlas (*.dset)',           'DSET'; ...
                    {'_scout'},  'Brainstorm scouts (*scout*.mat)', 'BST'; ...
                    {'.mri', '.fif', '.img', '.ima', '.nii', '.mgh', '.mgz', '.mnc', '.mni', '.gz', '_subjectimage'}, 'Volume mask or atlas (dilated, subject space)', 'MRI-MASK'; ...
                    {'.mri', '.fif', '.img', '.ima', '.nii', '.mgh', '.mgz', '.mnc', '.mni', '.gz'},                  'Volume mask or atlas (dilated, MNI space)',     'MRI-MASK-MNI'; ...
                    {'.mri', '.fif', '.img', '.ima', '.nii', '.mgh', '.mgz', '.mnc', '.mni', '.gz', '_subjectimage'}, 'Volume mask or atlas (no overlap, subject space)', 'MRI-MASK-NOOVERLAP'; ...
                    {'.mri', '.fif', '.img', '.ima', '.nii', '.mgh', '.mgz', '.mnc', '.mni', '.gz'},                  'Volume mask or atlas (no overlap, MNI space)',     'MRI-MASK-NOOVERLAP-MNI'; ...
                    };
            case 'resultsout'
                argout1 = {...
                    {'_sources'}, 'Brainstorm sources (*sources*.mat)',        'BST'; ...
                    {'.mat'},   'FieldTrip sources (*.mat)',                   'FT-SOURCES'; ...
                    {'.txt'},   'ASCII: Space-separated, fixed columns size (*.txt)', 'ASCII-SPC'; ...
                    {'.txt'},   'ASCII: Space-separated with header, fixed column size (*.txt)', 'ASCII-SPC-HDR'; ...
                    {'.tsv'},   'ASCII: Tab-separated (*.tsv)',                'ASCII-TSV'; ...
                    {'.tsv'},   'ASCII: Tab-separated with header (*.tsv)',    'ASCII-TSV-HDR'; ...
                    {'.tsv'},   'ASCII: Tab-separated with header transposed (*.tsv)',   'ASCII-TSV-HDR-TR'; ...
                    {'.csv'},   'ASCII: Comma-separated (*.csv)',              'ASCII-CSV'; ...
                    {'.csv'},   'ASCII: Comma-separated with header (*.csv)',  'ASCII-CSV-HDR'; ...
                    {'.csv'},   'ASCII: Comma-separated with header transposed (*.csv)', 'ASCII-CSV-HDR-TR'; ...
                    {'.xlsx'},  'Microsoft Excel (*.xlsx)',                    'EXCEL'; ...
                    {'.xlsx'},  'Microsoft Excel transposed (*.xlsx)',         'EXCEL-TR'; ...
                    };
            case 'timefreqout'
                argout1 = {...
                    {'_timefreq'}, 'Brainstorm structure (*timefreq*.mat)',    'BST'; ...
                    {'.mat'},   'FieldTrip freq (*.mat)',                      'FT-FREQ'; ...
                    {'.txt'},   'ASCII: Space-separated, fixed columns size (*.txt)', 'ASCII-SPC'; ...
                    {'.txt'},   'ASCII: Space-separated with header, fixed column size (*.txt)', 'ASCII-SPC-HDR'; ...
                    {'.tsv'},   'ASCII: Tab-separated (*.tsv)',                'ASCII-TSV'; ...
                    {'.tsv'},   'ASCII: Tab-separated with header (*.tsv)',    'ASCII-TSV-HDR'; ...
                    {'.tsv'},   'ASCII: Tab-separated with header transposed (*.tsv)',   'ASCII-TSV-HDR-TR'; ...
                    {'.csv'},   'ASCII: Comma-separated (*.csv)',              'ASCII-CSV'; ...
                    {'.csv'},   'ASCII: Comma-separated with header (*.csv)',  'ASCII-CSV-HDR'; ...
                    {'.csv'},   'ASCII: Comma-separated with header transposed (*.csv)', 'ASCII-CSV-HDR-TR'; ...
                    {'.xlsx'},  'Microsoft Excel (*.xlsx)',                    'EXCEL'; ...
                    {'.xlsx'},  'Microsoft Excel transposed (*.xlsx)',         'EXCEL-TR'; ...
                    };
            case 'matrixout'
                argout1 = {...
                    {'_matrix'}, 'Brainstorm structure (*matrix*.mat)',        'BST'; ...
                    {'.mat'},   'FieldTrip timelock (*.mat)',                  'FT-TIMELOCK'; ...
                    {'.txt'},   'ASCII: Space-separated, fixed columns size (*.txt)', 'ASCII-SPC'; ...
                    {'.txt'},   'ASCII: Space-separated with header, fixed column size (*.txt)', 'ASCII-SPC-HDR'; ...
                    {'.tsv'},   'ASCII: Tab-separated (*.tsv)',                'ASCII-TSV'; ...
                    {'.tsv'},   'ASCII: Tab-separated with header (*.tsv)',    'ASCII-TSV-HDR'; ...
                    {'.tsv'},   'ASCII: Tab-separated with header transposed (*.tsv)',   'ASCII-TSV-HDR-TR'; ...
                    {'.csv'},   'ASCII: Comma-separated (*.csv)',              'ASCII-CSV'; ...
                    {'.csv'},   'ASCII: Comma-separated with header (*.csv)',  'ASCII-CSV-HDR'; ...
                    {'.csv'},   'ASCII: Comma-separated with header transposed (*.csv)', 'ASCII-CSV-HDR-TR'; ...
                    {'.xlsx'},  'Microsoft Excel (*.xlsx)',                    'EXCEL'; ...
                    {'.xlsx'},  'Microsoft Excel transposed (*.xlsx)',         'EXCEL-TR'; ...
                    };
            case 'montagein'
                argout1 = {...
                    {'.sel'},     'MNE selection files (*.sel)',              'MNE'; ...
                    {'.mon'},     'Text montage files (*.mon)',               'MON'; ...
                    {'_montage'}, 'Brainstorm montage files (montage_*.mat)', 'BST';
                    {'.csv'},     'Comma-separated montage files (*.csv)',    'CSV'; ...
                    {'.xml'},     'Compumedics ProFusion montages (*.xml)',   'EEG-COMPUMEDICS-PFS'};
            case 'montageout'
                argout1 = {...
                    {'.sel'},     'MNE selection files (*.sel)',              'MNE'; ...
                    {'.mon'},     'Text montage files (*.mon)',               'MON'; ...
                    {'_montage'}, 'Brainstorm montage files (montage_*.mat)', 'BST'};
            case 'fibers'
                argout1 = {...
                    {'.trk'},    'TrackVis (*.trk)',                       'TRK'; ...
                    {'_fibers'}, 'Brainstorm fibers files (fibers_*.mat)', 'BST'};
        end
        

%% ===== FONTS =====
    case 'FigFont'
        if ispc
            argout1 = 8;
        else
            argout1 = 9;
        end
        InterfaceScaling = bst_get('InterfaceScaling');
        if (InterfaceScaling ~= 100)
            argout1 = argout1 * InterfaceScaling / 100;
        end
        
    case 'Font'
        % Default font size
        if (nargin < 2)
            if strncmp(computer,'MAC',3)
                fontSize = 12;
            else
                fontSize = 11;
            end
        % Font size in input
        else
            fontSize = varargin{2};
        end
        % Adjust for interface scaling
        fontSize = fontSize * bst_get('InterfaceScaling') / 100;
        
        % Font types
        fontTypes = {};
        if (nargin >= 3)
            fontTypes{end + 1} = varargin{3};
        end
        fontTypes{end + 1} = 'Arial';  % Default font
        fontTypes{end + 1} = 'Liberation Sans';  % Free Arial substitute
        
        % Check for cached font
        foundFont = 0;
        for iFont = 1 : length(fontTypes)
            strCache = strrep(sprintf('%s%d', fontTypes{iFont}, round(fontSize*100)), ' ', '_');
            if ~isempty(GlobalData) && isfield(GlobalData, 'Program') && isfield(GlobalData.Program, 'FontCache') && isfield(GlobalData.Program.FontCache, strCache)
                argout1 = GlobalData.Program.FontCache.(strCache);
                foundFont = 1;
                break;
            end
        end
            
        % If font not cached, find first supported font
        if ~foundFont
            ge = java.awt.GraphicsEnvironment.getLocalGraphicsEnvironment();
            allFonts = cell(ge.getAvailableFontFamilyNames());
            
            for iFont = 1 : length(fontTypes)
                if any(strcmp(fontTypes{iFont}, allFonts))
                    fontType = fontTypes{iFont};
                    foundFont = 1;
                    break;
                end
            end
            
            if ~foundFont
                fontType = 'SansSerif';  % If nothing else works.
            end
            
            strCache = strrep(sprintf('%s%d', fontType, round(fontSize*100)), ' ', '_');
            argout1 = java.awt.Font(fontType, java.awt.Font.PLAIN, fontSize);
            GlobalData.Program.FontCache.(strCache) = argout1;
        end
        
%% ==== PANEL CONTAINERS ====
    case 'PanelContainer'    
        % Get Brainstorm GUI context structure
        bst_GUI = GlobalData.Program.GUI;
        if (isempty(bst_GUI) || ~isfield(bst_GUI, 'panelContainers'))
            error('Brainstorm GUI is not yet initialized');
        end
        
        % Get ContainerName in argument
        if ((nargin >= 2) && (ischar(varargin{2})))
            ContainerName = varargin{2};
        % If no container name in argument : just display all the container names
        else
            disp('Registered panel containers :');
            for iContainer = 1:length(bst_GUI.panelContainers)
                disp(['    - ' bst_GUI.panelContainers(iContainer).name]);
            end
            return
        end

        % Look for containerName in all the registered panel containers
        iContainer = 1;
        found = 0;
        while (~found (iContainer <= length(bst_GUI.panelContainers)))
             if (strcmpi(ContainerName, bst_GUI.panelContainers(iContainer).name))
                 found = 1;
             else
                 iContainer = iContainer + 1;
             end
        end
        % If container is found : return it
        if (found)
            argout1 = bst_GUI.panelContainers(iContainer).jHandle;
        else
            % warning('Brainstorm:InvalidContainer', 'Container ''%s'' could not be found.', ContainerName);
        end

        
%% ==== PANELS ====
    case 'Panel'    
        % Get Brainstorm GUI context structure
        if (isempty(GlobalData) || isempty(GlobalData.Program.GUI) || ~isfield(GlobalData.Program.GUI, 'panels'))
            return
        end
        listPanels = GlobalData.Program.GUI.panels;
        % Get Panel in argument
        if ((nargin >= 2) && (ischar(varargin{2})))
            PanelName = varargin{2};
        % If no panel name in argument : just display all the panels names
        else
            disp('Registered panels :');
            for iContainer = 1:length(listPanels)
                disp(['    - ' get(listPanels(iContainer), 'name')]);
            end
            return
        end
        % Look for panelName in all the registered panels
        iPanel = find(strcmpi(PanelName, get(listPanels, 'name')), 1);
        if ~isempty(iPanel)
            argout1 = listPanels(iPanel);
            argout2   = iPanel;
        end
                
        
%% ==== PANEL CONTROLS ====
%  Calls : bst_get('PanelControls', PanelName)
    case 'PanelControls'
        % Get Panel name in argument
        if ((nargin >= 2) && (ischar(varargin{2})))
            PanelName = varargin{2};
        else
            error('Invalid call to bst_get()');
        end
        % Find BstPanel with this name
        bstPanel = bst_get('Panel', PanelName);
        % If panel was found : return its controls
        if ~isempty(bstPanel)
            argout1 = get(bstPanel, 'sControls');
        end
        
%% ===== NODES COPY =====
%  Calls : bst_get('Clipboard')
    case 'Clipboard'
        argout1 = GlobalData.Program.Clipboard.Nodes;
        argout2 = GlobalData.Program.Clipboard.isCut;
        
%% ==== DIRECTORIES ====
    case 'DirDefaultSubject'
        argout1 = '@default_subject';
    case 'DirDefaultStudy'
        argout1 = '@default_study';
    case 'DirAnalysisIntra'
        argout1 = '@intra';
    case 'DirAnalysisInter'
        argout1 = '@inter';
    case 'NormalizedSubjectName'
        argout1 = 'Group_analysis';
        
%% ==== OTHER ====
    case 'ResizeFunction'
        if (bst_get('MatlabVersion') <= 803)
            argout1 = 'ResizeFcn';
        else
            argout1 = 'SizeChangedFcn';
        end
    case 'groot'
        if (bst_get('MatlabVersion') <= 803)
            argout1 = 0;
        else
            argout1 = groot;
        end
    case 'JFrame'
        hFig = varargin{2};
        MatlabVersion = bst_get('MatlabVersion');
        jFrame = [];
        try
            if (MatlabVersion <= 705)
                jf = get(hFig, 'javaframe');
                jFrame = jf.fFigureClient.getWindow();
            elseif (MatlabVersion <= 712)
                jf = get(handle(hFig), 'javaframe');
                jFrame = jf.fFigureClient.getWindow();
            elseif (MatlabVersion <= 803)
                jf = get(handle(hFig), 'javaframe');
                jFrame = jf.fHG1Client.getWindow();
            elseif (MatlabVersion < 907)    % Matlab >= 2019b deprecated the JavaFrame property
                warning('off', 'MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
                jf = get(hFig, 'javaframe');
                warning('on', 'MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
                jFrame = jf.fHG2Client.getWindow();
            else
                disp('BST> Error: Matlab 2019b deprecated the JavaFrame property.');
            end
        catch
            disp('BST> Warning: Cannot get the JavaFrame property for the selected figure.');
        end
        argout1 = jFrame;
        
%% ==== ERROR ====
    otherwise
        error(sprintf('Invalid context : "%s"', contextName));
end
end




%% ==== HELPERS ====
% Return all the protocol studies that have a given file in its structures
% Possible field names: Result.DataFile, Result.FileName, Data.FileName, Channel.FileName
%
% USAGE:  [sFoundStudy, iFoundStudy, iItem] = findFileInStudies(fieldGroup, fieldName, fieldFile, iStudiesList)
%         [sFoundStudy, iFoundStudy, iItem] = findFileInStudies(fieldGroup, fieldName, fieldFile)
function [sFoundStudy, iFoundStudy, iItem] = findFileInStudies(fieldGroup, fieldName, fieldFile, iStudiesList)
    global GlobalData;
    sFoundStudy = [];
    iFoundStudy = [];
    iItem       = [];
    % If no file provided, return
    if isempty(fieldFile)
        return;
    end
    % Extract folder(s) of the file(s) we're looking for
    fieldParts = strfind(fieldFile, '|');
    if ~isempty(fieldParts)
        fieldParts(end+1) = length(fieldFile);
        fieldFolders = {};
        iLast = 1;
        for iPart = 1:length(fieldParts)
            folder = fileparts(fieldFile(iLast:fieldParts(iPart)-1));
            if ~isempty(folder)
                fieldFolders{end + 1} = folder;
            end
            iLast = fieldParts(iPart) + 1;
        end
    else
        fieldFolders = {fileparts(fieldFile)};
    end
    % Get protocol information
    ProtocolStudies = GlobalData.DataBase.ProtocolStudies(GlobalData.DataBase.iProtocol);
    % List studies to process
    if (nargin < 4) || isempty(iStudiesList)
        iStudiesList = [-2, -3, 1:length(ProtocolStudies.Study)];
    end
    
    % NORMAL STUDIES: Look for surface file in all the surfaces of all subjects
    for iStudy = iStudiesList
        % Get study
        switch (iStudy)
            case -2,    sStudy = ProtocolStudies.AnalysisStudy;
            case -3,    sStudy = ProtocolStudies.DefaultStudy;
            otherwise,  sStudy = ProtocolStudies.Study(iStudy);
        end
        % Check if field is available for the study
        if isempty(sStudy.(fieldGroup))
            continue;
        end
         % Check we are in the correct folder
        if ~any(file_compare(fieldFolders, fileparts(sStudy.FileName)))
            continue;
        end
        % Get list of files from study
        filesList = {sStudy.(fieldGroup).(fieldName)};
        if isempty(filesList)
            continue;
        end
        % Replace empty cells with empty strings
        iValidFiles = find(cellfun(@ischar, filesList));
        if isempty(iValidFiles)
            continue;
        end
        % Find target in this list
        iItem = find(file_compare(filesList(iValidFiles), fieldFile));
        if ~isempty(iItem)
            sFoundStudy  = sStudy;
            iFoundStudy  = iStudy;
            iItem = iValidFiles(iItem);
            return
        end
    end
end


%% ===== FILL MISSING FIELDS =====
function bstPref = FillMissingFields(PrefName, defPref)
    global GlobalData;
    if isfield(GlobalData, 'Preferences') && isfield(GlobalData.Preferences, PrefName)
        bstPref = GlobalData.Preferences.(PrefName);
        bstPref = struct_copy_fields(bstPref, defPref, 0);
    else
        bstPref = defPref;
    end
end







