function [clusrand] = clusterrandanalysis(cfg,varargin);

% CLUSTERRANDANALYSIS
%
% Use as
%   [clusrand] = clusterrandanalysis(cfg, data)
% where the configuration can contain
%
% 1. Options for data selection and averaging over selected dimensions
% --------------------------------------------------------------------
%
%   cfg.channel          = Nx1 cell-array with selection of channels (default = 'all'),
%                          see CHANNELSELECTION for details
% or,
%   cfg.channelcmb       = Mx2 cell-array with selection of channel pairs (default = {'all' 'all'}),
%                          see CHANNELCOMBINATION for details
%   cfg.latency          = [begin end] in seconds or 'all' (default = 'all')
%   cfg.frequency        = [begin end], can be 'all'       (default = 'all')
%   cfg.avgoverchan      = 'yes' or 'no'                   (default = 'no')
%   cfg.avgovertime      = 'yes' or 'no'                   (default = 'no')
%   cfg.avgoverfreq      = 'yes' or 'no'                   (default = 'no')
%
% 2. Statistics options
% ---------------------
%
%   cfg.statistic       = 'indepsamplesT' (independent samples T-statistic) 
%                       | 'indepsamplesregrT' (independent samples regression coefficient T-statistic)
%                       | 'indepsamplesZcoh' (independent samples Z-statistic for coherence)
%                       | 'indepsamplesTsqcrs' (independent samples T-square-statistic for the cross-spectrum)
%                       | 'depsamplesT' (dependent samples T-statistic)
%                       | 'actvsblT' (activation versus baseline T-statistic)
%                       | 'depsamplesregrT' (dependent samples regression coefficient T-statistic) 
%                       | 'indepsamplesF' (independent samples F-statistic)
%                       | 'depsamplesF' (dependent samples F-statistic)
%
%   cfg.alpha           = alpha level of the randomization test
%   cfg.clusterteststat = 'maxsum' (default) or 'orderedsums'
%   cfg.smallestcluster = smallest cluster size that is large enough to be considered (only
%                         relevant when clusterteststat='orderedsums')
%
% For every statistic, we now list (a) the constraints with respect to the 
% data that should be passed as an argument to this function (in varargin) and 
% (b) a number of fields that are relevant for a particular choice of statistic. 
%
% For the independent samples T-statistic:
% Two data sets must be passed. The first data set belongs to condition 1
% and the second one belongs to condition 2.
%   cfg.onetwo         = 'onesided_1<2' (H_0:condition_1<condition_2), 
%                        'onesided_2<1' (H_0:condition_2<condition_1), or
%                        'twosided' (H_0:condition_1=condition_2, the
%                        default).
%
% For the independent samples regression coefficient T-statistic:
% When this statistic is chosen, a predictor variable has to be specified.
% This can be done in two ways: (1) as an array in the configuration 
% field cfg.ext, or (2) in de data field data.ext. If the second way of 
% specifying .ext is used, and there are multiple data sets in the argument, 
% then every data set must contain a scalar field data.ext. 
%   cfg.onetwo         = 'onesided_B<0' (H_0: regr. coefficient B<0), 
%                        'onesided_B>0' (H_0: regr. coefficient B>0), or
%                        'twosided' (H_0:regr. coefficient B=0, the default).
%
% For the independent samples Z-statistic for coherence:
% When this statistic is chosen, the spatial dimension of the data matrix must 
% contain (a) the cross-spectrum for a number of channel combinations, and
% (b) the auto-spectrum (power) for the channels that are involved in the
% channel combinations. With this test statistic, one examines whether the
% two conditions differ with respect to coherence. 
%   cfg.onetwo         = 'onesided_1<2' (H_0:condition_1<condition_2), 
%                        'onesided_2<1' (H_0:condition_2<condition_1), or
%                        'twosided' (H_0:condition_1=condition_2, the
%                        default).
%
% For the independent samples T-square-statistic for the cross-spectrum:
% When this statistic is chosen, the data matrix must contain the
% cross-spectrum for a number of channel combinations. 
% With this test statistic, one examines whether the two conditions differ 
% with respect to the average cross-spectrum.
%
% For the dependent samples T-statistic:
% Two data sets must be passed. The first data set belongs to condition 1
% and the second one belongs to condition 2. The two data sets must have
% the same number of replications (subjects or trials) and the order of
% these replications must be the same in both.
%   cfg.onetwo         = 'onesided_1<2' (H_0:condition_1<condition_2), 
%                        'onesided_2<1' (H_0:condition_2<condition_1), or
%                        'twosided' (H_0:condition_1=condition_2, the
%                        default).
%
% For the activation versus baseline T-statistic:
% Two data sets must be passed. The first data set belongs to the
% activation period and the second data set to the baseline period. 
% The two data sets must have the same number of replications 
% (subjects or trials) and the order of these replications must be the
% same. The first data set must have a positive time axis (the
% activation condition) and the second data set a negative time axis (the
% baseline condition). In the statistical comparison, the time axis of the baseline period is
% shifted forward in time such that its first sample corresponds to
% the first sample of the activation period. The activation versus baseline
% T-statistic is a dependent samples T-statistic that compares every sample
% in the activation period with the time average of the baseline period.
%   cfg.onetwo         = 'onesided_act<bl' (H_0:act<bl), 'onesided_bl<act' (H_0:bl<act), or
%                        'twosided' (H_0:act=bl, the default).
%
% For the dependent samples regression coefficient T-statistic:
% When this statistic is chosen, a predictor variable has to be specified.
% This can be done in two ways: (1) as an array in the configuration 
% field cfg.ext, or (2) in de data field data.ext. Contrary to the
% INdependent samples regression coefficient T-statistic, its DEpendent
% counterpart requires that multiple data sets are passed in the argument;
% one data set for every level of the predictor variable.
%   cfg.onetwo         = 'onesided_B<0' (H_0: regr. coefficient B<0), 
%                        'onesided_B>0' (H_0: regr. coefficient B>0), or
%                        'twosided' (H_0:regr. coefficient B=0, the default).
%
% For the independent samples F-statistic:
% Two or more data sets must be passed. 
% 
% For the dependent samples F-statistic:
% Two or more data sets must be passed. The order of the conditions 
% corresponds to the order of the data sets. The data sets must have
% the same number of replications (subjects or trials) and the order of
% these replications must be the same in both.
%   cfg.contrastcoefs  = matrix of contrast coefficients determining the
%                        effect being tested. The number of columns of this
%                        matrix must be equal to the number of conditions (ncond). 
%                        The default is a matrix that specifies the
%                        main effect of the independent variable. This matrix
%                        has size [(ncond-1),ncond]. 
%   cfg.allowedperms   = matrix of permutations of conditions that all have the
%                        same probability (given the data) under the null 
%                        hypothesis. The default is all ncond! possible permutations.
%
% 3. Clustering options
% ---------------------
% 
%   cfg.makeclusters = 'yes' (default) | 'no'
%   cfg.alphathresh  = alpha level of the (channel,timepoint,frequency)-specific 
%                    test statistic that will be used for thresholding 
%   cfg.minnbchan    = minimum number of neighbouring channels in which a 
%                    particular time-frequency-specific t-test has to be significant in order
%                    for it to be included in the clustering algorithm (default=0).
%   cfg.chancmbgeom = 'refchan' (one reference channel; the neighborhood geometry of the channel 
%                    combinations is determined by the non-reference channels) or 'free' 
%                    (the neighborhood geometry of the channel combinations is determined by the 
%                    channel combinations)
%
% 4. Neighborhood geometry
% ------------------------
%
% You can specify the neighbours of each channel (1) by providing them as a structured 
% cell array, (2) by looking at the sensor positions that are present in
% the data, or (3) by loading the neighborhood geometry from a file
% (cfg.geomfile). The latter option can also be used to load the geometry
% of the channel combinations.
% This is done either with
%   cfg.neighbourdist  = distance, default is 4 cm
% or with the field
%   cfg.neighbours     = definition of neighbours for each channel, see NEIGHBCHANSELECTION
% which should be structured like this:
%        cfg.neighbours{1}.label = 'Fz';
%        cfg.neighbours{1}.neighblabel = {'Cz', 'F3', 'F3A', 'FzA', 'F4A', 'F4'};
%        cfg.neighbours{2}.label = 'Cz';
%        cfg.neighbours{2}.neighblabel = {'Fz', 'F4', 'RT', 'RTP', 'P4', 'Pz', 'P3', 'LTP', 'LT', 'F3'};
%        cfg.neighbours{3}.label = 'Pz';
%        cfg.neighbours{3}.neighblabel = {'Cz', 'P4', 'P4P', 'Oz', 'P3P', 'P3'};
%     etc.
%        (Note that a channel is not considered to be a neighbour of itself.)
%     or which can be the name of a file that contains the logical matrix in which the
%        channel (combination) geometry is specified, together with the
%        channel (combination) labels.
% or with the field
%   cfg.geomfile = 'filename';
%
% Electrode or gradiometer positions are obtained from the first dataset, or can be specified
% using either one of 
%   cfg.gradfile      = string, file containing the gradiometer definition
%   cfg.elecfile      = string, file containing the electrode definition
% or alternatively
%   cfg.grad          = structure with gradiometer definition
%   cfg.elec          = structure with electrode definition
%
% 5. Miscellaneous options
% ------------------------
%
%   cfg.nranddraws     = number of draws from the randomization distribution
%   cfg.randomseed     = 'yes' (default), 'no' or seed-number
%   cfg.savegeom       = 'yes', 'no' (default), save the neighborhood
%                      geometry (of channels or channel combinations) in the file cfg.geomfile.
%   cfg.mirrordata     = 'yes', 'no' (default), to be used when the data are
%                      elements of the cross-spectrum. If the results of
%                      the statistical test depends on the order of the
%                      channel pairs (e.g., if data contain the imaginary
%                      part of the cross-spectrum and the the test statistic is
%                      one-dimensional) then cfg.mirrordata should be
%                      'yes'. Note that only the lower or upper diagonal of the
%                      cross-spectral matrix must be passed as an argument (in data).

% This function depends on PREPARE_TIMEFREQ_DATA which has the following options:
% cfg.avgoverchan, documented
% cfg.avgoverfreq, documented
% cfg.avgovertime, documented
% cfg.channel, documented            
% cfg.channelcmb, documented     
% cfg.datarepresentation (set in CLUSTERRANDANALYSIS; if between = 'concatenated' if within = 'cell-array')
% cfg.frequency, documented
% cfg.latency, documented         
% cfg.precision
% cfg.previous
% cfg.version

% Copyright (C) 2005-2006, Eric Maris, NICI, University Nijmegen
%
% $Log: clusterrandanalysis.m,v $
% Revision 1.26  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.25  2008/03/05 10:46:35  roboos
% moved electrode reading functionality from read_fcdc_elec to read_sens, switched to the use of the new function
%
% Revision 1.24  2007/04/03 15:37:07  roboos
% renamed the checkinput function to checkdata
%
% Revision 1.23  2007/03/27 11:05:18  ingnie
% changed call to fixdimord in call to checkinput
%
% Revision 1.22  2006/10/19 15:32:55  roboos
% updated documentation
%
% Revision 1.21  2006/10/04 07:10:07  roboos
% updated documentation
%
% Revision 1.20  2006/06/20 16:29:11  ingnie
% updated documentation added default cfg.channelcmb
%
% Revision 1.19  2006/06/13 14:51:03  ingnie
% updated documentation to increase consistency in help of cfg options, added
% defaults cfg.channel ='all', cfg.latency = 'all'
%
% Revision 1.18  2006/04/10 16:33:46  ingnie
% updated documentation
%
% Revision 1.17  2006/04/06 13:05:00  erimar
% Added help wrt channel combinations.
%
% Revision 1.16  2006/02/28 11:56:23  erimar
% Added functionality for channel combination data: test statistics indepsamplesZcoh (for coherence differences) and
% indepsamplesTsqcrs (for cross-spectrum differences), routines for calculating the neighborhood geometry for channel
% combinations, and the option to save the neighborhood geometry on file (cfg.geomfile, cfg.savegeom).
%

% Revision 1.15  2006/02/23 10:28:16  roboos
% changed dimord strings for consistency, changed toi and foi into time and freq, added fixdimord where neccessary
%
% Revision 1.14  2006/02/01 12:26:00  roboos
% made all uses of dimord consistent with the common definition of data dimensions, see the fixdimord() function
%
% Revision 1.13  2005/12/06 13:18:08  erimar
% (1) Removed all functionality with respect to channel combinations
% (because it was not sufficiently tested for data sets of realistic
% sizes). (2) Improved and extended the handling of external predictor
% variables (for the independent and the dependent samples regression
% coefficient T-statistics).
%
% Revision 1.12  2005/11/28 11:15:09  erimar
% Improved help
%
% Revision 1.11  2005/08/10 15:35:41  roboos
% fixed cfg.elecfile (accidentally called gradfile), thanks to chrfor
%
% Revision 1.10  2005/08/05 12:19:03  erimar
% Correct an error in the handling of the .ext-field for
% cfg.indepsamplesregrT detected by Vladimir Litvak.
%
% Revision 1.9  2005/08/05 07:48:37  roboos
% added support for cfg.elec/grad/elecfile/gradfile
% cleaned up the help
% removed the "defaults set in xxx" help comments, since those subfunctions are not accessible for the end-user
%
% Revision 1.8  2005/04/22 07:44:41  roboos
% added/corrected copyrights, added a Log tag for CVS, converted to unix ascii
%

fieldtripdefs

warning off;
pack;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(cfg, 'neighbours'),          cfg.neighbours = [];            end;
if ~isfield(cfg, 'neighbourdist'),       cfg.neighbourdist = 4;          end;
if ~isfield(cfg, 'geomfile'),            cfg.geomfile = '';              end;
if ~isfield(cfg, 'chancmbgeom'),         cfg.chancmbgeom = 'refchan';    end;
if ~isfield(cfg, 'savegeom'),            cfg.savegeom = 'no';            end;
if ~isfield(cfg, 'randomseed'),          cfg.randomseed = 'yes';         end;
if ~isfield(cfg, 'mirrordata'),          cfg.mirrordata = 'no';          end;
if ~isfield(cfg, 'channel'),             cfg.channel = 'all';            end;
if ~isfield(cfg, 'channelcmb'),          cfg.channelcmb = {'all' 'all'}; end;
if ~isfield(cfg, 'latency'),             cfg.latency = 'all';            end;

% for backward compatibility with old data structures
for i=1:length(varargin)
  varargin{i} = checkdata(varargin{i});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform some checks.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine whether a beween or a within-replications design is requested.
if any(strcmp(cfg.statistic,{'indepsamplesT','indepsamplesregrT','indepsamplesZcoh','indepsamplesTsqcrs','indepsamplesF'}))
    between = 1;
end;
if any(strcmp(cfg.statistic,{'depsamplesregrT','depsamplesT','actvsblT','depsamplesF'}))
    between = 0;
end;
if ~(exist('between')==1)
    error('Unknown test statistic.');
end;

onedimtest = any(strcmp(cfg.statistic,{'indepsamplesT','indepsamplesregrT','indepsamplesZcoh','depsamplesregrT','depsamplesT','actvsblT'}));

Nvarargin = length(varargin);

if any(strcmp(cfg.statistic,{'indepsamplesT','indepsamplesZcoh','indepsamplesTsqcrs','depsamplesT','actvsblT'})) && (Nvarargin~=2)
    error('Number of data sets is not equal to 2, as is required by the test statistic.');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If required by cfg, add the cross- and the autospectrum to cfg.channelcmb.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(strcmp(cfg.statistic,{'indepsamplesZcoh'}))
    if ~isfield(cfg,'channelcmb')
        if isfield(varargin{1},'labelcmb')
            cfg.channelcmb=varargin{1}.labelcmb;
        else
            error('The first dataset does not contain cross-spectra, as is required with this test statistic.');
        end;
    else
        if isfield(varargin{1},'labelcmb')
            cfg.channelcmb=channelcombination(cfg.channelcmb,varargin{1}.label);
        else
            error('The first dataset does not contain cross-spectra, as is required with this test statistic.');
        end;
    end;
    if strcmp(cfg.statistic,'indepsamplesZcoh')
        uniqchanlabels = unique(cfg.channelcmb(:));
        nuniqchanlabels=length(uniqchanlabels);
        autospctrselvec=strcmp(cfg.channelcmb(:,1),cfg.channelcmb(:,2));
        crsspctrselvec=~autospctrselvec;
        autospctrchannelcmb=cat(2,uniqchanlabels,uniqchanlabels);
        cfg.channelcmb=cat(1,autospctrchannelcmb, cfg.channelcmb(crsspctrselvec,:));
    end;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocessing, specific for the activation-versus-baseline statistic.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(cfg.statistic,'actvsblT')
    poststim=find(varargin{1}.time>=0);
    if length(poststim)~=length(varargin{1}.time)
        error('Not all time points in the first data set are in the activation (post-stimulus) period.');
    end;
    prestim=find(varargin{2}.time<=0);
    if length(prestim)~=length(varargin{2}.time)
        error('Not all time points in the second data set are in the baseline (pre-stimulus) period.');
    end;
    if length(prestim)~=length(poststim)
        error('The time axes of the first (activation) and the second (baseline) data set are of different length.');
    end;
    varargin{2}.time=varargin{2}.time + min(varargin{1}.time) - min(varargin{2}.time);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data selection and averaging over selected dimensions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if between
    cfg.datarepresentation = 'concatenated';
else % a within-replication design    
    cfg.datarepresentation = 'cell-array';
end;

fprintf('Selecting and formatting the data.\n');
[cfg,data]=prepare_timefreq_data(cfg, varargin{1:Nvarargin});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Reshape the data in a format suitable for clusterrandstatistics.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if between
    [nrepl,lengthspatialdim,nfreq,ntime]=size(data.biol);
    data.biol=reshape(data.biol,[nrepl,1,lengthspatialdim,nfreq,ntime]);
else   % within-replication conditions
    data.biolcell=data.biol;
    nwcond=length(data.biolcell);
    [firstnrepl,lengthspatialdim,nfreq,ntime]=size(data.biolcell{1});
    data.biol=reshape(data.biolcell{1},[firstnrepl,1,lengthspatialdim,nfreq,ntime]);
    data.biolcell{1}=[];
    for condindx=2:nwcond
        nrepl=size(data.biolcell{condindx},1);
        if nrepl~=firstnrepl
            error('The number of replications in the different conditions are unequal. This is not allowed with a test statistic for a within-replications design.');
        end;
        data.biol=cat(2,data.biol,reshape(data.biolcell{condindx},[nrepl,1,lengthspatialdim,nfreq,ntime]));
        data.biolcell{condindx}=[];
    end;
end;
% add the dimension wcond to the dimord
s = strfind(data.dimord, 'chan');
data.dimord = [data.dimord(1:(s-1)) 'wcond_chan' data.dimord((s+4):end)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adding the external variable data.ext
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if between || strcmp(cfg.statistic,'depsamplesregrT')
    if Nvarargin>1
        data.ext=zeros(Nvarargin,1);
    else % Nvarargin=1
        if strcmp(cfg.statistic,'indepsamplesregrT')
            data.ext=zeros(nrepl,1);
        elseif strcmp(cfg.statistic,'depsamplesregrT')
            data.ext=zeros(nwcond,1);
        end;
    end;
    if strcmp(cfg.statistic,'indepsamplesregrT') | strcmp(cfg.statistic,'depsamplesregrT')
    % constructing the predictor variable for the dependent and the
    % independent samples regression T-statistic
        if isfield(cfg,'ext')
            if size(cfg.ext,1)==1
                data.ext=cfg.ext'; % take the transpose of the row vector
            else
                data.ext=cfg.ext;
            end;
        else
            if Nvarargin>1
                for argindx=1:Nvarargin
                    if isfield(varargin{argindx},'ext')
                        data.ext(argindx)=varargin{argindx}.ext;
                    else
                        error('Neither the configuration nor the data contain the field "ext", as is required by the test statistic.');
                    end;
                end;
            else
                if isfield(varargin{1},'ext')
                    if size(varargin{1}.ext,1)==1
                        data.ext=varargin{1}.ext';
                    else
                        data.ext=varargin{1}.ext;
                    end;
                else
                    error('Neither the configuration nor the data contain the field "ext", as is required by the test statistic.');
                end;
            end;
        end;
    else % if the statistic is not the independent nor the dependent samples regression T-statistic
        data.ext=data.design(:,3);
    end;
    if strcmp(cfg.statistic,'indepsamplesregrT') && (length(data.ext)~=size(data.biol,1))
        error('The number of entries in the predictor variable is not equal to the number of replications in the data.');
    end;
    if strcmp(cfg.statistic,'depsamplesregrT') && (length(data.ext)~=size(data.biol,2))
        error('The number of entries in the predictor variable is not equal to the number of replications in the data.');
    end;
end;    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocessing, specific for channel combination data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(data,'labelcmb')
    data.label=unique(data.labelcmb(:));
    if isfield(cfg,'mirrordata') && strcmp(cfg.mirrordata,'yes')
        reverseddatalabelcmb=cell(size(data.labelcmb,1),2);
        reverseddatalabelcmb=data.labelcmb(:,[2 1]);
        data.labelcmb=[data.labelcmb;reverseddatalabelcmb];
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the neighbourhood geometry.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Calculating the neighbourhood structure of the channels.\n');
[cfg,data] = getneighbgeometry(cfg,data,varargin{1});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run clusterrandstatistics.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[clusrand] = clusterrandstatistics(cfg, data);

warning on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN SUBFUNCTION GETNEIGHBGEOMETRY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cfg,data] = getneighbgeometry(cfg,data,firstdataset);
% get the neighbourhood geometry

singlechans = ~isfield(data,'labelcmb');

if  ~strcmp(cfg.savegeom,'yes') && ~isempty(cfg.geomfile) && exist(cfg.geomfile,'file')
    load(cfg.geomfile); % load the structure geom into memory
end;

if ~exist('geom') && isempty(cfg.neighbours)
    % Add the grad and elec fields.
    swapmemfile;  % initialization
    firstdataset = swapmemfile(firstdataset);
    if isfield(cfg, 'grad')
        fprintf('Obtaining the gradiometer configuration from the configuration.\n');
        data.grad = cfg.grad;
    elseif isfield(cfg, 'elec')
        fprintf('Obtaining the electrode configuration from the configuration.\n');
        data.elec = cfg.elec;
    elseif isfield(cfg, 'gradfile')
        fprintf('Obtaining the gradiometer configuration from a file.\n');
        data.grad = read_sens(cfg.gradfile);
    elseif isfield(cfg, 'elecfile')
        fprintf('Obtaining the electrode configuration from a file.\n');
        data.elec = read_sens(cfg.elecfile);
    elseif isfield(firstdataset, 'grad')
        fprintf('Obtaining the gradiometer configuration from the first dataset.\n');
        data.grad = firstdataset.grad;
    elseif isfield(firstdataset, 'elec')
        fprintf('Obtaining the electrode configuration from the first dataset.\n');
        data.elec = firstdataset.elec;
    end
    firstdataset = swapmemfile(firstdataset);
    if ~(isfield(data,'grad') | isfield(data,'elec'))
        error('Did not find gradiometer or electrode information.');
    end;
    cfg.neighbours=compneighbstructfromgradelec(data,cfg.neighbourdist);
end;

% At this point, either (1) cfg.neighbours contains the neighborhood
% structure, (2) cfg.neighbours does NOT contain the neighborhood
% structure, but geom exist, or (3) cfg.neighbours does NOT contain the neighborhood
% structure and geom does not exist.
    
% Find the non-reference channels if cfg.chancmbgeom='refchan'
if ~singlechans && (strcmp(cfg.chancmbgeom,'refchan') | strcmp(cfg.chancmbgeom,'free'))
    crsspctrselvec = find(~strcmp(data.labelcmb(:,1),data.labelcmb(:,2)));
    crsspctrlabelcmb = data.labelcmb(crsspctrselvec,:);
    if strcmp(cfg.chancmbgeom,'refchan')
        uniquelabels=unique(crsspctrlabelcmb(:));
        nuniqlabels=length(uniquelabels);
        npairsperchan=zeros(nuniqlabels,1);
        for indx=1:nuniqlabels
            npairsperchan(indx)=length(strmatch(uniquelabels{indx},crsspctrlabelcmb(:),'exact'));
        end;
        [sortednpairs,sorti]=sort(npairsperchan,'descend');
        % perform a check on the channel combinations to see if they conform to cfg.chancmbgeom='refchan'
        if sortednpairs(1)~=(nuniqlabels-1) && sortednpairs(2)~=1
            error('The channel combinations in the data do not conform to the constraints imposed by cfg.chancmbgeom=''refchan''.');
        end;
        refchanlabel=uniquelabels(sorti(1));
        nonrefchanlabels=cell(nuniqlabels-1,1);
        nonrefchanlabels(~strcmp(refchanlabel,crsspctrlabelcmb(:,1)))=crsspctrlabelcmb(~strcmp(refchanlabel,crsspctrlabelcmb(:,1)),1);
        nonrefchanlabels(~strcmp(refchanlabel,crsspctrlabelcmb(:,2)))=crsspctrlabelcmb(~strcmp(refchanlabel,crsspctrlabelcmb(:,2)),2);
    end;
end;
           
if exist('geom')
    if singlechans
        [checkvec,selvec]=match_str(data.label,geom.label);
        if length(checkvec)~=length(data.label)
            error('For some channels in the data, the neighbourhood geometry file does not contain information.');
        end;
        data.channeighbstructmat=geom.channeighbstructmat(selvec,selvec);
    else    % channel combination data
        if strcmp(cfg.chancmbgeom,'refchan')
            [checkvec,selvec]=match_str(nonrefchanlabels,geom.label);
            if length(checkvec)~=length(nonrefchanlabels)
                error('For some non-reference channels in the data, the neighbourhood geometry file does not contain information.');
            end;
            data.chancmbneighbstructmat=geom.channeighbstructmat(selvec,selvec);
        elseif strcmp(cfg.chancmbgeom,'free')
            datalabelcmb=cmb2label(data.labelcmb);
            geomlabelcmb=cmb2label(geom.labelcmb);
            [checkvec,selvec]=match_str(datalabelcmb,geomlabelcmb);
            if length(checkvec)~=length(datalabelcmb)
                error('For some channel combinations in the data, the neighbourhood geometry file does not contain information.');
            end;
            data.chancmbneighbstructmat=geom.chancmbneighbstructmat(selvec,selvec);
            if isfield(geom,'chancmbneighbselmat');
                data.chancmbneighbselmat=geom.chancmbneighbselmat(selvec,selvec);
            end;
        end;
    end;
    clear geom;
else  % the neighbourhood geometry structure geom does not exist
    % compute CHANNEIGHBSTRUCTMAT, which will be used in the clustering algorithm
    if singlechans
        data.channeighbstructmat = makechanneighbstructmat(cfg.neighbours,data.label);
        if strcmp(cfg.savegeom,'yes')
            geom.label = data.label;
            geom.channeighbstructmat = data.channeighbstructmat;
        end;
    else  % channel combination data
        if strcmp(cfg.chancmbgeom,'refchan')
            data.chancmbneighbstructmat = makechanneighbstructmat(cfg.neighbours,nonrefchanlabels);
            if strcmp(cfg.savegeom,'yes')
                geom.label = nonrefchanlabels;
                geom.channeighbstructmat = data.chancmbneighbstructmat;
            end;
        elseif strcmp(cfg.chancmbgeom,'free')
            % compute CHANNEIGHBSTRUCTMAT, which will be used in the clustering algorithm
            singlechanlabels = unique(data.labelcmb(:));
            data.channeighbstructmat = makechanneighbstructmat(cfg.neighbours,singlechanlabels);
            % compute CHANCMBNEIGHBSTRUCTMAT and CHANCMBNEIGHBSELMAT, which will be
            % used for clustering channel combinations.
            orderedchancmbs=strcmp(cfg.mirrordata,'yes');
            % orderedchancmbs is true if the result of the statistical test depends 
            % on the order of the elements of the cross-spectrum (as, for
            % example, when these are coherencies or the imaginary
            % parts of the coherencies). 
            original = true;
            % with original=false, we use a calculation that consumes less memory.
            [data.chancmbneighbstructmat data.chancmbneighbselmat] = makechancmbneighbmats(data.channeighbstructmat, ...
                crsspctrlabelcmb,singlechanlabels,orderedchancmbs,original);
            if strcmp(cfg.savegeom,'yes')
                geom.labelcmb = crsspctrlabelcmb;
                geom.chancmbneighbstructmat = data.chancmbneighbstructmat;
                geom.chancmbneighbselmat = data.chancmbneighbselmat;
            end;
        else
            error('Unsupported value for cfg.chancmbgeom.');
        end;
    end;
    if strcmp(cfg.savegeom,'yes') && ~isempty(cfg.geomfile)
        save(cfg.geomfile,'geom');
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN SUBFUNCTION COMPNEIGHBSTRUCTFROMGRADELEC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [neighbours]=compneighbstructfromgradelec(data,neighbourdist);
    
% compute the neighbourhood geometry from the gradiometer/electrode positions provided in the data    

if isfield(data, 'grad')
    sens = data.grad;
elseif isfield(data, 'elec')
    sens = data.elec;
else
    error('No gradiometer or electrode configuration present in the data');
end
nsensors = length(sens.label);

% compute the distance between all sensors
dist = zeros(nsensors,nsensors);
for i=1:nsensors
    dist(i,:) = sqrt(sum((sens.pnt(1:nsensors,:) - repmat(sens.pnt(i,:), nsensors, 1)).^2,2))';
end;

% find the neighbouring electrodes based on distance, later we have to restrict the neighbouring 
% electrodes to those actually selected in the dataset
channeighbstructmat = (dist<neighbourdist);

% electrode istelf is not a neighbour
channeighbstructmat=channeighbstructmat.*~eye(nsensors);

% construct a structured cell array with all neighbours
neighbours=cell(1,nsensors);
for i=1:nsensors
    neighbours{i}.label       = sens.label{i};
    neighbours{i}.neighblabel = sens.label(find(channeighbstructmat(i,:)));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN SUBFUNCTION MAKECHANNEIGHBSTRUCTMAT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [channeighbstructmat] = makechanneighbstructmat(neighbours,label);

% MAKECHANNEIGHBSTRUCTMAT makes the makes the matrix containing the channel
% neighbourhood structure.

nchan=length(label);
channeighbstructmat = logical(zeros(nchan,nchan));
if length(neighbours)==0
    error('Empty neighborhood structure ''cfg.neighbours''.');
end;
for chan=1:length(neighbours)
    [seld] = match_str(label, neighbours{chan}.label);
    [seln] = match_str(label, neighbours{chan}.neighblabel);
    if isempty(seld)
        % this channel was not present in the data
        continue;
    else
        % add the neighbours of this channel to the matrix
        channeighbstructmat(seld, seln) = 1;
    end
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN SUBFUNCTION MAKECHANCMBNEIGHBMATS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [chancmbneighbstructmat, chancmbneighbselmat] = makechancmbneighbmats(channeighbstructmat,labelcmb,label,orderedchancmbs,original);

% compute CHANCMBNEIGHBSTRUCTMAT and CHANCMBNEIGHBSELMAT, which will be
% used for clustering channel combinations.

nchan=length(label);
nchansqr=nchan^2;

% Construct an array labelcmbnrs from labelcmb by replacing the label pairs
% in labelcmb by numbers that correspond to the order of the labels in label. 
nchancmb = size(labelcmb,1);
labelcmbnrs=zeros(nchancmb,2);
for chanindx=1:nchan
    [chansel1] = match_str(labelcmb(:,1),label(chanindx));
    labelcmbnrs(chansel1,1)=chanindx;
    [chansel2] = match_str(labelcmb(:,2),label(chanindx));
    labelcmbnrs(chansel2,2)=chanindx;
end;
% Calculate the row and column indices (which are identical) of the
% channel combinations that are present in the data.
chancmbindcs=zeros(nchancmb,1);
for indx=1:nchancmb
    chancmbindcs(indx)=(labelcmbnrs(indx,1)-1)*nchan + labelcmbnrs(indx,2);
end;

% put all elements on the diagonal of CHANNEIGHBSTRUCTMAT equal to one
channeighbstructmat=channeighbstructmat | logical(eye(nchan));

% Compute CHANCMBNEIGHBSTRUCTMAT
% First compute the complete CHANCMBNEIGHBSTRUCTMAT (containing all 
% ORDERED channel combinations that can be formed with all channels present in 
% the data) and later select and reorder the channel combinations actually
% present in data). In the complete CHANCMBNEIGHBSTRUCTMAT, the row and the
% column pairs are ordered lexicographically.

chancmbneighbstructmat = false(nchansqr);
if original
	% two channel pairs are neighbours if their first elements are
	% neighbours
	chancmbneighbstructmat = logical(kron(channeighbstructmat,ones(nchan)));
	
	% or if their second elements are neighbours
	chancmbneighbstructmat = chancmbneighbstructmat | logical(kron(ones(nchan),channeighbstructmat));
else    %  version that consumes less memory
	for chanindx=1:nchan
        % two channel pairs are neighbours if their first elements are neighbours
        % or if their second elements are neighbours
        chancmbneighbstructmat(:,((chanindx-1)*nchan + 1):(chanindx*nchan)) = ...
            logical(kron(channeighbstructmat(:,chanindx),ones(nchan))) | ...
            logical(kron(ones(nchan,1),channeighbstructmat));
	end;
end;

if ~orderedchancmbs
    if original
		% or if the first element of the row channel pair is a neighbour of the
		% second element of the column channel pair
		chancmbneighbstructmat = chancmbneighbstructmat | logical(repmat(kron(channeighbstructmat,ones(nchan,1)), [1 nchan]));
		
		% or if the first element of the column channel pair is a neighbour of the
		% second element of the row channel pair
		chancmbneighbstructmat = chancmbneighbstructmat | logical(repmat(kron(channeighbstructmat,ones(1,nchan)),[nchan 1]));
    else
		for chanindx=1:nchan
            % two channel pairs are neighbours if 
            % the first element of the row channel pair is a neighbour of the
            % second element of the column channel pair
            % or if the first element of the column channel pair is a neighbour of the
            % second element of the row channel pair
            chancmbneighbstructmat(:,((chanindx-1)*nchan + 1):(chanindx*nchan)) = ...
                chancmbneighbstructmat(:,((chanindx-1)*nchan + 1):(chanindx*nchan)) | ... 
                logical(kron(channeighbstructmat,ones(nchan,1))) | ...
                logical(repmat(kron(channeighbstructmat(:,chanindx),ones(1,nchan)),[nchan 1]));
		end;
    end;
end;

% reorder and select the entries in chancmbneighbstructmat such that they correspond to labelcmb.
chancmbneighbstructmat = sparse(chancmbneighbstructmat);
chancmbneighbstructmat = chancmbneighbstructmat(chancmbindcs,chancmbindcs);

% compute CHANCMBNEIGHBSELMAT
% CHANCMBNEIGHBSELMAT identifies so-called parallel pairs. A channel pair
% is parallel if (a) all four sensors are different and (b) all elements (of the first 
% and the second pair) are neighbours of an element of the other pair. 
% if orderedpairs is true, then condition (b) is as follows: the first
% elements of the two pairs are neighbours, and the second elements of the 
% two pairs are neighbours.

% put all elements on the diagonal of CHANNEIGHBSTRUCTMAT equal to zero
channeighbstructmat = logical(channeighbstructmat.*(ones(nchan)-diag(ones(nchan,1))));

chancmbneighbselmat = false(nchansqr);
if orderedchancmbs
    if original
		% if the first element of the row pair is a neighbour of the
		% first element of the column pair
		chancmbneighbselmat = logical(kron(channeighbstructmat,ones(nchan)));
		% and the second element of the row pair is a neighbour of the
		% second element of the column pair
		chancmbneighbselmat = chancmbneighbselmat & logical(kron(ones(nchan),channeighbstructmat));
    else    %  version that consumes less memory 
		for chanindx=1:nchan
			% if the first element of the row pair is a neighbour of the
			% first element of the column pair
			% and the second element of the row pair is a neighbour of the
			% second element of the column pair
            chancmbneighbselmat(:,((chanindx-1)*nchan + 1):(chanindx*nchan)) = ...
				logical(kron(channeighbstructmat(:,chanindx),ones(nchan))) & logical(kron(ones(nchan,1),channeighbstructmat));
		end;
    end;
else  % unordered channel combinations
    if original
        % if the first element of the row pair is a neighbour of one of the
        % two elements of the column pair
        chancmbneighbselmat = logical(kron(channeighbstructmat,ones(nchan))) | logical(repmat(kron(channeighbstructmat,ones(nchan,1)), [1 nchan]));
		% and the second element of the row pair is a neighbour of one of the
		% two elements of the column pair
		chancmbneighbselmat = chancmbneighbselmat & (logical(kron(ones(nchan),channeighbstructmat)) | ... 
            logical(repmat(kron(channeighbstructmat,ones(1,nchan)), [nchan 1])));
    else    %  version that consumes less memory 
		for chanindx=1:nchan
			% if the first element of the row pair is a neighbour of one of the
			% two elements of the column pair
			% and the second element of the row pair is a neighbour of one of the
			% two elements of the column pair
            chancmbneighbselmat(:,((chanindx-1)*nchan + 1):(chanindx*nchan)) = ...
			(logical(kron(channeighbstructmat(:,chanindx),ones(nchan))) | logical(kron(channeighbstructmat,ones(nchan,1)))) ...
			& (logical(kron(ones(nchan,1),channeighbstructmat)) | ...
                logical(repmat(kron(channeighbstructmat(:,chanindx),ones(1,nchan)), [nchan 1])));
		end;
    end;
end;

if original
	% remove all pairs of channel combinations that have one channel in common.
	% common channel in the first elements of the two pairs.
	chancmbneighbselmat = chancmbneighbselmat & kron(~eye(nchan),ones(nchan));
	
	% common channel in the second elements of the two pairs.
	chancmbneighbselmat = chancmbneighbselmat & kron(ones(nchan),~eye(nchan));
	
	% common channel in the first element of the row pair and the second
	% element of the column pair.
	tempselmat=logical(zeros(nchansqr));
	tempselmat(:)= ~repmat([repmat([ones(nchan,1); zeros(nchansqr,1)],[(nchan-1) 1]); ones(nchan,1)],[nchan 1]);
	chancmbneighbselmat = chancmbneighbselmat & tempselmat;
	
	% common channel in the second element of the row pair and the first 
	% element of the column pair.
	chancmbneighbselmat = chancmbneighbselmat & tempselmat';
else
	noteye=~eye(nchan);
	tempselmat=logical(zeros(nchansqr,nchan));
	tempselmat(:)= ~[repmat([ones(nchan,1); zeros(nchansqr,1)],[(nchan-1) 1]); ones(nchan,1)];
	for chanindx=1:nchan
		% remove all pairs of channel combinations that have one channel in common.
		% common channel in the first elements of the two pairs.
		% common channel in the second elements of the two pairs.
		% common channel in the first element of the row pair and the second
		% element of the column pair.
        chancmbneighbselmat(:,((chanindx-1)*nchan + 1):(chanindx*nchan)) = ...
            chancmbneighbselmat(:,((chanindx-1)*nchan + 1):(chanindx*nchan)) ...
			& logical(kron(noteye(:,chanindx),ones(nchan))) ...
			& logical(kron(ones(nchan,1),noteye)) ...
            & tempselmat;
	end;
	for chanindx=1:nchan
		% remove all pairs of channel combinations that have one 
        % common channel in the second element of the row pair and the first 
		% element of the column pair.
        chancmbneighbselmat(((chanindx-1)*nchan + 1):(chanindx*nchan),:) = ...
            chancmbneighbselmat(((chanindx-1)*nchan + 1):(chanindx*nchan),:) & tempselmat';
	end;
end;

% reorder and select the entries in chancmbneighbselmat such that they correspond to labelcmb.
chancmbneighbselmat = sparse(chancmbneighbselmat);
chancmbneighbselmat = chancmbneighbselmat(chancmbindcs,chancmbindcs);

% put all elements below and on the diagonal equal to zero
nchancmbindcs = length(chancmbindcs);
for chancmbindx=1:nchancmbindcs
    selvec=[true(chancmbindx-1,1);false(nchancmbindcs-chancmbindx+1,1)];
    chancmbneighbstructmat(:,chancmbindx) = chancmbneighbstructmat(:,chancmbindx) & selvec;
    chancmbneighbselmat(:,chancmbindx) = chancmbneighbselmat(:,chancmbindx) & selvec;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION convert Nx2 cell array with channel combinations into Nx1 cell array with labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function label = cmb2label(labelcmb);
label = {};
for i=1:size(labelcmb)
  label{i,1} = [labelcmb{i,1} '&' labelcmb{i,2}];
end
