function [clusrand] = clusterrandstatistics(cfg,data);

% CLUSTERRANDSTATISTICS
%
% Use as
%   [clusrand] = clusterrandstatistics(cfg, data)
%
% The argument data is of the statsdata type. The statsdata type is a
% struct containing the required field biol and several optional fields: 
% label and labelcmb (one of these two has to be present), ext, 
% grad (for MEG), and elec (for EEG). The field biol contains a five-dimensional 
% array with dimord repl_wcond_chan(cmb)_freq_time containing the 
% biological data. If some dimension is not present in the data, it has to be represented as a
% singleton dimension in data.biol. For example, in a study with 120
% replications, 1 within-replication condition, 151 channels, no frequecies (time domain 
% data), and 600 time samples, the dimensionality of data.biol is [120x1x151x1x600]. 
%
% The the variable data contains a field label, containing the labels of the 
% channels, then the spatial dimension of data.biol (the third one) 
% corresponds to single channels (and not to channel combinations). 
% This information is used in combination with the 
% information in data.grad, data.elec, or cfg.neighbours to determine the 
% neighbourhood structure of the channels.
%
% If the variable data contains a field labelcmb, containing the labels of the 
% channel combinations, then the spatial dimension of data.biol (the third
% one) corresponds to channel combinations (and not to single channels). 
% The field data.labelcmb contains a cell array of size [nchancmb,2], with 
% nchancmb being the number of channel combinations. From data.labelcmb,
% this function computes data.label, containing the labels of the single
% channels.
%
% The optional field data.ext contains a one-dimensional array with length 
% equal to the number of replications. Statistical tests for regression and
% comparisons involving between-replication conditions (e.g., independent 
% samples T and F) require that data.ext specifies the independent variable
% with integers between 1 and the number of levels of the 
% independent variable. 
%
% The optional fields data.grad (for MEG) and data.elec (for EEG) contain the 
% gradiometer/electrode positions. The structure of these fields must be
% identical to the grad and elec fields produced by PREPROCESSING. 
%
% cfg.statistic      = 'indepsamplesT' (independent samples T-statistic) 
%                      | 'indepsamplesregrT' (independent samples regression coefficient T-statistic)
%                      | 'indepsamplesZcoh' (independent samples Z-statistic for coherence)
%                      | 'indepsamplesTsqcrs' (independent samples T-square-statistic for the cross-spectrum)
%                      | 'depsamplesregrT' (dependent samples regression coefficient T-statistic)
%                      | 'depsamplesT' (dependent samples T-statistic)
%                      | 'actvsblT' (activation versus baseline T-statistic)
%                      | 'indepsamplesF' (independent samples F-statistic)
%                      | 'depsamplesF' (dependent samples F-statistic)
%
% For every statistic, we now list (a) the constraints with respect to the 
% data structure and (b) a number of fields that are relevant for a particular choice
% of statistic. 
%
% For the independent samples T-statistic:
% When this statistic is chosen, data.ext must be column array containing
% the value 1 for replications that belong to the first condition and
% the value 2 for replications that belong to the second condition.
% cfg.onetwo         = 'onesided_1<2' (H_0:condition_1<condition_2), 
%                      'onesided_2<1' (H_0:condition_2<condition_1), or
%                      'twosided' (H_0:condition_1=condition_2, the
%                      default).
% 
% For the independent samples regression coefficient T-statistic:
% When this statistic is chosen, data.ext must be a column array containing
% the predictor variable.
% cfg.onetwo         = 'onesided_B<0' (H_0: regr. coefficient B<0), 
%                      'onesided_B>0' (H_0: regr. coefficient B>0), or
%                      'twosided' (H_0:regr. coefficient B=0, the default).
%
% For the independent samples Z-statistic for coherence:
% When this statistic is chosen, the spatial dimension of the data matrix must 
% contain (a) the cross-spectrum for a number of channel combinations, and
% (b) the auto-spectrum (power) for the channels that are involved in the
% channel combinations. With this test statistic, one examines whether the
% two condition differ with respect to coherence. 
% When this statistic is chosen, data.ext must be column vector containing
% the value 1 for replications that belong to the first condition and
% the value 2 for replications that belong to the second condition.
% cfg.onetwo         = 'onesided_1<2' (H_0:condition_1<condition_2), 
%                      'onesided_2<1' (H_0:condition_2<condition_1), or
%                      'twosided' (H_0:condition_1=condition_2, the
%                      default).
%
% For the independent samples T-square-statistic for the cross-spectrum:
% When this statistic is chosen, the data matrix must contain the
% cross-spectrum for a number of channel combinations. 
% With this test statistic, one examines whether the two conditions differ 
% with respect to the average cross-spectrum.
% When this statistic is chosen, data.ext must be column vector containing
% the value 1 for replications that belong to the first condition and
% the value 2 for replications that belong to the second condition.
%
% For the dependent samples T-statistic:
% When this statistic is chosen, the second dimension of data.biol (wcond)
% must have two levels. The first level contains the data of the first condition 
% and the second one the data of the second condition.
% cfg.onetwo         = 'onesided_1<2' (H_0:condition_1<condition_2), 
%                      'onesided_2<1' (H_0:condition_2<condition_1), or
%                      'twosided' (H_0:condition_1=condition_2, the
%                      default).
%
% For the activation versus baseline T-statistic:
% When this statistic is chosen, the second dimension of data.biol (wcond)
% must have two levels. The first level contains the data of the activation 
% period and the second level contains the data of the baseline. 
% cfg.onetwo         = 'onesided_act<bl' (H_0:act<bl), 'onesided_bl<act' (H_0:bl<act), or
%                      'twosided' (H_0:act=bl, the default).
%
% For the dependent samples regression coefficient T-statistic:
% When this statistic is chosen, data.ext must be a column array containing
% the predictor variable.
% cfg.onetwo         = 'onesided_B<0' (H_0: regr. coefficient B<0), 
%                      'onesided_B>0' (H_0: regr. coefficient B>0), or
%                      'twosided' (H_0:regr. coefficient B=0, the default).
%
% For the independent samples F-statistic:
% When this statistic is chosen, data.ext must be column array containing
% values between 1 and the number of conditions (ncond): 1 for the first condition, 
% 2 for the second, etc. 
% 
% For the dependent samples F-statistic:
% When this statistic is chosen, the second dimension of data.biol (wcond)
% must have as many levels as the number of conditions (ncond). The first level 
% contains the data of the first condition, the second one the data of the second 
% condition, etc.
% cfg.contrastcoefs  = matrix of contrast coefficients determining the
%                      effect being tested. The number of columns of this
%                      matrix has to be equal to the number of conditions. 
%                      The default is a matrix that specifies the
%                      main effect of the independent variable. This matrix
%                      has size [(ncond-1),ncond]. 
% cfg.allowedperms   = matrix of permutations of conditions that all have the
%                      same probability (given the data) under the null 
%                      hypothesis. The default is all ncond! possible permutations.
%
% We now list a number of fields that determine the way clustering is
% performed:
% cfg.alphathresh    = alpha level of the (channel(comb),frequency,timepoint)-specific 
%                      test statistic that will be used for thresholding 
% cfg.minnbchan      = minimum number of neighbouring channel(s) (combinations) in which a 
%                      particular time-frequency-specific t-test has to be significant in order for it 
%                      to be included in the clustering algorithm (default=0).
% cfg.alpha          = alpha level of the randomization test
% cfg.makeclusters = 'yes' (default) | 'no'
% cfg.clusterteststat = 'maxsum' (default) or 'orderedsums'
% cfg.smallestcluster= smallest cluster size that is large enough to be considered (only
%                      relevant when clusterteststat='orderedsums')
%
% cfg.nranddraws     = number of draws from the randomization distribution
% cfg.randomseed     = 'yes' (default), 'no' or seed-number

% STATISTICAL DESCRIPTION
% This function approximates (via the Monte Carlo method) the
% randomization distribution of the largest cluster level statistic (the n
% largest if cfg.clusterteststat=='orderedsums'). A cluster level statistic 
% is the sum of the statistics within a connected cluster
% of (channel,freq,time)-points. These connected clusters are found by
% the FINDCLUSTER function. The randomization distribution is for the 
% largest positive cluster level statistic. For t-statistics, because of symmetry, 
% this distribution is identical to the randomization distribution of minus the 
% smallest negative cluster level statistic.
% 
% The (channel,freq,time)-specific statistics are thresholded at cfg.alphathresh, 
% and the cluster level critical value(s) are computed such that, under the randomization
% distribution, with probability cfg.alpha, the observed cluster level statistics 
% are larger than these critical values. For t-statistics and cfg.onetwo set at 'twosided', 
% the positive (channel,freq,time)-specific t-statistics are thresholded at 
% cfg.alphathresh/2, and the cluster level critical value(s) are computed such that, 
% under the randomization distribution, with probability cfg.alpha, the observed cluster level
% statistics are more extreme than these critical values. 

% Copyright (C) 2005, Eric Maris, NICI, University of Nijmegen
%
% $Log: clusterrandstatistics.m,v $
% Revision 1.14  2007/09/23 14:12:44  erimar
% Correct a few bugs related to the collection of the results
% with respect to one-sided testing.
%
% Revision 1.13  2006/04/06 13:07:37  erimar
% Added ";" to suppress screen output.
%
% Revision 1.12  2006/03/10 08:22:41  erimar
% Corrected incorrect handling of mirror data (due to treating
% cfg.mirrordata as a boolean).
%
% Revision 1.11  2006/02/28 12:04:31  erimar
% Added the test statistic indepsamplesZcoh (for testing coherence differences).
%
% Revision 1.10  2005/12/06 13:20:58  erimar
% Changed handling of external predictor variables, consistent with the
% updated version of clusterrandanalysis.
%
% Revision 1.9  2005/11/28 11:21:53  erimar
% Corrected an error in the code that is specific for
% clusterteststat='orderedsums'.
%
% Revision 1.8  2005/08/05 12:20:28  erimar
% Correct and error in the calculation of the independent samples
% regression T-statistic detected by Vladimir Litvak.
%
% Revision 1.7  2005/07/07 14:38:29  erimar
% Replaced .*-operator by logical AND (&) operator for logical variables.
%
% Revision 1.6  2005/04/22 07:44:41  roboos
% added/corrected copyrights, added a Log tag for CVS, converted to unix ascii
%
% Revision 1.5  2005/04/20 11:59:01  erimar
% Added the dependent samples T-statistic for regression coefficients. Added
% functionality for mirror data. This functionality is essential for proper
% clustering of channel combination data (coherence and coherency) in the
% case the data consists of only half of the sensor combinations. Removed
% the old handling of NaNs (i.c., sum(data.biol(data.notnan))) and replaced
% it by the Matlab function nan_sum. Corrected errors in the calculation of
% the activation-versus-baseline T-statistic.
%
% Revision 1.1  2004/11/04 15:27:09  erimar
% Statistics engine. Used to be a part of the first version of clusterrandanalysis (not in CVS).
%

% Turn divideByZero warnings off.
warning off;

fprintf('Running the statistics engine.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the default configuration parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(cfg, 'onetwo'),             cfg.onetwo = 'twosided';                 end  
if ~isfield(cfg, 'alphathresh'),        cfg.alphathresh = 0.05;                  end
if ~isfield(cfg, 'minnbchan'),          cfg.minnbchan = 0;                       end
if ~isfield(cfg, 'alpha'),              cfg.alpha = 0.05;                        end  
if ~isfield(cfg, 'makeclusters'),       cfg.makeclusters = 'yes';                end
if ~isfield(cfg, 'clusterteststat'),    cfg.clusterteststat = 'maxsum';          end
if ~isfield(cfg, 'smallestcluster'),    cfg.smallestcluster = 0;                 end
if ~isfield(cfg, 'nranddraws'),         cfg.nranddraws = 100;                    end
if ~isfield(cfg, 'randomseed'),         cfg.randomseed = 'yes';                  end
if ~isfield(cfg, 'mirrordata'),         cfg.mirrordata = 'no';                   end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform some checks on data and cfg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(data, 'biol')
    error('You must specify the field data.biol containing the biological data.');
end;
if ~isfield(data, 'label') & ~isfield(data, 'labelcmb') 
    error('You must specify the field data.label or data.labelcmb containing the labels of the channels or channel combinations in data.biol.');
end;
if ~isfield(cfg, 'statistic')
    error('You must specify a statistic (in cfg.statistic).');
end;
if any(strcmp(cfg.statistic,{'indepsamplesT','indepsamplesZcoh','indepsamplesTsqcrs','indepsamplesF','indepsamplesregrT'}))
    if ~isfield(data, 'ext')
        error('You must specify an array that contains the independent variable (in data.ext).');
    end;
end;    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize some variables (not statistic-specific)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(strcmp(cfg.statistic,{'indepsamplesT','indepsamplesZcoh','depsamplesT','actvsblT','indepsamplesregrT','depsamplesregrT'}))
    onedimtest=1;
    if any(strcmp(cfg.onetwo,{'twosided' 'onesided_2<1' 'onesided_bl<act' 'onesided_B>0'}))
        critnegativetail = 1;
    else
        critnegativetail = 0;
    end;
else
    onedimtest=0;
    critnegativetail = 0;
end;

% make the channel combination neighbourhood structure matrix full
if isfield(data, 'chancmbneighbstructmat')
    data.chancmbneighbstructmat=full(data.chancmbneighbstructmat);
end;
% make the channel combination neighbourhood selection matrix full
if isfield(data, 'chancmbneighbselmat')
    data.chancmbneighbselmat=full(data.chancmbneighbselmat);
end;

% determine the length of the dimensions
nrepl=size(data.biol,1);
nwcond=size(data.biol,2);
nfreq=size(data.biol,4);
ntime=size(data.biol,5);
if isfield(data,'labelcmb')
    nchancmb=size(data.biol,3);
    lengthspatdim=nchancmb;
    singlechans=false;
else
    nchan=size(data.biol,3);
    lengthspatdim=nchan;
    singlechans=true;
end;

% remove the wcond dimension if it is a singleton
if nwcond==1
    data.biol=reshape(data.biol,[nrepl lengthspatdim nfreq ntime]);
end;

% initialize the random number generator.
if strcmp(cfg.randomseed, 'no')
   % do nothing
elseif strcmp(cfg.randomseed, 'yes')
   rand('state',sum(100*clock));
else
   % seed with the user-given value
   rand('state',cfg.randomseed);
end;

containsnans=any(isnan(data.biol(:)));
if containsnans
    data.notnan=~isnan(data.biol);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start with performing the statistic-specific preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Statistic-specific preprocessing (calculating critical values, initializing draws from the randomization distribution, ...).\n');

if any(strcmp(cfg.statistic,{'indepsamplesT','indepsamplesZcoh','indepsamplesTsqcrs','indepsamplesF','indepsamplesregrT'}))
    if containsnans
        nobs=reshape(sum(data.notnan,1),lengthspatdim,nfreq,ntime);
    else
        nobs=nrepl;
    end;

    if strcmp(cfg.statistic,'indepsamplesT')
        % preprocessing, specific for the independent samples T-statistic.
        % compute the critical value for the (channel,frequency,timepoint)-specific t-tests
        ncond = 2;
        df = nobs - ncond;
        if any(strcmp(cfg.onetwo,{'onesided_1<2','onesided_2<1'})) 
            crit = tinv(1-cfg.alphathresh,df);
        else
            crit = tinv(1-cfg.alphathresh/2,df);
        end;
        if strcmp(cfg.mirrordata,'yes') & containsnans
            crit = [crit;crit];
        end;
    end;
    
    if strcmp(cfg.statistic,'indepsamplesZcoh')
        % preprocessing, specific for the independent samples Z-statistic
        % for coherence.
        % compute the critical value for the (channel,frequency,timepoint)-specific t-tests
        ncond = 2;
        if any(strcmp(cfg.onetwo,{'onesided_1<2','onesided_2<1'})) 
            crit = norminv(1-cfg.alphathresh,0,1);
        else
            crit = norminv(1-cfg.alphathresh/2,0,1);
        end;
    end;
    
    if strcmp(cfg.statistic,'indepsamplesTsqcrs')
        % preprocessing, specific for the independent samples T-square-statistic for the cross-spectrum.
        % compute the critical value for the (channel,frequency,timepoint)-specific F-tests
        dfnum = 2;
        dfdenom = nobs - 3;
        crit = (nobs-2)*dfnum./dfdenom.*finv(1-cfg.alphathresh,dfnum,dfdenom);
    end;
 
    if any(strcmp(cfg.statistic,{'indepsamplesT','indepsamplesZcoh','indepsamplesTsqcrs'}))
        % construct the matrix randselect, containing the specification of the
        % nranddraws different draws according to the randomization scheme.
        randselect=zeros((cfg.nranddraws+1),nrepl);
        initrandselect=zeros(1,nrepl);
        replselc1 = find(data.ext==1);
        initrandselect(replselc1)=1;
        replselc2 = find(data.ext==2);
        initrandselect(replselc2)=2;
        randselect(1,:)=initrandselect;
        nreplc1=length(replselc1);
        nreplc2=length(replselc2);
        for i=2:(cfg.nranddraws+1)
            permreplsel=randperm(nrepl);
            randselect(i,permreplsel(1:nreplc1))=1;
            randselect(i,permreplsel((nreplc1+1):nrepl))=2;
        end;
    end;
    
    if strcmp(cfg.statistic,'indepsamplesT')
        % compute statistics that can be reused in the different random draws
        sumsq = reshape(nan_sum(data.biol.^2,1),lengthspatdim,nfreq,ntime);
        % the replication and the within-replication-conditions dimension dropped out
    end;

    if strcmp(cfg.statistic,'indepsamplesZcoh')
        % compute statistics that can be reused in the different random draws
        % determine the relation between the cross- and the autospectra
        powindcs=find(strcmp(data.labelcmb(:,1),data.labelcmb(:,2)));
        tmpbool=true(size(data.labelcmb,1),1);
        tmpbool(powindcs)=false;
        crspctrindcs=find(tmpbool);
        powlabels=data.labelcmb(powindcs,1);
        ncrspctr=length(crspctrindcs);
        powindcsforcrspctr=zeros(ncrspctr,2);
        for indx=1:ncrspctr
            powindcsforcrspctr(indx,1)=powindcs(strmatch(data.labelcmb{crspctrindcs(indx),1},powlabels,'exact'));
            powindcsforcrspctr(indx,2)=powindcs(strmatch(data.labelcmb{crspctrindcs(indx),2},powlabels,'exact'));
        end;
    end;
    
    if strcmp(cfg.statistic,'indepsamplesTsqcrs')
        % compute statistics that can be reused in the different random draws
        sumsqcp=zeros(2,2,lengthspatdim,nfreq,ntime);
        sumsqcp(1,1,:,:,:) = nan_sum(real(data.biol).^2,1); 
        sumsqcp(2,2,:,:,:) = nan_sum(imag(data.biol).^2,1);
        sumsqcp(1,2,:,:,:) = nan_sum(real(data.biol).*imag(data.biol),1);
        sumsqcp(2,1,:,:,:) = sumsqcp(1,2,:,:,:);
    end;
    
    if strcmp(cfg.statistic,'indepsamplesregrT')
        % preprocessing, specific for the regression coefficient T-statistic.
        % compute the critical value for the (channel,frequency,timepoint)-specific t-tests
        df = nobs - 2;
        if any(strcmp(cfg.onetwo,{'B<0','B>0'})) 
            crit = tinv(1-cfg.alphathresh,df);
        else
            crit = tinv(1-cfg.alphathresh/2,df);
        end;
  
        % construct the matrix randselect, containing the specification of the
        % nranddraws different draws according to the randomization scheme.
        randselect=zeros((cfg.nranddraws+1),nrepl);
        initrandselect=zeros(1,nrepl);
    
        randselect(1,:)=1:nrepl;
        for i=2:(cfg.nranddraws+1)
            randselect(i,:)=randperm(nrepl);
        end;
        
        % compute statistics that can be reused in the different random draws
        meanbiol = reshape(nan_sum(data.biol,1),lengthspatdim,nfreq,ntime)./nobs;
        expandext = repmat(data.ext,[1 lengthspatdim nfreq ntime]);
        if containsnans
            meanext = reshape(sum(data.notnan.*expandext, 1),lengthspatdim,nfreq,ntime)./nobs;
            varext = reshape(sum((data.notnan.*expandext).^2, 1),lengthspatdim,nfreq,ntime)./nobs - meanext.^2;
        else
            meanext = reshape(sum(expandext, 1),lengthspatdim,nfreq,ntime)./nobs;
            varext = reshape(sum((expandext).^2, 1),lengthspatdim,nfreq,ntime)./nobs - meanext.^2;
        end;
        sdext = sqrt(varext);
        sqrtnobs = sqrt(nobs);
        % the replication and the within-replication-conditions dimension dropped out
    end;    

    if strcmp(cfg.statistic,'indepsamplesF')
        % preprocessing, specific for the independent samples F-statistic.
        % compute the critical value for the (sensor,frequency,timepoint)-specific F-tests
        
        ncond = length(unique(data.ext));
        dfnum = ncond - 1;
        dfdenom = nobs - ncond;
        crit = finv(1-cfg.alphathresh,dfnum,dfdenom);
        
        % construct the matrix randselect, containing the specification of the
        % nranddraws different draws according to the randomization scheme.
        randselect=zeros((cfg.nranddraws+1),nrepl);
        initrandselect=zeros(1,nrepl);
        replsel=cell(1,ncond);
        for condindx=1:ncond
            replsel{condindx} = find(data.ext==condindx);
            initrandselect(replsel{condindx}) = condindx;
        end;
        randselect(1,:)=initrandselect;
        for i=2:(cfg.nranddraws+1)
            randselect(i,:)=zeros(1,nrepl);
            permreplsel=randperm(nrepl);
            shift=0;
            for condindx=1:ncond
                randselect(i,permreplsel((shift+1):(shift+length(replsel{condindx}))))=condindx;
                shift=shift+length(replsel{condindx});
            end;
        end;
        
        % compute statistics that can be reused in the different random draws
        sumsq = reshape(nan_sum(data.biol.^2,1),lengthspatdim,nfreq,ntime);
        overallmean = reshape(nan_sum(data.biol,1),lengthspatdim,nfreq,ntime)./nobs;
        % the replication and the within-replication-conditions dimension dropped out
    end;
    
end;

if any(strcmp(cfg.statistic,{'depsamplesT','depsamplesF','actvsblT','depsamplesregrT'}))
% preprocessing, specific for the dependent samples T- and F-statistic.
    if containsnans
    % guarantee that notnan(:,1,:,:,:) is identical to notnan(:,2,:,:,:), notnan(:,3,:,:,:), etc.
        mask=data.notnan(:,1,:,:,:);
        for condindx=2:nwcond
            mask=mask & data.notnan(:,condindx,:,:,:);
        end;
        for condindx=1:nwcond
            data.notnan(:,condindx,:,:,:)=mask & data.notnan(:,condindx,:,:,:);
        end;
        data.biol(~data.notnan)=NaN;
        nobs=reshape(sum(data.notnan(:,1,:,:,:),1),lengthspatdim,nfreq,ntime);
        % note that nobs is identical for all within-condition levels.
    else
        nobs=nrepl;
    end;

    if strcmp(cfg.statistic,'depsamplesT') | strcmp(cfg.statistic,'actvsblT')
        ncond=2;
        % compute the critical value for the (channel,frequency,timepoint)-specific t-tests
        df = nobs - 1;
        if sum(strcmp(cfg.onetwo,{'onesided_1<2' 'onesided_2<1' 'onesided_act<bl' 'onesided_bl<act'}))
            crit = tinv(1-cfg.alphathresh,df);
        else
            crit = tinv(1-cfg.alphathresh/2,df);
        end;
    end;

    if strcmp(cfg.statistic,'depsamplesregrT')
        ncond=nwcond;
        % compute the critical value for the (channel,frequency,timepoint)-specific t-tests
        df = nobs - 1;
        if sum(strcmp(cfg.onetwo,{'B<0','B>0'})) 
            crit = tinv(1-cfg.alphathresh,df);
        else
            crit = tinv(1-cfg.alphathresh/2,df);
        end;
    end;
    
    if strcmp(cfg.statistic,'depsamplesF')
        if ~isfield(cfg, 'contrastcoefs')
            cfg.contrastcoefs=[];
        end;
        if ~isfield(cfg, 'allowedperms')
            cfg.allowedperms=[];
        end;
        ncond=nwcond;
        if isempty(cfg.contrastcoefs)   % specify the default contrast coefficient matrix.
            ncontrasts = ncond-1;
            cfg.contrastcoefs = zeros(ncontrasts,ncond);
            cfg.contrastcoefs(:,1) = 1;
            for contrastindx=1:ncontrasts
                cfg.contrastcoefs(contrastindx,contrastindx+1)=-1;
            end;
        else
            ncontrasts = size(cfg.contrastcoefs,1);
        end;
        % compute the critical value for the (channel,frequency,timepoint)-specific F-tests
        dfnum = ncontrasts;
        dfdenom = nobs - ncontrasts;
        crit = ((nobs-1).*ncontrasts./(nobs-ncontrasts)).*finv(1-cfg.alphathresh,dfnum,dfdenom);
    end;
    % construct the matrix randselect, containing the specification of the
    % nranddraws different draws according to the randomization scheme.

    % the default randomization scheme involves that the ncond! possible 
    % condition orders are determined by independent multinomially distributed
    % events (thus, NO sampling with replacement, as in the case where the 
    % number of replications with a particular condition order is fixed).
    % Note that the actual condition order is not identical to the order
    % specified in randselect(drawindx,replindx,:). Instead, to save
    % memory, the random order of the n-th random draw is obtained by
    % applying the random reordering in randselect(n,replindx,:) to the
    % random order of the (n-1)-th random draw.
    % Note that randselect is always initialized according to the default 
    % randomization scheme. However, it will be overwritten if the user has
    % specified his own randomization scheme via the field
    % cfg.allowedperms, which is possible with the statistic depsamplesF.
    
    randselect = zeros(cfg.nranddraws+1,nrepl,ncond);
    for replindx=1:nrepl
        randselect(1,replindx,:)=[1:ncond];
    end;
    % first row corresponds to the observed assignment of condition orders to replications.    
    for drawindx=2:(cfg.nranddraws+1)
        for replindx=1:nrepl
            randselect(drawindx,replindx,:)=randperm(ncond);
        end;
    end;
    
    if strcmp(cfg.statistic,'depsamplesT')
        % compute the difference between the two conditions
        diff = squeeze(data.biol(:,1,:,:,:)-data.biol(:,2,:,:,:));
        % compute statistics that can be reused in the different random draws
        sumsqdiff=reshape(nan_sum(diff.^2,1),lengthspatdim,nfreq,ntime);
    end;
    
    if strcmp(cfg.statistic,'depsamplesregrT')
        % compute statistics that can be reused in the different random draws
        expandext = repmat(data.ext',[nrepl 1 lengthspatdim nfreq ntime]);
        centerfactor = reshape(nan_sum(data.biol,2),nrepl,lengthspatdim,nfreq,ntime)*mean(data.ext);
        centeredssqext = ncond*var(data.ext,1);
    end;
    
    if strcmp(cfg.statistic,'actvsblT')
        if containsnans
            actavg = reshape(nan_sum(data.biol(:,1,:,:,:),5)./sum(data.notnan(:,1,:,:,:),5),nrepl,lengthspatdim,nfreq);
            blavg = reshape(nan_sum(data.biol(:,2,:,:,:),5)./sum(data.notnan(:,2,:,:,:),5),nrepl,lengthspatdim,nfreq);
        else
            actavg = reshape(mean(data.biol(:,1,:,:,:),5),nrepl,lengthspatdim,nfreq);
            blavg = reshape(mean(data.biol(:,2,:,:,:),5),nrepl,lengthspatdim,nfreq);
        end;
        diff = zeros(nrepl,lengthspatdim,nfreq,ntime);
    end;
    
    if strcmp(cfg.statistic,'depsamplesF')
        if ~isempty(cfg.allowedperms)
            nallowedperms=size(cfg.allowedperms,1);
            for drawindx=2:(cfg.nranddraws+1)
                for replindx=1:nrepl
                    randint=ceil(rand*nallowedperms);
                    randselect(drawindx,replindx,:)=cfg.allowedperms(randint,:);
                end;
            end;
        end;
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ready with performing statistic-specific preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(cfg.clusterteststat, 'maxsum') | strcmp(cfg.makeclusters, 'no')
    maxsumstats = zeros(1,cfg.nranddraws);
elseif strcmp(cfg.clusterteststat,'orderedsums')
    maxsumstats = cell(1,cfg.nranddraws);
end;

for drawindx=1:(cfg.nranddraws+1)
    
    fprintf('randomization %d %d\n',(drawindx-1),cfg.nranddraws);
    
    if any(strcmp(cfg.statistic,{'indepsamplesT','indepsamplesZcoh','indepsamplesTsqcrs'}))
        % perform the randomization
        replselc1 = find(randselect(drawindx,:)==1);
        replselc2 = find(randselect(drawindx,:)==2);
        
        if containsnans
            % count the number of observations (non-NaNs) in the
            % replications replselc1. these counts are stored in an array with 
            % dimord chan_freq_time.
            nobsc1 = reshape(sum(data.notnan(replselc1,:,:,:,:),1),lengthspatdim,nfreq,ntime);
            % count the number of observations (non-NaNs) in the
            % replications replselc2.
            nobsc2 = reshape(sum(data.notnan(replselc2,:,:,:,:),1),lengthspatdim,nfreq,ntime);
        else
            nobsc1 = length(replselc1);
            nobsc2 = length(replselc2);
        end;

        % compute the means
        meanc1 = reshape(nan_sum(data.biol(replselc1,:,:,:,:),1),lengthspatdim,nfreq,ntime)./nobsc1;
        meanc2 = reshape(nan_sum(data.biol(replselc2,:,:,:,:),1),lengthspatdim,nfreq,ntime)./nobsc2;
        % note that meanc1 and meanc2 may contain NaNs (resulting from a 0/0 division), 
        % namely for those (chan,freq,time)-elements that contain NaNs for all replications
        % (in which case there are 0 observations for this element).
        
        if strcmp(cfg.statistic,'indepsamplesT')
            % compute the t-statistics
            if strcmp(cfg.mirrordata,'yes')
                stats = zeros(2*lengthspatdim,nfreq,ntime);
                stats(1:lengthspatdim,:,:)=sqrt(nobsc1.*nobsc2./nobs.*(nobs-2)).* ...
                    (meanc1-meanc2)./sqrt(sumsq-nobsc1.*(meanc1.^2)-nobsc2.*(meanc2.^2));
                stats((lengthspatdim+1):(2*lengthspatdim),:,:)=-stats(1:lengthspatdim,:,:);
            else
                stats = sqrt(nobsc1.*nobsc2./nobs.*(nobs-2)).* ...
                    (meanc1-meanc2)./sqrt(sumsq-nobsc1.*(meanc1.^2)-nobsc2.*(meanc2.^2));
            end;
        end;

        if strcmp(cfg.statistic,'indepsamplesZcoh')
            % compute the Z-statistics for coherence
            powc1=meanc1(powindcs,:,:);
            powc2=meanc2(powindcs,:,:);
            crspctrc1=meanc1(crspctrindcs,:,:);
            crspctrc2=meanc2(crspctrindcs,:,:);
            
            cohc1=abs(crspctrc1)./sqrt(powc1(powindcsforcrspctr(:,1),:,:).*powc1(powindcsforcrspctr(:,2),:,:));
            cohc2=abs(crspctrc2)./sqrt(powc2(powindcsforcrspctr(:,1),:,:).*powc2(powindcsforcrspctr(:,2),:,:));
            
            dofc1=reshape(sum(data.dof(replselc1,:,:),1),1,nfreq,ntime);
            dofc2=reshape(sum(data.dof(replselc2,:,:),1),1,nfreq,ntime);
            
            denomZ=sqrt(1./(repmat(dofc1,[ncrspctr,1])-2) + 1./(repmat(dofc2,[ncrspctr,1])-2));
            
            stats=(atanh(cohc1)-atanh(cohc2))./denomZ;
        end;
        
        if strcmp(cfg.statistic,'indepsamplesTsqcrs')
            % compute the T-square-statistics for the cross-spectrum
            decomplmeanc1=zeros(2,lengthspatdim,nfreq,ntime);
            decomplmeanc1(1,:,:,:)=real(meanc1);
            decomplmeanc1(2,:,:,:)=imag(meanc1);
            decomplmeanc2=zeros(2,lengthspatdim,nfreq,ntime);
            decomplmeanc2(1,:,:,:)=real(meanc2);
            decomplmeanc2(2,:,:,:)=imag(meanc2);
            covmat=zeros(2);
            diffmean=zeros(2,1);
            stats=zeros(lengthspatdim,nfreq,ntime);
            for chancmbindx=1:lengthspatdim
                for freqindx=1:nfreq
                    for timeindx=1:ntime
                        covmat=sumsqcp(:,:,chancmbindx,freqindx,timeindx)-nobsc1(chancmbindx,freqindx,timeindx)*...
                            decomplmeanc1(:,chancmbindx,freqindx,timeindx)*decomplmeanc1(:,chancmbindx,freqindx,timeindx)'...
                            -nobsc2(chancmbindx,freqindx,timeindx)*decomplmeanc2(:,chancmbindx,freqindx,timeindx)*...
                            decomplmeanc2(:,chancmbindx,freqindx,timeindx)';
                        covmat=covmat/(nobs(chancmbindx,freqindx,timeindx)-2);
                        diffmean=decomplmeanc1(:,chancmbindx,freqindx,timeindx)-decomplmeanc2(:,chancmbindx,freqindx,timeindx);
                        stats(chancmbindx,freqindx,timeindx)=nobsc1(chancmbindx,freqindx,timeindx)* ...
                            nobsc2(chancmbindx,freqindx,timeindx)/nobs(chancmbindx,freqindx,timeindx)*diffmean'*covmat^(-1)*diffmean;
                    end;
                end;
            end;
        end;
        
        if drawindx==1
            if strcmp(cfg.mirrordata,'yes')
                raweffect=zeros(2*lengthspatdim,nfreq,ntime);
                raweffect(1:lengthspatdim,:,:)=meanc1-meanc2;
                obsmeanc1(1:lengthspatdim,:,:)=meanc1;
                obsmeanc2(1:lengthspatdim,:,:)=meanc2;
                raweffect((lengthspatdim+1):(2*lengthspatdim),:,:)=-raweffect(1:lengthspatdim,:,:);
                obsmeanc1((lengthspatdim+1):(2*lengthspatdim),:,:)=-obsmeanc1(1:lengthspatdim,:,:);
                obsmeanc2((lengthspatdim+1):(2*lengthspatdim),:,:)=-obsmeanc2(1:lengthspatdim,:,:);
            else
                if strcmp(cfg.statistic,'indepsamplesZcoh')
                    raweffect=(cohc1-cohc2);
                    obsmeanc1=cohc1;
                    obsmeanc2=cohc2;
                else
                    raweffect=(meanc1-meanc2);
                    obsmeanc1=meanc1;
                    obsmeanc2=meanc2;
                end;
            end;
        end;
    end;

    if strcmp(cfg.statistic,'indepsamplesregrT')
        B=reshape(nan_sum(data.biol(randselect(drawindx,:),:,:,:).*expandext, 1),[lengthspatdim nfreq ntime])./nobs;
        B=B-meanbiol.*meanext;
        B=B./varext;
        
        sdreserror = reshape(nan_sum((data.biol(randselect(drawindx,:),:,:,:)-repmat(reshape(meanbiol,[1 lengthspatdim nfreq ntime]),[nrepl 1]) - ...
                repmat(reshape(B,[1 lengthspatdim nfreq ntime]),[nrepl 1]).* ...
                (expandext - repmat(reshape(meanext,[1 lengthspatdim nfreq ntime]),[nrepl 1]))).^2, 1), ...
                [lengthspatdim nfreq ntime]);
        sdreserror = sqrt(sdreserror./(nobs-2));

        stats = (B.*sdext.*sqrtnobs)./sdreserror;

        if drawindx==1
            raweffect=B;
        end;
    end;

    if strcmp(cfg.statistic,'depsamplesT')
        selflip = find(randselect(drawindx,:,1)==2);
        diff(selflip,:,:,:) = -diff(selflip,:,:,:);
        
        % compute the mean difference
        meandiff = reshape(nan_sum(diff,1),lengthspatdim,nfreq,ntime)./nobs;
        
        % compute the t-statistics
        stats = sqrt(nobs).*meandiff./sqrt((sumsqdiff-nobs.*meandiff.^2)./(nobs-1));

        if drawindx==1
            raweffect=meandiff;
        end;
    end;
    
    if strcmp(cfg.statistic,'actvsblT')
        selflip = find(randselect(drawindx,:,1)==2);
        selnoflip = find(randselect(drawindx,:,1)==1);
        
        for timeindx=1:ntime
            diff(selnoflip,:,:,timeindx)=reshape(data.biol(selnoflip,1,:,:,timeindx),length(selnoflip),lengthspatdim,nfreq,1)- ...
                                         reshape(blavg(selnoflip,:,:),length(selnoflip),lengthspatdim,nfreq,1);
            diff(selflip,:,:,timeindx)=reshape(data.biol(selflip,2,:,:,timeindx),length(selflip),lengthspatdim,nfreq,1)- ...
                                         reshape(actavg(selflip,:,:),length(selflip),lengthspatdim,nfreq,1);
        end;
            
        % compute the mean difference
        meandiff = reshape(nan_sum(diff,1),lengthspatdim,nfreq,ntime)./nobs;
        % compute the sum of the squared differences
        sumsqdiff=reshape(nan_sum(diff.^2,1),lengthspatdim,nfreq,ntime);
        
        % compute the t-statistics
        stats = sqrt(nobs).*meandiff./sqrt((sumsqdiff-nobs.*meandiff.^2)./(nobs-1));
        if drawindx==1
            raweffect=meandiff;
        end;
    end;

    if strcmp(cfg.statistic,'depsamplesregrT')
        for replindx=1:nrepl
            data.biol(replindx,:,:,:,:)=data.biol(replindx,randselect(drawindx,replindx,:),:,:,:);
        end;
        B=reshape(nan_sum(data.biol.*expandext,2),[nrepl lengthspatdim nfreq ntime]) - centerfactor;
        B=B./centeredssqext;
        meanB=reshape(nan_sum(B,1),[lengthspatdim nfreq ntime])./nobs;
        
        stats = (meanB.*sqrt(nobs))./reshape(std(B),[lengthspatdim nfreq ntime]);

        if drawindx==1
            raweffect=meanB;
        end;
    end;

    if strcmp(cfg.statistic,'indepsamplesF')
        % perform the randomization
        for condindx=1:ncond
            replsel{condindx} = find(randselect(drawindx,:)==condindx);
        end;
        
        if containsnans
            % count the number of observations (non-NaNs) in the
            % replications replsel{condindx}. these counts are stored in
            % arrays with dimord chan_freq_time.
            for condindx=1:ncond
                nobscell{condindx} = reshape(sum(data.notnan(replsel{condindx},:,:,:,:),1),lengthspatdim,nfreq,ntime);
            end;
        else
            for condindx=1:ncond
                nobscell{condindx} = length(replsel{condindx});
            end;
        end;

        % compute the group means
        for condindx=1:ncond
            meancell{condindx} = reshape(nan_sum(data.biol(replsel{condindx},:,:,:,:),1),lengthspatdim,nfreq,ntime)./nobscell{condindx};
        end;
        % remember that the NaNs in data.biol were replaced by zeros.
        % note that meancell{levelindx} may contain NaNs (resulting from a 0/0 division), 
        % namely for those (chan,freq,time)-elements that contain NaNs for all replications
        % (in which case there are 0 observations for this element).
        
        % compute the F-statistics
        ssbetween=zeros(lengthspatdim,nfreq,ntime);        
        sserror=sumsq;        
        for condindx=1:ncond
            ssbetween = ssbetween + nobscell{condindx}*(meancell{condindx}-overallmean).^2;
            sserror = sserror - nobscell{condindx}*meancell{condindx}.^2;
        end;
        stats = (ssbetween.*dfdenom)./(sserror.*dfnum);
    end;
    
    if strcmp(cfg.statistic,'depsamplesF')
        for replindx=1:nrepl
            data.biol(replindx,:,:,:,:)=data.biol(replindx,randselect(drawindx,replindx,:),:,:,:);
        end;
        contrasts=zeros(nrepl,ncontrasts,lengthspatdim,nfreq,ntime);
        for chanindx=1:lengthspatdim
            for freqindx=1:nfreq
                for timeindx=1:ntime
                    contrasts(:,:,chanindx,freqindx,timeindx)= ...
                        reshape(data.biol(:,:,chanindx,freqindx,timeindx),nrepl,ncond)*cfg.contrastcoefs';
                end;
            end;
        end;
        contrastavg=zeros(ncontrasts,lengthspatdim,nfreq,ntime);
        for contrastindx=1:ncontrasts
            contrastavg(contrastindx,:,:,:)=reshape(nan_sum(contrasts(:,contrastindx,:,:,:),1),lengthspatdim,nfreq,ntime)./nobs;
        end;
        covmat=zeros(ncontrasts,ncontrasts);
        if ~containsnans
            thisnobs = nobs;
        end;
        for chanindx=1:lengthspatdim
            for freqindx=1:nfreq
                for timeindx=1:ntime
                    if containsnans
                        thisnobs = nobs(chanindx,freqindx,timeindx);
                    end;
                    covmat=contrasts(:,:,chanindx,freqindx,timeindx)'*contrasts(:,:,chanindx,freqindx,timeindx) - ...
                        thisnobs*contrastavg(:,chanindx,freqindx,timeindx)*contrastavg(:,chanindx,freqindx,timeindx)';
                    covmat=covmat/(thisnobs-1);
                    stats(chanindx,freqindx,timeindx)=thisnobs*contrastavg(:,chanindx,freqindx,timeindx)'* ...
                        covmat^(-1)*contrastavg(:,chanindx,freqindx,timeindx);
                end;
            end;
        end;
    end;
    

    if strcmp(cfg.makeclusters, 'yes')
        % perform the clustering
        
        % positive statistics
        onoff = logical(stats>crit);
        if singlechans
            [posclusterslabelmat,nposclusters] = findcluster(onoff,data.channeighbstructmat,cfg.minnbchan);
        else
            if isfield(data,'chancmbneighbselmat')
                [posclusterslabelmat,nposclusters] = findcluster(onoff,data.chancmbneighbstructmat,data.chancmbneighbselmat,cfg.minnbchan);
            else
                [posclusterslabelmat,nposclusters] = findcluster(onoff,data.chancmbneighbstructmat,data.chancmbneighbstructmat,cfg.minnbchan);
            end;
        end;
        sumposstatspercluster=[];
        for clusterindx=1:nposclusters
            selstats = stats(posclusterslabelmat==clusterindx);
            if length(selstats)>=cfg.smallestcluster
                sumposstatspercluster(end+1)= sum(selstats(:));
            end;
        end;
        
        % negative statistics
        % For one-dimensional tests with a critical value in the
        % left tail (i.e., two sided tests and one-sided tests of the null
        % hypothesis 2<1) the cluster-level statistics are only computed
        % for the observed data (drawindx=1). For testing these negative
        % statistics, we will use the same reference distribution as for
        % the positive statistics, but mirrored around zero (positive
        % becomes negative).
        if (drawindx==1 & critnegativetail)
            onoff = logical(stats<=-crit);
            % perform the clustering
            if singlechans
                [negclusterslabelmat,nnegclusters] = findcluster(onoff,data.channeighbstructmat,cfg.minnbchan);
            else
                if isfield(data,'chancmbneighbselmat')
                    [negclusterslabelmat,nnegclusters] = findcluster(onoff,data.chancmbneighbstructmat,data.chancmbneighbselmat,cfg.minnbchan);
                else
                    [negclusterslabelmat,nnegclusters] = findcluster(onoff,data.chancmbneighbstructmat,data.chancmbneighbstructmat,cfg.minnbchan);
                end;
            end;
        end;
        
        if drawindx==1
            if strcmp(cfg.onetwo,'twosided') | ~critnegativetail
              obsposclusterslabelmat=posclusterslabelmat;
              obsnposclusters=nposclusters;
            else
              obsposclusterslabelmat=[];
              obsnposclusters=0;
            end;
            if critnegativetail
              obsnegclusterslabelmat=negclusterslabelmat;
              obsnnegclusters=nnegclusters;
            else
              obsnegclusterslabelmat=[];
              obsnnegclusters=0;
            end;
            obsstats=stats;
        else
            if strcmp(cfg.clusterteststat,'maxsum')
                if length(sumposstatspercluster)>0
                    maxsumstats(drawindx-1)=max(sumposstatspercluster);
                else
                    maxsumstats(drawindx-1)=-inf;
                end;
            elseif strcmp(cfg.clusterteststat,'orderedsums')
                if length(sumposstatspercluster)>0
                    % sort in descending order
                    maxsumstats{drawindx-1} = -sort(-sumposstatspercluster);
                end;
            end;        
        end;
    elseif strcmp(cfg.makeclusters, 'no')
        if drawindx==1
            obsstats=stats;
        else
            maxsumstats(drawindx-1)=max(stats(:));
        end;
    end;  % if strcmp(cfg.makeclusters, 'yes')
end;   % for drawindx=0:cfg.nranddraws


if strcmp(cfg.onetwo,'twosided')
    critalpha = cfg.alpha/2;
else
    critalpha = cfg.alpha;
end;
critvals = [];

if cfg.nranddraws>0
    if strcmp(cfg.makeclusters, 'yes')
        posclusters=[];
        negclusters=[];
        sorteer=[];
        if strcmp(cfg.clusterteststat, 'maxsum') 
            maxsumstats = sort(maxsumstats);
            critvals = maxsumstats(max([1,round((1-critalpha)*cfg.nranddraws)]));
            if onedimtest 
                if strcmp(cfg.onetwo,'twosided')
                    critvals = [-critvals,critvals];
                elseif critnegativetail
                    critvals = -critvals;
                end;
            end;
            
            for clusterindx=1:obsnposclusters
                selstats = obsstats(obsposclusterslabelmat==clusterindx);
                posclusters(end+1).sumstat = sum(selstats(:));
                posclusters(end  ).pval     = sum(maxsumstats>posclusters(end).sumstat)/cfg.nranddraws;
                posclusters(end  ).size     = length(selstats);
                sorteer(end+1)= posclusters(end).sumstat;
            end;    
            [dum,i]=sort(-sorteer);
            posclusters=posclusters(i);
            posclusterslabelmat=zeros(size(obsposclusterslabelmat));
            for clusterindx=1:obsnposclusters
                posclusterslabelmat(find(obsposclusterslabelmat==i(clusterindx))) = clusterindx;
            end;
            obsposclusterslabelmat=posclusterslabelmat;
            
            sorteer=[];
            for clusterindx=1:obsnnegclusters
                selstats = obsstats(obsnegclusterslabelmat==clusterindx);
                negclusters(end+1).sumstat = sum(selstats(:));
                % we now make use of the fact that the permutation
                % distribution for the negative cluster statistics is the
                % mirror image of the permutation distribution for the
                % positive cluster statistics
                negclusters(end  ).pval     = sum(maxsumstats>abs(negclusters(end).sumstat))/cfg.nranddraws;
                negclusters(end  ).size     = length(selstats);
                sorteer(end+1)= abs(negclusters(end).sumstat);
            end;    
            [dum,i]=sort(-sorteer);
            negclusters=negclusters(i);
            negclusterslabelmat=zeros(size(obsnegclusterslabelmat));
            for clusterindx=1:obsnnegclusters
                negclusterslabelmat(find(obsnegclusterslabelmat==i(clusterindx))) = clusterindx;
            end;
            obsnegclusterslabelmat=negclusterslabelmat;
            
        elseif  strcmp(cfg.clusterteststat, 'orderedsums') 
            % find the critical values for the first n order statistics of the summed
            % over the elements in a cluster) t-statistics, where n is a random
            % variable determined by cfg.smallestcluster.
            
            % store all data in the cell array maxsumstats in an array of
            % size [cfg.nranddraws maxnlargeclusters], where maxnlargeclusters
            % is the maximum (over the cfg.nranddraws draws) number of clusters
            maxnlargeclusters = 0;
            for drawindx=1:cfg.nranddraws
                if length(maxsumstats{drawindx})>maxnlargeclusters
                    maxnlargeclusters=length(maxsumstats{drawindx});
                end;
            end;
            % initialize
            orderstats=-inf*ones(cfg.nranddraws,maxnlargeclusters);
            for drawindx=1:cfg.nranddraws
                if length(maxsumstats{drawindx})>0
                    % remember that the elements in maxsumstats{drawindx} are
                    % sorted in descending order (the largest order statistics appear 
                    % in the first columns).
                    orderstats(drawindx,1:length(maxsumstats{drawindx})) = maxsumstats{drawindx};
                end;
            end;
            % sort the order statistics per column (rank ordered cluster)
            sortedorderstats=zeros(cfg.nranddraws,maxnlargeclusters);
            for clusterindx=1:maxnlargeclusters
                % sort in ascending order
                sortedorderstats(:,clusterindx)=sort(orderstats(:,clusterindx));
            end;
            % initializing the value of the marginal critical quantile
            margcritquantile=1;
            % compute the gap between the cumulative proportions (values of the empirical cumulative 
            % distribution function) of adjacent empirical quantiles
            gap=1/cfg.nranddraws;
            % next line to be used in combination with CEIL below
            margcritquantile=1+gap/2;
            
            margcritvals = [];
            actualalpha=0;
            while actualalpha<critalpha & margcritquantile>0 & maxnlargeclusters>0
                margcritquantile=margcritquantile-gap;
                margcritvals=sortedorderstats(max([1,ceil(margcritquantile*cfg.nranddraws)]),:);
                nlarger=0;
                for drawindx=1:cfg.nranddraws
                    if any(orderstats(drawindx,:)>margcritvals)
                        nlarger=nlarger+1;
                    end;
                end;
                actualalpha=nlarger/cfg.nranddraws;
            end;
            % find marginal critical values such that actualalpha is less than
            % critalpha
            if maxnlargeclusters>0
                margcritquantile=margcritquantile+gap;
                margcritvals=sortedorderstats(max([1,ceil(margcritquantile*cfg.nranddraws)]),:);
                critvals=margcritvals';
                if strcmp(cfg.onetwo,'twosided')
                    critvals=[-critvals,critvals];
                end;
                nlarger=0;
                for drawindx=1:cfg.nranddraws
                    if orderstats(drawindx,:)>margcritvals
                        nlarger=nlarger+1;
                    end;
                end;
                actualalpha=nlarger/cfg.nranddraws;
            end;
            
            % store the results in the clusters struct array.
            posclusterslabelmat=zeros(size(obsposclusterslabelmat));
            for clusterindx=1:obsnposclusters
                selstats = obsstats(obsposclusterslabelmat==clusterindx);
                if length(selstats)>=cfg.smallestcluster
                    posclusters(end+1).sumstat = sum(selstats(:));
                    posclusters(end  ).size     = length(selstats);
                    sorteer(end+1)= posclusters(end).sumstat;
                    clusterlabel=length(posclusters);
                    posclusterslabelmat(obsposclusterslabelmat==clusterindx)=clusterlabel;
                end;
            end;
            obsnlargeposclusters = length(posclusters);
            [dum,i]=sort(-sorteer);
            posclusters=posclusters(i);
            tempclusterslabelmat=posclusterslabelmat;
            posclusterslabelmat=zeros(size(obsposclusterslabelmat));
            for clusterindx=1:obsnlargeposclusters
                posclusterslabelmat(find(tempclusterslabelmat==i(clusterindx))) = clusterindx;
            end;
            obsposclusterslabelmat=posclusterslabelmat;
            
            % if needed, enlarge margcritvals to avoid problems later on.
            margcritvals = [margcritvals -inf*ones(1,(obsnlargeposclusters-length(margcritvals)))]; 
            for clusterindx=1:obsnlargeposclusters
                posclusters(clusterindx).margcritval=margcritvals(clusterindx);
                if posclusters(clusterindx).sumstat>posclusters(clusterindx).margcritval
                    posclusters(clusterindx).significant='yes';
                else
                    posclusters(clusterindx).significant='no';
                end;
            end;
            
            sorteer=[];
            negclusterslabelmat=zeros(size(obsnegclusterslabelmat));
            for clusterindx=1:obsnnegclusters
                selstats = obsstats(obsnegclusterslabelmat==clusterindx);
                if length(selstats)>=cfg.smallestcluster
                    negclusters(end+1).sumstat = sum(selstats(:));
                    negclusters(end  ).size     = length(selstats);
                    sorteer(end+1)= abs(negclusters(end).sumstat);
                    clusterlabel=length(posclusters);
                    negclusterslabelmat(obsnegclusterslabelmat==clusterindx)=clusterlabel;
                end
            end;    
            obsnlargenegclusters = length(negclusters);
            [dum,i]=sort(-sorteer);
            negclusters=negclusters(i);
            tempclusterslabelmat=negclusterslabelmat;
            negclusterslabelmat=zeros(size(obsnegclusterslabelmat));
            for clusterindx=1:obsnlargenegclusters
                negclusterslabelmat(find(tempclusterslabelmat==i(clusterindx))) = clusterindx;
            end;
            obsnegclusterslabelmat=negclusterslabelmat;
            
            % if needed, enlarge margcritvals to avoid problems later on.
            margcritvals = [margcritvals -inf*ones(1,(obsnlargenegclusters-length(margcritvals)))]; 
            for clusterindx=1:obsnlargenegclusters
                negclusters(clusterindx).margcritval=-margcritvals(clusterindx);
                if negclusters(clusterindx).sumstat<negclusters(clusterindx).margcritval
                    negclusters(clusterindx).significant='yes';
                else
                    negclusters(clusterindx).significant='no';
                end;
            end;
        end;  %  if strcmp(cfg.clusterteststat, ) 
    elseif strcmp(cfg.makeclusters, 'no')  % do not make clusters
        nsigelements=[];
        pvals=[];
        
        maxsumstats = sort(maxsumstats);
        uppercritval = maxsumstats(max([1,round((1-critalpha)*cfg.nranddraws)]));
        if onedimtest
            if strcmp(cfg.onetwo,'twosided')
                critvals = [-uppercritval,uppercritval];
                onoff = (obsstats>uppercritval) | (obsstats<-uppercritval);
            elseif critnegativetail
                critvals = -critvals;
                onoff = (obsstats<-uppercritval);
            else
                onoff = (obsstats>uppercritval);
            end;
        else
            onoff = (obsstats>uppercritval);
        end;
        nsigelements = sum(onoff(:));
        dims = size(obsstats);
        dims((length(dims)+1):3)=1;
        pvals = zeros(dims);
        for chanindx=1:dims(1)
            for freqindx=1:dims(2)
                for timeindx=1:dims(3)
                    if onedimtest
                        if strcmp(cfg.onetwo,'twosided')
                            pvals(chanindx,freqindx,timeindx)=sum(maxsumstats<abs(obsstats(chanindx,freqindx,timeindx)))/cfg.nranddraws;
                        elseif critnegativetail
                            pvals(chanindx,freqindx,timeindx)=sum(-maxsumstats<obsstats(chanindx,freqindx,timeindx))/cfg.nranddraws;
                        else
                            pvals(chanindx,freqindx,timeindx)=sum(maxsumstats>obsstats(chanindx,freqindx,timeindx))/cfg.nranddraws;
                        end;
                    else
                        pvals(chanindx,freqindx,timeindx)=sum(maxsumstats>stats(chanindx,freqindx,timeindx))/cfg.nranddraws;
                    end;
                end;
            end;
        end;
    end;
end;    %   if cfg.nranddraws>0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% collect the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clusrand.stats = squeeze(obsstats);
if exist('raweffect')==1
    clusrand.raweffect=squeeze(raweffect);
end;
if exist('obsmeanc1')==1
    clusrand.obsmeanc1=squeeze(obsmeanc1);
end;
if exist('obsmeanc2')==1
    clusrand.obsmeanc2=squeeze(obsmeanc2);
end;
if  strcmp(cfg.makeclusters, 'yes')
  clusrand.posclusters = posclusters;
  clusrand.negclusters = negclusters;
  clusrand.posclusterslabelmat    = obsposclusterslabelmat;
  clusrand.negclusterslabelmat    = obsnegclusterslabelmat;
else
  if cfg.nranddraws>0
    clusrand.pvals    = pvals;
    clusrand.nsigelements    = nsigelements;
  end;
end;

if strcmp(cfg.clusterteststat, 'orderedsums') 
  clusrand.nlargeposclusters = obsnlargeposclusters;
  clusrand.nlargenegclusters = obsnlargenegclusters;
  clusrand.margcritquantile = margcritquantile*100;
  clusrand.MonteCarloalpha = actualalpha;
end
if cfg.nranddraws>0
  clusrand.critvals = critvals;
end;
if isfield(data,'labelcmb')
    if strcmp(cfg.statistic,'indepsamplesZcoh')
        clusrand.labelcmb = data.labelcmb(crspctrindcs,:);
    else
        clusrand.labelcmb=data.labelcmb;
    end;
else
    if isfield(data,'label')
        clusrand.label=data.label;
    end;
end;
if isfield(data,'freq')
    clusrand.freq=data.freq;
end;
if isfield(data,'time')
    clusrand.time=data.time;
end;
clusrand.cfg  = cfg;			% remember the configuration details

% Turn warnings on.
warning on;


