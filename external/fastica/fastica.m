function [Out1, Out2, Out3] = fastica(mixedsig, varargin)
%FASTICA - Fast Independent Component Analysis
%
%   NOTICE!
% THIS IS NOT THE ORIGINAL VERSION BUT A MODIFIED VERSION: 
% this one does not compute the icasig if this is not necessary (otherwise
% there is no change, thus there is (should be) no change in behaviour
%   kaigoe  (kai.goergen@bccn-berlin.de)
%
% FastICA for Matlab 7.x and 6.x
% Version 2.5, October 19 2005
% Copyright (c) Hugo Gävert, Jarmo Hurri, Jaakko Särelä, and Aapo Hyvärinen.
%
% FASTICA(mixedsig) estimates the independent components from given
% multidimensional signals. Each row of matrix mixedsig is one
% observed signal.  FASTICA uses Hyvarinen's fixed-point algorithm,
% see http://www.cis.hut.fi/projects/ica/fastica/. Output from the
% function depends on the number output arguments:
%
% [icasig] = FASTICA (mixedsig); the rows of icasig contain the
% estimated independent components.
%
% [icasig, A, W] = FASTICA (mixedsig); outputs the estimated separating
% matrix W and the corresponding mixing matrix A.
%
% [A, W] = FASTICA (mixedsig); gives only the estimated mixing matrix
% A and the separating matrix W.
%
% Some optional arguments induce other output formats, see below.
%
% A graphical user interface for FASTICA can be launched by the
% command FASTICAG
%
% FASTICA can be called with numerous optional arguments. Optional
% arguments are given in parameter pairs, so that first argument is
% the name of the parameter and the next argument is the value for
% that parameter. Optional parameter pairs can be given in any order.
%
% OPTIONAL PARAMETERS:
%
% Parameter name        Values and description
%
%======================================================================
% --Basic parameters in fixed-point algorithm:
%
% 'approach'            (string) The decorrelation approach used. Can be
%                       symmetric ('symm'), i.e. estimate all the
%                       independent component in parallel, or
%                       deflation ('defl'), i.e. estimate independent
%                       component one-by-one like in projection pursuit.
%                       Default is 'defl'.
%
% 'numOfIC'             (integer) Number of independent components to
%                       be estimated. Default equals the dimension of data.
%
%======================================================================
% --Choosing the nonlinearity:
%
% 'g'                   (string) Chooses the nonlinearity g used in 
%                       the fixed-point algorithm. Possible values:
%
%                       Value of 'g':      Nonlinearity used:
%                       'pow3' (default)   g(u)=u^3
%                       'tanh'             g(u)=tanh(a1*u)
%                       'gauss             g(u)=u*exp(-a2*u^2/2)
%                       'skew'             g(u)=u^2
% 
% 'finetune'		(string) Chooses the nonlinearity g used when 
%                       fine-tuning. In addition to same values
%                       as for 'g', the possible value 'finetune' is:
%                       'off'              fine-tuning is disabled.
%
% 'a1'                  (number) Parameter a1 used when g='tanh'.
%                       Default is 1.
% 'a2'                  (number) Parameter a2 used when g='gaus'.
%                       Default is 1.
%
% 'mu'			(number) Step size. Default is 1.
%                       If the value of mu is other than 1, then the
%                       program will use the stabilized version of the
%                       algorithm (see also parameter 'stabilization').
%
%
% 'stabilization'       (string) Values 'on' or 'off'. Default 'off'. 
%                       This parameter controls wether the program uses
%                       the stabilized version of the algorithm or
%                       not. If the stabilization is on, then the value
%                       of mu can momentarily be halved if the program
%                       senses that the algorithm is stuck between two
%                       points (this is called a stroke). Also if there
%                       is no convergence before half of the maximum
%                       number of iterations has been reached then mu
%                       will be halved for the rest of the rounds.
% 
%======================================================================
% --Controlling convergence:
%
% 'epsilon'             (number) Stopping criterion. Default is 0.0001.
%
% 'maxNumIterations'    (integer) Maximum number of iterations.
%                       Default is 1000.
%
% 'maxFinetune'         (integer) Maximum number of iterations in 
%                       fine-tuning. Default 100.
%
% 'sampleSize'          (number) [0 - 1] Percentage of samples used in
%                       one iteration. Samples are chosen in random.
%                       Default is 1 (all samples).
%
% 'initGuess'           (matrix) Initial guess for A. Default is random.
%                       You can now do a "one more" like this: 
%                       [ica, A, W] = fastica(mix, 'numOfIC',3);
%                       [ica2, A2, W2] = fastica(mix, 'initGuess', A, 'numOfIC', 4);
%
%======================================================================
% --Graphics and text output:
%
% 'verbose'             (string) Either 'on' or 'off'. Default is
%                       'on': report progress of algorithm in text format.
%
% 'displayMode'         (string) Plot running estimates of independent
%                       components: 'signals', 'basis', 'filters' or
%                       'off'. Default is 'off'.
%
% 'displayInterval'     Number of iterations between plots.
%                       Default is 1 (plot after every iteration).
%
%======================================================================
% --Controlling reduction of dimension and whitening:
%
% Reduction of dimension is controlled by 'firstEig' and 'lastEig', or
% alternatively by 'interactivePCA'. 
%
% 'firstEig'            (integer) This and 'lastEig' specify the range for
%                       eigenvalues that are retained, 'firstEig' is
%                       the index of largest eigenvalue to be
%                       retained. Default is 1.
%
% 'lastEig'             (integer) This is the index of the last (smallest)
%                       eigenvalue to be retained. Default equals the
%                       dimension of data.
%
% 'interactivePCA'      (string) Either 'on' or 'off'. When set 'on', the
%                       eigenvalues are shown to the user and the
%                       range can be specified interactively. Default
%                       is 'off'. Can also be set to 'gui'. Then the user
%                       can use the same GUI that's in FASTICAG.
%
% If you already know the eigenvalue decomposition of the covariance
% matrix, you can avoid computing it again by giving it with the
% following options:
%
% 'pcaE'                (matrix) Eigenvectors
% 'pcaD'                (matrix) Eigenvalues
%
% If you already know the whitened data, you can give it directly to
% the algorithm using the following options:
%
% 'whiteSig'            (matrix) Whitened signal
% 'whiteMat'            (matrix) Whitening matrix
% 'dewhiteMat'          (matrix) dewhitening matrix
%
% If values for all the 'whiteSig', 'whiteSig' and 'dewhiteMat' are
% supplied, they will be used in computing the ICA. PCA and whitening
% are not performed. Though 'mixedsig' is not used in the main
% algorithm it still must be entered - some values are still
% calculated from it.
%
% Performing preprocessing only is possible by the option:
%
% 'only'                (string) Compute only PCA i.e. reduction of
%                       dimension ('pca') or only PCA plus whitening
%                       ('white'). Default is 'all': do ICA estimation
%                       as well.  This option changes the output
%                       format accordingly. For example: 
%
%                       [whitesig, WM, DWM] = FASTICA(mixedsig, 
%                       'only', 'white') 
%                       returns the whitened signals, the whitening matrix
%                       (WM) and the dewhitening matrix (DWM). (See also
%                       WHITENV.) In FastICA the whitening matrix performs
%                       whitening and the reduction of dimension. Dewhitening
%                       matrix is the pseudoinverse of whitening matrix.
%                        
%                       [E, D] = FASTICA(mixedsig, 'only', 'pca') 
%                       returns the eigenvector (E) and diagonal 
%                       eigenvalue (D) matrices  containing the 
%                       selected subspaces. 
%
%======================================================================
% EXAMPLES
%
%       [icasig] = FASTICA (mixedsig, 'approach', 'symm', 'g', 'tanh');
%               Do ICA with tanh nonlinearity and in parallel (like
%               maximum likelihood estimation for supergaussian data).
%
%       [icasig] = FASTICA (mixedsig, 'lastEig', 10, 'numOfIC', 3);
%               Reduce dimension to 10, and estimate only 3
%               independent components.
%
%       [icasig] = FASTICA (mixedsig, 'verbose', 'off', 'displayMode', 'off');
%               Don't output convergence reports and don't plot
%               independent components.
%
%
% A graphical user interface for FASTICA can be launched by the
% command FASTICAG
%
%   See also FASTICAG

% @(#)$Id$

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check some basic requirements of the data
if nargin == 0,
  error ('You must supply the mixed data as input argument.');
end

if length (size (mixedsig)) > 2,
  error ('Input data can not have more than two dimensions.');
end

if any (any (isnan (mixedsig))),
  error ('Input data contains NaN''s.');
end

if ~isa (mixedsig, 'double')
  fprintf ('Warning: converting input data into regular (double) precision.\n');
  mixedsig = double (mixedsig);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove the mean and check the data

[mixedsig, mixedmean] = remmean(mixedsig);

[Dim, NumOfSampl] = size(mixedsig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default values for optional parameters

% All
verbose           = 'on';

% Default values for 'pcamat' parameters
firstEig          = 1;
lastEig           = Dim;
interactivePCA    = 'off';

% Default values for 'fpica' parameters
approach          = 'defl';
numOfIC           = Dim;
g                 = 'pow3';
finetune          = 'off';
a1                = 1;
a2                = 1;
myy               = 1;
stabilization     = 'off';
epsilon           = 0.0001;
maxNumIterations  = 1000;
maxFinetune       = 5;
initState         = 'rand';
guess             = 0;
sampleSize        = 1;
displayMode       = 'off';
displayInterval   = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for fastICA - i.e. this file

b_verbose = 1;
jumpPCA = 0;
jumpWhitening = 0;
only = 3;
userNumOfIC = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the optional parameters

if (rem(length(varargin),2)==1)
  error('Optional parameters should always go by pairs');
else
  for i=1:2:(length(varargin)-1)
    if ~ischar (varargin{i}),
      error (['Unknown type of optional parameter name (parameter' ...
	      ' names must be strings).']);
    end
    % change the value of parameter
    switch lower (varargin{i})
     case 'stabilization'
      stabilization = lower (varargin{i+1});
     case 'maxfinetune'
      maxFinetune = varargin{i+1};
     case 'samplesize'
      sampleSize = varargin{i+1};
     case 'verbose'
      verbose = lower (varargin{i+1});
      % silence this program also
      if strcmp (verbose, 'off'), b_verbose = 0; end
     case 'firsteig'
      firstEig = varargin{i+1};
     case 'lasteig'
      lastEig = varargin{i+1};
     case 'interactivepca'
      interactivePCA = lower (varargin{i+1});
     case 'approach'
      approach = lower (varargin{i+1});
     case 'numofic'
      numOfIC = varargin{i+1};
      % User has supplied new value for numOfIC.
      % We'll use this information later on...
      userNumOfIC = 1;
     case 'g'
      g = lower (varargin{i+1});
     case 'finetune'
      finetune = lower (varargin{i+1});
     case 'a1'
      a1 = varargin{i+1};
     case 'a2'
      a2 = varargin{i+1};
     case {'mu', 'myy'}
      myy = varargin{i+1};
     case 'epsilon'
      epsilon = varargin{i+1};
     case 'maxnumiterations'
      maxNumIterations = varargin{i+1};
     case 'initguess'
      % no use setting 'guess' if the 'initState' is not set
      initState = 'guess';
      guess = varargin{i+1};
     case 'displaymode'
      displayMode = lower (varargin{i+1});
     case 'displayinterval'
      displayInterval = varargin{i+1};
     case 'pcae'
      % calculate if there are enought parameters to skip PCA
      jumpPCA = jumpPCA + 1;
      E = varargin{i+1};
     case 'pcad'
      % calculate if there are enought parameters to skip PCA
      jumpPCA = jumpPCA + 1;
      D = varargin{i+1};
     case 'whitesig'
      % calculate if there are enought parameters to skip PCA and whitening
      jumpWhitening = jumpWhitening + 1;
      whitesig = varargin{i+1};
     case 'whitemat'
      % calculate if there are enought parameters to skip PCA and whitening
      jumpWhitening = jumpWhitening + 1;
      whiteningMatrix = varargin{i+1};
     case 'dewhitemat'
      % calculate if there are enought parameters to skip PCA and whitening
      jumpWhitening = jumpWhitening + 1;
      dewhiteningMatrix = varargin{i+1};
     case 'only'
      % if the user only wants to calculate PCA or...
      switch lower (varargin{i+1})
       case 'pca'
	only = 1;
       case 'white'
	only = 2;
       case 'all'
	only = 3;
      end
      
     otherwise
      % Hmmm, something wrong with the parameter string
      error(['Unrecognized parameter: ''' varargin{i} '''']);
    end;
  end;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print information about data
if b_verbose
  fprintf('Number of signals: %d\n', Dim);
  fprintf('Number of samples: %d\n', NumOfSampl);
end

% Check if the data has been entered the wrong way,
% but warn only... it may be on purpose

if Dim > NumOfSampl
  if b_verbose
    fprintf('Warning: ');
    fprintf('The signal matrix may be oriented in the wrong way.\n');
    fprintf('In that case transpose the matrix.\n\n');
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating PCA

% We need the results of PCA for whitening, but if we don't
% need to do whitening... then we dont need PCA...
if jumpWhitening == 3
  if b_verbose,
    fprintf ('Whitened signal and corresponding matrises supplied.\n');
    fprintf ('PCA calculations not needed.\n');
  end;
else
  
  % OK, so first we need to calculate PCA
  % Check to see if we already have the PCA data
  if jumpPCA == 2,
    if b_verbose,
      fprintf ('Values for PCA calculations supplied.\n');
      fprintf ('PCA calculations not needed.\n');
    end;
  else
    % display notice if the user entered one, but not both, of E and D.
    if (jumpPCA > 0) & (b_verbose),
      fprintf ('You must suply all of these in order to jump PCA:\n');
      fprintf ('''pcaE'', ''pcaD''.\n');
    end;
    
    % Calculate PCA
    [E, D]=pcamat(mixedsig, firstEig, lastEig, interactivePCA, verbose);
  end
end

% skip the rest if user only wanted PCA
if only > 1
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Whitening the data
  
  % Check to see if the whitening is needed...
  if jumpWhitening == 3,
    if b_verbose,
      fprintf ('Whitening not needed.\n');
    end;
  else
    
    % Whitening is needed
    % display notice if the user entered some of the whitening info, but not all.
    if (jumpWhitening > 0) & (b_verbose),
      fprintf ('You must suply all of these in order to jump whitening:\n');
      fprintf ('''whiteSig'', ''whiteMat'', ''dewhiteMat''.\n');
    end;
    
    % Calculate the whitening
    [whitesig, whiteningMatrix, dewhiteningMatrix] = whitenv ...
						     (mixedsig, E, D, verbose);
  end
  
end % if only > 1

% skip the rest if user only wanted PCA and whitening
if only > 2
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Calculating the ICA
  
  % Check some parameters
  % The dimension of the data may have been reduced during PCA calculations.
  % The original dimension is calculated from the data by default, and the
  % number of IC is by default set to equal that dimension.
  
  Dim = size(whitesig, 1);
  
  % The number of IC's must be less or equal to the dimension of data
  if numOfIC > Dim
    numOfIC = Dim;
    % Show warning only if verbose = 'on' and user supplied a value for 'numOfIC'
    if (b_verbose & userNumOfIC)
      fprintf('Warning: estimating only %d independent components\n', numOfIC);
      fprintf('(Can''t estimate more independent components than dimension of data)\n');
    end
  end
  
  % Calculate the ICA with fixed point algorithm.
  [A, W] = fpica (whitesig,  whiteningMatrix, dewhiteningMatrix, approach, ...
		  numOfIC, g, finetune, a1, a2, myy, stabilization, epsilon, ...
		  maxNumIterations, maxFinetune, initState, guess, sampleSize, ...
		  displayMode, displayInterval, verbose);
  
  % Check for valid return
%   if ~isempty(W)   % ORIGINAL VERSION

% MODIFIED VERSION: CHECK ALSO IF ICASIG IS RETURNED AT ALL, before we
% comput it
% CHANGED BY kaigoe (kai.goergen@bccn-berlin.de)

    % modified version
  if ~isempty(W)  ...
          && nargout ~= 2  %% if nargout == 2, we return [W. A], and NOT ICASIG
    % Add the mean back in.
    if b_verbose
      fprintf('Adding the mean back to the data.\n');
    end
    icasig = W * mixedsig + (W * mixedmean) * ones(1, NumOfSampl);
    %icasig = W * mixedsig;
    if b_verbose & ...
	  (max(abs(W * mixedmean)) > 1e-9) & ...
	  (strcmp(displayMode,'signals') | strcmp(displayMode,'on'))
      fprintf('Note that the plots don''t have the mean added.\n');
    end
  else
    icasig = [];
  end

end % if only > 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The output depends on the number of output parameters
% and the 'only' parameter.

if only == 1    % only PCA
  Out1 = E;
  Out2 = D;
elseif only == 2  % only PCA & whitening
  if nargout == 2
    Out1 = whiteningMatrix;
    Out2 = dewhiteningMatrix;
  else
    Out1 = whitesig;
    Out2 = whiteningMatrix;
    Out3 = dewhiteningMatrix;
  end
else      % ICA
  if nargout == 2
    Out1 = A;
    Out2 = W;
  else
    Out1 = icasig;
    Out2 = A;
    Out3 = W;
  end
end
