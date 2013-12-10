function sR=icassoEst(mode,X,M,varargin)
%function sR=icassoEst(mode,X,M,['FastICAparamName1',value1,'FastICAparamName2',value2,...])
%
%PURPOSE
%
%To compute randomized ICA estimates M times from data X. Output of
%this function (sR) is called 'Icasso result structure' (see
%icassoStruct). sR keeps the on all the methods, parameters, and
%results in the Icasso procedure.  
%
%EXAMPLES OF BASIC USAGE
%
%   sR=icassoEst('randinit', X, 30); 
%
%estimates ICA for 30 times on data matrix X using Icasso
%default parameters for FastICA: symmetrical approach, kurtosis as
%contrast function. In maximum 100 iterations are used for
%estimating ICA in each round. Randomizes only initial conditions.   
%
%   sR=icassoEst('both', X, 30, 'g', 'tanh', 'approach', 'defl');
%
%estimates ICA for 15 times on data matrix X using 'tanh' as the
%contrast function and the deflatory approach in FastICA. Applies
%both bootstrapping the data and randomizing initial conditions. 
%
%INPUT 
%
% mode (string) 'randinit' | 'bootstrap | 'both' 
% X    (dxN matrix) data where d=dimension, N=number of vectors 
% M    (scalar) number of randomizations (estimation cycles) 
%
%Optional input arguments are given as argument identifier - value
%pairs: 'identifier1', value1, 'identifier2', value2,... 
%(case insensitive)  
%
%FastICA parameters apply here (see function fastica)  
%Default: 'approach', 'symm', 'g', 'pow3', 'maxNumIterations', 100 
%
%OUTPUTS 
%
% sR (struct) Icasso result data structure
%
%DETAILS 
%
%Meaning of different choices for input arg. 'mode'
% 'randinit': different random initial condition each time.
% 'bootstrap': the same initial cond. each time, but data is
%   bootstrapped. The initial condition can be explicitly
%   specified using FastICA parameter 'initGuess'.
% 'both': use both data bootstrapping and randomization of 
%    initial condition.
%
%FASTICA PARAMETERS See function 'fastica' in FastICA toolbox for
%more information. Note that the following FastICA parameters
%cannot be used: 
% 
% In all modes ('randinit','bootstrap','both'): 
%   using 'interactivePCA','sampleSize', 'displayMode', 'displayInterval',
%   and 'only' are not allowed for obvious reasons. In addition, 
% in modes 'randinit' and 'both': 
%   using 'initGuess' is not allowed since initial guess is
%   randomized, and
% in modes 'bootstrap' and 'both': 
%   using 'whiteMat', 'dewhiteMat', and 'whiteSig' are not allowed
%   since they need to be computed for each bootstrap sample  
%   individually.
%
%ESTIMATE INDEXING CONVENTION: when function icassoEst is run
%each estimate gets a unique, integer label in order of
%appearance. The same order and indexing is used throughout the
%Icasso software. In many functions, one can pick a subset of
%estimates sR by giving vector whose elements refers to this unique
%label.    
%
%SEE ALSO
% icasso
% fastica
% icassoStruct
% icassoExp
% icassoGet
% icassoShow
% icassoResult
%
%When icassoEst is accomplished, use icassoExp to obtain clustering
%results and to store them in sR After this, the results can be
%examined visually using icassoShow. Results and other information
%an be finally retrieved also by functions icassoResult and icassoGet. 

%COPYRIGHT NOTICE
%This function is a part of Icasso software library
%Copyright (C) 2003-2005 Johan Himberg
%
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

% ver 1.1 johan 210704

% Set the Icasso struct

sR=icassoStruct(X); 

% Check compulsatory input arguments
if nargin<3,
  error('At least three input args. required');   
end

%% Check mode.
mode=lower(mode);
switch mode
 case {'randinit','bootstrap','both'}
  ;
 otherwise
  error(['Randomization mode must be ''randinit'', ''bootstrap'' or' ...
	 ' ''both''.']);  
end

sR.mode=mode;

%% Set some values
num_of_args=length(varargin);
whitening='not done';

%% Default values for some FastICA options
fasticaoptions={'g','pow3','approach','symm',...
		'maxNumIterations',100};

%% flag: initial conditions given (default: not given)
isInitGuess=0;

%% flag: elements for whitening given (default: not given) 
isWhitesig=0; isWhitemat=0; isDewhitemat=0;

%% Check varargin & set defaults
fasticaoptions=processvarargin(varargin,fasticaoptions);

num_of_args=length(fasticaoptions);


%% Check fasticaoptions:
i=1;
while i<num_of_args,
  switch fasticaoptions{i}
   case {'approach','firstEig','lastEig','numOfIC','finetune','mu','g','a1','a2',...
	 'stabilization','epsilon','maxNumIterations','maxFinetune','verbose',...
	 'pcaE','pcaD'}
    ; % these are ok
      
   %% Get explicit whitening if given & update flags; note that the
   %% arguments are dropped away from fasticaoptions to avoid
   %% duplicate storage in Icasso result struct & in FastICA input
   
   case 'whiteSig'
    w=fasticaoptions{i+1};
    isWhitesig=1;
    fasticaoptions(i:i+1)=[]; i=i-2; num_of_args=num_of_args-2;
   case 'whiteMat'
    White=fasticaoptions{i+1};
    isWhitemat=1;       
    fasticaoptions(i:i+1)=[]; i=i-2; num_of_args=num_of_args-2;
   case 'dewhiteMat'
    deWhite=fasticaoptions{i+1};
    isDewhitemat=1;
    fasticaoptions(i:i+1)=[]; i=i-2; num_of_args=num_of_args-2;
    
   case {'sampleSize','displayMode','displayInterval','only','interactivePCA'}
    error(['You are not allowed to set FastICA option ''' fasticaoptions{i} ''' in Icasso.']);
    % initGuess depends on mode 
   case 'initGuess'
    switch mode
     case {'randinit','both'}
      error(['FastICA option ''initGuess'' cannot be used in sampling mode ''' ...
	     mode '''.']);
     case 'bootstrap'
      isInitGuess=1;
     otherwise
      error('Internal error!?');
    end
   otherwise
    error(['Doesn''t recognize FastICA option ''' fasticaoptions{i} '''.']);
  end
  % add counter
  i=i+2;
end

%% Whitening:

%% Check if some of whitening arguments have been given: 
if (isWhitesig | isWhitemat | isDewhitemat),
  %% both/bootstrap use each time different whitening... better to
  %% give error
  
  if (strcmp(mode,'bootstrap') | strcmp(mode,'both')),
    error(['FastICA options ''whiteSig'',''whiteMat'',''dewhiteMat'' cannot be' ...
	   ' used in modes ''bootstrap'' and ''both''.']);
  end
  
  %% FastICA expects that all of the three arguments are given (see
  %help fastica): if not, error

  if isWhitesig & isWhitemat & isDewhitemat,
    disp('Using user specified whitening.');
  else
    error(['To prewhiten, each of ''whiteSig'',''whiteMat'',''dewhiteMat'' have to' ...
	   ' be given (see help fastica)']);
  end
else
  % compute whitening for original data
  [w,White,deWhite]=fastica(X,'only','white',fasticaoptions{:});
end

% store whitening for original data:
sR.whiteningMatrix=White;
sR.dewhiteningMatrix=deWhite;


% Icasso uses the same random initial condition for every sample in
% 'bootstrap'. It has to be computed if not given!!

if strcmp(mode,'bootstrap') & ~isInitGuess,
  warning(sprintf('\n\n%s\n\n',['Initial guess not given for mode ''bootstrap'': I will' ...
	   ' set a (fixed) random initial condition']));

  % Randomize init conditions and add to fastica options: this
  % keeps it fixed in every estimation round.
  
  fasticaoptions{end+1}='initGuess';
  fasticaoptions{end+1}=rand(size(White'))-.5;
end

% store options (except whitening which is
% stored separately)

sR.fasticaoptions=fasticaoptions;

%% Compute N times FastICA
k=0; index=[];
for i=1:M,
  %clc; 
  fprintf('\n\n%s\n\n',['Randomization using FastICA: Round ' num2str(i) '/' ...
		    num2str(M)]);
  
  switch mode
   case 'randinit'
    % data is fixed; 
    X_=X;
   case {'bootstrap','both'}
    % Bootstrap and compute whitening for _bootstrapped_ data
    X_=bootstrap(X);
    [w,White,deWhite]=fastica(X_,'only','white',fasticaoptions{:});
   otherwise
    error('Internal error?!');
  end
  
  % Estimate FastICA set displayMode off
  [dummy,A_,W_]=fastica(X_,fasticaoptions{:},...
			'whiteMat',White,'dewhiteMat',deWhite,'whiteSig',w,...
			'sampleSize',1,'displayMode','off');
  
  % Store results if any
  n=size(A_,2);
  if n>0, 
    k=k+1;
    sR.index(end+1:end+n,:)=[repmat(k,n,1), [1:n]'];
    sR.A{k}=A_; sR.W{k}=W_; 
  end
end

function X=bootstrap(X)

N=size(X,2);
index=round(rand(N,1)*N+.5);
X=X(:,index);