function [gp, indA] = passgp(gp, x, y, varargin)
%PASSGP   Optimize active set and hyperparameters of PASS-GP 
%
%  Description
%    [GP, INDA] = PASSGP(GP, X, Y, OPTIONS)
%    Returns GP structure with hyperparameters optimized according to
%    PASS-GP routine (Henao & Winther, 2012) and active set indices 
%    INDA for X and Y
%
%   OPTIONS is optional parameter-value pair
%      npass - Number of passes the algorithm takes over the whole training
%              data set
%      ninit - Initial active set size
%      nsub  - Number of subsets we process at each pass (in how many
%              parts we divide the whole data)
%      pinc  - Predictive density threshold for inclusion in active set.
%      pdel  - LOO-predictive density threshold for deletion from active
%              set
%      pexc  - Exchange proportion for fixed PASS-GP
%      opt   - Options structure for optimizer
%      fixed - Whether we use fixed size of active set or not. Default
%               'off'
%      display - Whether to display additional info or not. Default 'off'
%      optimn - Whether to optimize always or only after every nth
%               deletion/addition to active set. Default 1 (every time). If
%               given e.g. value 3, optimizes after every 3rd
%               addition/deletion.
%
%  See also
%    GP_SET, LIK_*
%
%  Reference:
%    Ricardo Henao & Ole Winther (2012). Preditive active set selection
%    methods for Gaussian processes. Neurocomputing 80 (2012), 10-18.

% Copyright (c) 2013 Ville Tolvanen

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'PASSGP';
  ip.addRequired('gp',@(x) isstruct(x) || iscell(x));
  ip.addRequired('x', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
  ip.addRequired('y', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
  ip.addParamValue('npass', 3, @(x) isscalar(x) && x > 0)
  ip.addParamValue('ninit', 100, @(x) isscalar(x) && x > 0)
  ip.addParamValue('nsub', 10, @(x) isscalar(x) && x > 0)
  ip.addParamValue('pinc', 0.5, @(x) isscalar(x) && x > 0)
  ip.addParamValue('pdel', 0.99, @(x) isscalar(x) && x > 0)
  ip.addParamValue('pexc', 0.1, @(x) isscalar(x) && x > 0)
  ip.addParamValue('opt', [], @(x) isstruct(x))
  ip.addParamValue('fixed', 'off', @(x) ismember(x,{'on','off'}))
  ip.addParamValue('display', 'off', @(x) ismember(x,{'on','off'}))
  ip.addParamValue('optimn', 1', @(x) isscalar(x) && x > 0 && rem(10*x,2)==0)
  ip.addParamValue('z', [], @(x) isreal(x) && all(isfinite(x(:))))
  ip.parse(gp, x, y, varargin{:});
  opt=ip.Results.opt;
  fixed=ip.Results.fixed;
  options.z = ip.Results.z;
  npass=ip.Results.npass;
  ninit=ip.Results.ninit;
  nsub=ip.Results.nsub;
  pinc=ip.Results.pinc;
  pdel=ip.Results.pdel;
  pexc=ip.Results.pexc;
  display=ip.Results.display;
  optimn=ip.Results.optimn;
  
  if isequal(display,'on')
    display=1;
  else
    display=0;
  end
  if isequal(fixed, 'on')
    fixed=1;
  else
    fixed=0;
  end  
  
  if ninit > size(x,1)
    error('Initial active set must be subset of original data');
  end
  if nsub > (size(x,1) - ninit)
    error('nsub must be lower than size(x,1) - ninit');
  end
  
  [n,nin]=size(x);
  
  % Initial active set
  
  indA=sort(randperm(n, ninit),'ascend');
  % Inclusions/deletions per iteration for fixed pass-gp
  nexc=floor(ninit*pexc);
  iter=optimn-1;
  for i=1:npass
    if display
      fprintf('Pass %d / %d.\n', i, npass)
    end
    [tmp,indSub]=cvit(n, nsub, floor(10*rand(1)));
    for j=1:nsub      
      iter=iter+1;
      inds=indSub{j};
      % Remove indices that are already in active set
      inds(ismember(inds,indA))=[];
      
      if iter==optimn
        % Optimize hyperparameters
        gp = gp_optim(gp, x(indA,:), y(indA), 'opt', opt, options);
        iter=0;
      end
      
      % Calculate weights for active set inputs (loo predictive densities)
      [tmp,tmp,lpyt]=gp_loopred(gp,x(indA,:),y(indA), options);
      
      % Remove active set indices according to removal rule      
      if ~fixed
        indA(find(exp(lpyt)>pdel))=[];
      else
        [tmp,ii]=sort(lpyt, 'descend');
        indA(ii(1:nexc))=[];
      end
      % Calculate weights for inputs not in active set (predictive density)
      [tmp,tmp,lpyt]=gp_pred(gp, x(indA,:), y(indA), x(inds,:), 'yt', y(inds), options);
      
      % Add indices to active set according to addition rule
      if ~fixed
        ind=find(exp(lpyt)<pinc)';
      else
        [tmp,ii]=sort(lpyt, 'ascend');
        ind=ii(1:nexc);
      end
      indA=[indA inds(ind)];
        
    end
    
      
  end
  
end

