function [b, iCv, iC] = gp_looprep(gp, x, y)
%GP_LOOPREP  Calculates help terms needed in Gaussian GP LOO
%
%  Description
%    [b, iCv] = gp_looprep(gp, x, y, varargin) takes GP structure,
%    training input and output and returns following terms
%      b = C\y
%      iCv = diag(inv(C))
%    [b, iCv, iC] = gp_looprep(gp, x, y, varargin) returns additionally
%      iC = inv(C)
%    In case of CS, FIC, PIC and CS+FIC, these terms are computed
%    using appropriate efficient sparse matrix and inverse lemma
%    computations.
% 
  
% Copyright (c) 2012 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.
  
  
  switch gp.type
    case 'FULL' 
      % FULL GP (and compact support GP)
      [K, C] = gp_trcov(gp, x);
      
      if issparse(C)
        if nargout==3
          iC = full(spinv(C));        % sparse inverse
          iCv = full(diag(iC));       % the diagonal of sparse inverse
        else
          iCv = full(diag(spinv(C))); % the diagonal of sparse inverse
        end
        [LD, notpositivedefinite] = ldlchol(C);
        if notpositivedefinite
          b = NaN;
          iCv = NaN;
          iC = NaN;
          return
        end
        b = ldlsolve(LD,y);
      else
        if nargout==3
          iC = inv(C);                % full inverse
          iCv = diag(inv(C));         % the diagonal of full inverse
        else
          iCv = diag(inv(C));         % the diagonal of full inverse
        end
        b = C\y;
      end
      
    case 'FIC' 
      % FIC
      % Use inverse lemma for FIC low rank covariance matrix approximation
      % Code adapated from gp_pred
      
      u = gp.X_u;
      m = size(u,1);
      % Turn the inducing vector on right direction
      if size(u,2) ~= size(x,2)
        u=u';
      end
      [Kv_ff, Cv_ff] = gp_trvar(gp, x);  % 1 x f  vector
      K_fu = gp_cov(gp, x, u);   % f x u
      K_uu = gp_trcov(gp, u);     % u x u, noiseles covariance K_uu
      [Luu, notpositivedefinite] = chol(K_uu,'lower');
      if notpositivedefinite
        b = NaN;
        iCv = NaN;
        iC = NaN;
        return
      end
      B=Luu\(K_fu');
      Qv_ff=sum(B.^2)';
      Lav = Cv_ff-Qv_ff;   % 1 x f, Vector of diagonal elements
      
      % iLaKfu = diag(inv(Lav))*K_fu = inv(La)*K_fu
      iLaKfu = zeros(size(K_fu));  % f x u,
      n=size(x,1);
      for i=1:n
        iLaKfu(i,:) = K_fu(i,:)./Lav(i);  % f x u
      end
      A = K_uu+K_fu'*iLaKfu;
      A = (A+A')./2;
      L = iLaKfu/chol(A);
      
      if nargout==3
        iC = diag(1./Lav) - L*L'; % Cv = inv(C);
        iCv = diag(iC); %           iCv = diag(inv(C));
      else
        iCv = 1./Lav - sum(L.^2,2); % iCv = diag(inv(C));
      end
      b = y./Lav - L*(L'*y);      % b = C\y;
      
    case {'PIC' 'PIC_BLOCK'}
      % PIC
      % Use inverse lemma for PIC low rank covariance matrix approximation
      % Code adapated from gp_pred
      
      u = gp.X_u;
      ind = gp.tr_index;
      if size(u,2) ~= size(x,2)
        % Turn the inducing vector on right direction
        u=u';
      end

      % Calculate some help matrices
      [Kv_ff, Cv_ff] = gp_trvar(gp, x);  % 1 x f  vector
      K_fu = gp_cov(gp, x, u);         % f x u
      K_uu = gp_trcov(gp, u);    % u x u, noiseles covariance K_uu
      Luu = chol(K_uu)';

      % Evaluate the Lambda (La) for specific model
      % Q_ff = K_fu*inv(K_uu)*K_fu'
      % Here we need only the diag(Q_ff), which is evaluated below
      B=Luu\K_fu';
      iLaKfu = zeros(size(K_fu));  % f x u
      for i=1:length(ind)
        Qbl_ff = B(:,ind{i})'*B(:,ind{i});
        [Kbl_ff, Cbl_ff] = gp_trcov(gp, x(ind{i},:));
        La{i} = Cbl_ff - Qbl_ff;
        iLaKfu(ind{i},:) = La{i}\K_fu(ind{i},:);    
      end
      A = K_uu+K_fu'*iLaKfu;
      A = (A+A')./2;            % Ensure symmetry
      L = iLaKfu/chol(A);

      % From this on evaluate the prediction
      % See Snelson and Ghahramani (2007) for details
      n=size(y,1);
      iCv=zeros(n,1);
      if nargout==3
        iC=zeros(n,n);
      end
      for i=1:length(ind)
        if nargout==3
          iC(ind{i},ind{i}) = inv(La{i});
        else
          iCv(ind{i},:) = diag(inv(La{i}));
        end
        b(ind{i},:) = La{i}\y(ind{i},:);
      end
      if nargout==3
        iC = iC - L*L';          % iC  = inv(C);
        iCv = diag(iC);          % iCv = diag(inv(C));
      else
        iCv = iCv - sum(L.^2,2); % iCv = diag(inv(C));
      end
      b = b - L*(L'*y);        % b = C\y;

    case 'CS+FIC'
      % CS+FIC
      % Use inverse lemma for CS+FIC
      % Code adapated from gp_pred
      
      u = gp.X_u;
      if size(u,2) ~= size(x,2)
        % Turn the inducing vector on right direction
        u=u';
      end
      
      n = size(x,1);
      m = size(u,1);
      ncf = length(gp.cf);
      
      % Indexes to all non-compact support and compact support covariances.
      cf1 = [];
      cf2 = [];

      % Loop through all covariance functions
      for i1 = 1:ncf        
        if ~isfield(gp.cf{i1},'cs') 
          % Non-CS covariances
          cf1 = [cf1 i1];
        else
          % CS-covariances
          cf2 = [cf2 i1];           
        end
      end
      
      % First evaluate needed covariance matrices
      % v defines that parameter is a vector
      [Kv_ff, Cv_ff] = gp_trvar(gp, x, cf1); % f x 1  vector    
      K_fu = gp_cov(gp, x, u, cf1);          % f x u
      K_uu = gp_trcov(gp, u, cf1);    % u x u, noiseles covariance K_uu
      K_uu = (K_uu+K_uu')./2;         % ensure the symmetry of K_uu

      Luu  = chol(K_uu)';
      % Evaluate the Lambda (La)
      % Q_ff = K_fu*inv(K_uu)*K_fu'
      B=Luu\(K_fu');       % u x f
      Qv_ff=sum(B.^2)';
      Lav = Cv_ff-Qv_ff;   % f x 1, Vector of diagonal elements

      K_cs = gp_trcov(gp,x,cf2);
      La = sparse(1:n,1:n,Lav,n,n) + K_cs;

      iLaKfu = La\K_fu;
      A = K_uu+K_fu'*iLaKfu;
      A = (A+A')./2;     % Ensure symmetry
      L = iLaKfu/chol(A);
   
      if nargout==3
        iC = inv(La) - L*L'; % iCv = diag(inv(C));
        iCv = diag(iC)     ; % iCv = diag(inv(C));
      else
        iCv = diag(inv(La)) - sum(L.^2,2); % iCv = diag(inv(C));
      end
      b = La\y - L*(L'*y);         % b = C\y;
      
    case 'SSGP' % SSGP
      error('GP_LOOPRED is not yet implemented for SSGP and Gaussian likelihood')
      
    otherwise
      error('Unknown type of Gaussian process')
  end

end
