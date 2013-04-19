function [RB RAR] = mean_jpredf(gp,x,xt,K_nf,L,Ksy,latent_method,S)
%MEAN_PREDF  Calculates help terms needed in prediction with mean function
%
%  Description
%    [RB, RAR] = MEAN_PREDF(gp,x,xt,K_nf,L,Ksy,latent_method,S)
%
%    Gaussian likelihood:
%      gp              - a GP structure
%      x               - training inputs
%      xt              - test inputs
%      K_nf            - covariance matrix K(x,xt)
%      L               - chol (K(x,x) + sigmaI)
%      Ksy             - L'\(L\y) 
%      S               - [] (empty)
%      latent_method   - gaussian
%
%    EP:
%      gp              - a GP structure
%      x               - training inputs
%      xt              - test inputs
%      K_nf            - covariance matrix K(x,xt)
%      L               - inv(K(x,x) + S^-1)*S^-1
%      Ksy             - inv(K + S^-1)*S^-1*nutilde  
%      S               - diag(tautilde)
%      latent_method   - EP
%
%    Laplace:
%      NOT IMPLEMENTED YET
%
%    Returns the help term RB for the posterior predicative mean
%    and help term RAR for posterior predicative variance.
%
%      RB  = R'*Beta = R'*(inv(B1)*B2) 
%      RAR = R'*inv(B1)*R

%    The vague prior functionalities commented. Uncommenting vague
%    prior rows here doesn't make vague prior compatibible with
%    other functions.

%        See GPstuff doc and (Rasmussen and Williams 2006) page 28 for further
%        explaining.

% Copyright (c) 2010 Tuomas Nikoskinen
% Copyright (c) 2011 Jarno Vanhatalo

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.
  
    % prior assumption for weights, w ~ N(b,B) 
    % b_m = prior mean for weights, B_m prior covariance matrix for weights
    [H,b,B,Hs]=mean_prep(gp,x,xt);
    
    if isequal(latent_method,'gaussian')
        if ~isempty(L)
            if issparse(L)
                KsK = ldlsolve(L,K_nf);                  % inv(C)*K(x,xt)
                KsH = ldlsolve(L,H');                    % inv(C)*H'
            else
                KsK = L'\(L\K_nf);                       % inv(C)*K(x,xt)
                KsH = L'\(L\H');                         % inv(C)*H'
            end
        else
            [nh mh]=size(H);
            KsK=zeros(length(x),length(xt));
            KsH=zeros(length(x),nh);
        end
    elseif isequal(latent_method,'EP')
        KsK = L*(S*K_nf);                        % inv(K + S^-1)*S^-1*(S*K(x,xt)) 
        KsH = L*(S*H');                          % inv(K + S^-1)*S^-1*(S*H')
    elseif isequal(latent_method,'laplace') 
        error('latent method = laplace not implemented yet')
    else
        error('no correct latent method specified')
    end
        
    R = Hs - H*KsK;

%     if gp.mf{1}.p.vague==0        % non-vague prior
        invB = B\eye(size(B));
        B1 = invB + H*KsH;
        B2 = H*Ksy + invB*b;
%     else                          % vague prior
%         B1 = H*KsH;
%         B2 = H*Ksy;
%     end
    
    Beta = B1\B2;
    invAR=B1\R;
    RARapu=R'*invAR;          % Calculate only the necessary.
    
    
    RB  = R'*Beta;                % For predictive mean
    %RAR = sum(RARapu,2);          % For predictive variance
    RAR = RARapu;          % For predictive covariance
    
end