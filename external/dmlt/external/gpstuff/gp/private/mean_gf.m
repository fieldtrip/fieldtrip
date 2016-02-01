function [dMNM trA HinvC] = mean_gf(gp,x,Ky,invKy,DKff,Stildesqroot,y,latent_method)
%MEAN_GF  Calculates help terms needed in gradient calculation with
%         mean function
%
%  Description
%    [dMNM trA dyKy dyCy trAv]
%       = MEAN_GF(gp,x,C,invC,DKff,Stildesqroot,y,latent_method)
%    Gaussian likelihood:
%      gp              - a GP structure
%      x               - training inputs
%      Ky              - cov. matrix K(x,x) + sigma*I
%      invKy           - inv(Ky)
%      DKff            - d Ky / d th, (th = hyperparameters)
%      Stildesqroot    - [] (empty)
%      y               - noisy latent values
%      latent_method   - gaussian
%
%    EP:
%      gp              - a GP structure
%      x               - training inputs
%      Ky              - cov. matrix K(x,x) + sigma*I
%      invKy           - inv(Ky + S^-1)*S^-1
%      DKff            - d Ky / d th, (th = hyperparameters)
%      Stildesqroot    - sqrt( diag(tautilde) )
%      y               - S*mutilde
%      latent_method   - EP
%
%    Laplace:
%      NOT IMPLEMENTED YET
%
%    Returns the help terms dMNM and trA
%
%      dMNM = d M'*inv(N)*M / d th
%      trA  = d log|A| / dth
%      dyKy = d y'*Ky*y/ d th
%      dyCy = d y'*C*y / d th
%      trAv = d log|Av|/ d th

%    The vague prior functionalities commented. The function should
%    return help terms dyKy dyCy trAv with vague prior. 
%    Uncommenting vague prior rows here doesn't make vague prior
%    compatibible with other functions.

%    See GPstuff doc and (Rasmussen and Williams 2006) for further
%    details.

% Copyright (c) 2010 Tuomas Nikoskinen
% Copyright (c) 2011 Jarno Vanhatalo

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.


dMNM = cell(1,length(DKff));
trA  = cell(1,length(DKff));
%         dyKy = cell(1,length(DKff));
%         dyCy = cell(1,length(DKff));      % with vauge prior
%         trAv = cell(1,length(DKff));

% prior assumption for weights, w ~ N(b,B)
% b_m = prior mean for weights, B_m prior covariance matrix for weights
[H,b_m,B_m]=mean_prep(gp,x,[]);

% help arguments
if issparse(Ky)
    % in case of CS covariance function invKy contains the LDL Cholesky decomposition of the covariance function
    KH = ldlsolve(invKy, H');   
    HinvC = KH';
else
    HinvC = H*invKy;
    N = Ky + H'*B_m*H;
end

%         if gp.mf{1}.p.vague==0   % non-vague prior

% help arguments that don't depend on DKff; non-vague p
if isequal(latent_method,'gaussian')
    % help arguments with gaussian latent method
    M = H'*b_m-y;
    
    LB = chol(B_m);
    LA = chol(LB\(LB'\eye(size(B_m))) + HinvC*H');
    invAt=LA\(LA'\eye(size(LA)));
    
    if issparse(Ky)
        invNM = ldlsolve(invKy, M) - KH*(LA\(LA'\(KH'*M)));
    else
        invNM = N\M;
    end
    
    dMNM = invNM;
    trA = invAt;
        
elseif isequal(latent_method,'EP')
    % help arguments with EP latent method
    S=Stildesqroot.^2;
    M = S*H'*b_m-y;                                     % M is now (S*H'b_m - S*mutilde)
    HKH = HinvC*S*H';                                   % inv(Ky + S^-1)*S^-1*S*H'
    A = B_m\eye(size(B_m)) + HKH;
    invAt=A\eye(size(A));
    invAt=invAt';
    B_h = eye(size(N)) + Stildesqroot*N*Stildesqroot;
    L_m=chol(B_h,'lower');
    zz=Stildesqroot*(L_m'\(L_m\(Stildesqroot*N)));
    invN = eye(size(zz))-zz;                            % inv(Ky + H'*B_m*H + S^-1)*S^-1
    
    % Calculate the arguments which are to be returned
    for i2=1:length(DKff)
        dA = -1*HinvC*S*DKff{i2}*(invKy*S*H');          % d A / d th
        trA{i2} = sum(invAt(:).*dA(:));                 % d log(|A|) / dth
        dMNM{i2} = M'*(S^(-1)*invN*S*DKff{i2}*invN*M);  % with EP the d M'*N*M / d th
    end
    
    
elseif isequal(latent_method,'laplace')
    error('latent method = laplace not implemented yet')
else
    error('no correct latent method specified')
end


%         else  % vague prior

%             if isequal(latent_method,'gaussian')
%                 % help arguments that don't depend on DKff; vague p
%                 HKH = HinvC*H';
%                 A     = HKH;
%                 AH    = A\H;
%                 invAt = A\eye(size(A));
%                 invAt = invAt';
%                 G     = H'*AH*invKy*y;
%                 b     = invKy*y;
%
%                 for i2 = 1:length(DKff)
%                     % help arguments that depend on DKff; vague p
%                     dyKy{i2} = b'*(DKff{i2}*b);            % d y'*Kyâ?»*y / d th
%                     dA  = -1*HinvC*DKff{i2}*HinvC';        % d A / d th
%                     trAv{i2} = sum(invAt(:).*dA(:));       % d log(|A|)/dth = trace(inv(A) * dA/dth)
%                     P   = invKy*DKff{i2}*invKy;
%
%                     dyCy1 = y'*P*G;
%                     dyCy3 = -G'*P*G;
%                     dyCy{i2} = 2*dyCy1 + dyCy3;          % d y'*C*y /d th
%                 end
%             else
%                 error('vague prior only for gaussian latent method at the moment')
%             end
%         end


end