%Example:
% config = gmmreg_load_config('./fish_full.ini');
% config.motion = 'grbf';
% config.init_param = zeros(25,2);
% [fp,fm] = gmmreg_cpd(config);
% DisplayPoints(fm,config.scene,2);


function [param, model] = gmmreg_cpd(config)
%%=====================================================================
%% $Author: bjian $
%% $Date: 2008/12/07 00:43:34 $
%% $Revision: 1.2 $
%%=====================================================================

% todo: use the statgetargs() in statistics toolbox to process parameter name/value pairs
% Set up shared variables with OUTFUN
if nargin<1
    error('Usage: gmmreg_cpd(config)');
end
[n,d] = size(config.model); % number of points in model set
if (d~=2)&&(d~=3)
    error('The current program only deals with 2D or 3D point sets.');
end

tic
model = config.model;
scene = config.scene;
ctrl_pts = config.ctrl_pts;
sigma = config.init_sigma;
anneal_rate = config.anneal_rate;
outliers = config.outliers;
lambda = config.lambda;
max_iter = config.max_iter;
max_em_iter = config.max_em_iter;

tol = config.tol;
EMtol = config.emtol;

[n,d] = size(ctrl_pts);
[m,d] = size(model);


% Rescaling and shifting to the origin
[model, centroid, scale] = cpd_normalize(model);
[ctrl_pts, centroid, scale] = cpd_normalize(ctrl_pts);
[scene, centroid, scale] = cpd_normalize(scene);
model0 = model;


switch lower(config.motion)
    case 'tps'
        K = tps_compute_kernel(ctrl_pts, ctrl_pts);
        Pn = [ones(n,1) ctrl_pts];
        PP = null(Pn');  % or use qr(Pn)
        kernel = PP'*K*PP;
        U = tps_compute_kernel(model, ctrl_pts);
        Pm = [ones(m,1) model];
        [q,r]   = qr(Pm);
        Q1      = q(:, 1:d+1);
        Q2      = q(:, d+2:m);
        R       = r(1:d+1,1:d+1);
        TB = U*PP;
        QQ = Q2*Q2';
        A = inv(TB'*QQ*TB + lambda*kernel)*TB'*QQ;
        basis = [Pm U*PP];
    case 'grbf'
        beta = config.beta;
        basis = cpd_G(model,ctrl_pts,beta);
        kernel = cpd_G(ctrl_pts,ctrl_pts,beta);
        %A = inv(basis'*basis+lambda*sigma*sigma*kernel)*basis';
        A = basis'*basis+lambda*sigma*sigma*kernel;
    otherwise
        error('Unknown motion model');
end % end of switch

param = config.init_param;  % it should always be of size n*d
model = model + basis*param;

%it_total = 1;
%flag_stop = 0;

iter=0; E=1; ntol=tol+10;

while (iter < max_iter) && (ntol > tol)
    EMiter=0; EMtol=tol+10;  % repeat at each termperature.
    model_old = model;
    while (EMiter < max_em_iter) && (EMtol > tol)
        %disp(sprintf('EMiter=%d',EMiter));
        %disp(sprintf('E=%f',E));
        E_old = E;
        % E-step: Given model&scene, update posterior probability matrix P.
        [P,Eu] = cpd_P(model, scene, sigma, outliers);
        % M-step: Given correspondence, solve warp parameter.
        %
        switch lower(config.motion)
            case 'tps'
                tps = param(d+2:end,:);
                E = Eu + lambda/2*trace(tps'*kernel*tps); % CPD energy function.
                s=sum(P,2);
                %P=P./repmat(s,1,m); % normalization such that each row sums to 1
                motion = P * scene - model0;
                tps = A*motion;
                affine = inv(R)*Q1'*(motion-TB*tps);
                param = [affine; tps];
            case 'grbf'
                E = Eu + lambda/2*trace(param'*kernel*param); % CPD energy function.
                %param = A\(basis'*motion);
                dP=spdiags(sum(P,2),0,m,m); % precompute diag(P)

                % with ctrl_pts
                param =(basis'*dP*basis+lambda*sigma^2*kernel)\(basis'*(P*scene-dP*model0));
                % when ctrl_pts is same as model0
                % param =(dP*basis+lambda*sigma^2*eye(m))\(P*scene-dP*model0);
        end
        % update model
        model = model0 + basis*param;
        EMtol = norm((E_old-E)/E_old);
        EMiter = EMiter + 1;
    end  % end of iteration/perT
    % Anneal
    sigma = sigma * anneal_rate;
    iter = iter + 1;
    ntol = norm(model_old - model);
end % end of annealing.
toc
model = cpd_denormalize(model, centroid, scale);


end % end of function



function [P, E] = cpd_P(x, y, sigma, outliers)

if nargin<3, error('cpd_P.m error! Not enough input parameters.'); end;
if ~exist('outliers','var') || isempty(outliers), outliers = 0; end;

k=-2*sigma^2;
[n, d]=size(x);[m, d]=size(y);

P=repmat(x,[1 1 m])-permute(repmat(y,[1 1 n]),[3 2 1]);

P=squeeze(sum(P.^2,2));
P=P/k;
P=exp(P);

% compute column sums -> s
if outliers
  Pn=outliers*(-k*pi)^(0.5*d)*ones(1,m);
  s=sum([P;Pn]);
else
  s=sum(P);
end

if nnz(s)==numel(s)
    E=-sum(log(s));  % log-likelihood
    P=P./repmat(s,n,1); % normalization such that each column sums to 1
    s=sum(P,2);
    P=P./repmat(s,1,m); % normalization such that each row sums to 1
else
    P=[];E=[];
end

end

function G=cpd_G(x,y,beta)

if nargin<3, error('cpd_G.m error! Not enough input parameters.'); end;

k=-2*beta^2;
[n, d]=size(x); [m, d]=size(y);

G=repmat(x,[1 1 m])-permute(repmat(y,[1 1 n]),[3 2 1]);
G=squeeze(sum(G.^2,2));
G=G/k;
G=exp(G);

end


function  [X, centroid, scale] = cpd_normalize(x)

[n, d]=size(x);
centroid=mean(x);
x=x-repmat(centroid,n,1);
scale=sqrt(sum(sum(x.^2,2))/n);
X=x/scale;

end


function x =cpd_denormalize(X, centroid, scale)

[m, d]=size(X);
x=X*scale;         % scale
x=x+repmat(centroid,m,1); % move
end

