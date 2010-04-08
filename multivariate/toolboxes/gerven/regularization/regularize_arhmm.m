function [path,diagnostics] = regularize_arhmm(cfg,data)
% REGULARIZE_ARHMM learns the regularization path for autoregressive HMM
% models
%
% [path,diagnostics] = regularize_arhmm(cfg,data)
%
% cfg contains the configuration options:
%  
% cfg.maxtime = 1e5; % maximum number of seconds to run the algorithm
% cfg.maxiter = 1e4; % maximum number of iterations
% cfg.cvar = 1; the column for the class (or regression) variable
% cfg.nclasses; the number of classes (automatically determined if unspecified)
% cfg.epsilon = 1; % step size
% cfg.tolerance = 1e-5; % error tolerance
% cfg.p = 2; % L1/LP norm
% cfg.groups; % grouping of the parameter vector model (default: no grouping)
%
% data contains the examples X variables X time (2 slice) data for the model
%
% we assume all continuous features!
%
% Copyright (C) 2008, Marcel van Gerven
%
% $Log: regularize_arhmm.m,v $
% Revision 1.3  2008/03/20 14:38:47  marvger
% regular update
%
% Revision 1.2  2008/03/17 12:47:41  marvger
% regular update
%
% Revision 1.1  2008/03/10 09:12:29  marvger
% add arhmm regularization code
%
% Revision 1.2  2008/03/01 16:38:28  marvger
% removed somes files
%
% Revision 1.1.1.1  2008/02/27 14:42:55  roboos
% Marcel van Gerven, version 27 Feb 2008
%

%% initialization

cfg.recompute = true; % recompute the lambdas

if ~isfield(cfg,'cvar'), cfg.cvar = 1; end
 
if ~isfield(cfg,'nclasses')

    if any((mod(data(:,cfg.cvar),1)))
        cfg.nclasses = 1;
    else
        cfg.nclasses = length(unique(data(:,cfg.cvar,1)));
    end
    
end

if ~isfield(cfg,'p'), cfg.p = 2; end

if ~(cfg.p == 2 || cfg.p == Inf), error('cfg.p can only be 2 or Inf'); end

% determine features
cfg.features = [1:(cfg.cvar-1) (cfg.cvar+1):size(data,2)]; 
cfg.nfeatures = length(cfg.features);

%% prepare objective function

% space for mean, standard deviation, and linear terms per feature
model = rand(cfg.nfeatures,cfg.nclasses*(cfg.nfeatures+2)); 

%% data is split up according to classes for faster processing

cfg.data = cell(1,cfg.nclasses);
for c=1:cfg.nclasses, cfg.data{c} = data(data(:,cfg.cvar) == c,:,:); end

%% specify groups

% index: <= cfg.nfeatures is mean/sigma >cfg.nfeatures is beta
cfg.groups = cell(1,((cfg.nfeatures+1)*cfg.nfeatures));

% groups for mu and sigma
for m=1:cfg.nfeatures 
    x = zeros(size(model)); x(m,1:(2*cfg.nclasses)) = 1;
    cfg.groups{m} = find(x)';
end

idx = cfg.nfeatures+1;
for i=1:cfg.nfeatures
    for j=1:cfg.nfeatures
        x = zeros(size(model)); x(i,((j+1)*cfg.nclasses+1):((j+2)*cfg.nclasses)) = 1;
        cfg.groups{idx} = find(x)';
        idx = idx + 1;
    end
end

% everything is regularized
cfg.isregularized = true(1,length(cfg.groups));

%% specify functions

cfg.loss = @arhmm_loss;
cfg.loss_grad = @arhmm_loss_grad;
cfg.reg = @arhmm_reg;
cfg.reg_grad = @arhmm_reg_grad;

cfg.stability = @arhmm_stability;

% keep track of total loss
cfg.optimize = @arhmm_optimize;

%% initialize means and standard deviations

tdiff = Inf;
while tdiff > 0
    
    tdiff = 0;
    for i=1:cfg.nfeatures

        params = model(i,:);
        mu = params(1:cfg.nclasses);
        beta = reshape(params(((2*cfg.nclasses)+1):end),[cfg.nclasses cfg.nfeatures]);

        % reevaluate mean
        tmodel = model(i,1:cfg.nclasses);
        model(i,1:cfg.nclasses) = 0;
        for c=1:cfg.nclasses
            model(i,1:cfg.nclasses) = model(i,1:cfg.nclasses) + sum(cfg.data{c}(:,cfg.features(i),1) - cfg.data{c}(:,cfg.features,2)*beta(c,:)');
        end
        model(i,1:cfg.nclasses) = model(i,1:cfg.nclasses)./size(data,1);
        tdiff = tdiff + sum(abs(tmodel - model(i,1:cfg.nclasses)));
                
        % reevaluate sigma
        tmodel = model(i,cfg.nclasses + (1:cfg.nclasses));
        model(i,cfg.nclasses + (1:cfg.nclasses)) = 0;
        for c=1:cfg.nclasses
            model(i,cfg.nclasses + (1:cfg.nclasses)) = model(i,cfg.nclasses + (1:cfg.nclasses)) + sum((cfg.data{c}(:,cfg.features(i),1) - cfg.data{c}(:,cfg.features,2)*beta(c,:)' - mu(c)).^2);
        end
        model(i,cfg.nclasses + (1:cfg.nclasses)) = model(i,cfg.nclasses + (1:cfg.nclasses))./size(data,1);
        tdiff = tdiff + sum(abs(tmodel - model(i,cfg.nclasses + (1:cfg.nclasses))));
        
        % reevaluate all betas
        for j=1:cfg.nfeatures

            tmodel = model(i,(j+1)*cfg.nclasses + (1:cfg.nclasses));
            model(i,(j+1)*cfg.nclasses + (1:cfg.nclasses)) = 0;
            alpha = beta; alpha(:,j) = 0;
            for c=1:cfg.nclasses
                model(i,(j+1)*cfg.nclasses + (1:cfg.nclasses)) = model(i,(j+1) * cfg.nclasses + (1:cfg.nclasses)) + sum((cfg.data{c}(:,cfg.features(i),1) - cfg.data{c}(:,cfg.features,2)*alpha(c,:)' - mu(c)).*cfg.data{c}(:,cfg.features(j),1));
            end
            model(i,(j+1)*cfg.nclasses + (1:cfg.nclasses)) = model(i,(j+1)*cfg.nclasses + (1:cfg.nclasses)) ./ (size(data,1)*sum(data(:,cfg.features(j),1).^2));
            tdiff = tdiff + sum(abs(tmodel - model(i,(j+1)*cfg.nclasses + (1:cfg.nclasses))));
        end

    end
end


%% run algorithm

global LOSS;

LOSS = 0;
for f=1:cfg.nfeatures, LOSS = LOSS + cfg.loss(cfg,model,data,f); end

cfg.iloss = LOSS; % specify initial loss

[path,diagnostics] = regularize(cfg,model,data);


function loss = arhmm_loss(cfg,model,data,gg)

% set to feature
g = mod(gg,cfg.nfeatures); if ~g, g = cfg.nfeatures; end

mu = model(g,1:cfg.nclasses)'; sigma = model(g,(cfg.nclasses+1):(2*cfg.nclasses))';
beta = reshape(model(g,(2*cfg.nclasses+1):end),[cfg.nclasses cfg.nfeatures]);

% return infinite loss if sigma becomes negative
if any(sigma <= 0)
    loss = Inf;
    return
end

% compute mean plus linear component
bmu = mu(data(:,cfg.cvar,2))+sum(data(:,cfg.features,1).*beta(data(:,cfg.cvar,2),:),2);

loss = sum((data(:,cfg.features(g),2) - bmu).^2 ./ (2 * sigma(data(:,cfg.cvar)).^2) + log(sigma(data(:,cfg.cvar,2)).*sqrt(2*pi)));


function grad = arhmm_loss_grad(cfg,model,data,gg)
% gradient of the loss term w.r.t. group gg

% set to feature
g = mod(gg,cfg.nfeatures); if ~g, g = cfg.nfeatures; end

mu = model(g,1:cfg.nclasses); sigma = model(g,(cfg.nclasses+1):(2*cfg.nclasses));
beta = reshape(model(g,(2*cfg.nclasses+1):end),[cfg.nclasses cfg.nfeatures]);

if gg > cfg.nfeatures % beta
    
    gin = ceil((gg - cfg.nfeatures)/cfg.nfeatures); % input feature for betas

    grad = zeros(1,cfg.nclasses);
    
    for c=1:cfg.nclasses
        grad(c) = sum(((mu(c)+cfg.data{c}(:,cfg.features,1)*beta(c,:)') - cfg.data{c}(:,cfg.features(g),2))'*cfg.data{c}(:,cfg.features(gin),1) ./ sigma(c).^2);
    end
        
else % mu and sigma
    
    gradmu = zeros(cfg.nclasses,1);
    gradsigma = zeros(cfg.nclasses,1);

    for c=1:cfg.nclasses
        
        gradmu(c) = sum(((mu(c)+cfg.data{c}(:,cfg.features,1)*beta(c,:)') - cfg.data{c}(:,cfg.features(g),2)) ./ sigma(c).^2);      
        gradsigma(c) = sum( - (cfg.data{c}(:,cfg.features(g),2) - (mu(c)+cfg.data{c}(:,cfg.features,1)*beta(c,:)')).^2 ./ sigma(c).^3 + 1/sigma(c));

    end

    grad = [gradmu; gradsigma]';
end


function reg = arhmm_reg(cfg,model,gg)

g = mod(gg,cfg.nfeatures); if ~g, g = cfg.nfeatures; end

if cfg.isregularized(gg) && cfg.selected(gg)

    if gg > cfg.nfeatures % betas

        beta = model(cfg.groups{gg});

        reg = norm(beta - mean(beta),cfg.p);

    else % mu and sigma

        mu = model(g,1:cfg.nclasses);
        sigma = model(g,(cfg.nclasses+1):(2*cfg.nclasses));

        reg = norm([mu - mean(mu) sigma - mean(sigma)],cfg.p);
    end
else
    reg = 0;
end

function grad = arhmm_reg_grad(cfg,model,gg)
% gradient of the regularization term w.r.t. group g

g = mod(gg,cfg.nfeatures); if ~g, g = cfg.nfeatures; end

if gg > cfg.nfeatures % betas
      
    if cfg.isregularized(gg)
        
        bmbetak = (model(cfg.groups{gg}) - mean(model(cfg.groups{gg})));
        
        % f is the denominator in the partial derivative of the regularization term
        f = sum(bmbetak.^cfg.p).^((cfg.p-1)/cfg.p);

        % compute second component of partial derivative
        betasgn = abs(bmbetak).^(cfg.p-1) .* sign(bmbetak);        
        betasum = sum(betasgn) / cfg.nclasses;

        % deal with betas
        grad = (betasgn - betasum) ./ f;
        
    else
        grad = zeros(1,cfg.nclasses);
    end
else

    if cfg.isregularized(gg)
        
        mmuk = (model(g,1:cfg.nclasses) - mean(model(g,1:cfg.nclasses)));
        msigmak = (model(g,(cfg.nclasses+1):(2*cfg.nclasses)) - mean(model(g,(cfg.nclasses+1):(2*cfg.nclasses))));
 
        % f is the denominator in the partial derivative of the
        % regularization term
        f = (sum(mmuk.^cfg.p) + sum(msigmak.^cfg.p)).^((cfg.p-1)/cfg.p);

        % compute second component of partial derivative
        musgn = abs(mmuk).^(cfg.p-1) .* sign(mmuk);
        musum = sum(musgn) / cfg.nclasses;
        
        sigmasgn = abs(msigmak).^(cfg.p-1) .* sign(msigmak);
        sigmasum = sum(sigmasgn) / cfg.nclasses;

        % deal with means
        grad = [(musgn - musum) ./ f (sigmasgn - sigmasum) ./f];

    else
        grad = zeros(1,2*cfg.nclasses);
    end
        
end

% prevent problems at the constrained solution
grad(isnan(grad)) = 0;


function lambdas = arhmm_stability(cfg,model,data)

lambdas = zeros(1,cfg.ngroups);

if cfg.p == 2

    for g=1:cfg.ngroups

        % we recompute the stability conditions each time, so we need to
        % set the current group to the constrained minimum

        tmp = model;
        if g > cfg.nfeatures % reset beta

            i = ceil((g - cfg.nfeatures)/cfg.nfeatures); % the target feature
            j = mod(g - cfg.nfeatures,cfg.nfeatures); % the source feature
            if j == 0, j = cfg.nfeatures; end
            
            params = tmp(i,:);
            mu = params(1:cfg.nclasses);
            sigma = params(cfg.nclasses+(1:cfg.nclasses));
            beta = reshape(params(((2*cfg.nclasses)+1):end),[cfg.nclasses cfg.nfeatures]);

            % reevaluate beta

            tmp(i,(j+1)*cfg.nclasses + (1:cfg.nclasses)) = 0;
            alpha = beta; alpha(:,j) = 0;
            for c=1:cfg.nclasses
                tmp(i,(j+1)*cfg.nclasses + (1:cfg.nclasses)) = tmp(i,(j+1)*cfg.nclasses + (1:cfg.nclasses)) + sum(sigma(c).^(-2)*(cfg.data{c}(:,cfg.features(i),1) - cfg.data{c}(:,cfg.features,2)*alpha(c,:)' - mu(c)).*cfg.data{c}(:,cfg.features(j),1));
            end
            tmp(i,(j+1)*cfg.nclasses + (1:cfg.nclasses)) = tmp(i,(j+1)*cfg.nclasses + (1:cfg.nclasses)) ./ (size(data,1)*sum(sigma(c).^(-2)*data(:,cfg.features(j),1).^2));

        else % reset mu and sigma

            i = g;

            tdiff = Inf;
            while tdiff > 0

                tdiff = 0;

                params = tmp(i,:);
                mu = params(1:cfg.nclasses);
                beta = reshape(params(((2*cfg.nclasses)+1):end),[cfg.nclasses cfg.nfeatures]);

                % reevaluate mean; note that sigma must still be
                % constrained in this case
                tmodel = tmp(i,1:cfg.nclasses);
                tmp(i,1:cfg.nclasses) = 0;
                for c=1:cfg.nclasses
                    tmp(i,1:cfg.nclasses) = tmp(i,1:cfg.nclasses) + sum(cfg.data{c}(:,cfg.features(i),1) - cfg.data{c}(:,cfg.features,2)*beta(c,:)');
                end
                tmp(i,1:cfg.nclasses) = tmp(i,1:cfg.nclasses)./size(data,1);
                tdiff = tdiff + sum(abs(tmodel - tmp(i,1:cfg.nclasses)));

                % reevaluate sigma; note that mu must still be constrained
                % in this case
                tmodel = tmp(i,cfg.nclasses + (1:cfg.nclasses));
                tmp(i,cfg.nclasses + (1:cfg.nclasses)) = 0;
                for c=1:cfg.nclasses
                    tmp(i,cfg.nclasses + (1:cfg.nclasses)) = tmp(i,cfg.nclasses + (1:cfg.nclasses)) + sum((cfg.data{c}(:,cfg.features(i),1) - cfg.data{c}(:,cfg.features,2)*beta(c,:)' - mu(c)).^2);
                end
                tmp(i,cfg.nclasses + (1:cfg.nclasses)) = tmp(i,cfg.nclasses + (1:cfg.nclasses))./size(data,1);
                tdiff = tdiff + sum(abs(tmodel - tmp(i,cfg.nclasses + (1:cfg.nclasses))));

            end
        end
                
        grad = cfg.loss_grad(cfg,tmp,data,g);
        
        lambdas(g) = norm(grad,2);

    end

else
    error('exact expression for p=%d is not available',cfg.p);
end

function [gmodel,gloss,greg] = arhmm_optimize(cfg,model,data,g,regularized)

global LOSS;

% note that greg is the regularization term for group g only

% the backward step is only performed for the groups that are nonzero

% compute loss for the old model
oldreg = cfg.reg(cfg,model,g);
oldloss = cfg.loss(cfg,model,data,g);

% gradient of the loss term plus regularization term
if ~regularized
    grad = cfg.epsilon * cfg.loss_grad(cfg,model,data,g);
else
    grad = cfg.epsilon * (cfg.loss_grad(cfg,model,data,g) + cfg.lambda .* cfg.reg_grad(cfg,model,g));
end

% return immediately in case of gradient problems
if any(isnan(grad)) 
    gmodel = model;
    gloss = LOSS;
    greg = oldreg;
    return
end

% compute loss for the new model
newmodel = update_model(cfg,model,-grad,g);

reg = cfg.reg(cfg,newmodel,g);
loss = cfg.loss(cfg,newmodel,data,g);

% update

oldgamma = oldloss;
gamma = loss;
if regularized
    oldgamma = oldgamma + cfg.lambda * oldreg;
    gamma = gamma + cfg.lambda * reg;
end

if oldgamma >= gamma

    gmodel = newmodel;
    LOSS = LOSS - oldloss + loss;
    greg = reg;

else
    
    if regularized
        ofun = @(x) (lnsrch_backward(cfg,model,data,g,x));
    else
        ofun = @(x) (lnsrch_forward(cfg,model,data,g,x));
    end

    dw = lnsrch(model(cfg.groups{g}), oldgamma, grad, -0.1 * grad, ofun);

    newmodel = update_model(cfg,model,dw,g);
    
    reg = cfg.reg(cfg,newmodel,g);
    loss = cfg.loss(cfg,newmodel,data,g);
    
    gamma = loss;
    if regularized
        gamma = gamma + cfg.lambda * reg;
    end
    
    if oldgamma >= gamma % accept

        gmodel = newmodel;
        LOSS = LOSS - oldloss + loss;
        greg = reg;
    else
        gmodel = model;
        greg = oldreg;
    end
end

gloss = LOSS;

function fval = lnsrch_forward(cfg,model,data,g,x)
% linesearch objective for forward step

newmodel = update_model(cfg,model,x,g);
fval = cfg.loss(cfg,newmodel,data,g);

function fval = lnsrch_backward(cfg,model,data,g,x)
% linesearch objective for backward step

newmodel = update_model(cfg,model,x,g);
fval = cfg.loss(cfg,newmodel,data,g) + cfg.lambda * cfg.reg(cfg,newmodel,g);

function newmodel = update_model(cfg,model,update,g)
% the function to update the model

newmodel = model;
newmodel(cfg.groups{g}) = model(cfg.groups{g}) + update;


