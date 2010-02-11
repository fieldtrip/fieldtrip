function [path,diagnostics] = regularize_nb(cfg,data)
% REGULARIZE_NB returns the regularization path for naive Bayes
%
% [path,diagnostics] = regularize_nb(cfg,data)
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
% data contains the examples X variables data for the model
%
% we assume all continuous features!
%
% Copyright (C) 2008, Marcel van Gerven
%
% $Log: regularize_nb.m,v $
% Revision 1.4  2008/04/16 12:22:18  marvger
% regular update
%
% Revision 1.3  2008/04/09 10:49:48  marvger
% changed activeset=cfg.nonzero to activeset=selected
%
% Revision 1.2  2008/03/01 16:38:28  marvger
% removed somes files
%
% Revision 1.1.1.1  2008/02/27 14:42:55  roboos
% Marcel van Gerven, version 27 Feb 2008
%


%% initialization

cfg.recompute = false; % do not recompute the lambdas

if ~isfield(cfg,'cvar'), cfg.cvar = 1; end
 
if ~isfield(cfg,'nclasses')

    if any((mod(data(:,cfg.cvar),1)))
        cfg.nclasses = 1;
    else
        cfg.nclasses = length(unique(data(:,cfg.cvar)));
    end
    
end

if ~isfield(cfg,'p'), cfg.p = 2; end

if ~(cfg.p == 2 || cfg.p == Inf)
    error('cfg.p can only be 2 or Inf');
end

% determine features
cfg.features = [1:(cfg.cvar-1) (cfg.cvar+1):size(data,2)]; 
cfg.nfeatures = length(cfg.features);

%% prepare objective function

model = zeros(cfg.nclasses*2,cfg.nfeatures); % space for mean and variance per feature

%% initialize means and standard deviations

for m=1:cfg.nfeatures

    model(1:cfg.nclasses,m) = mean(data(:,cfg.features(m))) * ones(cfg.nclasses,1);
    model((1:cfg.nclasses) + cfg.nclasses,m) = std(data(:,cfg.features(m))) * ones(cfg.nclasses,1);
end

%% data is split up according to classes for faster processing

cfg.data = cell(1,cfg.nclasses);
for c=1:cfg.nclasses
    cfg.data{c} = data(data(:,cfg.cvar) == c,:);
end

%% specify groups

cfg.groups = cell(1,length(cfg.features));
for m=1:cfg.nfeatures

    cfg.groups{m} = (1:(2*cfg.nclasses)) + (m-1)*2*cfg.nclasses;    
end

% everything is regularized
cfg.isregularized = true(1,length(cfg.groups));

%% specify functions

cfg.loss = @nb_loss;
cfg.loss_grad = @nb_loss_grad;
cfg.reg = @nb_reg;
cfg.reg_grad = @nb_reg_grad;

cfg.stability = @nb_stability;

% keep track of total loss
cfg.optimize = @nb_optimize;

%% run algorithm

global LOSS;

LOSS = 0;
for f=1:cfg.nfeatures
    LOSS = LOSS + cfg.loss(cfg,model,data,f);
end

% specify initial loss
cfg.iloss = LOSS;

[path,diagnostics] = regularize(cfg,model,data);


function loss = nb_loss(cfg,model,data,g)
% loss term for naive bayes

mu = model(1:cfg.nclasses,g);
sigma = model((1:cfg.nclasses) + cfg.nclasses,g);

% return infinite loss if sigma becomes negative
if any(sigma <= 0)
    loss = Inf;
    return
end

loss = sum(((data(:,cfg.features(g)) - mu(data(:,cfg.cvar))).^2) ./ (2 * sigma(data(:,cfg.cvar)).^2) + log(sqrt(2*pi)*sigma(data(:,cfg.cvar))));


function grad = nb_loss_grad(cfg,model,data,g)
% gradient of the loss term w.r.t. group g

gradmu = zeros(cfg.nclasses,1);
gradsigma = zeros(cfg.nclasses,1);

mu = model(1:cfg.nclasses,g);
sigma = model((1:cfg.nclasses) + cfg.nclasses,g);

for c=1:cfg.nclasses
    
    gradmu(c) = sum((mu(c) - cfg.data{c}(:,cfg.features(g))) ./ sigma(c).^2);

    gradsigma(c) = sum(- (cfg.data{c}(:,cfg.features(g)) - mu(c)).^2 ./ sigma(c).^3 + 1/sigma(c));

end

grad = [gradmu; gradsigma]';



function reg = nb_reg(cfg,model,g)
% regularization term for naive bayes

if cfg.isregularized(g) && cfg.selected(g)

    mu = model(cfg.groups{g}(1:(end/2)));
    sigma = model(cfg.groups{g}(((end/2)+1):end));

    reg = norm([mu - mean(mu) sigma - mean(sigma)],cfg.p);
else
    reg = 0;
end

function grad = nb_reg_grad(cfg,model,g)
% gradient of the regularization term w.r.t. group g

if cfg.isregularized(g)

    gradmu = zeros(cfg.nclasses,1);
    gradsigma = zeros(cfg.nclasses,1);

    mu = model(1:cfg.nclasses,g);
    sigma = model((1:cfg.nclasses) + cfg.nclasses,g);

    % f is the denominator in the partial derivative of the regularization term

    f = 0;
    muk = mean(mu);
    musigmak = mean(sigma);
    for c=1:cfg.nclasses
        f = f + (mu(c) - muk).^cfg.p;
        f = f + (sigma(c) - musigmak).^cfg.p;
    end
    f = f.^((cfg.p-1)/cfg.p);

    if f

        % compute second component of partial derivative

        musgn = zeros(1,cfg.nclasses);
        for c=1:cfg.nclasses
            musgn(c) = abs(mu(c) - muk).^(cfg.p-1) * sign(mu(c) - muk);
        end
        musum = sum(musgn) / cfg.nclasses;

        sigmasgn = zeros(1,cfg.nclasses);
        for c=1:cfg.nclasses
            sigmasgn(c) = abs(sigma(c) - musigmak).^(cfg.p-1) * sign(sigma(c) - musigmak);
        end
        sigmasum = sum(sigmasgn) / cfg.nclasses;

        % deal with means

        for c=1:cfg.nclasses
            gradmu(c) = gradmu(c) + (musgn(c) - musum) / f;
        end

        % deal with sigmas

        for c=1:cfg.nclasses
            gradsigma(c) = gradsigma(c) + (sigmasgn(c) - sigmasum) / f;
        end

        grad = [gradmu; gradsigma]';

    else
        grad = zeros(1,2*cfg.nclasses); 
    end

else
   grad = zeros(1,2*cfg.nclasses); 
end

function lambdas = nb_stability(cfg,model,data)
% stability condition for naive bayes

lambdas = zeros(1,length(cfg.groups));

if cfg.p == 2

    for m=1:cfg.nfeatures

        grad = cfg.loss_grad(cfg,model,data,m);
        
        lambdas(m) = norm(grad,2);

    end

elseif cfg.p == Inf % CHECK THIS; ALSO Inf in other computations

  
     for m=1:cfg.nfeatures

        grad = cfg.loss_grad(cfg,model,data,m);
    
        % subtract median for p=Inf
        grad(1:(end/2)) = grad(1:(end/2)) - median(grad(1:(end/2)));
        grad(((end/2)+1):end) = grad(((end/2)+1):end) - median(grad(((end/2)+1):end));
        
        lambdas(m) = sum(abs(grad));

     end

else
    error('exact expression for p=%d is not available',cfg.p);
end


function [gmodel,gloss] = nb_forward_step(cfg,model,data,g)

global LOSS;

% compute loss for the old model
oldloss = cfg.loss(cfg,model,data,g);

grad = cfg.epsilon * cfg.loss_grad(cfg,model,data,g);

% compute loss for the new model
newmodel = update_model(cfg,model,-grad,g);

% loss may depend just on group g
loss = cfg.loss(cfg,newmodel,data,g);

% update
if oldloss < loss || isnan(loss) % do a line search

    ofun = @(x) (lnsrch_forward(cfg,model,data,g,x));

    dw = lnsrch(model(cfg.groups{g}), oldloss, grad, -0.1 * grad, ofun);

    newmodel = update_model(cfg,model,dw,g);
    
    loss = cfg.loss(cfg,newmodel,data,g);
    if oldloss < loss || isnan(loss)

        gmodel = model;
    else
        gmodel = newmodel;
        gloss = LOSS - oldloss + loss;
    end

else
    gmodel = newmodel;
    gloss = LOSS - oldloss + loss;
end



function [gmodel,gloss,greg] = nb_optimize(cfg,model,data,g,regularized)

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
    
    if regularized && cfg.lambda ~= 0
        ofun = @(x) (lnsrch_backward(cfg,model,data,g,x));
    else
        ofun = @(x) (lnsrch_forward(cfg,model,data,g,x));
    end

    m = model(cfg.groups{g});
    dw = lnsrch(m(:)', oldgamma, grad, -0.1 * grad, ofun);

    newmodel = update_model(cfg,model,dw,g);
    
    reg = cfg.reg(cfg,newmodel,g);
    loss = cfg.loss(cfg,newmodel,data,g);
    
    gamma = loss;
    if regularized && cfg.lambda ~= 0
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
m = model(cfg.groups{g});
newmodel(cfg.groups{g}) = m(:)' + update;


