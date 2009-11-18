function [path,diagnostics] = regularize_lr_transfer(cfg,data)
% REGULARIZE_LR_TRANSFER returns the regularization path for logistic regression
% where we use a transfer learning approach that couples tasks (i.e.,
% subjects or sessions).
%
% [path,diagnostics] = regularize_lr_transfer(cfg,data)
%
% cfg contains the configuration options:
%  
% cfg.maxtime = 1e5; % maximum number of seconds to run the algorithm
% cfg.maxgroup = 1e4; % maximum number of used groups
% cfg.cvar = 1; the column for the class (or regression) variable
% cfg.nclasses; the number of classes (automatically determined if unspecified)
% cfg.tolerance = 1e-5; % error tolerance
% cfg.p = 2; % L1/LP norm
% cfg.groups; % grouping of the parameter vector model; groups are
% specified as class - feature combinations (linear indices are converted)
%
% at the moment we only allow a standard grouping for tranfer learning!
%
% data is a cell array that contains the examples X variables data for the
% model for each task
%
% Copyright (C) 2008, Marcel van Gerven
%
% $Log: regularize_lr_transfer.m,v $
% Revision 1.2  2008/02/29 12:28:39  marvger
% added unregularized logistic regression
%
% Revision 1.1.1.1  2008/02/27 14:42:55  roboos
% Marcel van Gerven, version 27 Feb 2008
%


%% initialization

if ~iscell(data)
    error('expecting task-specific data as a cell array');
end

% lambdas need to be recomputed in each iteration
cfg.recompute = true;

if ~isfield(cfg,'cvar'), cfg.cvar = 1; end

if ~isfield(cfg,'nclasses')

    if any((mod(data{1}(:,cfg.cvar),1)))
        cfg.nclasses = 1;
    else
        % the first task must have all classes
        cfg.nclasses = length(unique(data{1}(:,cfg.cvar)));
    end
    
end

if ~isfield(cfg,'p'), cfg.p = 2; end

% change data representation

cfg.ntasks = length(data);
cfg.nexamples = cell(1,cfg.ntasks);

dat = cell(1,cfg.ntasks);
for s=1:cfg.ntasks
    cfg.nexamples{s} = size(data{s},1);
    dat{s}.input = [data{s}(:,[1:(cfg.cvar-1) (cfg.cvar+1):size(data{s},2)]) ones(cfg.nexamples{s},1)];
    dat{s}.targets = data{s}(:,cfg.cvar);
end

cfg.nfeatures = size(dat{1}.input,2);

%% grouping of weight entries

if ~isfield(cfg,'groups')
        
    if cfg.nclasses == 2 % only regularize first class
        cfg.nclasses = 1; % temporarily change to one class
    end

    if cfg.p == 1 % no grouping whatsoever
     
        cfg.groups = cell(1,cfg.ntasks * cfg.nclasses * cfg.nfeatures);
        cfg.isregularized = true(1,cfg.ntasks * cfg.nfeatures * cfg.nclasses);

        i = 1;
        for f=1:(cfg.nfeatures - 1)
            for c=1:cfg.nclasses
                for t=1:cfg.ntasks
                    cfg.groups{i} = [t c f];
                    i = i + 1;
                end
            end

        end

        for c=1:cfg.nclasses
            for t=1:cfg.ntasks
                cfg.groups{i} = [t c cfg.nfeatures];

                % bias term should never be regularized
                cfg.isregularized(i) = false;

                i = i + 1;
            end
        end
        
    else % group over tasks

        cfg.groups = cell(1,cfg.nclasses * cfg.nfeatures);
        cfg.isregularized = true(1,cfg.nfeatures*cfg.nclasses);

        i = 1;
        for f=1:(cfg.nfeatures - 1)
            for c=1:cfg.nclasses
                cfg.groups{i} = [];
                for t=1:cfg.ntasks
                    cfg.groups{i} = [cfg.groups{i}; t c f];
                end
                i = i + 1;
            end

        end

        for c=1:cfg.nclasses
            for t=1:cfg.ntasks
                cfg.groups{i} = [t c cfg.nfeatures];

                % bias term should never be regularized
                cfg.isregularized(i) = false;

                i = i + 1;
            end
        end
    end
            
    if cfg.nclasses == 1 % this must have been 2
        cfg.nclasses = 2;
    end
    
end

cfg.ngroups = length(cfg.groups);

cfg.gtasks = cell(1,cfg.ngroups);
cfg.gfeatures = cell(1,cfg.ngroups);
cfg.gclasses = cell(1,cfg.ngroups);

for g=1:cfg.ngroups % class - feature indexing
    
     % group indices for each group and each task
    cfg.gtasks{g} = cfg.groups{g}(:,1);
    
    % group indices for each group and each task
    cfg.gclasses{g} = cfg.groups{g}(:,2);
    
    % feature indices for each group and each task
    cfg.gfeatures{g} = cfg.groups{g}(:,3);
    
    % change back to linear indices for each group and each task
    cfg.groups{g} = sub2ind([cfg.ntasks cfg.nclasses cfg.nfeatures],cfg.gtasks{g},cfg.gclasses{g},cfg.gfeatures{g})';

    cfg.gtasks{g} = transpose(cfg.gtasks{g});
    
end

% keep track of the size of each group per task
cfg.gngroup = cellfun(@numel,cfg.gfeatures);

%% prepare objective function

model = zeros(cfg.ntasks,cfg.nclasses,cfg.nfeatures);

%% specify functions

cfg.loss = @lr_loss;
cfg.loss_grad = @lr_loss_grad;
cfg.reg = @lr_reg;
cfg.reg_grad = @lr_reg_grad;
cfg.optimize = @lr_optimize;
cfg.stability = @lr_stability;

%% some additional changes in the representation

cfg.ptargets = cell(1,cfg.ntasks);
for s=1:cfg.ntasks
    classidxs = (1:cfg.nexamples{s})' + (dat{s}.targets - 1) .* cfg.nexamples{s};
    cfg.ptargets{s} = zeros(cfg.nexamples{s},cfg.nclasses);
    cfg.ptargets{s}(classidxs) = 1;
end

% indices of targets per example
cfg.targets = cell(1,cfg.ntasks);
for s=1:cfg.ntasks
    cfg.targets{s} = (1:cfg.nexamples{s})' + (dat{s}.targets - 1) * cfg.nexamples{s};
end

% compute exp_r and softmax

global EXPR;
global SOFTMAX;

EXPR = cell(1,cfg.ntasks);
SOFTMAX = cell(1,cfg.ntasks);
for s=1:cfg.ntasks
    EXPR{s} = exp(dat{s}.input * squeeze(model(s,:,:))');
    z = sum(EXPR{s}, 2);
    SOFTMAX{s} = EXPR{s} ./ z(:, ones(size(EXPR{s},2), 1));
end

cfg.iloss = cfg.loss(cfg,SOFTMAX);

%% run algorithm

[path,diagnostics] = regularize(cfg,model,dat);

%% function redefinition

function [gmodel,gloss,greg] = lr_optimize(cfg,model,data,g,regularized)

global EXPR;
global SOFTMAX;

% compute loss for the old model
oldreg = cfg.reg(cfg,model,g);
oldloss = cfg.loss(cfg,SOFTMAX);

if ~regularized
    grad = cfg.epsilon * cfg.loss_grad(cfg,SOFTMAX,data,g);
else
    grad = cfg.epsilon * (cfg.loss_grad(cfg,SOFTMAX,data,g) + cfg.lambda .* cfg.reg_grad(cfg,model,g));
end

if any(isnan(grad)) % return immediately in case of problems
    gmodel = model;
    gloss = oldloss;
    greg = oldreg;
    return
end

% compute loss for the new model
newmodel = lr_update_model(cfg,model,-grad,g);

% recompute global variables
[expr,softmax] = recompute_globals(cfg,EXPR,data,-grad,g);

reg = cfg.reg(cfg,newmodel,g);
loss = cfg.loss(cfg,softmax);

% update

oldgamma = oldloss;
gamma = loss;
if regularized
    oldgamma = oldgamma + cfg.lambda * oldreg;
    gamma = gamma + cfg.lambda * reg;
end

if  oldgamma >= gamma && ~isnan(gamma)

    EXPR = expr;
    SOFTMAX = softmax;
    gmodel = newmodel;
    gloss = loss;
    greg = reg;

else

    if regularized
        dw = lnsrch(model(cfg.groups{g}), oldgamma, grad, -0.5 * grad, @(x) (lr_lnsrch_backward(cfg,model,data,g,x)));
    else
        dw = lnsrch(model(cfg.groups{g}), oldgamma, grad, -0.5 * grad, @(x) (lr_lnsrch_forward(cfg,model,data,g,x)));
    end
    
    newmodel = lr_update_model(cfg,model,dw,g);
    [expr,softmax] = recompute_globals(cfg,EXPR,data,dw,g);

    reg = cfg.reg(cfg,newmodel,g);
    loss = cfg.loss(cfg,softmax);
    
    gamma = loss;
    if regularized
        gamma = gamma + cfg.lambda * reg;
    end

    if oldgamma >= gamma % accept

        EXPR = expr;
        SOFTMAX = softmax;
        gmodel = newmodel;
        gloss = loss;
        greg = reg;
    else
        
        dw =  cfg.oldmodel(cfg.groups{g}) - model(cfg.groups{g});
        
        gmodel = lr_update_model(cfg,model,dw,g);
        [expr,softmax] = recompute_globals(cfg,EXPR,data,dw,g);
        
        greg = cfg.reg(cfg,gmodel,g);
        
        gloss = cfg.loss(cfg,softmax);
                
    end
end

function fval = lr_lnsrch_forward(cfg,model,data,g,x)
% linesearch objective for forward step

global EXPR;

[expr,softmax] = recompute_globals(cfg,EXPR,data,x,g);

fval = cfg.loss(cfg,softmax);

function fval = lr_lnsrch_backward(cfg,model,data,g,x)
% linesearch objective for backward step

global EXPR;

newmodel = lr_update_model(cfg,model,x,g);
[expr,softmax] = recompute_globals(cfg,EXPR,data,x,g);

fval = cfg.loss(cfg,softmax) + cfg.lambda * cfg.reg(cfg,newmodel,g);

function newmodel = lr_update_model(cfg,model,update,g)
% the function to update the model

newmodel = model;
newmodel(cfg.groups{g}) = model(cfg.groups{g}) + update;


%% computational speedup

function [newexpr,newsoftmax] = recompute_globals(cfg,expr,data,x,g)

global SOFTMAX;

if any(x) % if x is nonzero

    newexpr = expr;
    newsoftmax = SOFTMAX;
    
    for i=1:cfg.gngroup(g)
        s = cfg.gtasks{g}(i);
        
        newexpr{s}(:,cfg.gclasses{g}(i)) = expr{s}(:,cfg.gclasses{g}(i)) .* exp(data{s}.input(:,cfg.gfeatures{g}(i)) * x(i)); 
        z = sum(newexpr{s}, 2);
        newsoftmax{s} = newexpr{s} ./ z(:, ones(size(newexpr{s},2), 1));
    end
    
else
    newexpr = expr;
    newsoftmax = SOFTMAX;
end

%% subroutines

function loss = lr_loss(cfg,softmax)
% loss term for logistic regression

loss = 0;
for s=1:cfg.ntasks
    loss = loss - sum(log(softmax{s}(cfg.targets{s})));
end

function grad = lr_loss_grad(cfg,softmax,data,g)
% gradient of the loss term w.r.t. group g

% this code assumes no grouping within a task!

grad = zeros(1,cfg.gngroup(g));
for i=1:cfg.gngroup(g)
    s = cfg.gtasks{g}(i);
    grad(i) = - sum(data{s}.input(:,cfg.gfeatures{g}(i)) .* (cfg.ptargets{s}(:,cfg.gclasses{g}(i)) - softmax{s}(:,cfg.gclasses{g}(i))));
end

function reg = lr_reg(cfg,model,g)
% regularization term for logistic regression

reg = 0;
if cfg.isregularized(g) && cfg.selected(g)
    reg = norm(model(cfg.groups{g}),cfg.p);
end

function grad = lr_reg_grad(cfg,model,g)
% gradient of the regularization term w.r.t. group g

if cfg.isregularized(g) && any(model(cfg.groups{g}))
    grad = (sign(model(cfg.groups{g})) .* (abs(model(cfg.groups{g})) ./ norm(model(cfg.groups{g}),cfg.p)).^(cfg.p-1));
    grad(isnan(grad)) = 0;
else
    grad = zeros(1,cfg.gngroup(g));
end

function lambdas = lr_stability(cfg,model,data)
% stability condition for logistic regression

global EXPR; % from the original model
global SOFTMAX; 

lambdas = zeros(1,length(cfg.groups));

% only consider groups that have not yet been considered when computing new
% values of lambda
for g=1:cfg.ngroups

    % set this group to zero
        
    if any(cfg.groups{g})
    
        x = zeros(1,length(cfg.groups{g}));
        
        [expr,softmax] = recompute_globals(cfg,EXPR,data,x,g);
        
    else
        softmax = SOFTMAX;
    end

    grad = lr_loss_grad(cfg,softmax,data,g);
             
    if cfg.p == Inf
        lambdas(g) = norm(grad,1);
    elseif cfg.p == 1
        lambdas(g) = norm(grad,Inf);
    else
        lambdas(g) = norm(grad,cfg.p/(cfg.p-1));
    end
end
