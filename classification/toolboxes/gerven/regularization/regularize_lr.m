function [path,diagnostics] = regularize_lr(cfg,data)
% REGULARIZE_LR returns the regularization path for logistic regression
%
% [path,diagnostics] = regularize_lr(cfg,data)
%
% cfg contains the configuration options:
%  
% cfg.maxtime = 1e5; % maximum number of seconds to run the algorithm
% cfg.maxgroup = 1e4; % maximum number of used groups
% cfg.cvar = 1; the column for the class (or regression) variable
% cfg.nclasses; the number of classes (automatically determined if unspecified)
% cfg.tolerance = 1e-5; % error tolerance
% cfg.p = 1; % L1 norm
% cfg.groups; % grouping of the parameter vector model; groups are
% specified as class - feature combinations (linear indices are converted)
%
% data contains the examples X variables data for the model
%
% To Do:
% - rethink cfg.overlap
%
% Copyright (C) 2008, Marcel van Gerven
%
% $Log: regularize_lr.m,v $
% Revision 1.4  2008/03/03 16:32:21  marvger
% changed lr_reg_grad computation exception: partial reason for problems
% when p becomes large
%
% Revision 1.3  2008/02/28 08:07:26  marvger
% small step size for forward regularization
%
% Revision 1.2  2008/02/27 17:49:14  marvger
% removed some backup files
%
% Revision 1.1.1.1  2008/02/27 14:42:55  roboos
% Marcel van Gerven, version 27 Feb 2008
%


%% initialization

% lambdas need to be recomputed in each iteration
cfg.recompute = true;

if ~isfield(cfg,'cvar'), cfg.cvar = 1; end
 
if ~isfield(cfg,'nclasses')

    if any((mod(data(:,cfg.cvar),1)))
        cfg.nclasses = 1;
    else
        cfg.nclasses = length(unique(data(:,cfg.cvar)));
    end
    
end

% assume standard L1 regularization by default!
if ~isfield(cfg,'p'), cfg.p = 1; end

% change data representation

cfg.nexamples = size(data,1);
dat.input = [data(:,[1:(cfg.cvar-1) (cfg.cvar+1):size(data,2)]) ones(cfg.nexamples,1)];
dat.targets = data(:,cfg.cvar);

cfg.nfeatures = size(dat.input,2);

%% grouping of weight entries

if ~isfield(cfg,'groups')

    % in case of L1 regularization each weight has its own group
    if cfg.p == 1
    
        % for a two class problem we only update the
        % weights associated with class 1
        if cfg.nclasses == 2

            cfg.groups = cell(1,cfg.nfeatures);
            
            for f=1:cfg.nfeatures
                cfg.groups{f} = [1 f];
            end
            
            % bias term should not be regularized
            cfg.isregularized = true(1,cfg.nfeatures); 
            cfg.isregularized(cfg.nfeatures) = false;
            
        else
    
            % [c f] indexing
            cfg.groups = cell(1,cfg.nclasses * cfg.nfeatures);

            cfg.isregularized = true(1,cfg.nfeatures*cfg.nclasses); 
            
            i = 1;
            for f=1:(cfg.nfeatures - 1)
                for c=1:cfg.nclasses

                    cfg.groups{i} = [c f];
                    i = i + 1;
                end
                
            end
            
            for c=1:cfg.nclasses

                cfg.groups{i} = [c cfg.nfeatures];
                
                % bias term should not be regularized
                cfg.isregularized(i) = false;
                
                i = i + 1;
            end
            
        end

    else % by default we group the classes

        cfg.groups = cell(1,cfg.nfeatures);

        i = 1;
        for f=1:cfg.nfeatures

            cfg.groups{i} = [(1:cfg.nclasses)' f * ones(cfg.nclasses,1)];

            i = i + 1;
        end

        % bias term should not be regularized
        cfg.isregularized = true(1,cfg.nfeatures);
        cfg.isregularized(cfg.nfeatures) = false;

    end
end

% convert linear indices to [c f] representation
cfg.ngroups = length(cfg.groups);
for g=1:cfg.ngroups

    if size(cfg.groups{g},2) == 1 % linear indexing using a column vector

        [gclasses gfeatures] = ind2sub([cfg.nclasses cfg.nfeatures], cfg.groups{g});
        cfg.groups{g} = [gclasses' gfeatures'];
    end
end

% we precompute the feature, class, and linear indices for each group
cfg.gfeatures = cell(1,cfg.ngroups);
cfg.gclasses = cell(1,cfg.ngroups);

for g=1:cfg.ngroups % class - feature indexing

    % group indices for each group and each task
    cfg.gclasses{g} = cfg.groups{g}(:,1);
    
    % feature indices for each group and each task
    cfg.gfeatures{g} = cfg.groups{g}(:,2);
    
    % change back to linear indices for each group and each task
    cfg.groups{g} = sub2ind([cfg.nclasses cfg.nfeatures],cfg.gclasses{g},cfg.gfeatures{g})';
end

% keep track of the size of each group
cfg.gngroup = cellfun(@numel,cfg.groups);
    
if isfield(cfg,'overlap') && cfg.overlap

    % make bit matrix
    bm = false(cfg.ngroups,cfg.ngroups);

    for g=1:cfg.ngroups
        bm(g,cfg.groups{g}) = 1;
    end

    % we assume a group can occur only once per group
    cfg.groupcount = cfg.groups;
    for g=1:cfg.ngroups
        cfg.groupcount{g} = sum(bm(:,cfg.groups{g}));        
    end

    cfg.dependence = cell(1,length(cfg.groups));
    for g=1:cfg.ngroups        
        cfg.dependence{g} = find(any(bm(:,cfg.groups{g}),2));
    end
   
%     % keep track of how many times group elements are represented in all groups
%     cfg.groupcount = cfg.groups;    
%     for i=1:length(cfg.groups)
%         
%         for j=1:length(cfg.groupcount{i})
%             
%             g = cfg.groupcount{i}(j);
%             cfg.groupcount{i}(j) = 0;
%             
%             for k=1:length(cfg.groups)
%                 cfg.groupcount{i}(j) = cfg.groupcount{i}(j) + sum(cfg.groups{k} == g);
%             end
%         end
%     end           
%     
%     % keep track of dependence
%     cfg.dependence = cell(1,length(cfg.groups));
%     for g=1:length(cfg.groups)
%         for j=1:length(cfg.groups)
%        
%             if ~isempty(intersect(cfg.groups{g},cfg.groups{j}))
%                 cfg.dependence{g} = [cfg.dependence{g} j];
%             end
%             
%         end
%     end

end

%% prepare objective function

model = sparse(cfg.nclasses,cfg.nfeatures); 

%% specify functions

cfg.loss = @lr_loss;
cfg.loss_grad = @lr_loss_grad;
cfg.reg = @lr_reg;
cfg.reg_grad = @lr_reg_grad;
cfg.optimize = @lr_optimize;
cfg.stability = @lr_stability;

%% some additional changes in the representation

classidxs = (1:cfg.nexamples)' + (dat.targets - 1) .* cfg.nexamples;
cfg.ptargets = zeros(cfg.nexamples,cfg.nclasses);
cfg.ptargets(classidxs) = 1;

% indices of targets per example
cfg.targets = (1:cfg.nexamples)' + (dat.targets - 1) * cfg.nexamples;

% compute exp_r and softmax

global EXPR;
global SOFTMAX;

EXPR = exp(dat.input * model');
z = sum(EXPR, 2);
SOFTMAX = EXPR ./ z(:, ones(size(EXPR,2), 1));

%% run algorithm

cfg.iloss = cfg.loss(cfg,SOFTMAX);

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
    oldgamma = oldgamma + cfg.lambda * sum(oldreg);
    gamma = gamma + cfg.lambda * sum(reg);
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
    if regularized, gamma = gamma + cfg.lambda * sum(reg); end

    if oldgamma >= gamma && any(dw) % accept

        EXPR = expr;
        SOFTMAX = softmax;
        gmodel = newmodel;
        gloss = loss;
        greg = reg;
        
    else
        
        % in this case the update has failed;
        % we revert to the situation prior to the forward step

        dw =  cfg.oldmodel(cfg.groups{g}) - model(cfg.groups{g});

        gmodel = lr_update_model(cfg,model,dw,g);
        [expr,softmax] = recompute_globals(cfg,EXPR,data,dw,g);

        EXPR = expr;
        SOFTMAX = softmax;

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

reg = cfg.reg(cfg,newmodel,g);

fval = cfg.loss(cfg,softmax) + cfg.lambda * sum(reg);

function newmodel = lr_update_model(cfg,model,update,g)
% the function to update the model

newmodel = model;
newmodel(cfg.groups{g}) = model(cfg.groups{g}) + update;

%% computational speedup

function [newexpr,newsoftmax] = recompute_globals(cfg,expr,data,x,g)

global SOFTMAX;

if any(x)
    
    if cfg.gngroup(g) == 1 % just one class and feature; speed up computation
    
        newexpr = expr;
        newexpr(:,cfg.gclasses{g}) = expr(:,cfg.gclasses{g}) .* exp(data.input(:,cfg.gfeatures{g}) * x);
        z = sum(newexpr, 2);
        newsoftmax = newexpr ./ z(:, ones(cfg.nclasses, 1));
        
    else
        
        sw = sparse(cfg.nclasses,cfg.nfeatures);
        sw(cfg.groups{g}) = x;

        d = data.input * sw';
        dnonzero = (d ~= 0);

        newexpr = expr;
        newexpr(dnonzero) = newexpr(dnonzero) .* exp(d(dnonzero));
        z = sum(newexpr, 2);
        newsoftmax = newexpr ./ z(:, ones(cfg.nclasses, 1));
    end
else

    newexpr = expr;
    newsoftmax = SOFTMAX;
end

%% subroutines

function loss = lr_loss(cfg,softmax)
% loss term for logistic regression

loss = - sum(log(softmax(cfg.targets)));

function grad = lr_loss_grad(cfg,softmax,data,g)
% gradient of the loss term w.r.t. group g

grad = - sum(data.input(:,cfg.gfeatures{g}) .* (cfg.ptargets(:,cfg.gclasses{g}) - softmax(:,cfg.gclasses{g})));

function reg = lr_reg(cfg,model,g)
% regularization term for logistic regression

if cfg.overlap ~= 1 
    
    if cfg.isregularized(g) && cfg.selected(g)

        reg = norm(model(cfg.groups{g}),cfg.p);

    else
        reg = 0;
    end

else % if groups are overlapping then we recompute the full regularization
    
%     reg = zeros(1,cfg.ngroups);
%     for g=1:cfg.ngroups
%         if cfg.isregularized(g) %&& cfg.selected(g) % YES OR NO?????????????????????????????? 
%             
%             % bij yes gaan we opeens sterker regularizeren als hij erin komt
%             % bij no gaan we pas neighbours selecteren als iedereen actief
%             % wordt
%             
%             reg(g) = norm(model(cfg.groups{g}),cfg.p);
%         end
%     end

    reg = zeros(1,length(cfg.dependence{g}));
    for j=1:length(cfg.dependence{g});
        
        if cfg.isregularized(j)            
            reg(j) = norm(model(cfg.groups{cfg.dependence{g}(j)}),cfg.p);
        end
    end
end

function grad = lr_reg_grad(cfg,model,g)
% gradient of the regularization term w.r.t. group g

if cfg.isregularized(g)

    % special handling of p = inf
    if isinf(cfg.p)

        grad = zeros(1,numel(cfg.groups{g}));
        ag = abs(model(cfg.groups{g}));
        [a,b] = max(ag);
        grad(b) = sign(ag(b));

    else

        nrm = norm(model(cfg.groups{g}),cfg.p);

        if nrm
            grad = sign(model(cfg.groups{g})) .* (abs(model(cfg.groups{g})) ./ nrm).^(cfg.p-1);
        else
            grad = zeros(1,numel(cfg.groups{g})); %sign(model(cfg.groups{g})) .* abs(model(cfg.groups{g})).^(cfg.p-1);
        end
    end
else
    grad = zeros(1,numel(cfg.groups{g}));
end

function lambdas = lr_stability(cfg,model,data)
% stability condition for logistic regression

global EXPR; % from the original model
global SOFTMAX; 

lambdas = Inf*ones(1,length(cfg.groups));

thegroups = 1:cfg.ngroups;

% only consider groups that have not yet been considered when computing new values of lambda
if cfg.pursuit
    
    % disregard unregularized
    thegroups = thegroups(cfg.isregularized); 
else
    
    thegroups = thegroups(~cfg.selected & cfg.isregularized); 
end

for g=thegroups
    
    % set this group to zero
        
    % THIS SHOULD NOT HAPPEN (lambdas only decrease?)
%    if any(model(cfg.groups{g}))    
%        [expr,softmax] = recompute_globals(cfg,EXPR,data,zeros(1,length(cfg.groups{g})),g);        
%    else
%        softmax = SOFTMAX;
%    end    
    
    %grad = lr_loss_grad(cfg,softmax,data,g);    
    grad = - sum(data.input(:,cfg.gfeatures{g}) .* (cfg.ptargets(:,cfg.gclasses{g}) - SOFTMAX(:,cfg.gclasses{g})));
    
    % correction terms for overlapping groups
    if cfg.overlap
        grad = grad ./ cfg.groupcount{g};
    end
    
    if cfg.p == Inf        
        lambdas(g) = norm(grad,1);
    elseif cfg.p == 1        
        lambdas(g) = norm(grad,Inf);
    else
        lambdas(g) = norm(grad,cfg.p/(cfg.p-1));
    end
        
end
