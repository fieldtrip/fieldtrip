function [path,diagnostics] = regularize(cfg,model,data)
% REGULARIZE returns the regularization path at those values of the
% regularization parameter for which the obtained solutions change.
%
% [path,diagnostics] = regularize(cfg,model,data)
%
% cfg contains the configuration options:
%  
% cfg.lambda = Inf; % can be used to impose a starting value for the
% regularization parameter
% cfg.pursuit = false; % reconsider all features in each iteration when computing lambdas
% cfg.gridsize = 0; % define the minimum step in lambda
% cfg.recompute = true; % recompute all lambdas in each iteration (needed
% whenever the stability condition changes)
% cfg.maxtime = 1e5; % maximum number of seconds to run the algorithm
% cfg.maxgroup = 1e4; % maximum number of used groups
% cfg.cvar = 1; the column for the class (or regression) variable
% cfg.nclasses; the number of classes (automatically determined if unspecified)
% cfg.tolerance = 1e-3; % error tolerance
% cfg.p = 1; % L1/LP norm
% cfg.groups; % grouping of the parameter vector model (default: no grouping)
% cfg.epsilon = 1e-2; % step size < 1 may significantly speed up convergence
% cfg.isregularized = 1:ngroups % define which groups need to be
% regularized
% cfg.overlap = false; overlapping groups requires recomputing all
%   regularization terms; if cfg.overlap > 1 then we only look at the
%   regularization term when accepting the update step
%
% model contains the vector of initial parameters
% data contains the examples X variables data for the model
%
% only the loss term, regularization term, loss gradient, regularization
% gradient, and stability conditions need to be defined.
% forward and backward steps may be redefined in order to speed up certain computations
%
% Remarks:
% - this code requires that the gradient at zero equals zero in order to
% force a forward step that is subsequently minimized by a line search
%
% Copyright (C) 2008, Marcel van Gerven
%
% $Log: regularize.m,v $
% Revision 1.9  2008/04/09 10:49:48  marvger
% changed activeset=cfg.nonzero to activeset=selected
%
% Revision 1.8  2008/03/20 14:38:47  marvger
% regular update
%
% Revision 1.7  2008/03/07 12:22:22  marvger
% major changes w.r.t. ar-hmm models
%
% Revision 1.6  2008/03/03 16:32:21  marvger
% changed lr_reg_grad computation exception: partial reason for problems
% when p becomes large
%
% Revision 1.5  2008/03/01 16:38:28  marvger
% removed somes files
%
% Revision 1.4  2008/02/29 12:28:39  marvger
% added unregularized logistic regression
%
% Revision 1.3  2008/02/28 08:07:26  marvger
% small step size for forward regularization
%
% Revision 1.1.1.1  2008/02/27 14:42:55  roboos
% Marcel van Gerven, version 27 Feb 2008
%


%% initialization

% maximum number of optimization steps per iteration
maxupdate = 1e4; 

if ~isfield(cfg,'epsilon'), cfg.epsilon = 1e-2; end
epsilon = cfg.epsilon; % used to keep track of original epsilon

if ~isfield(cfg,'pursuit'), cfg.pursuit = false; end

if ~isfield(cfg,'unregularized'), cfg.unregularized = false; end
   
% step size in lambdas
if ~isfield(cfg,'gridsize'), cfg.gridsize = 0; end

% recompute lambdas in each iteration
if ~isfield(cfg,'recompute'), cfg.recompute = true; end

if ~isfield(cfg,'maxtime'), cfg.maxtime = 1e5; end

if ~isfield(cfg,'maxgroup'), cfg.maxgroup = 1e4; end

if ~isfield(cfg,'cvar'), cfg.cvar = 1; end
 
if ~isfield(cfg,'nclasses')

    if any((mod(data(:,cfg.cvar),1)))
        cfg.nclasses = 1;
    else
        cfg.nclasses = length(unique(data(:,cfg.cvar))); % assumes 1,2,3,..
    end
    
end

if ~isfield(cfg,'tolerance'), cfg.tolerance = 1e-3; end

if ~isfield(cfg,'p'), cfg.p = 1; end

if ~isfield(cfg,'groups')
    cfg.groups = num2cell(1:numel(model));
end
cfg.ngroups = length(cfg.groups);

% which groups are regularized
if ~isfield(cfg,'isregularized'), cfg.isregularized = true(1,cfg.ngroups); end

% current lambda
if ~isfield(cfg,'lambda'), cfg.lambda = Inf; end

% iterate lambda over an interval of values instead of 
% determining lambda
if ~isscalar(cfg.lambda)
    
    cfg.interval = cfg.lambda;
    cfg.lambda = Inf;

else
    cfg.interval = false;
end


% assume no overlapping groups
if ~isfield(cfg,'overlap'), cfg.overlap = false; end

%% we may override the forward and backward steps

if ~isfield(cfg,'optimize'), cfg.optimize = @optimize; end

%% minimize the unregularized features

tic

fiter=1; biter=0; 

% initial loss
if ~isfield(cfg,'iloss') 
    loss = cfg.loss(cfg,model,data);
else
    loss = cfg.iloss;
end

%% main loop

diagnostics = [];

diagnostics.gamma = [];
diagnostics.activeset = cell(1,cfg.ngroups);

cfg.nonzero = 0; 
path = cell(1,cfg.ngroups);

if cfg.recompute
    diagnostics.lambdahist = sparse(cfg.ngroups,cfg.ngroups);
    diagnostics.lambdas = sparse(1,cfg.ngroups);
end

cfg.selected = false(1,cfg.ngroups);

% first iteration updates the unregularized variables
if any(cfg.interval)
    lambdas = zeros(1,cfg.ngroups);
else
    % lambdas are computed based on initial lambda
    lambdas = cfg.stability(cfg,model,data);
end

lambdas(~cfg.isregularized) = Inf; % unregularized features are always unstable

% regularization term per group
regs = zeros(1,cfg.ngroups); 

% in the first iteration the unregularized features are minimized
while toc < cfg.maxtime
    
    % these features have been selected
    % minimization ensures that features move from active->inactive or
    % inactive->active
    
    % old selected is used to push actives back to inactive; other selected
    % are the ones which satisfy the stability condition and the
    % unregularized groups
    cfg.selected = cfg.selected | (lambdas >= cfg.lambda & lambdas < Inf) | ~cfg.isregularized;
    
%     if cfg.overlap
%         sel = 1:length(cfg.selected);
%     else
%         sel = find(cfg.selected);
%     end

    if cfg.overlap        
       dep = cfg.dependence(cfg.selected);
       cfg.selected(unique(cat(1,dep{:}))) = 1;       
    end

    sel = find(cfg.selected);

    if cfg.recompute
        diagnostics.lambdahist(fiter,:) = lambdas;
        diagnostics.lambdas(fiter) = cfg.lambda;
    end

    % compute gamma
    oldgamma = loss + max(cfg.lambda * sum(regs),0);
   
    % take several backward steps
    cfg.oldmodel = model;
    for i=1:maxupdate % use a maximum of 1000 updates!

        % iterate over all active groups in the backward step
        for cur=sel
            
            biter = biter + 1;
           
            if cfg.isregularized(cur)

                % regularized group
                [model,loss,reg] = cfg.optimize(cfg,model,data,cur,true);
                
                cfg.oldmodel = model;
                
                if isscalar(reg)
                    regs(cur) = reg; % only reg for current group returned
                else
                    regs(cfg.dependence{cur}) = reg; % multiple regs returned (overlapping groups)
                end

            else

                % unregularized group
                [model,loss] = cfg.optimize(cfg,model,data,cur,false);
                cfg.oldmodel = model;        
            end

        end
                                    
        % set to zero when lambda=Inf and reg=0
        regterm = max(cfg.lambda * sum(regs),0);         

        % compute new gamma
        gamma = loss + regterm;
               
        % we stop the optimization if the update drops below the tolerance
        update = oldgamma - gamma; 
                
        if update <= cfg.tolerance || isnan(update)
             break;
        end
        
        oldgamma = gamma;        

    end
    
    if i==maxupdate
        warning('maximum number of updates reached; proceeding to next group');
    end

    %% check the nonzero status of all groups
    cfg.nonzero = 0;
    for g = sel

        if iscell(model) % take transfer learning formulation into account
            
             for s=1:length(model)
               if  any(model{s}(cfg.groups{g}))
                   cfg.nonzero = cfg.nonzero + 1;
                   break
               end
            end
                        
        else
            if  any(model(cfg.groups{g}))
                cfg.nonzero = cfg.nonzero + 1;
            end
        end
    end

    diagnostics.gamma = [diagnostics.gamma gamma];
    
    fprintf('%d %d %g %g %g %g %d\n',fiter, biter,cfg.lambda,loss,regterm,gamma,cfg.nonzero);

    % save all distinct models
    path{fiter} = model;
   
    diagnostics.activeset{fiter} = cfg.selected;
    
    % lambda = 0 evaluated so stop
    if cfg.lambda == 0, break; end

    % maximum number of features reached
    if cfg.nonzero >= cfg.maxgroup, break; end
    
    % compute all lambdas
    if cfg.recompute || fiter == 1
               
        if any(cfg.interval)

            % set all lambdas to this lambda
            % this grid based optimization also prevents having to compute
            % the stability conditions
            if fiter <= length(cfg.interval)
                lambdas = cfg.interval(fiter)*ones(1,cfg.ngroups);
            end
        else
            lambdas = cfg.stability(cfg,model,data);
        end

        % unregularized features are always unstable
        lambdas(~cfg.isregularized) = Inf;
        
        if ~cfg.recompute, diagnostics.lambdas = lambdas; end
        
    end    

    % find next smallest lambda
    if any(cfg.interval)
        
        if fiter > length(cfg.interval)
            cfg.lambda = []; % none left
        else
            cfg.lambda = cfg.interval(fiter);
        end
    else

        if cfg.pursuit % not preferred since it is slow

            % pursuit reconsiders all recomputed lambdas in each run
            cfg.lambda = max(lambdas(lambdas < (cfg.lambda - cfg.gridsize)));

        else

            % we only consider lambdas for features that have not been included in the past
            cfg.lambda = max(lambdas(lambdas < (cfg.lambda - cfg.gridsize) & ~cfg.selected));
        end
    end
    
    % no smaller lambda found so return
    if isempty(cfg.lambda)
      if cfg.unregularized
        cfg.lambda = 0;
        diagnostics.lambdas = [diagnostics.lambdas 0];
      else
        break;
      end
    end 

    fiter = fiter + 1;
end

% last model was not specified
if toc >= cfg.maxtime, fiter = fiter-1; end

path = path(1:fiter); % truncate
diagnostics.activeset = diagnostics.activeset(1:fiter);
diagnostics.selected = sel;
diagnostics.toc = toc;

if cfg.recompute
    diagnostics.lambdahist = diagnostics.lambdahist(1:fiter,:);
end

function [gmodel,gloss,greg] = optimize(cfg,model,data,g,regularized)

% compute loss for the old model
oldloss = cfg.loss(cfg,model,data);
oldreg = cfg.reg(cfg,model,g);

% gradient of the loss term plus regularization term
if ~regularized %|| ~cfg.selected(g) % unregularized
    
    grad = cfg.epsilon * cfg.loss_grad(cfg,model,data,g);
    
else % regularized
    
    grad = cfg.epsilon * (cfg.loss_grad(cfg,model,data,g) + cfg.lambda .* cfg.reg_grad(cfg,model,g));
end

if any(isnan(grad)) % return immediately
    
    gmodel = model;
    gloss = oldloss;
    greg = oldreg;
    return
end

% compute loss for the new model
newmodel = update_model(cfg,model,-grad,g);

loss = cfg.loss(cfg,newmodel,data,g);
reg = cfg.reg(cfg,newmodel,g);

% update

oldgamma = oldloss;
gamma = loss;
if regularized
    oldgamma = oldgamma + cfg.lambda * oldreg;
    gamma = gamma + cfg.lambda * reg;
end


if oldgamma >= gamma && ~isnan(gamma)

    gmodel = newmodel;
    gloss = loss;
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
    loss = cfg.loss(cfg,newmodel,data);
    
    gamma = loss;
    if regularized
       gamma = gamma + cfg.lambda * reg;
    end
    
    if oldgamma >= gamma % accept

        gmodel = newmodel;
        gloss = loss;
        greg = reg;
    else

        % in this case the update has failed;
        % we revert to the situation prior to the forward step

        dw =  cfg.oldmodel(cfg.groups{g}) - model(cfg.groups{g});
        
        gmodel = lr_update_model(cfg,model,dw,g);
        [expr,softmax] = recompute_globals(cfg,EXPR,data,dw,g);

        if cfg.overlap
            greg = cfg.reg(cfg,gmodel); % recompute all terms
        else
            greg = cfg.reg(cfg,gmodel,g);
        end

        gloss = cfg.loss(cfg,softmax);
        
%         gmodel = model;
%         gloss = oldloss;
%         greg = oldreg;
    end
end


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



