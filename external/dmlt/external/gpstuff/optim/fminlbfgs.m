function [x,fval,exitflag,output,grad]=fminlbfgs(funfcn,x_init,optim)
%FMINLBFGS finds a local minimum of a function of several variables.
%          This optimizer is developed for image registration
%          methods with large amounts of unknown variables.
%
%  Description
%    Optimization methods supported:
%      - Quasi Newton Broyden-Fletcher-Goldfarb-Shanno (BFGS)
%      - Limited memory BFGS (L-BFGS)
%      - Steepest Gradient Descent optimization.
%
%    [X,FVAL,EXITFLAG,OUTPUT,GRAD] = FMINLBFGS(FUN,X0,OPTIONS)
%
%    Inputs,
%      FUN      Function handle or string which is minimized,
%               returning an error value and optional the error
%               gradient.
%      X0       Initial values of unknowns can be a scalar, vector
%               or matrix
%      OPTIONS  Structure with optimizer options, made by a struct or
%               optimset. (optimset does not support all input options)
%    Outputs,
%      X        The found location (values) which minimize the function.
%      FVAL     The minimum found
%      EXITFLAG Gives value, which explains why the minimizer stopped
%      OUTPUT   Structure with all important output values and parameters
%      GRAD     The gradient at this location
%
%  Extended description of inputs/outputs
%    OPTIONS,
%      GoalsExactAchieve - If set to 0, a line search method is
%                          used which uses a few function calls to
%                          do a good line search. When set to 1 a
%                          normal line search method with Wolfe
%                          conditions is used (default).
%      GradConstr        - Set this variable to 'on' if gradient
%                          calls are cpu-expensive. If 'off' more
%                          gradient calls are used and less
%                          function calls (default).
%      HessUpdate        - If set to 'bfgs', Broyden-Fletcher-Goldfarb-Shanno
%                          optimization is used (default), when the
%                          number of unknowns is larger then 3000
%                          the function will switch to Limited
%                          memory BFGS, or if you set it to
%                          'lbfgs'. When set to 'steepdesc',
%                          steepest decent optimization is used.
%      StoreN            - Number of iterations used to approximate
%                          the Hessian, in L-BFGS, 20 is default. A
%                          lower value may work better with non
%                          smooth functions, because than the
%                          Hessian is only valid for a specific
%                          position. A higher value is recommend
%                          with quadratic equations.
%      GradObj           - Set to 'on' if gradient available otherwise
%                          finite difference is used.
%      Display           - Level of display. 'off' displays no output;
%                          'plot' displays all linesearch results
%                          in figures. 'iter' displays output at
%                          each iteration; 'final' displays just
%                          the final output; 'notify' displays
%                          output only if the function does not
%                          converge;
%      TolX              - Termination tolerance on x, default 1e-6.
%      TolFun            - Termination tolerance on the function value,
%                          default 1e-6.
%      MaxIter           - Maximum number of iterations allowed, default 400.
%      MaxFunEvals       - Maximum number of function evaluations allowed,
%                          default 100 times the amount of unknowns.
%      DiffMaxChange     - Maximum stepsize used for finite difference
%                          gradients.
%      DiffMinChange     - Minimum stepsize used for finite difference
%                          gradients.
%      OutputFcn         - User-defined function that an optimization
%                          function calls at each iteration.
%      rho               - Wolfe condition on gradient (c1 on Wikipedia),
%                          default 0.01.
%      sigma             - Wolfe condition on gradient (c2 on Wikipedia),
%                          default 0.9.
%      tau1              - Bracket expansion if stepsize becomes larger,
%                          default 3.
%      tau2              - Left bracket reduction used in section phase,
%                          default 0.1.
%      tau3              - Right bracket reduction used in section phase,
%                          default 0.5.
%    FUN,
%      The speed of this optimizer can be improved by also
%      providing the gradient at X. Write the FUN function as
%      follows function [f,g]=FUN(X)
%       f , value calculation at X;
%       if ( nargout > 1 )
%         g , gradient calculation at X;
%       end
%    EXITFLAG,
%      Possible values of EXITFLAG, and the corresponding exit
%      conditions are
%        1, 'Change in the objective function value was less than TolFun.';
%        2, 'Change in x was smaller than the specified tolerance TolX.';
%        3, 'Magnitude of gradient smaller than the specified tolerance';
%        4, 'Boundary fminimum reached.';
%        0, 'Number of iterations exceeded options.MaxIter or
%            number of function evaluations exceeded
%            options.FunEvals.';
%       -1, 'Algorithm was terminated by the output function.';
%       -2, 'Line search cannot find an acceptable point along the
%            current search direction';
%
%   Examples
%       options = optimset('GradObj','on');
%       X = fminlbfgs(@myfun,2,options)
%
%       % where myfun is a MATLAB function such as:
%       function [f,g] = myfun(x)
%       f = sin(x) + 3;
%       if ( nargout > 1 ), g = cos(x); end
%
%  See also OPTIMSET, FMINSEARCH, FMINBND, FMINCON, FMINUNC, @, INLINE.
%
%  Function is written by D.Kroon University of Twente (Updated Nov. 2010)

% 2010-10-29 Aki Vehtari : GradConstr is 'off' by default.
% 2011-9-28 Ville Tolvanen : Reduce step size until function returns finite
%                            value
% 2012-10-10 Aki Vehtari : Improved robustness if function returns NaN or Inf
%                          Spell check and some beautification

%Copyright (c) 2009, Dirk-Jan Kroon
%All rights reserved.
%
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are
%met:
%
%    * Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in
%      the documentation and/or other materials provided with the distribution
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%POSSIBILITY OF SUCH DAMAGE.

% Read Optimization Parameters
defaultopt = struct('Display','final',...
                    'HessUpdate','bfgs',...
                    'GoalsExactAchieve',1,...
                    'GradConstr','off',  ...
                    'TolX',1e-6,...
                    'TolFun',1e-6,...
                    'GradObj','off',...
                    'MaxIter',400,...
                    'MaxFunEvals',100*numel(x_init)-1,...
                    'DiffMaxChange',1e-1,...
                    'DiffMinChange',1e-8,...
                    'OutputFcn',[], ...
                    'rho',0.0100,...
                    'sigma',0.900,...
                    'tau1',3,...
                    'tau2', 0.1,...
                    'tau3', 0.5,...
                    'StoreN',20);

if (~exist('optim','var'))
  optim=defaultopt;
else
  f = fieldnames(defaultopt);
  for i=1:length(f),
    if (~isfield(optim,f{i})||(isempty(optim.(f{i})))), optim.(f{i})=defaultopt.(f{i}); end
  end
end

% Initialize the data structure
data.fval=0;
data.gradient=0;
data.fOld=[];
data.xsizes=size(x_init);
data.numberOfVariables = numel(x_init);
data.xInitial = x_init(:);
data.alpha=1;
data.xOld=data.xInitial;
data.iteration=0;
data.funcCount=0;
data.gradCount=0;
data.exitflag=[];
data.nStored=0;
data.timeTotal=tic;
  data.timeExtern=0;
  % Switch to L-BFGS in case of more than 3000 unknown variables
  if(optim.HessUpdate(1)=='b')
    if(data.numberOfVariables<3000),
      optim.HessUpdate='bfgs';
    else
      optim.HessUpdate='lbfgs';
    end
  end

  if(optim.HessUpdate(1)=='l')
    succes=false;
    while(~succes)
      try
        data.deltaX=zeros(data.numberOfVariables,optim.StoreN);
        data.deltaG=zeros(data.numberOfVariables,optim.StoreN);
        data.saveD=zeros(data.numberOfVariables,optim.StoreN);
        succes=true;
      catch ME
        warning('fminlbfgs:memory','Decreasing StoreN value because out of memory');
        succes=false;
        data.deltaX=[]; data.deltaG=[]; data.saveD=[];
        optim.StoreN=optim.StoreN-1;
        if(optim.StoreN<1)
          rethrow(ME);
        end
      end
    end
  end

  exitflag=[];

  % Display column headers
  if(strcmp(optim.Display,'iter'))
    disp(' Iteration  Func-count   Grad-count         f(x)         Step-size');
  end

  % Calculate the initial error and gradient
  data.initialStepLength=1;
  [data,fval,grad]=gradient_function(data.xInitial,funfcn, data, optim);
  data.gradient=grad;
  data.dir = -data.gradient;
  data.fInitial = fval;
  data.fPrimeInitial= data.gradient'*data.dir(:);
  data.fOld=data.fInitial;
  data.xOld=data.xInitial;
  data.gOld=data.gradient;

  gNorm = norm(data.gradient,Inf);  % Norm of gradient
  data.initialStepLength = min(1/gNorm,5);

  % Show the current iteration
  if(strcmp(optim.Display,'iter'))
    s=sprintf(' %5.0f       %5.0f       %5.0f       %13.6g    ',data.iteration,data.funcCount,data.gradCount,data.fInitial); disp(s);
  end

  % Hessian intialization
  if(optim.HessUpdate(1)=='b')
    data.Hessian=eye(data.numberOfVariables);
  end

  % Call output function
  if(call_output_function(data,optim,'init')), exitflag=-1; end

  % Start Minimizing
  while(true)
    % Update number of iterations
    data.iteration=data.iteration+1;

    % Set current lineSearch parameters
    data.TolFunLnS = eps(max(1,abs(data.fInitial )));
    data.fminimum = data.fInitial - 1e16*(1+abs(data.fInitial));

    % Make arrays to store linesearch results
    data.storefx=[]; data.storepx=[]; data.storex=[]; data.storegx=[];

    % If option display plot, than start new figure
    if(optim.Display(1)=='p'), figure, hold on; end

    % Find a good step size in the direction of the gradient: Linesearch
    if(optim.GoalsExactAchieve==1)
      data=linesearch(funfcn, data,optim);
    else
      data=linesearch_simple(funfcn, data, optim);
    end

    % Make linesearch plot
    if(optim.Display(1)=='p');
      plot(data.storex,data.storefx,'r*');
      plot(data.storex,data.storefx,'b');

      alpha_test= linspace(min(data.storex(:))/3, max(data.storex(:))*1.3, 10);
      falpha_test=zeros(1,length(alpha_test));
      for i=1:length(alpha_test)
        [data,falpha_test(i)]=gradient_function(data.xInitial(:)+alpha_test(i)*data.dir(:),funfcn, data, optim);
      end
      plot(alpha_test,falpha_test,'g');
      plot(data.alpha,data.f_alpha,'go','MarkerSize',8);
    end

    % Check if exitflag is set
    if(~isempty(data.exitflag)),
      exitflag=data.exitflag;
      data.xInitial=data.xOld;
      data.fInitial=data.fOld;
      data.gradient=data.gOld;
      break,
    end;

    % Update x with the alpha step
    data.xInitial = data.xInitial + data.alpha*data.dir;

    % Set the current error and gradient
    data.fInitial =  data.f_alpha;
    data.gradient = data.grad;

    % Set initial step-length to 1
    data.initialStepLength = 1;


    gNorm = norm(data.gradient,Inf);  % Norm of gradient

    % Set exit flags
    if(gNorm <optim.TolFun), exitflag=1; end
    if(max(abs(data.xOld-data.xInitial)) <optim.TolX), exitflag=2; end
    if(data.iteration>=optim.MaxIter), exitflag=0; end

    % Check if exitflag is set
    if(~isempty(exitflag)), break, end;

    % Update the inverse Hessian matrix
    if(optim.HessUpdate(1)~='s')
      % Do the Quasi-Newton Hessian update.
      data = updateQuasiNewtonMatrix_LBFGS(data,optim);
    else
      data.dir = -data.gradient;
    end

    % Derivative of direction
    data.fPrimeInitial= data.gradient'*data.dir(:);

    % Call output function
    if(call_output_function(data,optim,'iter')), exitflag=-1; end

    % Show the current iteration
    if(strcmp(optim.Display(1),'i')||strcmp(optim.Display(1),'p'))
      s=sprintf(' %5.0f       %5.0f       %5.0f       %13.6g   %13.6g',data.iteration,data.funcCount,data.gradCount,data.fInitial,data.alpha); disp(s);
    end

    % Keep the variables for next iteration
    data.fOld=data.fInitial;
    data.xOld=data.xInitial;
    data.gOld=data.gradient;
  end
  % Set output parameters
  fval=data.fInitial;
  grad=data.gradient;
  x = data.xInitial;

  % Reshape x to original shape
  x=reshape(x,data.xsizes);

  % Call output function
  if(call_output_function(data,optim,'done')), exitflag=-1; end

  % Make exist output structure
  if(optim.HessUpdate(1)=='b'), output.algorithm='Broyden-Fletcher-Goldfarb-Shanno (BFGS)';
  elseif(optim.HessUpdate(1)=='l'), output.algorithm='limited memory BFGS (L-BFGS)';
  else output.algorithm='Steepest Gradient Descent';
  end
  output.message=getexitmessage(exitflag);
  output.iteration = data.iteration;
  output.funccount = data.funcCount;
  output.fval = data.fInitial;
  output.stepsize = data.alpha;
  output.directionalderivative = data.fPrimeInitial;
  output.gradient = reshape(data.gradient, data.xsizes);
  output.searchdirection = data.dir;
  output.timeTotal=toc(data.timeTotal);
  output.timeExtern=data.timeExtern;
  output.timeIntern=output.timeTotal-output.timeExtern;
  % Display final results
  if(~strcmp(optim.Display,'off'))
    disp(' Optimizer Results')
    disp(['  Algorithm Used: ' output.algorithm]);
    disp(['  Exit message : ' output.message]);
    disp(['  Iterations : '  int2str(data.iteration)]);
    disp(['  Function Count : ' int2str(data.funcCount)]);
    disp(['  Minimum found : ' num2str(fval)]);
    disp(['  Intern Time : ' num2str(output.timeIntern) ' seconds']);
    disp(['  Total Time : ' num2str(output.timeTotal) ' seconds']);
  end

function message=getexitmessage(exitflag)
switch(exitflag)
  case 1, message='Change in the objective function value was less than TolFun.';
  case 2, message='Change in x was smaller than the specified tolerance TolX.';
  case 3, message='Magnitude of gradient smaller than the specified tolerance';
  case 4, message='Boundary fminimum reached.';
  case 0, message='Number of iterations exceeded options.MaxIter or number of function evaluations exceeded options.FunEvals.';
  case -1, message='Algorithm was terminated by the output function.';
  case -2, message='Line search cannot find an acceptable point along the current search direction';
  otherwise, message='Undefined exit code';
end


function stopt=call_output_function(data,optim,where)
stopt=false;
if(~isempty(optim.OutputFcn))
  output.iteration = data.iteration;
  output.funccount = data.funcCount;
  output.fval = data.fInitial;
  output.stepsize = data.alpha;
  output.directionalderivative = data.fPrimeInitial;
  output.gradient = reshape(data.gradient, data.xsizes);
  output.searchdirection = data.dir;
  stopt=feval(optim.OutputFcn,reshape(data.xInitial,data.xsizes),output,where);
end


function data=linesearch_simple(funfcn, data, optim)
% Find a bracket of acceptable points
data = bracketingPhase_simple(funfcn, data, optim);

if (data.bracket_exitflag  == 2)
  % BracketingPhase found a bracket containing acceptable points;
  % now find acceptable point within bracket
  data = sectioningPhase_simple(funfcn, data, optim);
  data.exitflag = data.section_exitflag;
else
  % Already acceptable point found or MaxFunEvals reached
  data.exitflag = data.bracket_exitflag;
end

function data = bracketingPhase_simple(funfcn, data,optim)
% Number of iterations
itw=0;

% Point with smaller value, initial
data.beta=0;
data.f_beta=data.fInitial;
data.fPrime_beta=data.fPrimeInitial;

% Initial step is equal to alpha of previous step.
alpha = data.initialStepLength;

% Going uphill?
hill=false;

% Search for brackets
while(true)
  % Calculate the error registration gradient
  if isequal(optim.GradConstr,'on')
    [data,f_alpha]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data, optim);
    while isnan(f_alpha) || isinf(f_alpha)
      % Go to smaller stepsize
      alpha=alpha*optim.tau3;
      % Set hill variable
      hill=true;
      [data,f_alpha]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data,optim);
    end
    fPrime_alpha=nan;
    grad=nan;
  else
    [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data,optim);
    while isnan(f_alpha) || isinf(f_alpha) || any(isnan(grad)) || any(isinf(grad))
      % Go to smaller stepsize
      alpha=alpha*optim.tau3;
      % Set hill variable
      hill=true;
      [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data,optim);
    end
    fPrime_alpha = grad'*data.dir(:);
  end

  % Store values linesearch
  data.storefx=[data.storefx f_alpha];
  data.storepx=[data.storepx fPrime_alpha];
  data.storex=[data.storex alpha];
  data.storegx=[data.storegx grad(:)];

  % Update step value
  if(data.f_beta<f_alpha),
    % Go to smaller stepsize
    alpha=alpha*optim.tau3;

    % Set hill variable
    hill=true;
  else
    % Save current minimum point
    data.beta=alpha; data.f_beta=f_alpha; data.fPrime_beta=fPrime_alpha; data.grad=grad;
    if(~hill)
      alpha=alpha*optim.tau1;
    end
  end

  % Update number of loop iterations
  itw=itw+1;

  if(itw>(log(optim.TolFun)/log(optim.tau3))),
    % No new optimum found, linesearch failed.
    data.bracket_exitflag=-2; break;
  end

  if(data.beta>0&&hill)
    % Get the brackets around minimum point
    % Pick bracket A from stored trials
    [t,i]=sort(data.storex,'ascend');
    storefx=data.storefx(i);storepx=data.storepx(i); storex=data.storex(i);
    [t,i]=find(storex>data.beta,1);
    if(isempty(i)), [t,i]=find(storex==data.beta,1); end
    alpha=storex(i); f_alpha=storefx(i); fPrime_alpha=storepx(i);

    % Pick bracket B from stored trials
    [t,i]=sort(data.storex,'descend');
    storefx=data.storefx(i);storepx=data.storepx(i); storex=data.storex(i);
    [t,i]=find(storex<data.beta,1);
    if(isempty(i)), [t,i]=find(storex==data.beta,1); end
    beta=storex(i); f_beta=storefx(i); fPrime_beta=storepx(i);

    % Calculate derivatives if not already calculated
    if isequal(optim.GradConstr,'on')
      gstep=data.initialStepLength/1e6;
      if(gstep>optim.DiffMaxChange), gstep=optim.DiffMaxChange; end
      if(gstep<optim.DiffMinChange), gstep=optim.DiffMinChange; end
      [data,f_alpha2]=gradient_function(data.xInitial(:)+(alpha+gstep)*data.dir(:),funfcn, data, optim);
      [data,f_beta2]=gradient_function(data.xInitial(:)+(beta+gstep)*data.dir(:),funfcn, data, optim);
      fPrime_alpha=(f_alpha2-f_alpha)/gstep;
      fPrime_beta=(f_beta2-f_beta)/gstep;
    end

    % Set the brackets A and B
    data.a=alpha; data.f_a=f_alpha; data.fPrime_a=fPrime_alpha;
    data.b=beta; data.f_b=f_beta; data.fPrime_b=fPrime_beta;

    % Finished bracketing phase
    data.bracket_exitflag  = 2; return
  end

  % Reached max function evaluations
  if(data.funcCount>=optim.MaxFunEvals), data.bracket_exitflag=0; return; end
end


function data = sectioningPhase_simple(funfcn, data, optim)
% Get the brackets
brcktEndpntA=data.a; brcktEndpntB=data.b;

% Calculate minimum between brackets
[alpha,f_alpha_estimated] = pickAlphaWithinInterval(brcktEndpntA,brcktEndpntB,data.a,data.b,data.f_a,data.fPrime_a,data.f_b,data.fPrime_b,optim);
if(isfield(data,'beta')&&(data.f_beta<f_alpha_estimated)), alpha=data.beta; end


[t,i]=find(data.storex==alpha,1);
if((~isempty(i))&&(~isnan(data.storegx(i))))
  f_alpha=data.storefx(i); grad=data.storegx(:,i);
else
  % Calculate the error and gradient for the next minimizer iteration
  [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data,optim);
  if(isfield(data,'beta')&&(data.f_beta<f_alpha)),
    alpha=data.beta;
    if((~isempty(i))&&(~isnan(data.storegx(i))))
      f_alpha=data.storefx(i); grad=data.storegx(:,i);
    else
      [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data,optim);
    end
  end
end

% Store values linesearch
data.storefx=[data.storefx f_alpha]; data.storex=[data.storex alpha];

fPrime_alpha = grad'*data.dir(:);
data.alpha=alpha;
data.fPrime_alpha= fPrime_alpha;
data.f_alpha= f_alpha;
data.grad=grad;

% Set the exit flag to success
data.section_exitflag=[];


function data=linesearch(funfcn, data, optim)

% Find a bracket of acceptable points
data = bracketingPhase(funfcn, data,optim);
%if abs(data.a-data.b)<eps
%  data = bracketingPhase_simple(funfcn, data,optim);
%end

if (data.bracket_exitflag  == 2)
  % BracketingPhase found a bracket containing acceptable points;
  % now find acceptable point within bracket
  data = sectioningPhase(funfcn, data, optim);
  data.exitflag = data.section_exitflag;
else
  % Already acceptable point found or MaxFunEvals reached
  data.exitflag = data.bracket_exitflag;
end

function data = sectioningPhase(funfcn, data, optim)
%
% sectioningPhase finds an acceptable point alpha within a given bracket [a,b]
% containing acceptable points. Notice that funcCount counts the total number of
% function evaluations including those of the bracketing phase.

while(true)

  % Pick alpha in reduced bracket
  brcktEndpntA = data.a + min(optim.tau2,optim.sigma)*(data.b - data.a);
  brcktEndpntB = data.b - optim.tau3*(data.b - data.a);

  % Find global minimizer in bracket [brcktEndpntA,brcktEndpntB] of 3rd-degree
  % polynomial that interpolates f() and f'() at "a" and at "b".
  alpha = pickAlphaWithinInterval(brcktEndpntA,brcktEndpntB,data.a,data.b,data.f_a,data.fPrime_a,data.f_b,data.fPrime_b,optim);

  % No acceptable point could be found
  if (abs( (alpha - data.a)*data.fPrime_a ) <= data.TolFunLnS), data.section_exitflag = -2; return; end

  % Calculate value (and gradient if no extra time cost) of current alpha
  if ~isequal(optim.GradConstr,'on')
    [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data, optim);
    fPrime_alpha = grad'*data.dir(:);
  else
    gstep=data.initialStepLength/1e6;
    if(gstep>optim.DiffMaxChange), gstep=optim.DiffMaxChange; end
    if(gstep<optim.DiffMinChange), gstep=optim.DiffMinChange; end
    [data,f_alpha]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data,optim);
    [data,f_alpha2]=gradient_function(data.xInitial(:)+(alpha+gstep)*data.dir(:),funfcn, data, optim);
    fPrime_alpha=(f_alpha2-f_alpha)/gstep;
  end

  % Store values linesearch
  data.storefx=[data.storefx f_alpha]; data.storex=[data.storex alpha];

  % Store current bracket position of A
  aPrev = data.a;
  f_aPrev = data.f_a;
  fPrime_aPrev = data.fPrime_a;

  % Update the current brackets
  if ((f_alpha > data.fInitial + alpha*optim.rho*data.fPrimeInitial) || (f_alpha >= data.f_a))
    % Update bracket B to current alpha
    data.b = alpha; data.f_b = f_alpha; data.fPrime_b = fPrime_alpha;
  else
    % Wolfe conditions, if true then acceptable point found
    if (abs(fPrime_alpha) <= -optim.sigma*data.fPrimeInitial),
      if isequal(optim.GradConstr,'on')
        % Gradient was not yet calculated because of time costs
        [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data, optim);
        fPrime_alpha = grad'*data.dir(:);
      end
      % Store the found alpha values
      data.alpha=alpha; data.fPrime_alpha= fPrime_alpha; data.f_alpha= f_alpha;
      data.grad=grad;
      data.section_exitflag = []; return,
    end

    % Update bracket A
    data.a = alpha; data.f_a = f_alpha;  data.fPrime_a = fPrime_alpha;

    if (data.b - data.a)*fPrime_alpha >= 0
      % B becomes old bracket A;
      data.b = aPrev; data.f_b = f_aPrev;  data.fPrime_b = fPrime_aPrev;
    end
  end

  % No acceptable point could be found
  if (abs(data.b-data.a) < eps)
    if f_alpha<data.fInitial;
      % however, point with a smaller function value found
      data.alpha=alpha; data.fPrime_alpha= fPrime_alpha; data.f_alpha= f_alpha;
      data.grad=grad;
      data.section_exitflag = []; return
    else
      % No acceptable point could be found
      data.section_exitflag = -2; return
    end
  end

  % maxFunEvals reached
  if(data.funcCount >optim.MaxFunEvals), data.section_exitflag = -1; return, end
end

function data = bracketingPhase(funfcn, data, optim)
% bracketingPhase finds a bracket [a,b] that contains acceptable
% points; a bracket is the same as a closed interval, except that a
% > b is allowed.
%
% The outputs f_a and fPrime_a are the values of the function and
% the derivative evaluated at the bracket endpoint 'a'. Similar
% notation applies to the endpoint 'b'.

% Parameters of bracket A
data.a = [];
data.f_a = [];
data.fPrime_a = [];

% Parameters of bracket B
data.b = [];
data.f_b = [];
data.fPrime_b = [];

% First trial alpha is user-supplied
% f_alpha will contain f(alpha) for all trial points alpha
% fPrime_alpha will contain f'(alpha) for all trial points alpha
alpha = data.initialStepLength;
f_alpha = data.fInitial;
fPrime_alpha = data.fPrimeInitial;

% Set maximum value of alpha (determined by fminimum)
alphaMax = (data.fminimum - data.fInitial)/(optim.rho*data.fPrimeInitial);
alphaPrev = 0;

here_be_dragons=false;

while(true)
  % Evaluate f(alpha) and f'(alpha)
  fPrev = f_alpha;
  fPrimePrev = fPrime_alpha;

  % Calculate value (and gradient if no extra time cost) of current alpha
  if ~isequal(optim.GradConstr,'on')
    [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data, optim);
    while isnan(f_alpha) || isinf(f_alpha) || any(isnan(grad)) || any(isinf(grad))
      % NaN or Inf encountered, switch to safe mode
      here_be_dragons=true;
      alphaMax=alpha;
      alpha = alphaPrev+0.25*(alpha-alphaPrev);
      [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data, optim);
    end
    fPrime_alpha = grad'*data.dir(:);
  else
    gstep=data.initialStepLength/1e6;
    if(gstep>optim.DiffMaxChange), gstep=optim.DiffMaxChange; end
    if(gstep<optim.DiffMinChange), gstep=optim.DiffMinChange; end
    [data,f_alpha]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data, optim);
    while isnan(f_alpha) || isinf(f_alpha)
      % NaN or Inf encountered, switch to safe mode
      here_be_dragons=true;
      alphaMax=alpha;
      alpha = alphaPrev+0.25*(alpha-alphaPrev);
      [data,f_alpha]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data, optim);
    end
    [data,f_alpha2]=gradient_function(data.xInitial(:)+(alpha+gstep)*data.dir(:),funfcn, data, optim);
    fPrime_alpha=(f_alpha2-f_alpha)/gstep;
  end

  % Store values linesearch
  data.storefx=[data.storefx f_alpha]; data.storex=[data.storex alpha];

  % Terminate if f < fminimum
  if (f_alpha <= data.fminimum), data.bracket_exitflag = 4; return; end

  % Bracket located - case 0 (near NaN or Inf switch to safe solution)
  if here_be_dragons && (f_alpha >= fPrev)
    % a smaller function value was found on previous step
    data.a = 0; data.f_a = data.fInitial;  data.fPrime_a = data.fPrimeInitial;
    data.b = alpha; data.f_b = f_alpha;  data.fPrime_b = fPrime_alpha;
    % Finished bracketing phase
    data.bracket_exitflag  = 2; return
  end

  % Bracket located - case 1 (Wolfe conditions)
  if (f_alpha > (data.fInitial + alpha*optim.rho*data.fPrimeInitial)) || (f_alpha >= fPrev)
    % Set the bracket values
    data.a = alphaPrev; data.f_a = fPrev;  data.fPrime_a = fPrimePrev;
    data.b = alpha; data.f_b = f_alpha;  data.fPrime_b = fPrime_alpha;
    % Finished bracketing phase
    data.bracket_exitflag  = 2; return
  end

  % Acceptable step-length found
  if (abs(fPrime_alpha) <= -optim.sigma*data.fPrimeInitial),
    if isequal(optim.GradConstr,'on')
      % Gradient was not yet calculated because of time costs
      [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data, optim);
      fPrime_alpha = grad'*data.dir(:);
    end
    % Store the found alpha values
    data.alpha=alpha;
    data.fPrime_alpha= fPrime_alpha; data.f_alpha= f_alpha; data.grad=grad;
    % Finished bracketing phase, and no need to call sectioning phase
    data.bracket_exitflag = [];  return
  end

  % Bracket located - case 2
  if (fPrime_alpha >= 0)
    % Set the bracket values
    data.a = alpha; data.f_a = f_alpha;  data.fPrime_a = fPrime_alpha;
    data.b = alphaPrev; data.f_b = fPrev; data.fPrime_b = fPrimePrev;
    % Finished bracketing phase
    data.bracket_exitflag  = 2; return
  end

  % Update alpha
  if (2*alpha - alphaPrev < alphaMax )
    brcktEndpntA = 2*alpha-alphaPrev;
    brcktEndpntB = min(alphaMax,alpha+optim.tau1*(alpha-alphaPrev));
    % Find global minimizer in bracket [brcktEndpntA,brcktEndpntB] of 3rd-degree polynomial
    % that interpolates f() and f'() at alphaPrev and at alpha
    alphaNew = pickAlphaWithinInterval(brcktEndpntA,brcktEndpntB,alphaPrev,alpha,fPrev, ...
                                       fPrimePrev,f_alpha,fPrime_alpha,optim);
    alphaPrev = alpha;
    alpha = alphaNew;
  else
    alpha = alphaMax;
  end

  % maxFunEvals reached
  if(data.funcCount >optim.MaxFunEvals), data.bracket_exitflag = -1; return, end
end

function [alpha,f_alpha]= pickAlphaWithinInterval(brcktEndpntA,brcktEndpntB,alpha1,alpha2,f1,fPrime1,f2,fPrime2,optim)
% finds a global minimizer alpha within the bracket [brcktEndpntA,brcktEndpntB] of the cubic polynomial
% that interpolates f() and f'() at alpha1 and alpha2. Here f(alpha1) = f1, f'(alpha1) = fPrime1,
% f(alpha2) = f2, f'(alpha2) = fPrime2.

% determines the coefficients of the cubic polynomial with c(alpha1) = f1,
% c'(alpha1) = fPrime1, c(alpha2) = f2, c'(alpha2) = fPrime2.
coeff = [(fPrime1+fPrime2)*(alpha2-alpha1)-2*(f2-f1) ...
         3*(f2-f1)-(2*fPrime1+fPrime2)*(alpha2-alpha1) (alpha2-alpha1)*fPrime1 f1];

% Convert bounds to the z-space
lowerBound = (brcktEndpntA - alpha1)/(alpha2 - alpha1);
upperBound = (brcktEndpntB - alpha1)/(alpha2 - alpha1);

% Swap if lower bound is higher than the upper bound
if (lowerBound  > upperBound), t=upperBound; upperBound=lowerBound; lowerBound=t; end

% Find minima and maxima from the roots of the derivative of the polynomial.
sPoints = roots([3*coeff(1) 2*coeff(2) coeff(3)]);

% Remove imaginary and points outside range

sPoints(imag(sPoints)~=0)=[];
sPoints(sPoints<lowerBound)=[]; sPoints(sPoints>upperBound)=[];

% Make vector with all possible solutions
sPoints=[lowerBound sPoints(:)' upperBound];

% Select the global minimum point
[f_alpha,index]=min(polyval(coeff,sPoints)); z=sPoints(index);

% Add the offset and scale back from [0..1] to the alpha domain
alpha = alpha1 + z*(alpha2 - alpha1);

% Show polynomial search
if(optim.Display(1)=='p');
  vPoints=polyval(coeff,sPoints);
  plot(sPoints*(alpha2 - alpha1)+alpha1,vPoints,'co');
  plot([sPoints(1) sPoints(end)]*(alpha2 - alpha1)+alpha1,[vPoints(1) vPoints(end)],'c*');
  xPoints=linspace(lowerBound/3, upperBound*1.3, 50);
  vPoints=polyval(coeff,xPoints);
  plot(xPoints*(alpha2 - alpha1)+alpha1,vPoints,'c');
end


function [data,fval,grad]=gradient_function(x,funfcn, data, optim)
% Call the error function for error (and gradient)
if ( nargout <3 )
  timem=tic;
    fval=funfcn(reshape(x,data.xsizes));
    data.timeExtern=data.timeExtern+toc(timem);
    data.funcCount=data.funcCount+1;
else
  if(strcmp(optim.GradObj,'on'))
    timem=tic;
      [fval, grad]=feval(funfcn,reshape(x,data.xsizes));
      data.timeExtern=data.timeExtern+toc(timem);
      data.funcCount=data.funcCount+1;
      data.gradCount=data.gradCount+1;
  else
    % Calculate gradient with forward difference if not provided by the function
    grad=zeros(length(x),1);
    fval=funfcn(reshape(x,data.xsizes));
    gstep=data.initialStepLength/1e6;
    if(gstep>optim.DiffMaxChange), gstep=optim.DiffMaxChange; end
    if(gstep<optim.DiffMinChange), gstep=optim.DiffMinChange; end
    for i=1:length(x),
      x_temp=x; x_temp(i)=x_temp(i)+gstep;
      timem=tic;
        [fval_g]=feval(funfcn,reshape(x_temp,data.xsizes)); data.funcCount=data.funcCount+1;
        data.timeExtern=data.timeExtern+toc(timem);
        grad(i)=(fval_g-fval)/gstep;
    end
  end
  grad=grad(:);
end

function data = updateQuasiNewtonMatrix_LBFGS(data,optim)
% updates the quasi-Newton matrix that approximates the inverse to the Hessian.
% Two methods are support BFGS and L-BFGS, in L-BFGS the hessian is not
% constructed or stored.
% Calculate position, and gradient difference between the
% iterations
deltaX=data.alpha* data.dir;
deltaG=data.gradient-data.gOld;

if ((deltaX'*deltaG) >= sqrt(eps)*max( eps,norm(deltaX)*norm(deltaG) ))

  if(optim.HessUpdate(1)=='b')
    % Default BFGS as described by Nocedal
    p_k = 1 / (deltaG'*deltaX);
    Vk = eye(data.numberOfVariables) - p_k*deltaG*deltaX';
    % Set Hessian
    data.Hessian = Vk'*data.Hessian *Vk + p_k * deltaX*deltaX';
    % Set new Direction
    data.dir = -data.Hessian*data.gradient;
  else
    % L-BFGS with scaling as described by Nocedal

    % Update a list with the history of deltaX and deltaG
    data.deltaX(:,2:optim.StoreN)=data.deltaX(:,1:optim.StoreN-1); data.deltaX(:,1)=deltaX;
    data.deltaG(:,2:optim.StoreN)=data.deltaG(:,1:optim.StoreN-1); data.deltaG(:,1)=deltaG;

    data.nStored=data.nStored+1; if(data.nStored>optim.StoreN), data.nStored=optim.StoreN; end

    % Initialize variables
    a=zeros(1,data.nStored);
    p=zeros(1,data.nStored);

    q = data.gradient;
    for i=1:data.nStored
      p(i)= 1 / (data.deltaG(:,i)'*data.deltaX(:,i));
      a(i) = p(i)* data.deltaX(:,i)' * q;
      q = q - a(i) * data.deltaG(:,i);
    end
    % Scaling of initial Hessian (identity matrix)
    p_k = data.deltaG(:,1)'*data.deltaX(:,1) / sum(data.deltaG(:,1).^2);

    % Make r = - Hessian * gradient
    r = p_k * q;
    for i=data.nStored:-1:1,
      b = p(i) * data.deltaG(:,i)' * r;
      r = r + data.deltaX(:,i)*(a(i)-b);
    end

    % Set new direction
    data.dir = -r;
  end
end
