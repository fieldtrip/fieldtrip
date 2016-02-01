function glmnetPlot( x, xvar, label )

%--------------------------------------------------------------------------
% glmnetPlot.m: plot coefficients from a "glmnet" object
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%    Produces a coefficient profile plot fo the coefficient paths for a
%    fitted "glmnet" object.
%
% USAGE: 
%    glmnetPlot(fit);
%    glmnetPlot(fit, xvar);
%    glmnetPlot(fit, xvar, label);
%
% INPUT ARGUMENTS:
% x           fitted "glmnet" model.
% xvar        What is on the X-axis. "norm" plots against the L1-norm of
%             the coefficients, "lambda" against the log-lambda sequence,
%             and "dev" against the percent deviance explained.
% label       if TRUE, label the curves with variable sequence numbers.
%
% DETAILS:
%    A coefficient profile plot is produced. If x is a multinomial model, a
%    coefficient plot is produced for each class.
%
% LICENSE: GPL-2
%
% DATE: 14 Jul 2009
%
% AUTHORS:
%    Algorithm was designed by Jerome Friedman, Trevor Hastie and Rob Tibshirani 
%    Fortran code was written by Jerome Friedman 
%    R wrapper (from which the MATLAB wrapper was adapted) was written by Trevor Hasite
%    MATLAB wrapper was written and maintained by Hui Jiang, jiangh@stanford.edu 
%    Department of Statistics, Stanford University, Stanford, California, USA.
%
% REFERENCES:
%    Friedman, J., Hastie, T. and Tibshirani, R. (2009)
%    Regularization Paths for Generalized Linear Models via Coordinate Descent.
%    Journal of Statistical Software, 33(1), 2010
%
% SEE ALSO:
%    glmnet, glmnetSet, glmnetPrint, glmnetPredict and glmnetCoef methods.
% 
% EXAMPLES:
%    x=randn(100,20);
%    y=randn(100,1);
%    g2=randsample(2,100,true);
%    g4=randsample(4,100,true);
%    fit1=glmnet(x,y);
%    glmnetPlot(fit1);
%    glmnetPlot(fit1, 'lambda', true);
%    fit3=glmnet(x,g4,'multinomial');
%    glmnetPlot(fit3);
%
% DEVELOPMENT: 14 Jul 2009: Original version of glmnet.m written.

if nargin < 2
    xvar = 'norm';
end

if nargin < 3
    label = false;
end

if strcmp(x.class,'multnet')
    beta=x.beta;
    if strcmp(xvar,'norm')
        norm = 0;
        for i=1:length(beta);
            which = nonzeroCoef(beta{i});
            beta{i} = beta{i}(which,:);
            norm = norm + sum(abs(beta{i}),1);
        end
    else 
        norm = 0;
    end
    dfmat=x.dfmat;
    ncl=size(dfmat,1);
    for i=1:ncl
        plotCoef(beta{i},norm,x.lambda,dfmat(i,:),x.dev,label,xvar,'',sprintf('Coefficients: Class %d', i));
    end
else 
    plotCoef(x.beta,[],x.lambda,x.df,x.dev,label,xvar,'','Coefficients');
end

%----------------------------------------------------------------
% End function glmnetPlot
%----------------------------------------------------------------

function plotCoef(beta,norm,lambda,df,dev,label,xvar,xlab,ylab)

which = nonzeroCoef(beta);
beta = beta(which,:);
if strcmp(xvar, 'norm')    
    if isempty(norm)
        index = sum(abs(beta),1);
    else
        index = norm;
    end
    iname = 'L1 Norm';
elseif strcmp(xvar, 'lambda')
    index=log(lambda);
    iname='Log Lambda';
elseif strcmp(xvar, 'dev')
    index=dev;
    iname='Fraction Deviance Explained';
end

if isempty(xlab)
    xlab = iname;
end

plot(index,transpose(beta));
xlabel(xlab);
ylabel(ylab);

%----------------------------------------------------------------
% End private function plotCoef
%----------------------------------------------------------------
