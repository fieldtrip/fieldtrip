function diff = derivativeCheck(funObj,x,order,useComplex,varargin)
% diff = derivativeCheck(funObj,x,order,useComplex,varargin)

if nargin < 3
    order = 1;
    if nargin < 4
        useComplex = 0;
    end
end

if order > 1
    [f,g,H] = funObj(x,varargin{:});
else
    [f,g] = funObj(x,varargin{:});
end

fprintf('Checking Gradient:\n');
[f2,g2] = autoGrad(x,useComplex,funObj,varargin{:});

fprintf('Max difference between user and numerical gradient: %f\n',max(abs(g-g2)));
if max(abs(g-g2)) > 1e-4
    fprintf('User NumDif:\n');
    [g g2]
    diff = abs(g-g2)
    pause
end

if order > 1
    fprintf('Check Hessian:\n');
    [f2,g2,H2] = autoHess(x,useComplex,funObj,varargin{:});
    
    fprintf('Max difference between user and numerical hessian: %f\n',max(abs(H(:)-H2(:))));
    if max(abs(H(:)-H2(:))) > 1e-4
        H
        H2
        diff = abs(H-H2)
        pause;
    end
end