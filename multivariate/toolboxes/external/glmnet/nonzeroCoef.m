function result = nonzeroCoef(beta,bystep)

if nargin < 2
    bystep = false;
end

result = abs(beta)>0;    
if ~bystep
    result = any(result,2);
end

%-------------------------------------------------------------
% End private function nonzeroCoef
%-------------------------------------------------------------
