function [pnew,x,f] = lnsrch(xold, fold, g, p, ofun)
% LNSRCH performs a line search to approximate a minimum of the objective function
% ofun by going along the direction p. xold and fold denote the old
% function value, and g the gradient. pnew, x, and f denote the updated step
% p, the new point in the direction of p, and the new value of the objective 
% functions. Described in Numerical Recipes in C pp. 385.
%
% Copyright (c) 2007, Marcel van Gerven
% F.C. Donders Centre for Cognitive Neuroimaging, Nijmegen, NL
%
% Notes:
%
% - the objective function expects the change in x and not xnew!
%

n = length(p);

thesum = sqrt(sum(p.^2));

if thesum > 100, p = p .* (100/thesum); end % 100 is maximum allowed step length

slope = sum(g .* p);

test = 0;

for i=1:n
   
    temp = abs(p(i)) / max(abs(xold(i)),1);
    
    if temp > test, test = temp; end
end

alamin = 1e-7/test; % convergence criterion on delta x
alam = 1; % full (Newton) step

while true
    
    pnew = alam * p; % new step    

    x = xold +  pnew; % make the step
       
    if alam < alamin
        
        x = xold;
        pnew = zeros(size(pnew));
        f = fold;
        
        %  we reach this place if our update becomes *very* small
        % and choose to do nothing in this case
        return;

    else

        % note that it is the change which acts as input to the objective
        % function!
        f = ofun(pnew);

        if f <= fold; % + 1e-4 * alam * slope % sufficient function decrease

            % this is where we want to be

            if isinf(f)
                error('infinite loss found in lnsrch');
            end

            return;

        else % backtrack

            if alam == 1 % first cycle

                tmplam = - slope / (2*(f-fold-slope));

            else % subsequent backtrack

                rhs1 = f - fold - alam * slope;
                rhs2 = f2 - fold2 - alam2  * slope;

                a = ((rhs1 / alam.^2) - (rhs2 / alam2.^2)) / (alam - alam2);
                b = ((-alam2*rhs1/alam.^2) + (alam*rhs2/alam2.^2)) / (alam - alam2);

                if a == 0
                    tmplam = -slope / (2*b);
                else
                    disc = b*b-3*a*slope;

                    if disc < 0
                        error('roundoff error in lnsrch');
                    else
                        tmplam = (-b+sqrt(disc))/(3*a);
                    end
                end

                if tmplam > 0.5*alam, tmplam = 0.5*alam; end

            end
        end
    end
    
    alam2 = alam;
    f2 = f;
    fold2 = fold;
    alam = max(tmplam,0.1*alam);
    
end
    
            
               
    
