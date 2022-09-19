function value= rouderS(g,y,X,grouping,continuousIx,options)
% The S(g) function of Eq 9 in Rouder et al.
% g = Matrix of g values, each row is a dimension(effect), each column a value
% that we're integrating over.
% y = Data values
% X = design matrix, without a constant term.
% OUTPUT
% value = Value of s [nrObservations 1]



[nrDims,nrPriorValues] = size(g);
[N,nrEffects] = size(X);

if ~isempty(continuousIx)
    contX = X(:,continuousIx);
    contGMatrix = inv(contX'*contX/N);
end
% NOTE It would be nice to evaluate this on a GPU but the builtin gpuArray
% cannot compute all, diag, det, etc, on the GPU using the gpuArray/arrayfun
% Just putting all variables on the gpuArray (and not using arrayfun) 
% actually slows things down as a lot of copying takes place...

nrObservations = size(X,1);
one = ones(nrObservations,1);
P0 = 1./nrObservations*(one*one');
yTilde = (eye(nrObservations)-P0)*y;
XTilde = (eye(nrObservations)-P0)*X;
value= nan(1,nrPriorValues);
%parfor (i=1:nrPriorValues,options.nrWorkers)
for i=1:nrPriorValues%,options.nrWorkers)
    if all(g(:,i)==0)
        value(i)=0;
    else
        G=[];
        for grp= 1:nrDims
            nrInDim = numel(grouping{grp});            
            if all(ismember(grouping{grp},continuousIx))
                % Continuous covariate
                thisG = g(grp,i).*contGMatrix;
            else
                thisG = g(grp,i)*eye(nrInDim);
            end
            G = blkdiag(G,thisG);
        end        
        Vg = XTilde'*XTilde + inv(G);
        yBar = one'*y/nrObservations;
        preFactor= 1./(sqrt(det(G))*sqrt(det(Vg)));
        numerator =    y'*y-nrObservations*yBar^2;
        denominator = ((yTilde'*yTilde) -yTilde'*XTilde*(Vg\XTilde'*yTilde));
        value(i)= preFactor*(numerator/denominator).^((nrObservations-1)/2);
    end
end
end


