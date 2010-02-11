
function [A,B]=platt_sigmoidest(svm_out,labels)

% PLATT_SIGMOIDEST function is based on Lin et al.'s implementation of Platt's 
% probabilistic output for SVM 
%
% INPUTS
%       svm_out - vector of outputs of the SVM on a validation set
%       labels  - validation labels
%
% OUTPUT
%       A, B    - coefficients for the posterior probability estimate
%                 P(y=1|x) = 1/(1+exp(A*f(x)+B)), where f(x) is the SVM output
%
% REFERENCES
% Lin, H.T. et al. (2003). A note on Platt's probabilistic outputs for 
%                          support vector machines (Technical Report). 

% Pawel Herman, 2009


class_labels = sort(unique(labels));
if length(class_labels) ~=2
    error('It should be a two-class problem');
elseif ~all(class_labels==[-1 1]')
    % labels are not [-1 1] so they need to be adjusted and svm_out has to be rescaled
    newlabels = zeros(size(labels));
    newlabels(labels==class_labels(1))= -1;
    newlabels(labels==class_labels(2))=  1;
    labels = newlabels;
    svm_out = 2 * (svm_out - 0.5*sum(class_labels))./(class_labels(2)-class_labels(1));  
end

% n_neg    - number of negative points
% n_pos    - number of positive points
n_neg = sum(labels<0);
n_pos = sum(labels>0);

%parameter setting
niter   = 100;     %Maximum number of iterations
minstep = 1e-10;   %Minimum step taken in line search
sigma   = 1e-12;   %Set to any value sufficient to make H' = H + sigma I always positive definite

%initialisation: target support in array t, initial function value in fval
hiTarget=(n_pos+1.0)/(n_pos+2.0);  loTarget=1/(n_neg+2.0);
t = (labels>=0)*hiTarget + (labels<0)*loTarget;
A = 0.0;
B = log((n_neg+1.0)/(n_pos+1.0));

%  t_i    * fApB_i + log(1.0 + exp(-fApB_i)) for fApB_i >= 0
% (t_i-1) * fApB_i + log(1.0 + exp(fApB_i))  for fApB_i < 0
fApB  = svm_out*A+B;
logF  = log(1.0+exp(-abs(fApB)));
t_aux = t - (fApB < 0);
fval  = sum(t_aux .* fApB + logF);

svm_out_sq = svm_out.*svm_out;

for it = 1:niter
   %update gradient and Hessian 
   expF = exp(-abs(fApB));
   oneexpFinv = (1.0+expF).^(-1);
   d2  = expF.*oneexpFinv.*oneexpFinv;
   d1  = t - max(expF, (fApB<0)).*oneexpFinv;
   h11 = sigma + sum(svm_out_sq.*d2);
   h22 = sigma + sum(d2);
   h21 = sum(svm_out.*d2);
   g1  = sum(svm_out.*d1);
   g2  = sum(d1);
   
   if (abs(g1)<1e-5 && abs(g2)<1e-5) %stopping criteria
      break;
   end;
   
   detinv=(h11*h22-h21*h21).^(-1);
   dA=-(h22*g1-h21*g2) .* detinv;
   dB=-(-h21*g1+h11*g2) .* detinv; %modified Newton direction
   gd=g1*dA+g2*dB;
   stepsize=1;
   
   while (stepsize >= minstep) %line search
      newA=A+stepsize*dA;
      newB=B+stepsize*dB;
      
      fApB  = svm_out*newA+newB;
      logF  = log(1.0+exp(-abs(fApB)));
      t_aux = t - (fApB < 0);
      newf  = sum(t_aux .* fApB + logF);  %update fval, compare with earlier initialisation
      
      if (newf < fval + 0.0001*stepsize*gd) %Check sufficient decrease
         A=newA; B=newB; fval=newf;
         break
      else
         stepsize = stepsize / 2.0;
      end;
      
      if (stepsize < minstep) %Line search fails
         disp('Line search fails'); 
         break
      end
   end
end