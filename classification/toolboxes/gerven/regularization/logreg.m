function [fv,grad] = logreg(w,data,targets,ptargets,nclasses)
% This function is used in a call to third-party software not included in
% the released distribution.
%
% returns the function value and gradient given a weight matrix w for
% logistic regression with training data: input and targets

w = reshape(w,nclasses,size(data,2));

softmaxes = exp(data * w');
softmaxes = softmaxes./repmat(sum(softmaxes,2),[1 nclasses]);
fv = - sum(log(softmaxes(targets)));
grad = - (ptargets-softmaxes)'*data(:,1:size(w,2));
grad = grad(:);


% OBSOLETE CODE
% 
% w = reshape(w,nclasses,size(data,2));
% 
% softmaxes = exp(data * w');
% softmaxes = softmaxes./repmat(sum(softmaxes,2),[1 nclasses]);
%     
% % compute function value; i.e., training set loss
% fv = - sum(log(softmaxes(targets)));
% 
% % compute gradient
% grad = zeros(size(w));
% 
% 
% 
% for i=1:size(w,1)
%     for j=1:size(w,2)
%         grad(i,j) = - sum(data(:,j) .* (ptargets(:,i) - softmaxes(:,i)));
%     end
% end
% 
% grad = grad(:);



