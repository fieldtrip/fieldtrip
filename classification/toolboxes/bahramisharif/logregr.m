function [fv,grad] = logregr(w,data,targets,ptargets,nclasses,lambda)

w = reshape(w,nclasses,size(data,2)); 
w2=w(:,1:(end-1));
regularization_term = abs(0.5*(w2(:)'*w2(:) ));
diff_term = w2 ;
diff_term = [diff_term zeros(nclasses,1)];

softmaxes = exp(data * w');
softmaxes = softmaxes./repmat(sum(softmaxes,2),[1 nclasses]);
    
% compute function value; i.e., training set loss
logsm=-log(softmaxes(targets));
logsm=sum(logsm);
fv = logsm+(lambda*regularization_term); 
grad = - (ptargets-softmaxes)'*data(:,1:size(w,2));
grad = grad(:)+(lambda.*diff_term(:)); 
