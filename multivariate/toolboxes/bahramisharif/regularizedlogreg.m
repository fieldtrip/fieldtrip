function [fv,grad] = regularizedlogreg(w,data,targets,ptargets,nclasses,lambda,divnum)
%   Copyright (c) 2008, Marcel van Gerven, Ali Bahramisharif

w = reshape(w,nclasses,size(data,2)); 

%%regularization of variation
w2=w(:,1:(end-1));%remove the bias term
w2res = reshape(w2,fix(nclasses*size(w2,2)/divnum),divnum);%classes*labels*freqnum,time divisions
mmm = mean(w2res,2);%average over segments
regularization_term = abs(0.5*(w2(:)'*w2(:) - divnum*(mmm(:)'*mmm(:))));
diff_term = w2 - reshape(repmat(mmm,1,divnum),nclasses,size(w2,2));
diff_term = [diff_term zeros(nclasses,1)];

%%loss term for LR
softmaxes = exp(data * w');
softmaxes = softmaxes./repmat(sum(softmaxes,2),[1 nclasses]);
    
logsm=-log(softmaxes(targets));
logsm=sum(logsm);

fv = logsm+(lambda*regularization_term); 
grad = - (ptargets-softmaxes)'*data(:,1:size(w,2));
grad = grad(:)+(lambda.*diff_term(:)); 
