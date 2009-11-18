function [fv,grad] = logregn(w,data,targets,ptargets,nclasses,lambda,divnum)
w = reshape(w,nclasses,size(data,2)); 
w2=w(:,1:(end-1));
w2res = reshape(w2,fix(nclasses*size(w2,2)/divnum),divnum);%labels,freqnum,time divisions
nd1=diff(w2res,1,2);
nd=nd1.*nd1;
regularization_term=sum(nd(:));
nd1(:,1+size(nd1,2))=zeros(size(nd1,1),1);
nd1=-nd1-nd1;
nd2=zeros(size(nd1));
nd2(:,2:end)=-nd1(:,1:end-1);
diff_term=nd2+nd1;
diff_term=[reshape(diff_term,nclasses,size(w2,2)) zeros(nclasses,1)];

softmaxes = exp(data * w');
softmaxes = softmaxes./repmat(sum(softmaxes,2),[1 nclasses]);
    
% compute function value; i.e., training set loss
logsm=-log(softmaxes(targets));
logsm=sum(logsm);
fv = logsm+(lambda*regularization_term); 
% disp([logsm regularization_term])
grad = - (ptargets-softmaxes)'*data(:,1:size(w,2));
grad = grad(:)+(lambda.*diff_term(:)); 
