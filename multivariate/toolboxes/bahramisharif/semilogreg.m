function [fv,g] = semilogreg(w,data,targets,ptargets,pidx,psz,nclasses,unlabl,nu,lambda)
% semisupervised logistic regression objective function

  nparts = length(data);

  weight=cell(1,nparts);
  for i=1:nparts
    weight{i} = w(pidx{i});
    weight{i} = reshape(weight{i},[nclasses psz(i)]);
  end

  fv=0;
  grad=cell(1,nparts);
  
  for i=1:nparts
    
    w1=weight{i};
    
    softmaxes = exp(data{i} * w1');
    softmaxes = bsxfun(@rdivide,softmaxes,sum(softmaxes,2));
  
    fv = fv - sum(log(softmaxes(targets)))+lambda*sum(w1(:).^2);
    
    ggrad = - (ptargets-softmaxes)'*data{i}(:,1:size(w1,2))+2*lambda*w1;
    grad{i}=ggrad(:);

  end
  
  if nu > 0 && ~isempty(unlabl)

    for i=1:length(unlabl)
      
      softmaxes1 = (unlabl{i} * weight{i}');

      for j=(i+1):length(unlabl)
        
        softmaxes2 = (unlabl{j} * weight{j}');

        S = (softmaxes1-softmaxes2);
        
        fv = fv + nu*sum(S(:).^2);
        
        gradd1 = 2*S'*unlabl{i};
        gradd2 = -2*S'*unlabl{j};
        
        grad{i} = grad{i}+nu*gradd1(:);
        grad{j} = grad{j}+nu*gradd2(:);
        
      end
    end
    
  end
  
  g=[];
  for i=1:length(grad)
    g = [g; grad{i}];
  end
  
end