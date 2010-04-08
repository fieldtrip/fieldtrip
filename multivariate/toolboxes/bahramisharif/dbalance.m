function [newdata,newdesign] = dbalance(data,design)

   fprintf('balancing data\n');
     sz = size(data);
     data = reshape(data,[sz(1) prod(sz(2:end))]);
               
    nclasses = max(design(:,1));
    maxsmp = 0;
    for j=1:nclasses
        maxsmp = max(maxsmp,sum(design(:,1) == j));
    end
    newdata = zeros(nclasses*maxsmp,size(data,2));
    for j=1:nclasses
        cdata = data(design(:,1) == j,:);
        if  size(cdata,1) ~= maxsmp

            % sample with replacement
            newdata(((j-1)*maxsmp+1):(j*maxsmp),:) = cdata(ceil(size(cdata,1)*rand(maxsmp,1)),:);
        else
            newdata(((j-1)*maxsmp+1):(j*maxsmp),:) = cdata;
        end
    end
    newdesign=ones(nclasses*maxsmp,1);
    for j=2:nclasses
        newdesign(((j-1)*maxsmp+1):(j*maxsmp)) = j;
    end
end