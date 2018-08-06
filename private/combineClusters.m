function cluster = combineClusters(labelmat, spatdimneighbstructmat, total)

% COMBINECLUSTERS is a helper function for FINDCLUSTER. It searches for
% adjacent clusters in neighbouring channels and combines them.

% A mex-file is available for this function, as it can take quite long.

replaceby=1:total;
spatdimlength = size(labelmat, 1);
for spatdimlev=1:spatdimlength
  neighbours=find(spatdimneighbstructmat(spatdimlev,:));
  for nbindx=neighbours
    indx = find((labelmat(spatdimlev,:)~=0) & (labelmat(nbindx,:)~=0));
    for i=1:length(indx)
      a = labelmat(spatdimlev, indx(i));
      b = labelmat(nbindx, indx(i));
      if replaceby(a)==replaceby(b)
        % do nothing
        continue;
      elseif replaceby(a)<replaceby(b)
        % replace all entries with content replaceby(b) by replaceby(a).
        replaceby(find(replaceby==replaceby(b))) = replaceby(a); 
      elseif replaceby(b)<replaceby(a)
        % replace all entries with content replaceby(a) by replaceby(b).
        replaceby(find(replaceby==replaceby(a))) = replaceby(b); 
      end
    end
  end
end

% renumber all the clusters
num = 0;
cluster = zeros(size(labelmat));
for uniquelabel=unique(replaceby(:))'
  num = num+1;
  cluster(ismember(labelmat(:),find(replaceby==uniquelabel))) = num;
end

end
