function sub = ind2subv(sz,index)
%IND2SUBV   Subscript vector from linear index.

csz = cumprod(sz(:)');
prev_csz = [1 csz(1:end-1)];
index = index(:) - 1;
sub = rem(repmat(index,1,length(sz)),repmat(csz,length(index),1));
sub = floor(sub ./ repmat(prev_csz,length(index),1))+1;
