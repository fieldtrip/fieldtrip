function M = projectLinf1BlockMatrix(M,indices,lambda)

nGroups = length(indices);

for j=1:nGroups
   for i=1:nGroups
      x = M(indices{i},indices{j});

      if isempty(x), continue; end;

      [m,n] = size(x); x = x(:);
      if i == j
         y = sign(x).*min(abs(x),lambda(i+nGroups*(j-1)));
      else
         y = sign(x).*projectRandom2C(abs(x),lambda(i+nGroups*(j-1)));
      end

      M(indices{i},indices{j}) = reshape(y,m,n);
   end
end

