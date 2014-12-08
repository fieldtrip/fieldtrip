function se = strel_bol(r)

% STREL_BOL constructs a 3D sphere with the specified radius
% that can be used as structural element in 3D image processing
%
% See STREL, IMERODE, IMDILATE (image processing toolbox)

dim = [2*r+1, 2*r+1, 2*r+1];
se  = zeros(dim);
for i=1:dim(1)
  for j=1:dim(2)
    for k=1:dim(3)
      x = i-1-r;
      y = j-1-r;
      z = k-1-r;
      if norm([x y z])<=r
        se(i,j,k) = 1;
      end
    end
  end
end

