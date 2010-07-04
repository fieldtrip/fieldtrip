function v = makemodel(p,ng)
%MAKEMODEL Helper function to make model matrix
%    P = max model term size, NG = number of grouping variables
% or
%    P = vector of term codes

% We want a matrix with one row per term.  Each row has a 1 for
% the variables participating in the term.  There will be rows
% summing up to 1, 2, ..., p.

if numel(p)==1
   % Create model matrix from a scalar max order value
   vgen = 1:ng;
   v = eye(ng);                      % linear terms
   for j=2:min(p,ng)
      c = nchoosek(vgen,j);          % generate column #'s with 1's
      nrows = size(c,1);             % generate row #'s
      r = repmat((1:nrows)',1,j);    %    and make it conform with c
      m = zeros(nrows,ng);           % create a matrix to get new rows
      m(r(:)+nrows*(c(:)-1)) = 1;    % fill in 1's
      v = [v; m];                    % append rows
   end
else
   % Create model matrix from terms encoded as bit patterms
   nterms = length(p);
   v = zeros(nterms,ng);
   for j=1:nterms
      tm = p(j);
      while(tm)
         % Get last-numbered effect remaining
         lne = 1 + floor(log2(tm));
         tm = bitset(tm, lne, 0);
         v(j,lne) = 1;
      end
   end
end
