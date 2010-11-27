function [U,pivcol,nonpivcol] = ref(A,tol)
%REF    The (R)educed (E)chelon (F)orm of A.
%       Computes the reduced row echelon form of A using partial pivoting.
%
%       Formats:   U = ref(A)
%                  U = ref(A,tol)     User specifies tolerance to determine
%                                     when an entry should be zero.
%                  [U,pivcol] = ref(A)    Also lists the pivot columns of A.
%                  [U,pivcol,nonpivcol] = ref(A)   And the nonpivot columns.

%Written by David Lay, University of Maryland, College Park, 6/20/94
%Based on the original rref(A) program written by Cleve B. Moler in 1985.
%       Version 12/15/96

[m,n] = size(A);
tiny = max(m,n)*eps*max(1,norm(A,'inf'))*10;    %default tolerance for zeros
if (nargin==2), tiny = tol; end                %reset tolerance, if specified 
pivcol = [];
nonpivcol = [];
U = A;
i = 1;                                  %row index
j = 1;                                  %column index
while (i <= m) & (j <= n)
   [x,k] = max(abs(U(i:m,j))); p = k+i-1; %value and row index of next pivot.
   if (x <= tiny)                       % This column is negligible.
      U(i:m,j) = zeros(m-i+1,1);        % So clean up the entries.
      nonpivcol = [nonpivcol j];
      j = j + 1;                        % Pivot row p must be recalculated.
   else                     
      U([i p],j:n) = U([p i],j:n);      % Swap the ith and pth rows.
      U(i,j:n) = U(i,j:n)/U(i,j);       % Divide pivot row by the pivot.
      for k = [1:i-1 i+1:m]             % Replacement operations on other rows.
         U(k,j:n) = U(k,j:n) - U(k,j)*U(i,j:n);
      end
      pivcol = [pivcol j];
      i = i + 1;
      j = j + 1;
   end
end
nonpivcol = [nonpivcol j:n];

