function [col] = column_variation(row,n,N)

% COLUMN_VARIATION is an auxiliary function used by MYVARIATION for
% producing all possible combinations of input parameters
%
% SEE ALSO
% myvariation.m

% Pawel Herman, 2008

col=[];

for i=1:n
  blok(i,:)=row;
end

J = round(N/(n*length(row)));
for k=1:J
    for j=1:length(row)
        col = [ col ; blok(:,mod(j,length(row)+1))];     
    end
end

return
