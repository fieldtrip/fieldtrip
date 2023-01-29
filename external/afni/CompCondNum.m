function [err,CondNum] = CompCondNum (N_tasks, N_basis, Task, N_TR, polort)
%
%   [err,] = CondNum ()
%
%Purpose:
%
%
%
%Input Parameters:
%
%
%
%Output Parameters:
%   err : 0 No Problem
%       : 1  Problems
%
%
%
%Key Terms:
%
%More Info :
%
%
%
%
%     Author : Gang Chen
%     Date : Wed Aug 20 15:46:08 EDT 2003
%     SSCC/NIMH/ National Institutes of Health, Bethesda Maryland


%Define the function name for easy referencing
FuncName = 'CondNum';

%Debug Flag
DBG = 1;

%initailize return variables
err = 1;
%N_tasks=4;
%N_basis = 1;

fprintf (1,'\n');
for (i=1:1:N_tasks),
   for (j=1:1:N_basis),
      fprintf(2, 'Loading %s\n', Task(i).Name(j).str);
      X = load (Task(i).Name(j).str);
      %length(X)
      Max(i) = max(X);
      if (i==1 & j==1),
         Xm = X(:);
         stmp = sprintf('%s ', Task(i).Name(j).str);
      else
         Xm = [Xm X(:)];
         stmp = sprintf('%s %s', stmp, Task(i).Name(j).str);
      end

   end
%     t = [0:1:(length(X)-1)];
%      plot (t, Xm);  title (stmp, 'Interpreter', 'none');
end
one = ones(N_TR,1);

%slope = 0;
for (i = 1:N_TR),  %create the coloumn(s) for the drifting
   for (j = 1:1:polort),
      slope(i, j) = (i-1)^j;
	end	
end
%design matrix
Dsgn = [one slope Xm];
%check condition number of the design matrix
%M = Dsgn'*Dsgn;
CondNum = cond(Dsgn'*Dsgn);

err = 0;
return;

