function [err,tasktest] = TaskTest (N_runs, polort, N_basis, N_tasks, Task, N_base)
%
%   [err,] = ContrTest ()
%
%Purpose:
%
% Generate matrix for the -glt option in 3dDevolve
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
%     Date : Mon Jul 21 18:25:31 EDT 2003
%     SSCC/NIMH/ National Institutes of Health, Bethesda Maryland


%Define the function name for easy referencing
FuncName = 'TaskTest';

%Debug Flag
DBG = 1;

%initailize return variables
err = 1;

%Create a matrix so that 3dDeconvolve uses for regressor test.

mtrx = zeros(N_tasks, (polort+1)*N_runs+N_basis*N_tasks+N_base);    %for each run, there are (polort+1) paramters in the baseline model
fprintf(1, 'Your task vectors for -glt options are:\n');

for (i = 1:1:N_tasks),
   fprintf(1,'No. %i task vector = ', i);
	tasktest(i).name = sprintf('%s_task.1D',Task(i).Label);
   fid = fopen(tasktest(i).name, 'w');
	for (j = 1:1:N_basis),
%		mtrx(i,(i-1)*N_basis+j+(polort+1)*N_runs) = 1;  %Assign one 1 for each basis function of every task
      mtrx((i-1)*N_basis+j,(i-1)*N_basis+j+(polort+1)*N_runs) = 1;
		for (k = 1:1:(polort+1)*N_runs+N_basis*N_tasks+N_base),
		   fprintf(1, '%g	', mtrx((i-1)*N_basis+j,k));
		   fprintf(fid, '%g ', mtrx((i-1)*N_basis+j,k));
	   end
      fprintf(fid, '\n');
	end
%	for (k = 1:1:(polort+1)*N_runs+N_basis*N_tasks),
%		fprintf(1, '%g	', mtrx(i,k));
%		fprintf(fid, '%g ', mtrx(i,k));
%	end
	fprintf(1, '\n');
	fclose(fid);	
end	

err = 0;
return;
