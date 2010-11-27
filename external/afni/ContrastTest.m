function [err,contrfile] = ContrastTest (N_runs, polort, N_basis, Task, N_tasks, N_contr, El, N_base)
%
%   [err,] = ContrastTest ()
%
%Purpose:
%   Construct the contrast test matrix used in -glt option of 3dDeconvolve.
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
%     Date : Tue Jul 29 13:01:05 EDT 2003
%     SSCC/NIMH/ National Institutes of Health, Bethesda Maryland


%Define the function name for easy referencing
FuncName = 'ContrastTest';

%Debug Flag
DBG = 1;

%initailize return variables
err = 1;

%Create a matrix so that 3dDeconvolve uses for regressor test.

length = (polort+1)*N_runs+N_base;
for (iT = 1:1:N_tasks),
   for (j = 1:1:N_basis)    ,
	   length = length + (Task(iT).BasisOpt_struct(j).maxlag - Task(iT).BasisOpt_struct(j).minlag + 1);
	end	
end

%mtrx = zeros(N_contr, (polort+1)*N_runs+N_basis*N_tasks);  %for each run, there are (polort+1) paramters in the baseline model

mtrx = zeros(N_contr, length);
fprintf(1, '\n\nYour contrast vectors for -glt options are:\n');

for (i = 1:1:N_contr),
	for (j = 1:1:El(i).cnt),
	   for (k = 1:1:N_tasks),
		   if El(i).order(j) == k,
			   for (m=1:1:N_basis),
				   for (p = 1:1:Task(k).BasisOpt_struct(m).maxlag - Task(k).BasisOpt_struct(m).minlag + 1)
%			         mtrx(i,m+(k-1)*N_basis+(polort+1)*N_runs) = El(i).sign(j);  
                  mtrx(i,p+(m-1)*(Task(k).BasisOpt_struct(m).maxlag - Task(k).BasisOpt_struct(m).minlag + 1)...
					      +(k-1)*(Task(k).BasisOpt_struct(m).maxlag - Task(k).BasisOpt_struct(m).minlag + 1)*N_basis...
					      +(polort+1)*N_runs) = El(i).sign(j);    %Don't ask me how I worked this out. It is so tedious and nasty 
							                                        %to housekeep the indices.
							                                        %Don't even remember how I figured it out now
					end
				end
	      end
		end  
	end
	contrfile(i).name = sprintf('%s_contr.1D',El(i).labels);
   fid = fopen(contrfile(i).name, 'w');
	fprintf(1,'No. %i contrast vector = ', i);
	for (n = 1:1:length),
		fprintf(1, '%g	', mtrx(i,n));
		fprintf(fid, '%g ', mtrx(i,n));
	end
	fprintf(1, '\n');
	fclose(fid);	
	fprintf(1, '\nThe vector for No. %i contrast is stored in file %s. \n\n', i, contrfile(i).name);
end	

err = 0;
return;
