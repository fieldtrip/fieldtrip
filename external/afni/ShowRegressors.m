function [err] = ShowRegressors (Task)
%
%   [err] = ShowRegressors ()
%
%Purpose:
%   
%   Gang
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
%     Author : Ziad Saad
%     Date : Fri Jul 18 15:24:43 EDT 2003
%     SSCC/NIMH/ National Institutes of Health, Bethesda Maryland


%Define the function name for easy referencing
FuncName = 'ShowRegressors';

%Debug Flag
DBG = 1;

%initailize return variables
err = 1;

figure(1); clf
N_tasks = length(Task);
%ln_style=['b' 'g' 'r' 'c' 'm' 'y' 'k' 'b-'];
m=colormap;
nm = rand(size(m,1),1); [nm, is] = sort(nm);
m(:,1) = m(is,1);
nm = rand(size(m,1),1); [nm, is] = sort(nm);
m(:,2) = m(is,2);
nm = rand(size(m,1),1); [nm, is] = sort(nm);
m(:,3) = m(is,3);

if (N_tasks),
   for (i=1:1:N_tasks),
      subplot (N_tasks,1,i); 
      for (j=1:1:size(Task(i).StimReg,1)),
%         plot (Task(i).StimReg', ln_style(j)); hold on
         plot (Task(i).StimReg(j,:)', 'color', m(j,:)); hold on
      end
      ylabel(sprintf('%s', Task(i).StimTime));
   end
end

figure(2);clf
if (N_tasks),
   for (i=1:1:N_tasks),
%      subplot (N_tasks,1,i); 
      for (j=1:1:size(Task(i).StimReg,1)),
		   plot (Task(i).StimReg(j,:)', 'color', m(i*j,:)); hold on
%         plot (Task(i).StimReg', ln_style(i)); hold on
      end
   end
end
title('All Regressors Together');

drawnow;

err = 0;
return;

