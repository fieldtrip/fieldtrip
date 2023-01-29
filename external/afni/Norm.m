function [err, out_fn] = Norm (InBrik, N_runs, run, run_norm, run_mask, foutid_log)
%
%   [err, norm_fn] = Norm(InBrik, N_runs, run)
%
%Purpose:
%
%   Normalize the original dataset so that the beta values in 3dDeconvolve would be automatically percent signal change.
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
%     Date : Mon Oct 20 11:26:38 EDT 2003
%     SSCC/NIMH/ National Institutes of Health, Bethesda Maryland


%Define the function name for easy referencing
FuncName = 'Norm';

%Debug Flag
DBG = 1;

%initailize return variables
err = 1;

[InBrikPrefix, InBrikView, InBrikExt] =  AfniPrefix (InBrik);
rm_cmmd = 'rm -f ';
	
%Originally 3dClipLevel was used to generate a cliplevel for each run,
%then a step function chops off those voxel outside of the brain.
%But due to the concern of losing the sensitive amygdala and orbitofrontal cortex,
%it might be better to run 3dAutomask, which is roughly equivalent to
%3dClipLevel but would protect amygdala and orbitofrontal cortex regions
%from being masked.

%create command for executing 3dAutomask
if (run_mask.do == 1),
   automask_fn = sprintf('%s_msk', InBrikPrefix);
   automask = sprintf('3dAutomask -dilate %i -prefix %s %s%s', run_mask.N_dilate, automask_fn, InBrikPrefix, InBrikView);
   fprintf(2, '\nRunning: %s\n', automask);
   fprintf(foutid_log, '%s\n', automask);
   [s,w] = system(automask);
   if (s),
      fprintf(2,'Error during 3dAutomask: \n\n%s\n\n',w);
      while (1); fprintf(2,'Halted: Ctrl+c to exit'); pause; end
   end
	rm_cmmd = sprintf('%s %s%s%s', rm_cmmd, automask_fn, InBrikView, InBrikExt);
	
	%If no normalization, do masking
   if (run_norm == 0),
	   out_fn = sprintf('%s_mask', InBrikPrefix);
      mask_signal = sprintf('3dcalc -fscale -a  %s%s  \\\n', InBrikPrefix, InBrikView);
 	   mask_signal = sprintf('%s    -b %s%s \\\n', mask_signal, automask_fn, InBrikView);
	   mask_signal = sprintf('%s    -expr  "a*b" \\\n', mask_signal);
      mask_signal = sprintf('%s    -prefix %s', mask_signal, out_fn);
	   fprintf(2, '\nRunning: %s\n', mask_signal);
	   fprintf(foutid_log, '%s\n', mask_signal);
	   [s,w] = system(mask_signal);
      if (s),
         fprintf(2,'Error during 3dcalc: \n\n%s\n\n',w);
         while (1); fprintf(2,'Halted: Ctrl+c to exit'); pause; end
      end
   end
	rm_cmmd = sprintf('%s %s%s%s', rm_cmmd, automask_fn, InBrikView, InBrikExt);
end	
	
if (run_norm == 1),   %for normalization
   out_fn = sprintf('%s_norm', InBrikPrefix);
   concat = sprintf('3dTcat -prefix %s ', out_fn);
	
   for (i=1:1:N_runs),

   %store each mean dataset in a temporary file, and delete it before leaving this function.
      mean_fn = sprintf('%s_mean%i', InBrikPrefix, i);
	   get_mean = sprintf('3dTstat -mean -prefix %s %s%s''[%i..%i]''', mean_fn, InBrikPrefix, InBrikView, run(i).first, run(i).last);
      fprintf(2, '\nRunning: %s\n', get_mean);
	   fprintf(foutid_log, '%s\n', get_mean);
	   [s, w] = system(get_mean);
	   if (s),
         fprintf(2,'Error at running 3dTstat: \n\n\n%s\n', w);
         while (1); fprintf(2,'Halted: Ctrl+c to exit');pause; end
      end
		
	%normalizing
%	norm_signal = sprintf('3dcalc -datum float -a  %s%s''[%i..%i]''  \\\n', InBrikPrefix, InBrikView, run(i).first, run(i).last);
	   norm_signal = sprintf('3dcalc -fscale -a  %s%s''[%i..%i]''  \\\n', InBrikPrefix, InBrikView, run(i).first, run(i).last);
   	norm_signal = sprintf('%s    -b %s%s \\\n', norm_signal, mean_fn, InBrikView);
		if (run_mask.do == 1),
		   norm_signal = sprintf('%s    -c %s%s \\\n', norm_signal, automask_fn,InBrikView);
		%	norm_signal = sprintf('%s    -expr  "(a/b*100)*step(b-clip_level)" \\\n', norm_signal);
   	   norm_signal = sprintf('%s    -expr  "(a/b*100)*c" \\\n', norm_signal);
		else norm_signal = sprintf('%s    -expr  "a/b*100" \\\n', norm_signal);
		end
	   norm_out_fn = sprintf('%s_norm%i', InBrikPrefix, i);
   	norm_signal = sprintf('%s    -prefix %s', norm_signal, norm_out_fn);
		fprintf(2, '\nRunning: %s\n', norm_signal);
	   fprintf(foutid_log, '%s\n', norm_signal);
	   [s,w] = system(norm_signal);
      if (s),
         fprintf(2,'Error during 3dcalc: \n\n%s\n\n',w);
         while (1); fprintf(2,'Halted: Ctrl+c to exit'); pause; end
      end
	
	%Concatenate the above files together
	   concat = sprintf('%s %s%s%s', concat, norm_out_fn, InBrikView, InBrikExt);
	
	%create remove command for deleting those individual percent signal change and mean file for each run.
	   rm_cmmd = sprintf('%s %s%s%s %s%s%s', rm_cmmd, norm_out_fn, InBrikView, InBrikExt, mean_fn, InBrikView, InBrikExt);	

   end

%Concatenate all runs back together
   fprintf(2, '\nRunning: %s\n', concat);
   fprintf(foutid_log, '%s\n', concat);
   [s,w] = system(concat);
   if (s),
      fprintf(2,'Error during running 3dTcat: \n\n%s\n\n',w);
      while (1); fprintf(2,'Halted: Ctrl+c to exit'); pause; end
   end	
end

%Don't forget to remove those temporary files above.
fprintf(2, '\nRemoving: %s\n', rm_cmmd);
fprintf(foutid_log, '%s\n', rm_cmmd);
[s,w] = system(rm_cmmd);
if (s),
    fprintf(2,'Error during rm: \n\n%s\n\n',w);
    while (1); fprintf(2,'Halted: Ctrl+c to exit'); pause; end
end

err = 0;
return;
