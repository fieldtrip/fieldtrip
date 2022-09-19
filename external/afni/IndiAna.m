%script IndiAna
%
%
%
%Purpose:
%   Run 3dDeconvolve in an interactive user-friendly way.
%
%Input:
%   Major ones are done through keyboard, while rare ones such as -censor , -nforist, etc. are left for user to implement.
%
%Output:
%	Wave functions from waver;
%  Various output files from 3dDeconvolve
%
%Key Terms:
%
%More Info :
%
%     Author : Gang Chen (with help from Ziad Saad)
%     Date : Fri Jul 18 13:58:44 EDT 2003
%     SSCC/NIMH/ National Institutes of Health, Bethesda MD 20892

clear all;

FuncName = 'Decon';

fprintf(1,'\nOK, let''s roll ...\n\n');

%parameters common to all regressors

%EPI time seris brick

flg1 = 0;    %flag for file existence
while flg1 == 0,
   InBrik = input('Input dataset (full filename with suffix BRIK): ','s');
   fid = fopen (InBrik,'r');
   if (fid == -1),
	   flg1 = 0; fprintf(2,'Error: File %s does not exist. Please try it again. \n', InBrik);
   else flg1 = 1; fclose (fid);
   end
end

%Use Ziad's AfniPrefix.m, which uses RemoveExtension.m
[InBrikPrefix, InBrikView, InBrikExt] =  AfniPrefix (InBrik);

flg2 = 0;
%Repition time                 %TR in seconds
while flg2 == 0,
   TR = input('TR (second): ');
	if (isnumeric(TR) == 0 | isempty(TR)),
	   flg2 = 0; fprintf(2,'Error: the input is not a number. Please try it again.\n');
	else flg2 = 1;
   end
end

flg1 = 0;    %flag for censor file
while flg1 == 0,
   Censor = input('\nNeed to censor any time point(s)? (0 - No; 1 - Yes):  ');
	if (Censor ~= 1 & Censor ~= 0),
	   flg1 = 0; fprintf(2,'Error: the input is not a valid number. Please try it again.\n');
	else flg1 = 1; end
end

if (Censor == 1),
   flg1 = 0;
	while flg1 == 0,
	   Censor_fn = input('Censor file name: ','s');
      fid2 = fopen (Censor_fn,'r');
      if (fid2 == -1),
	      flg1 = 0; fprintf(2,'Error: File %s does not exist. Please try it again. \n', Censor_fn);
      else flg1 = 1; fclose (fid2); end
   end
end

flg3=0;
%Number of tasks/conditions/regressors to be tested
while flg3 == 0,
   N_tasks = input('Number of stimuli/tasks/conditions/regressors: ');
	if (isnumeric(N_tasks) == 0 | isempty(N_tasks)),
	   flg3 = 0; fprintf(2,'Error: the input is not a number. Please try it again.\n');
	else flg3 = 1;
   end
end

flg4 = 0;
%Number of runs concatenated
while flg4 == 0,
   N_runs = input('Number of concatenated runs: ');
	if (isnumeric(N_runs) == 0 | isempty(N_runs)),
	   flg4 = 0; fprintf(2,'Error: the input is not a number. Please try it again.\n');
	else flg4 = 1;
   end
end
	
if (N_runs == 1),	
   run(1).first = 0;
end
	
%Get the concatenation points and create a concatenation file
if (N_runs > 1),
   fprintf('After concatenation, the 1st run usually starts at 0.\n');
	Concat_fn = sprintf('%s_concat.1D', InBrikPrefix);
   Concat_id = fopen(Concat_fn,'w');
   if (Concat_id < 0),
      fprintf(2,'Error in creating concatenation file: Failed to write %s to disk.\n', Concat_fn);
      while (1); fprintf(2,'Halted: Ctrl+c to exit');pause; end
   end
%   run(1).first = 0;
   for (i=1:1:N_runs),
	   fprintf('Number %i run starts ', i);
		flg5 = 0;
		while flg5 == 0,
   	   run(i).first = input('at (TR):  ');
			if (isnumeric(run(i).first) == 0 | isempty(run(i).first)),
			   flg5 = 0; fprintf(2,'Error: the input is not a number. Please try it again.\n');
	      else flg5 = 1;
         end
      end
		flg5 = 0;
%	   fprintf('Number %i run ends ', i);
%	   run(i).last = input('at (TR):  ');
%		run(i+1).first = run(i).last + 1;
		fprintf(Concat_id, '%i\n', run(i).first);
	end
end

%Numer of TR's for the whole dataset including concatenated ones
flg6 = 0;
while flg6 == 0,
   N_TR = input('Number of TRs for the whole dataset (run 3dinfo on the dataset if you are not sure): ');
	if (isnumeric(N_TR) == 0 | isempty(N_TR)),
	   flg6 = 0; fprintf(2,'Error: the input is not a number. Please try it again.\n');
	else flg6 = 1;
   end
end

for (i=1:1:N_runs),
	if (i<N_runs), run(i).last = run(i+1).first - 1;
	else run(i).last = N_TR-1; end
end

flg7 = 0;
while flg7 == 0,
   fprintf('\nNormalization would automatically covert the regressor coefficients into percent signal change.');
   run_norm = input('\nWant to normalize the dataset before running deconvolution/regression? (0 - No; 1 - Yes): ');
	if (run_norm ~= 0 & run_norm ~= 1),
	   flg7 = 0; fprintf(2,'Error: the input has to be 0 or 1. Please try it again.\n');
	else flg7 = 1;
   end
end

flg7 = 0;
while flg7 == 0,
   fprintf('\n3dAutomask uses 3dClipLevel algorithm to find clipping level.');
   run_mask.do = input('\nWant to mask the outside of brain? (0 - No; 1 - Yes): ');
	if (run_mask.do ~= 0 & run_mask.do ~= 1),
	   flg7 = 0; fprintf(2,'Error: the input has to be 0 or 1. Please try it again.\n');
	else flg7 = 1;
	end
end

if (run_mask.do == 1),
   flg7 = 0;
   while flg7 == 0,
	   run_mask.N_dilate = input('How many times do you want to dilate the mask outwards? (default - 5): ');
	   if isempty(run_mask.N_dilate), run_mask.N_dilate = 5; flg7 = 1;
		   elseif  (isnumeric(run_mask.N_dilate) == 0),
	      flg7 = 0; fprintf(2,'Error: the input is not a number. Please try it again.\n');
	      else flg7 = 1;
      end
	end	
end

flg8 = 0;
%polort = 2;             %Read the 3dDeconvolve manual ...
while flg8 == 0,
   polort = input('\nOrder of polynormial fitting for baseline (-1 - no basline; 0 - constant; 1 - linear; 2 - quadratic; ...): ');
	if (isnumeric(polort) == 0 | isempty(polort)),
	   flg8 = 0; fprintf(2,'Error: the input is not a number. Please try it again.\n');
	else flg8 = 1;
   end
end

%Analysis type
fprintf('\nYou can run 2 types of analysis: \n\n');
fprintf('(1) Deconvolution + Regression: generating impulse response functions and condition components; \n');
fprintf('(2) Regression: breaking the signal down into various regressor components. \n');

flg9 = 0;
while flg9 == 0,
   Run_Type = input('\nWhat type of analysis? 1 - Deconvolution; 2 - Regression: ');
	if (Run_Type ~= 1 & Run_Type ~= 2),
	   flg9 = 0; fprintf(2,'Error: the input is not a number. Please try it again.\n');
	else flg9 = 1;
   end
end

fprintf('\n');
if (Run_Type == 1),     %for deconvolution
   N_basis = 1;
   for (iT=1:1:N_tasks),
      flg10 = 0;
	   while flg10 == 0,
	      fprintf('Minimum lags for number %i stimulus ',iT);
%         Task(iT).MinLags = input('is: ');
			Task(iT).BasisOpt_struct(1).minlag = input('is: ');
%         if (isnumeric(Task(iT).MinLags) == 0 | isempty(Task(iT).MinLags)),
			if (isnumeric(Task(iT).BasisOpt_struct(1).minlag) == 0 | isempty(Task(iT).BasisOpt_struct(1).minlag)),
	         flg10 = 0; fprintf(2,'Error: the input is not a number. Please try it again. \n');
         else flg10 = 1;
         end
		end	
      flg10 = 0;
	   while flg10 == 0,
	      fprintf('Maximum lags for number %i stimulus ',iT);
%         Task(iT).MaxLags = input('is: ');
			Task(iT).BasisOpt_struct(1).maxlag = input('is: ');
%         if (isnumeric(Task(iT).MaxLags) == 0 | isempty(Task(iT).MaxLags)),
         if (isnumeric(Task(iT).BasisOpt_struct(1).maxlag) == 0 | isempty(Task(iT).BasisOpt_struct(1).maxlag)),
	         flg10 = 0; fprintf(2,'Error: the input is not a number. Please try it again. \n');
         else flg10 = 1;
         end
		end	
   end	
end

if (Run_Type == 2),     %for regression
   fprintf('\nAvailable basis function options are: \n\n');
   fprintf('GAM  -- Gamma variate \n');
   fprintf('tent -- tent function \n');
   fprintf('MGH  -- used by Massachusetts General Hospital group \n');
	fprintf('SPM1 -- double canonical  Gamma variates \n\n');
   fprintf('SPM2 -- mixture of Gamma variates with time derivative \n\n');

   BasisFunc = input('Basis function type (GAM, tent, MGH, SPM1, SPM2, ...): ', 's');

%There is a scale factor 1/1.25^2 with MGH, but it is OK with waver, and -peak 1 should take care of it.
   if (strcmp (lower(BasisFunc),'gam') | strcmp (lower(BasisFunc),'mgh')),
	
	   flg10 = 0;
		while flg10 == 0,
		   fprintf('\nThe default power paratmeter is: 8.6 - GAM; 2 - MGH.\n');
         BasisOpt.gamb = input('Input the Power parameter if different from the default value, otherwise hit Enter: ');
			if isempty(BasisOpt.gamb),
			   if (strcmp (lower(BasisFunc),'gam')),  BasisOpt.gamb = 8.6; end
				if (strcmp (lower(BasisFunc),'mgh')),  BasisOpt.gamb = 2.0; end
				flg10 = 1;
			elseif isnumeric(BasisOpt.gamb) == 0,
			   flg10 = 0; fprintf(2,'Error: the input is not a number. Please try it again.\n');
	      else flg10 = 1;
         end
      end
		
		flg10 = 0;
	   while flg10 == 0,
		   fprintf('\nThe default scaling paratmeter is: 0.547 - GAM; 1.25 - MGH.\n');
         BasisOpt.gamc = input('Input the scaling parameter if different from the default value, otherwise hit Enter:');
			if isempty(BasisOpt.gamc),
			   if (strcmp (lower(BasisFunc),'gam')), BasisOpt.gamc = 0.547; end
				if (strcmp (lower(BasisFunc),'mgh')), BasisOpt.gamc = 1.25; end
				flg10 = 1;
			elseif isnumeric(BasisOpt.gamc) == 0,
			   flg10 = 0;
		      fprintf(2,'Error: the input is not a number. Please try it again.\n');
	      else flg10 = 1;
         end
      end
		
		flg10 = 0;
	   while flg10 == 0,
		   fprintf('\nThe default delay paratmeter is: 0 - GAM; 2.25 - MGH.\n');		
         BasisOpt.gamd = input('Input the delay parameter if different from the default value, otherwise hit Enter:');
			if isempty(BasisOpt.gamd),
			   if (strcmp (lower(BasisFunc),'gam')), BasisOpt.gamd = 0.0; end
				if (strcmp (lower(BasisFunc),'mgh')), BasisOpt.gamd = 2.25; end
				flg10 = 1;
			elseif isnumeric(BasisOpt.gamd) == 0,
			   flg10 = 0; fprintf(2,'Error: the input is not a number. Please try it again.\n');
	      else
            flg10 = 1;
         end
      end
			
      flg10 = 0;
		fprintf('\nCurrently in each condition minlag and maxlag are set to be the same for all basis functions.\n');
		fprintf('The default for both is set to 0. Hit Enter if you want to set them 0.\n');
	   while flg10 == 0,
	      BasisOpt.minlag = input('Minimum lags (default: 0): ');   %Right now same lags for all regressors and basis functions
			if isempty(BasisOpt.minlag),
            BasisOpt.minlag = 0;
				flg10 = 1;
			elseif isnumeric(BasisOpt.minlag) == 0,
			   flg10 = 0; fprintf(2,'Error: the input is not a number. Please try it again.\n');
	      else flg10 = 1;
         end
      end
		
		flg10 = 0;
	   while flg10 == 0,
	      BasisOpt.maxlag = input('Maximum lags (default: 0): ');   %Right now same lags for all regressors and basis functions
			if isempty(BasisOpt.maxlag),
            BasisOpt.minlag = 0;
				flg10 = 1;
			elseif isnumeric(BasisOpt.maxlag) == 0,
			   flg10 = 0; fprintf(2,'Error: the input is not a number. Please try it again.\n');
	      else flg10 = 1;
         end
      end		
	
   end
		
   flg10 = 0;
   while flg10 == 0,
	   if (strcmp (lower(BasisFunc),'gam') | strcmp (lower(BasisFunc),'mgh') | strcmp (lower(BasisFunc),'spm1')), N_basis = 1;
		   elseif (strcmp (lower(BasisFunc),'spm2')), N_basis = 2;
         else N_basis = input('Number of basis functions: ');
		end
		if (isnumeric(N_basis) == 0 | isempty(N_basis)),
		   flg10 = 0; fprintf(2,'Error: the input is not a number. Please try it again.\n');
	   else flg10 = 1;
      end
   end

   %Total time span of all basis functions (i.e. duration of IRF) in seconds
   BasisOpt.tSpan = 0;
   if ((N_basis > 1) & (strcmp (lower(BasisFunc),'tent'))),
	   flg10 = 0;
		while flg10 == 0,
         BasisOpt.tSpan = input('Total time span of all basis functions in this regressor (sec): ');
			if (isnumeric(BasisOpt.tSpan) == 0 | isempty(BasisOpt.tSpan)),
		      flg10 = 0; fprintf(2,'Error: the input is not a number. Please try it again.\n');
	      else flg10 = 1;
         end
      end
   end	

   if (N_basis>1),
%Recon_dt = 0.1;         %Sampling period for reconstruction IRFs
      flg10 = 0;
		while flg10 == 0,
         Recon_dt = input('Sampling period for reconstruction IRFs, i.e., 0.1, (second): ');
			if (isnumeric(Recon_dt) == 0 | isempty(Recon_dt)),
		      flg10 = 0; fprintf(2,'Error: the input is not a number. Please try it again.\n');
	      else flg10 = 1;
         end
      end
   end
end

flg10 = 0;
while flg10 == 0,
   fprintf('\nThere are 2 kinds of stimulus file:\n');
	fprintf('1. Regular timing indexed with 0s and 1s (TR-locked));\n');
	fprintf('2. Irregular timing (not synchronized with TR).\n');
   Stim_Type = input('\nType of stimulus time files (1 - Regular; 2 - Irregular): ');
	if (Stim_Type ~= 1 & Stim_Type ~= 2),
	   flg10 = 0; fprintf(2,'Error: the input has to be 0 or 1. Please try it again.\n');
	else flg10 = 1;
	end
end

if (Stim_Type == 2),   %Irregular timing
   flg10 = 0;
   while flg10 == 0,
	   fprintf('Specify the precision of your timing files. For example, if time points reach the \n');
		fprintf('100th decimal, set it as 0.01 seconds.\n');
      WAV_dt = input('Irregular timing precision (second): ');
		if (isnumeric(WAV_dt) == 0 | isempty(WAV_dt)),
		   flg10 = 0; fprintf(2,'Error: the input is not a number. Please try it again.\n');
	   else flg10 = 1;
      end
   end
end

%Obtain stimulus files
for (iT=1:1:N_tasks),
   flg10 = 0;
	while flg10 == 0,
	   fprintf('Number %i stimulus timing file ',iT);
      Task(iT).StimTime = input('is: ','s');
		fid = fopen (Task(iT).StimTime,'r');
      if (fid == -1),
	      flg10 = 0; fprintf(2,'Error: File %s does not exist. Please try it again. \n', Task(iT).StimTime);
      else flg10 = 1; fclose (fid);
      end
   end
   flg10 = 0;
	while flg10 == 0,
		fprintf('Number %i stimulus label ',iT);
      Task(iT).Label = input('is: ','s');
		if (ischar(Task(iT).Label) == 0 | isempty(Task(iT).Label)),
	      flg10 = 0; fprintf(2,'Error: the input is not a string. Please try it again.\n');
   	else flg10 = 1;
      end
   end
end	

fprintf('\nSometimes covariates are considered in the modeling process. For example, if there were head motion');
fprintf('\nduring the scanning, the variation caused by the motion is better accounted for by adding those');
fprintf('\nparameters from 3dvolreg into the analysis. These regressors will be included as part of the baseline');
fprintf('\nso that the finalthe final full F is interpreted as baseline + pure noise versus baseline + real signal + noise.');

flg10 = 0;
while flg10 == 0,
   base = input('\n\nDo you want to add some covariates? (1 - yes; 0 - no) ');
	if (base ~= 0 & base ~= 1),
	   flg10 = 0; fprintf(2,'Error: the input is not a number. Please try it again.\n');
	   else flg10 = 1;
	end
end

N_base = 0;
if (base == 1),
   flg10 = 0;
   while flg10 == 0,
      N_base = input('\nHow many covariates? ');
		if (isnumeric(N_base) == 0 | isempty(N_base)),
	      flg10 = 0; fprintf(2,'Error: the input is not a number. Please try it again.\n');
   	else flg10 = 1;
      end
   end
	if (N_base > 0),
      for (i = 1:N_base),
	      flg10 = 0;
         while flg10 == 0,
			   fprintf('Number %i covariate ', i);
            stimbase(i).label = input('name is: ', 's');
			   if (ischar(stimbase(i).label) == 0 | isempty(stimbase(i).label)),
	            flg10 = 0; fprintf(2,'Error: the input is not a string. Please try it again.\n');
      	   else flg10 = 1;
            end
         end
			flg10 = 0;
         while flg10 == 0,		
   	      fprintf('Number %i covariate 1D file (or file with column specified) ', i);
            stimbase(i).file = input('expression is: ', 's');
	   		if (ischar(stimbase(i).file) == 0 | isempty(stimbase(i).file)),
	            flg10 = 0; fprintf(2,'Error: the input is not a string. Please try it again.\n');
   	      else flg10 = 1;
            end
         end
		end % for (i = 1:N_base)
	end % if (N_baser > 0)	
end % if (base == 1)	

%Get information for contrast tests from the user

flg10 = 0;
while flg10 == 0,
   contrast = input('\nDo you want to run contrast test? (1 - yes; 0 - no) ');
	if (contrast ~= 0 & contrast ~= 1),
	   flg10 = 0; fprintf(2,'Error: the input is not a number. Please try it again.\n');
	   else flg10 = 1;
	end
end	
	
if (contrast == 1),
   flg10 = 0;
   while flg10 == 0,
      N_contr = input('How many contrast tests do you want to run? ');
	   if (isnumeric(N_contr) == 0 | isempty(N_contr)),
	      flg10 = 0; fprintf(2,'Error: the input is not a number. Please try it again.\n');
   	else flg10 = 1;
      end
   end
		
   contr = [];
	if (N_contr > 0),
      for (i = 1:N_contr),
	      flg10 = 0;
         while flg10 == 0,
            fprintf('Number %i contrast ', i);
            contr(i).label = input('name is: ', 's');
			   if (ischar(contr(i).label) == 0 | isempty(contr(i).label)),
	            flg10 = 0; fprintf(2,'Error: the input is not a string. Please try it again.\n');
      	   else flg10 = 1;
            end
         end
	      flg10 = 0;
         while flg10 == 0,		
   	      fprintf('Number %i contrast ', i);
            contr(i).expr = input('expression is: ', 's');
	   		if (ischar(contr(i).expr) == 0 | isempty(contr(i).expr)),
	            flg10 = 0; fprintf(2,'Error: the input is not a string. Please try it again.\n');
   	      else flg10 = 1;
            end
         end

         if ((contr(i).expr(1) ~= '+') | (contr(i).expr(1) ~= '-')),
            contr(i).expr =['+' contr(i).expr];	   %add + in front of the first label if it does not have one
         end
      end
	end
end
if (contrast == 0), N_contr = 0; end

%Others = input('Any options not listed above (hit Enter if none): ', 's');

%All input from the user is done here.

FoutName_log = sprintf('%s_CommandLog', InBrikPrefix);
foutid_log = fopen(FoutName_log,'w');
if (foutid_log < 0),
   fprintf(2,'Error %s: Failed to write %s to disk.\n',FuncName, FoutName_log);
   while (1); fprintf(2,'Halted: Ctrl+c to exit');pause; end
end

% Start normalization or masking

if (run_norm == 1 | run_mask.do == 1),
  [err, preprc_fn] = Norm(InBrik, N_runs, run, run_norm, run_mask, foutid_log);
   if (err),
      fprintf(2,'Error %s: Failed in normalization.\n', preprc_fn);
      while (1); fprintf(2,'Halted: Ctrl+c to exit'); pause; end
   end
	FoutName_log = sprintf('%s_CommandLog', preprc_fn);
	else FoutName_log = sprintf('%s_CommandLog', InBrikPrefix);
end

if (contrast == 1),   %If the user requests for contrast tests
   for (i = 1:1:N_contr),
      cnt = 0;   								%number of stimulus types
	   El(i).labels = contr(i).label;
   	El(i).sign = [];
	   El(i).expr = '';
   	El(i).strln = [];
	   strlnth = 0;   %string length for each label
      for (j = 1:1:length(contr(i).expr)),
         switch contr(i).expr(j)
	         case '+'
		         cnt = cnt + 1;
				   El(i).strln = [El(i).strln strlnth];   %put the length of previous stimulus label into this array. The 1st one is 0
   				El(i).sign = [El(i).sign 1];           %put the sign of this stimulus label
   				strlnth = 0;
   		   case '-'
	   	      cnt = cnt + 1;
		   		El(i).strln = [El(i).strln, strlnth];  %put the length of previous stimulus label into this array. The 1st one is 0
			   	El(i).sign = [El(i).sign -1];          %put the sign of this stimulus label
   				strlnth = 0;
   	      otherwise if (contr(i).expr(j) ~= ' '),   %spaces are excluded
	   		   strlnth = strlnth+1;
		         El(i).expr = sprintf('%s%s',El(i).expr, contr(i).expr(j));
			   end	
   	   end
      end
   	El(i).strln = [El(i).strln strlnth];  %add the last which was not done in the loop
	   El(i).strln = El(i).strln(2:cnt+1);   %remove the first number which is 0
   	El(i).cnt=cnt;
	   El(i).order = [];
   	first = 1;
	   for (k = 1:1:cnt),
		   next = first+El(i).strln(k)-1;
   	   for (iT=1:1:N_tasks)
	         if strcmp(El(i).expr(first:next), Task(iT).Label),
		   	   El(i).order = [El(i).order iT];
   	         end
	   	end
		   first = next+1;
   	end
	%Make sure if there is any misspelling of the stimulus labels.
      if (size(El(i).order) ~= cnt),
	      fprintf(2,'Error: Misspelled stimulus label %s.\n', contr(i).expr);
         while (1); fprintf(2,'Halted: Ctrl+c to exit');pause; end
      end		
   end
end

%some more setup

fprintf (1,'\n');
if (Run_Type == 2),   %for regression
%create regressors using waver
   if (strcmp (lower(BasisFunc),'spm1') | strcmp(lower(BasisFunc),'spm2')),
      for (iT = 1:1:N_tasks),
         fprintf (1,'Creating regressor(s) for task %d...\n', iT);
	      Task(iT).StimReg = zeros(N_basis, N_TR);
         [err, Task(iT).BasisOpt_struct] = WaverBasisOption (BasisFunc, N_basis, BasisOpt);
         if (err),
            fprintf(2,'Error %s: Failed in WaverBasisOption.\n', FuncName);
            while (1); fprintf(2,'Halted: Ctrl+c to exit');pause; end
         end
      %	p(1) - delay of response (relative to onset)	   6        -gamd  in waver -GAM
      %	p(2) - delay of undershoot (relative to onset)    16
      %	p(3) - dispersion of response			   1
      %	p(4) - dispersion of undershoot			   1   X
      %	p(5) - ratio of response to undershoot		   6   X
      %	p(6) - onset (seconds)				   0     -gamd
      %	p(7) - length of kernel (seconds)		  32

%      p      = [6 16 1 1 6 0 32];   %default parameters
%	   [bf1 p]         = spm_hrf(TR,p);      %get gamma hrf from SPM function spm_hrf
%	   dp     = 1;
%	   p(6)   = p(6) + dp;    %for calculating time derivative
%	   bf2    = (bf1 - spm_hrf(TR,p))/dp;
%		D      = (bf(:,1) - spm_hrf(TR,p))/dp;   %time derivative
%		bf     = [bf D(:)];     %combine the two
%   fid = fopen('%s', Task(iT).Name(iB).str, 'w');
	      vec = load(Task(iT).StimTime);
   		Noext = RemoveExtension(Task(iT).StimTime, '.1D|.txt');
         Task(iT).Name(1).str = sprintf('%s_resp1.1D', Noext);
		   if (strcmp(lower(BasisFunc),'spm2')), Task(iT).Name(2).str = sprintf('%s_resp2.1D', Noext); end
         if (Stim_Type == 1),    %regular timing
	   	   [err, bf] = spm_bf (TR);
		   	if (err),
               fprintf(2,'Error running spm_bf: Failed in running spm_bf. \n');
               while (1); fprintf(2,'Halted: Ctrl+c to exit');pause; end
            end
	         reg1 = conv(bf(:,1), vec);     %use the convolution function in Matlab
		      reg2 = conv(bf(:,2), vec);
      	elseif (Stim_Type == 2), 	 %Irregular timing
   			[err, bf] = spm_bf (WAV_dt);
	   		if (err),
               fprintf(2,'Error running spm_bf: Failed in running spm_bf.\n');
               while (1); fprintf(2,'Halted: Ctrl+c to exit');pause; end
            end
			   irr_vec = zeros(N_TR*TR/WAV_dt, 1);    %Make finer stimulus vector
   			for (i = 1:1:length(vec)),
	   		   irr_vec(round(vec(i)/WAV_dt)+1) = 1;   %Only at those moments (i.e., 1.5 sec) it is set to 1.
		   	end
			   reg_f1 = conv(bf(:,1), irr_vec);    %use the convolution function in Matlab
   		   reg_f2 = conv(bf(:,2), irr_vec);
			
	   		reg1 = reg_f1(1:TR/WAV_dt:N_TR*TR/WAV_dt);   %Upsample those at TR's
		   	reg2 = reg_f2(1:TR/WAV_dt:N_TR*TR/WAV_dt);
			
      	end			
         [err1, Task(iT).Name(1).str] = wryte3(reg1(1:N_TR), Task(iT).Name(1).str);   %Ziad's routine for writing a matrix into a file
		   if (strcmp(lower(BasisFunc),'spm2')), [err2, Task(iT).Name(2).str] = wryte3(reg2(1:N_TR), Task(iT).Name(2).str);
			else err2 = 0; end
	      if (err1 | err2),
            fprintf(2,'Error %s: failed to write into %s\n',Task(iT).Name(1).str);
            while (1); fprintf(2,'Halted: Ctrl+c to exit');pause; end
         end
	   	Task(iT).StimReg(1,:) = load(Task(iT).Name(1).str);  %this is for calculating condition number
		   if (strcmp(lower(BasisFunc),'spm2')), Task(iT).StimReg(2,:) = load(Task(iT).Name(2).str); end
   	end

   else    %for other basis function options
      for (iT = 1:1:N_tasks),
         fprintf (1,'Creating regressor(s) for task %d...\n', iT);
         Task(iT).StimReg = zeros(N_basis, N_TR);
         [err, Task(iT).BasisOpt_struct] = WaverBasisOption (BasisFunc, N_basis, BasisOpt);
         if (err),
            fprintf(2,'Error %s: Failed in WaverBasisOption.\n', FuncName);
            while (1); fprintf(2,'Halted: Ctrl+c to exit');pause; end
         end
         for (iB = 1:1:N_basis),
            if (err),
               fprintf(2,'Error %s:Failed in WaverBasisOption\n', FuncName);
               while (1); fprintf(2,'Halted: Ctrl+c to exit');pause; end
            end
            Noext = RemoveExtension(Task(iT).StimTime, '.1D|.txt');
            Task(iT).Name(iB).str = sprintf('%s_resp%d.1D', Noext, iB);

      %form the command for waver, Probably better to put in a modular function later?
		      if (Stim_Type == 2),   %Irregular timing
				   if (strcmp (lower(BasisFunc),'gam') | strcmp (lower(BasisFunc),'mgh')),
                  comm_Waver = sprintf('waver -dt %g %s -gamb %g -gamc %g -gamd %g -peak 1 -numout %d -tstim `cat %s` > %s',...
                     TR, Task(iT).BasisOpt_struct(iB).opt, Task(iT).BasisOpt_struct.power, Task(iT).BasisOpt_struct.scale, ...
							Task(iT).BasisOpt_struct.delay, N_TR, Task(iT).StimTime, Task(iT).Name(iB).str); 							
		         elseif (strcmp (lower(BasisFunc),'tent')),
					   comm_Waver = sprintf('waver -dt %g %s -peak 1 -numout %d -tstim `cat %s` > %s',...
                     TR, Task(iT).BasisOpt_struct(iB).opt, N_TR, Task(iT).StimTime, Task(iT).Name(iB).str);
					end
				end		
      		if (Stim_Type == 1),    %regular timing
				   if (strcmp (lower(BasisFunc),'gam') | strcmp (lower(BasisFunc),'mgh')),
                  comm_Waver = sprintf('waver -dt %g %s -gamb %g -gamc %g -gamd %g -peak 1 -numout %d -input %s > %s',...
                     TR, Task(iT).BasisOpt_struct(iB).opt, Task(iT).BasisOpt_struct.power, Task(iT).BasisOpt_struct.scale, ...
							Task(iT).BasisOpt_struct.delay, N_TR, Task(iT).StimTime, Task(iT).Name(iB).str);
					elseif (strcmp (lower(BasisFunc),'tent')),
					   comm_Waver = sprintf('waver -dt %g %s -peak 1 -numout %d -input %s > %s',...
                     TR, Task(iT).BasisOpt_struct(iB).opt, N_TR, Task(iT).StimTime, Task(iT).Name(iB).str);
					end		
      		end
      %first log the waver command
		      fprintf(2, 'Running: %s\n', comm_Waver);
            fprintf(foutid_log, '%s\n', comm_Waver);
            [s,w] = system(comm_Waver);
            if (s),
               fprintf(2,'Error %s: calling waver.\n\n\n%s\n',FuncName, w);
               while (1); fprintf(2,'Halted: Ctrl+c to exit');pause; end
            end		
            Task(iT).StimReg(iB,:) = load(Task(iT).Name(iB).str);   %this is for calculating condition number
         end
      end
   end
   comm_Waver = ''; %useless from this point on

%Now display the results
   ShowRegressors(Task);

   Butt = questdlg(sprintf('Check the regressors now ...'),...
            'Continue ?', 'Yes', 'No', 'No');
   switch Butt,
      case 'No',
         while (1); fprintf(2,'Halted: Ctrl+c to exit'); pause; end
      case 'Yes',
         Continue = 1;
   end


%fprintf(1,'Check the regressors now, and hit Enter key if you want to continue ...\n');
%pause;

%Perhaps test the quality of their regressors by calculating the condition number

   [err,CondNum] = CompCondNum (N_tasks, N_basis, Task, N_TR, polort);
   fprintf(2, '\nCondition Number of design matrix is: %g\n. Continue?',CondNum);
   Butt = questdlg(sprintf('Condition number of design matrix is %g', CondNum),...
            'Continue ?', 'Yes', 'No', 'No');
   switch Butt,
      case 'No',
         while (1); fprintf(2,'Halted: Ctrl+c to exit'); pause; end
      case 'Yes',
         Continue = 1;
   end

end

%Now create the 3dDeconvolve thing

Mandatory_opt1 = '';

if (run_norm == 1 | run_mask.do == 1),
   Mandatory_opt1 = sprintf('%s   -input %s%s.BRIK \\\n   -num_stimts %d \\\n   -polort %d \\\n',...
                         Mandatory_opt1, preprc_fn, InBrikView, N_tasks*N_basis+N_base, polort);
else Mandatory_opt1 = sprintf('%s   -input %s%s.BRIK \\\n   -num_stimts %d \\\n   -polort %d \\\n',...
                         Mandatory_opt1, InBrikPrefix, InBrikView, N_tasks*N_basis+N_base, polort);
end
if (N_runs > 1),
   Mandatory_opt1 =sprintf('%s   -concat %s \\\n', Mandatory_opt1, Concat_fn);
end
if (Censor == 1),
   Mandatory_opt1 =sprintf('%s   -censor %s \\\n', Mandatory_opt1, Censor_fn);
end	
									
Mandatory_opt2 = '';
Mandatory_opt3 = '';
Mandatory_opt4 = '';
for (iT=1:1:N_tasks),
   if (Run_Type == 2),   %for regression
      for (iB=1:1:N_basis),
         Mandatory_opt2 = sprintf('%s   -stim_file %d %s \\\n',...
                     Mandatory_opt2, (iT-1)*N_basis+iB, Task(iT).Name(iB).str);
         Mandatory_opt3 = sprintf('%s   -stim_label %d %s \\\n',...
                     Mandatory_opt3, (iT-1)*N_basis+iB, sprintf('%s_b%d', Task(iT).Label, iB));
         Mandatory_opt4 = sprintf('%s   -stim_minlag %d %d -stim_maxlag %d %d \\\n',...
                     Mandatory_opt4, (iT-1)*N_basis+iB, Task(iT).BasisOpt_struct(iB).minlag,...
                     (iT-1)*N_basis+iB, Task(iT).BasisOpt_struct(iB).maxlag);
		end					
	elseif (Run_Type == 1),	%for deconvolution				
	   Mandatory_opt2 = sprintf('%s   -stim_file %d %s \\\n',...
                     Mandatory_opt2, iT, Task(iT).StimTime);
      Mandatory_opt3 = sprintf('%s   -stim_label %d %s \\\n',...
                     Mandatory_opt3, iT, sprintf('%s', Task(iT).Label));
      Mandatory_opt4 = sprintf('%s   -stim_minlag %d %d -stim_maxlag %d %d \\\n',...
                     Mandatory_opt4, iT, Task(iT).BasisOpt_struct(1).minlag, iT, Task(iT).BasisOpt_struct(1).maxlag);	   		
   end
end

if (N_base > 0), % for covariates (with stim_base)
	for (ii = N_base),		
			Mandatory_opt2 = sprintf('%s   -stim_file %d %s \\\n',...
                     Mandatory_opt2, N_tasks*N_basis+ii, stimbase(ii).file);
         Mandatory_opt3 = sprintf('%s   -stim_base %d -stim_label %d %s \\\n',...
                     Mandatory_opt3, N_tasks*N_basis+ii, N_tasks*N_basis+ii, stimbase(ii).label);
	end	
end

%create options for task effect test with genereal linear test
%here we assume the user implicitly wants to test all tasks
%first generate glt files by calling function TaskTest.m

% Only for regression(Run_Type == 2); No need for deconvolution (Run_Type == 1)
% since 3dDeconvolve is supposed to automatically take care of that?

gltopt1 = '';
if (Run_Type == 2),   %for regression
   if (N_basis ~= 1);
      [err,taskfile] = TaskTest(N_runs, polort, N_basis, N_tasks, Task, N_base);
		gltopt1 = sprintf('%s   -num_glt %i\\\n', gltopt1, N_tasks+N_contr);
      for (i = 1:1:N_tasks),		
%         gltopt1 = sprintf('%s   -glt 1 %s   -glt_label %i %s\\\n', gltopt1, taskfile(i).name, i, Task(i).Label);
			gltopt1 = sprintf('%s   -glt %i %s   -glt_label %i %s\\\n', gltopt1, N_basis, taskfile(i).name, i, Task(i).Label);
			% Here a matrix with N_basis rows is set for testing each task
      end
	end	
end

%create glt options for contrast tests

gltopt2 = '';
if (contrast)
   [err,contrfile] = ContrastTest (N_runs, polort, N_basis, Task, N_tasks, N_contr, El, N_base);
	if (err),
      fprintf(2,'Error %s: Failed in ContrastTest.\n', contrfile);
      while (1); fprintf(2,'Halted: Ctrl+c to exit'); pause; end
   end
	if (N_basis == 1), gltopt2 = sprintf('%s   -num_glt %i\\\n', gltopt2, N_contr); end
	% The reason for adding the above line is that with N_bais ~= 1 there has been an option of
	% -num_glt already implemented in gltopt1 above to get task effects
	if (Run_Type == 1),   %only for deconvolution	
      for (i = 1:1:N_contr),
         gltopt2 = sprintf('%s   -glt 1 %s   -glt_label %i %s\\\n', gltopt2, contrfile(i).name, i, contr(i).label);
      end
	elseif (Run_Type == 2),   %only for regression	
	   if (N_basis ~= 1);
	      for (i = 1:1:N_contr),
            gltopt2 = sprintf('%s   -glt 1 %s   -glt_label %i %s\\\n', gltopt2, contrfile(i).name, i+N_tasks, contr(i).label);
			end
		else for (i = 1:1:N_contr),
            gltopt2 = sprintf('%s   -glt 1 %s   -glt_label %i %s\\\n', gltopt2, contrfile(i).name, i, contr(i).label);
			end		
      end
	end
end

fitsprefix = sprintf('%s_fits', InBrikPrefix); statprefix = sprintf('%s_stat', InBrikPrefix);
FitOpts = sprintf('   -fitts %s\\\n   -fout -tout -full_first -bucket %s\\\n',...
                  fitsprefix, statprefix);
comm_3dDecon = sprintf ('3dDeconvolve %s %s %s %s %s %s %s',...
                 Mandatory_opt1, Mandatory_opt2,  Mandatory_opt3, Mandatory_opt4, gltopt1, gltopt2, FitOpts);

%log the command
fprintf(foutid_log, '%s\n', comm_3dDecon);

%let's run the 3dDeconvolve command
%check for output prefix existence:

Skip = 0;
if (PrefixStatus (sprintf('%s%s', statprefix, InBrikView)) ~= 1 | PrefixStatus (sprintf('%s%s', fitsprefix, InBrikView)) ~= 1),
   Butt = questdlg(sprintf('Output dataset %s%s and/or %s%s exists?', statprefix, InBrikView, fitsprefix, InBrikView),...
            'Skip ?', 'Yes', 'No', 'No');
   switch Butt,
      case 'No',
         while (1); fprintf(2,'Halted: Ctrl+c to exit'); pause; end
      case 'Yes',
         Skip = 1;
   end
end

tic,
if (~Skip),
   fprintf (1,'Running: %s\n',comm_3dDecon);
   [s,w] = system(comm_3dDecon);
   if (s),
      fprintf(2,'Error %s: calling 3dDeconvolve.\n\n%s\n\n',FuncName,w);
      while (1); fprintf(2,'Halted: Ctrl+c to exit'); pause; end
   end
end
t=toc;
fprintf(1, '3dDeconvolve done in %f minutes...\n\n', t/60);
	
if (N_basis>1),	
%Load the header of the output stats brik

   fprintf(1,'Recontructing each hemodynamic response function. ');
   [err, InfoStat] = BrikInfo(sprintf('%s%s.HEAD',statprefix, InBrikView));
%   [err, Task(iT).indx, SubLabel] = WhichSubBricks(InfoStat, 'Coef', Task(iT).Label);
   if (err),
      fprintf(2,'Error %s: Failed in BrikInfo.\n', FuncName);
      while (1); fprintf(2,'Halted: Ctrl+c to exit'); pause; end
   end

   %form the 3dcalc command for each task type
   for (iT=1:1:N_tasks),
      [err, Task(iT).indx, SubLabel] = WhichSubBricks(InfoStat, 'Coef', Task(iT).Label);
		if (err),
         fprintf(2,'Error %s: Failed in WhichSubBricks.\n', FuncName);
         while (1); fprintf(2,'Halted: Ctrl+c to exit'); pause; end
      end
      if (isempty(Task(iT).indx)),
         fprintf(2,'Error %s: Failed to find subbricks.\n', FuncName);
         while (1); fprintf(2,'Halted: Ctrl+c to exit'); pause; end
      end
      if (length(Task(iT).indx) ~= N_basis),
         fprintf(2,'Error %s: length(Task(iT).indx) ~= N_basis (%d, %d).\n', FuncName, length(Task(iT).indx), N_basis);
         while (1); fprintf(2,'Halted: Ctrl+c to exit'); pause; end
      end
		
		input_str='';
      expr_str='';
		
      if (strcmp (lower(BasisFunc),'spm1') | strcmp (lower(BasisFunc),'spm2')),
	      for (iB=1:1:N_basis),
            input_str = sprintf('%s   -%c ''%s%s[%d]'' \\\n',...
                           input_str, 96+iB, statprefix, InBrikView, Task(iT).indx(iB) - 1);
         end
	      for (iB=1:1:N_basis),
            input_str = sprintf('%s   -%c %s \\\n',...
                           input_str, 96+N_basis+ iB, Task(iT).Name(i).str);
            expr_str = sprintf('%s%c*%s+',...
                          expr_str, 96+iB, 96+N_basis+ iB);
         end
		
      %remove trailing +
         expr_str = expr_str(1:length(expr_str)-1);

      %create the junk.1D file
      %   fidout = fopen('junk.1D', 'w');
      %   for (i=1:1:BasisOpt.tSpan./Recon_dt+1),
      %     fprintf(fidout,'%g\n', i);
      %   end
      %   fclose(fidout);
         %augment input_str with other options
         input_str = sprintf('%s   -dt %g -datum float\\\n', input_str, Recon_dt);
		
		
      else
	      for (iB=1:1:N_basis),
            input_str = sprintf('%s   -%c ''%s%s[%d]'' \\\n',...
                           input_str, 96+iB, statprefix, InBrikView, Task(iT).indx(iB) - 1);
            expr_str = sprintf('%s%c*%s+',...
                          expr_str, 96+iB,Task(iT).BasisOpt_struct(iB).BareOpt);
         end

         %remove trailing +
         expr_str = expr_str(1:length(expr_str)-1);

         %create the junk.1D file
         fidout = fopen('junk.1D', 'w');
         for (i=1:1:BasisOpt.tSpan./Recon_dt+1),
            fprintf(fidout,'%g\n', i);
         end
         fclose(fidout);
         %augment input_str with junk.1D
         input_str = sprintf('%s   -%c junk.1D \\\n   -dt %g -datum float\\\n', input_str, 96+N_basis+1, Recon_dt);
		end
		

      %prefix of output for this task
      Task(iT).IRFprefix = sprintf('%s_irf', Task(iT).Label);

      %for the 3dcalc command
      comm_3dcalc = sprintf('3dcalc %s   -prefix %s \\\n   -expr ''%s''\n',...
                  input_str, Task(iT).IRFprefix, expr_str);

      %execute the command
      %add command to log file
      fprintf(foutid_log, '%s\n', comm_3dcalc);
      fprintf (1,'%s: Now running 3dcalc command (see %s)...\n', FuncName, FoutName_log);
      [s,w] = system(comm_3dcalc);
      if (s),
         fprintf(2,'Error %s: calling 3dcalc.\n\n%s\n\n',FuncName,w);
         while (1); fprintf(2,'Halted: Ctrl+c to exit'); pause; end
      end

      %for the refit
      comm_refit = sprintf('3drefit -epan %s%s', Task(iT).IRFprefix, InBrikView);
      %add command to log file
      fprintf(foutid_log, '%s\n', comm_refit);
      fprintf (1,'%s: Now running 3drefit command (see %s)...\n', FuncName, FoutName_log);
      [s,w] = system(comm_refit);
      if (s),
         fprintf(2,'Error %s: calling 3dcalc.\n\n%s\n\n',FuncName,w);
         while (1); fprintf(2,'Halted: Ctrl+c to exit!'); pause; end
      end
   end
end
%done, close log file
fclose(foutid_log);
%fprintf(1,'I am so glad that you have reached this point. Congratulations! Open AFNI and check the results ...\n');

fprintf(1,'%s: Done.\n\n', FuncName);
