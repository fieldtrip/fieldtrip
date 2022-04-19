function GroupAna
%GroupAna.m
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
diary('diary');    % save all subsequent command window input and most of the resulting command window output to be appended
                   % to file 'diary'

fprintf('\nWelcome to GroupAna, AFNI Group Analysis Package!');
fprintf('\n-----------------------------------------------------------\n');

fprintf('\nVersion 1.0.1, Nov. 23, 2005');
fprintf('\nAuthor: Gang Chen');
fprintf('\nSSCC/NIMH/ National Institutes of Health, Bethesda MD 20892');
fprintf('\n-----------------------------------------------------------\n');

%Define the function name for easy referencing
FuncName = 'GroupAna.m';

%Debug Flag
DBG = 1;

%initailize return variables
err = 1;

clear all;

% tolerance for numerical 0
tol = 1.0e-4;

fprintf('\nPlease read the following carefully about group analysis setup:\n');
fprintf('\n\t1. If the resolution of your EPI data is near n millimeters, during Talairach conversion use\n');
fprintf('\t  "the command adwarp -dxyz n" to greatly reduce runtime without sacrificing accuracy.\n');
fprintf('\n\t2. We strongly suggest that factor names be labeled with short (2 or 3 capital letters) names\n');
fprintf('\t   so that subbrik labels can be shown on AFNI viewer.\n');
fprintf('\n\t3. With nesting, arrange your design in such a way that the last factor is nested within the 1st factor.\n');
%fprintf('\n\t4. The program can only handle balanced design now. \n');
fprintf('\n\t4. Each input file should include only one subbrik. We suggest files be \n');
fprintf('\t   named by reflecting the hierarchy of the experiment design.\n');
fprintf('\n\t5. Currently all of the following terms are modeled: main effects and applicable interactions in various orders, .\n');
%fprintf('\n\t6. One covariate is currently allowed in the analysis, which should be in the format of one-column text file.\n');
%fprintf('\t   The column length has to be the same as the total number of input files.\n');


% Grouop analysis for Volume or Surface data?
flg = 0;
while flg == 0,
   data_type = input('\nGroups analysis for volume or surface data (0 - volume; 1 - surface 1D output; 2 - surface NIML output;)? ');
	if (data_type ~= 0 & data_type ~= 1 & data_type ~= 2),
	   flg = 0; fprintf(2,'Error: wrong input! Please try it again.\n');
	else flg = 1;
   end
end

% If for surface data, acquire number of nodes
flg = 0;
if (data_type == 1),
   fprintf(1, '\nInput files have to be in 1D format, and program @Purify_1D can be used to extract each regressor');
	fprintf(1, '\ncoefficient column to a 1D file. Type @Purify_1D for usage ...\n');
	Frame_N = input('\nWhich column corresponds to regressor coefficient in the 1D files? (1, 2, 3, ...) ');	
	while flg == 0,
	   node_n = input('\nHow many number of nodes in the surface data? ');
	   if (isnumeric(node_n) == 0 | isempty(node_n)),
	      flg = 0; fprintf(2,'Error: the input is not a number. Please try it again.\n');
	   else flg = 1; end
   end
elseif (data_type == 2),
   fprintf(1,...
'\nEven though the output is in NIML, the input files still have to be in 1D format.');
	Frame_N = input('\nWhich column corresponds to regressor coefficient in the 1D files? (1, 2, 3, ...) ');	
	while flg == 0,
	   node_file = input('Provide a 1D file containing node indices, or enter the total number of nodes if you have full datasets: ','s');
      node_n=str2num(node_file);
      if (isempty(node_n) & ~filexist(node_file)),
         flg = 0; fprintf(2,'Error: the input is not a number, or file %s not found. Please try it again.\n');
	   else flg = 1; end
   end
else Opt.Frames = 1;	 % In the case of volumetric data, it's supposed to have only ONE subbrik!
end

flg = 0;
while flg == 0,
   NF = input('\nHow many factors? ');
	if (isnumeric(NF) == 0 | isempty(NF)),
	   flg = 0; fprintf(2,'Error: the input is not a number. Please try it again.\n');
	else flg = 1;
   end
end

switch NF
   case 1,
	fprintf('\nAvailable design types: \n');
	fprintf('\n\tType 1: one factor with fixed effect.\n');
	
   case 2,
   fprintf('\nAvailabe design types:\n');
   fprintf('\n\tType 1: Factorial (crossed) design AXB - both factors are fixed.\n');
   fprintf('\n\tType 2: Factorial (crossed) design AXB - factor A is fixed while B is random. If B is subject, it is');
	fprintf('\n\t        usually called 1-way design with A varying within subject. Notice: It is inappropriate to');
	fprintf('\n\t        run any constrasts including mean estimates and differences for factor B with this design type.\n');
   fprintf('\n\tType 3: Factorial (crossed) design AXB - both factors are random.\n');	

   case 3,
   fprintf('\nAvailabe design types:\n');
	fprintf('\n\tType 1: Factorial (crossed) design AXBXC - all factors are fixed.\n');
   fprintf('\n\tType 2: Factorial (crossed) design AXBXC - factors A and B are fixed while C is random. If C is subject, it is');
	fprintf('\n\t        usually called 2-way design with A and B varying within subject. Notice: It is inappropriate to');
	fprintf('\n\t        run any constrast tests including mean estimates and differences for factor C with this design type.\n');
	fprintf('\n\tType 3: Mixed design BXC(A)              - A and B are fixed while C (usually subject) is random and nested within A.\n');
	fprintf('\n\tType 4: Mixed design BXC(A)              - Fixed factor C is nested within fixed factor A while B (usually subject) is random.\n');

   case 4,
   fprintf('\nAvailabe design types:\n');
   fprintf('\n\tType 1: Factorial (crossed) design AXBXCXD - all 4 factors are fixed.\n');
   fprintf('\n\tType 2: Factorial (crossed) design AXBXCXD - only factor D is random. If D is subject it is also');
   fprintf('\n\t        called 3-way design with all 3 factors A, B, and C varying within subject.');
   fprintf('\n\tType 3: Mixed design BXCXD(A)- only the nested (4th) factor D (usually subject) is random.');
   fprintf('\n\t        Also called 3-way design with factors B and C varying within subject and factor A between subjects.');
   fprintf('\n\tType 4: Mixed design BXCXD(A)- D is nested within A, but only the 3rd factor (usually subject)');
   fprintf('\n\t        is random.');
   fprintf('\n\tType 5: Mixed design CXD(AXB) - only the nested (4th) factor D (usually subject) is random,');
   fprintf('\n\t        but factor D is nested within both factors A and B. If D is subject it is also called 3-way');
   fprintf('\n\t        design with factor C varying within-subject and factors A and B between-subjects.\n');
%	fprintf('\nNotice: This is NOT an exhaustive list of design types for 4-way ANOVA. Other types might be implemented upon request.\n');
	
	case 5,
   fprintf('\nAvailabe design types:\n');
   fprintf('\n\tType 1: Factorial (crossed) design AXBXCXDXE - all 5 factors are fixed.\n');
   fprintf('\n\tType 2: Factorial (crossed) design AXBXCXDXE - only factor E is random. If E is subject it is also');
   fprintf('\n\t        called 4-way design with all 4 factors A, B, C and D varying within subject.\n');
   fprintf('\n\tType 3: Mixed design BXCXDXE(A) - only the nested (5th) factor E (usually subject) is random.');
   fprintf('\n\t        Also called 4-way design with factors B, C and D varying within subject and factor A between subjects.\n');
	fprintf('\n\tType 4: Mixed design BXCXDXE(A) - the 5th factor E is nested within factor A, but factor D (usually subject)');
   fprintf('\n\t        is random.\n');
	
end

flg = 0;
while flg == 0,
   dsgn = input('\nChoose design type (1, 2, 3, 4, 5): ');
   if (isnumeric(dsgn) == 0 | isempty(dsgn)),
	   flg = 0; fprintf(2,'Error: input is not a number. Please try it again.\n');
	   else flg = 1;
   end
end

%flg = 0;
%while flg == 0,
%   slices = input('How many slices along the Z axis (run 3dinfo on one of the input files to find out)? ');
%	if (isnumeric(slices) == 0 | isempty(slices)),
%	   flg = 0; fprintf(2,'Error: the input is not a number. Please try it again.\n');
%	else flg = 1;
%   end
%end

% Balanced design?
flg = 0;
while flg == 0,

	unbalanced.yes = 1 - input('\nIs the design balanced? (1 - Yes; 0 - No) ');
   if (unbalanced.yes ~= 0 & unbalanced.yes ~= 1),
 	   flg = 0; fprintf(2,'Error: inapproriate input. Please try it again.\n');
 	else flg = 1;
 	end
 	
	if (unbalanced.yes == 1),	
      fprintf('\nThe following two kinds of unbalanced designs are currently allowed:\n')
		fprintf('\n(1) All factors are fixed - ');
		fprintf('\n\t 1-way ANOVA; 2-way ANOVA type 1: AXB; and 3-way ANOVA type 1: AXBXC');
	   fprintf('\n\n(2) When a random factor (subject) is nested within another factor A,');
   	fprintf('\n\t each level of factor A contains a unique and unequal number of subjects - ');
	   fprintf('\n\t 3-way ANOVA type 3: BXC(A); 4-way ANOVA types 3: BXCXD(A); ');
		fprintf('\n\t and 4-way ANOVA type 5: CXD(AXB).')
	   if (input('\n\nDoes your unbalanced design belong to either of the above types? (1 - Yes; 0 - No) ') == 0);
		   while (1); fprintf(2,'Halted: Ctrl+c to exit'); pause; end
		end	
	end

%  fprintf('\nTwo kinds of unbalanced design are considered:\n');
% 	fprintf('\nType 1: When subject is a random factor and nested within another factor (3-way ANOVA type 3, 4-way ANOVA types 3 and 5),');
% 	fprintf('\n\t there can be different number of subjects for each level of factor A;\n');
% 	fprintf('\nType 2: Not type 1, but there are different sample sizes for some cells (factor level combinations).');
% 	fprintf('\n\t For example, add somethint here later...\n');
%	 	
% 	unbalanced.type = input('\nChoose design type (1, 2): ');
%    if (unbalanced.type ~= 1 & unbalanced.type ~= 2),
% 	   flg = 0; fprintf(2,'Error: inapproriate input. Please try it again.\n');
% 	   else flg = 1;
%    end
end

for (ii = 1:1:5)
   FL(ii).N_level = 0;  % initialization for ContrVec.m so that their values are available for lower-way ANOVA
end	


%Get the number of levels for each factor and also the total number of input files

flg = 0;
ntot = 1;

if (unbalanced.yes == 1),

   if ((NF == 1 | NF == 2 | NF == 3 | NF == 4) & dsgn == 1),  % Basically for 1,2,3,4-way ANCOVA: This loop is the same as balanced. Maybe also for 4-way ANCOVA?
      for (i=1:1:NF),
	      fprintf('\nLabel for No. %i ', i);
		   FL(i).expr = input('factor: ', 's');     % Factor Label i
         fprintf(2,'How many levels does factor %c (%s) ', 64+i, FL(i).expr);
	      FL(i).N_level = input('have? ');
	      if (isnumeric(FL(i).N_level) == 0 | isempty(FL(i).N_level)),
	         flg = 0; fprintf(2,'Error: the input is not a number. Please try again.\n');
	      else flg = 1; end
	      for (j=1:1:FL(i).N_level),
	         fprintf('Label for No. %i level of factor %c (%s)', j, 64+i, FL(i).expr);
		      FL(i).level(j).expr = input(' is: ', 's');
	      end
         sz(i) = FL(i).N_level;    % number of levels of factor i
         ntot = ntot * FL(i).N_level;  %total number of combinations
      end
   end

   if (((NF == 3 | NF == 4) & dsgn == 3)),	
%   if (unbalanced.type == 1),
      for (i=1:1:(NF-1)),
         fprintf('\nLabel for No. %i ', i);
         FL(i).expr = input('factor: ', 's');     % Factor Label i	
	      fprintf(2,'How many levels does factor %c (%s) ', 64+i, FL(i).expr);
 	      FL(i).N_level = input('have? ');
 	      if (isnumeric(FL(i).N_level) == 0 | isempty(FL(i).N_level)),
 	         flg = 0; fprintf(2,'Error: the input is not a number. Please try again.\n');
 	      else flg = 1;
         end
 	      for (j=1:1:FL(i).N_level),
 	         fprintf('Label for No. %i level of factor %c (%s)', j, 64+i, FL(i).expr);
 		      FL(i).level(j).expr = input(' is: ', 's');
 	      end
         sz(i) = FL(i).N_level;    % number of levels of factor i
         ntot = ntot * FL(i).N_level;  %total number of combinations
 	   end  % for (i=1:1:(NF-1))
	
	   FL(NF).N_level = 0;		
	   fprintf(2, '\nLabel for No. %i ', NF);
	   FL(NF).expr = input('factor: ', 's');     % Label this unbalanced factor
%		fprintf(2, '\nNote: Since this is a nested design the label for levels (usuall subject names) of factor No. %i (%c - %s)', NF, 64+NF, FL(NF).expr);
%		fprintf(2, '\nhas to be DIFFERENT for each level of factor %c (%s)!!!\n\n', 64+1, FL(1).expr);
 	
 	   flag = 0;
		while flag == 0,
		combine = [];
		for (i = 1:1:FL(1).N_level),
 	      fprintf(2,'How many levels does factor %c (%s) corresponding to level %i (%s) of factor %c (%s) ', ...
 		      64+NF, FL(NF).expr, i, FL(1).level(i).expr, 64+1, FL(1).expr);
 		   FL(NF).UL(i).N_level = input('have? ');   % unbalanced levels for this factor
 		   for (j=1:1:FL(NF).UL(i).N_level),
 		      fprintf('Label for No. %i level of factor %c (%s) in group %i of factor %c (%s)', j, 64+NF, FL(NF).expr, i, 64+1, FL(1).expr);
 		      FL(NF).UL(i).n(j).expr = input(' is: ', 's');
				combine = [combine {FL(NF).UL(i).n(j).expr}];  % Concatenate them to make a cell array
 		   end	
% 		   FL(NF).N_level = FL(NF).N_level + FL(NF).UL(i).N_level;
% 		   FL(NF).N_level = max([FL(NF).UL(:).N_level]);   % This is for positioning those contrast columns in the design matrix in PreProc.m
 	   end 		
		
		if (length(unique(combine)) == max([FL(NF).UL(:).N_level])),  % if same labels are used across groups
		   FL(NF).N_level = max([FL(NF).UL(:).N_level]); flag = 1;    % design matrix is built based on the longest group length
			ntot = ntot * sum([FL(NF).UL(:).N_level]) / FL(1).N_level;
		elseif (length(unique(combine)) == sum([FL(NF).UL(:).N_level])),  % if different labels are used across groups
		   FL(NF).N_level = sum([FL(NF).UL(:).N_level]); flag = 1;    % design matrix is built based on the total length of the groups
			ntot = ntot * FL(NF).N_level / FL(1).N_level;
		else 	
		   flag = 0;
			fprintf(2, '\nError: There is some overlap among the labels of factor %c (%s) across the groups of factor %c (%s)!\n', ...
			   64+NF, FL(NF).expr, 64+1, FL(1).expr);
	   end			
		end % while flag == 0
		
   end   % if (((NF == 3 | NF == 4) & dsgn == 3))
	
	if (NF == 4 & dsgn == 5),   % CXD(AXB)
      for (i=1:1:(NF-1)),   % Get information for factors A, B, and C
         fprintf('\nLabel for No. %i ', i);
         FL(i).expr = input('factor: ', 's');     % Factor Label i	
	      fprintf(2,'How many levels does factor %c (%s) ', 64+i, FL(i).expr);
 	      FL(i).N_level = input('have? ');
 	      if (isnumeric(FL(i).N_level) == 0 | isempty(FL(i).N_level)),
 	         flg = 0; fprintf(2,'Error: the input is not a number. Please try again.\n');
 	      else flg = 1;
         end
 	      for (j=1:1:FL(i).N_level),
 	         fprintf('Label for No. %i level of factor %c (%s)', j, 64+i, FL(i).expr);
 		      FL(i).level(j).expr = input(' is: ', 's');
 	      end
         sz(i) = FL(i).N_level;    % number of levels of factor i
         ntot = ntot * FL(i).N_level;  %total number of combinations
 	   end  % for (i=1:1:(NF-1))
	
	   FL(NF).N_level = 0;		% Info for factor D, which is nested with both A and B
	   fprintf(2, '\nLabel for No. %i ', NF);
	   FL(NF).expr = input('factor: ', 's');     % Label this unbalanced factor
 	
 	   flag = 0;
		while flag == 0,
		combine = [];
		for (ii = 1:1:FL(1).N_level),  % factor A
		for (jj = 1:1:FL(2).N_level), % factor B
 	      fprintf(2,'How many levels does factor %c (%s) corresponding to level %i (%s) of factor %c (%s) and level %i (%s) of factor %c (%s) ', ...
 		      64+NF, FL(NF).expr, ii, FL(1).level(ii).expr, 64+1, FL(1).expr, jj, FL(2).level(jj).expr, 64+2, FL(2).expr);
 		   FL(NF).UL(ii, jj).N_level = input('have? ');   % unbalanced levels for this factor
 		   for (kk=1:1:FL(NF).UL(ii, jj).N_level),
 		      fprintf('Label for No. %i level of factor %c (%s) in group %i of factor %c (%s) and group %i of factor %c (%s)', ...
				   kk, 64+NF, FL(NF).expr, ii, 64+1, FL(1).expr, jj, 64+2, FL(2).expr);
 		      FL(NF).UL(ii, jj).n(kk).expr = input(' is: ', 's');
				combine = [combine {FL(NF).UL(ii, jj).n(kk).expr}];  % Concatenate them to make a cell array
 		   end	
 	   end  % jj
		end  % ii		
		
		if (length(unique(combine)) == max([FL(NF).UL.N_level])),   % if same labels are used across groups
		   FL(NF).N_level = max([FL(NF).UL.N_level]); flag = 1;     % design matrix is built based on the longest group length
			ntot = ntot * sum([FL(NF).UL.N_level]) / (FL(1).N_level * FL(2).N_level);  % total number of input files
%		elseif (length(unique(combine)) == sum([FL(NF).UL.N_level])),  % if different labels are used across groups
%		   FL(NF).N_level = sum([FL(NF).UL.N_level]); flag = 1;    % design matrix is built based on the total length of the groups
%			ntot = ntot * FL(NF).N_level / (FL(1).N_level * FL(2).N_level);
%		else 	
%		   flag = 0;
%			fprintf(2, '\nError: There is some overlap among the labels of factor %c (%s) across the groups of factor %c (%s) and factor %c (%s)!\n', ...
%			   64+NF, FL(NF).expr, 64+1, FL(1).expr, 64+2, FL(2).expr);
		else 	
		   flag = 0;
			fprintf(2, '\nError: Currently the labels for different groups have to be the same for this design.'); 						
	   end			
		end % while flag == 0
		
   end   % (NF == 4 & dsgn == 5)
	

else  % Balanced designs

   for (i=1:1:NF),
      fprintf('\nLabel for No. %i ', i);
      FL(i).expr = input('factor: ', 's');     % Factor Label i
      fprintf(2,'How many levels does factor %c (%s) ', 64+i, FL(i).expr);
	   FL(i).N_level = input('have? ');
  	   if (isnumeric(FL(i).N_level) == 0 | isempty(FL(i).N_level)),
	      flg = 0; fprintf(2,'Error: the input is not a number. Please try again.\n');
	   else flg = 1;
      end
	   for (j=1:1:FL(i).N_level),
	      fprintf('Label for No. %i level of factor %c (%s)', j, 64+i, FL(i).expr);
		   FL(i).level(j).expr = input(' is: ', 's');
	   end
      sz(i) = FL(i).N_level;    % number of levels of factor i
      ntot = ntot * FL(i).N_level;  %total number of combinations
   end

end  % if (unbalanced)

if ~((NF == 1 | NF == 2 | NF == 3 | NF == 4) & dsgn == 1 & unbalanced.yes == 1),
fprintf(2,'\nEnter the sample size (number of observations) per combination');
FL(NF+1).N_level =  input(': ');  % Should it be changed to a different name instead of FL???
ntot = ntot * FL(NF+1).N_level;    % total number of factor combinations including repeats
sz(NF+1) = FL(NF+1).N_level;
end

% Running ANCOVA?
flg = 0;
%cov.label = [];
%while flg == 0,
%   cov.do = input('\nAny covariate (concomitant variable)? (1 - Yes; 0 - No) ');
%   if (cov.do ~= 0 & cov.do ~= 1),
%	   flg = 0; fprintf(2,'Error: inapproriate input. Please try it again.\n');
%	else flg = 1;
%	end
%end

cov.do = 0;

%Only allow one covariate now!
if (cov.do),

   FL(NF+1).N_level = 1;  % This is for PreProc.m to get the correct 'shift' position

   if (NF == 2 & dsgn == 3),  % 2-way ANOVA of A(R)XB(R) possible in FMRI?
	   fprintf('\nThe 2-way ANCOVA for this design is currently NOT available.\n');
		while (1); fprintf(2,'Halted: Ctrl+c to exit'); pause; end
   end

   flg = 0;
	cov.label = input('Label for the covariate is: ', 's');
	fprintf('The 1D file for the covariate has to be in the format of one-column text.\n');
	fprintf('And the column should have exactly the same number and order of those input files.\n');	
   while flg == 0,
   cov.FN = input('\nConvariate file name: ', 's');
   fid0 = fopen(cov.FN,'r');
   if (fid0 == -1),
      flg = 0; fprintf(2,'Error: File %s does not exist. Please try it again. \n', cov.FN);
   else flg = 1;
	   [cov.vec, count] = fscanf(fid0, '%f');
	   if (count ~= ntot & ~((NF ==1 | NF == 2 | NF == 3) & dsgn == 1)), % Check length of the 1D file
	      fprintf(2, '\nError: The column length of the covariate has to equal to the total number of input files!\n');
			fprintf(2,'Halted: Ctrl+c to exit'); pause;
	   end
		cov.vec = cov.vec - mean(cov.vec);
	end
	end   % while
end

%flg = 0;
%while flg == 0,
%   fprintf('\nIs input in the format of one-subbrik-one-file or one-file-per-subject?\n');
%	file_format = input('(1 - Single brik; 0 - Multiple subbriks) ');
%   if (file_format ~= 0 & file_format ~= 1),
%	   flg = 0; fprintf(2,'Error: inapproriate input. Please try it again.\n');
%	else flg = 1;
%	end
%end

fprintf('\nAll input files are supposed to contain only one subbrik.\n');
file_format = 1;  % Ignore that file_format = 0 now since it leads to too much trouble

if (file_format == 0),   % unused now!
   flg = 0;
   file_num = input('Number of input files: ');
	file_SB = input('Number of subbriks for analysis in each file: ');
	if (isnumeric(file_num) == 0 | isempty(file_num) | isnumeric(file_SB) == 0 | isempty(file_SB)),
	   flg = 0; fprintf(2,'Error: the input is not a number. Please try it again.\n');
	else flg = 1; end	
   if (file_num*file_SB ~= ntot),
	   fprintf(2, '\nError: The number of files and subbriks are not consistent!\n');
 		fprintf(2,'Halted: Ctrl+c to exit'); pause;
	end	
end   % if (file_format == 0)

if (file_format == 1 & ~((NF == 1 | NF == 2 | NF == 3 | NF == 4) & dsgn == 1)),
   fprintf(2,'\nThere should be totally %i input files. \n', ntot);	
	corr = input('Correct? (1 - Yes; 0 - No) ');
	if (corr == 0),
      fprintf(2,'Error: Check the inconsistency. \n');
		while (1); fprintf(2,'Halted: Ctrl+c to exit'); pause; end
   end	
end

%generate the subcripts and obtain all input files
GP = cell(NF, ntot);    %creat a cell array to reflect the structure of all the combinations

if (file_format == 1),

   if (unbalanced.yes == 1),    % Try 4-way design 3 first
	
%	if ((NF == 1 | NF == 2 | NF == 3 | NF == 4) & dsgn == 1),  % Meant for 1,2,3,4-way ANCOVA with unequal sample size
	if (dsgn == 1),
	   FI = 0; % File index
		flg = 0;
		ntot = input('Total number of input files: ');
		if (isnumeric(ntot) == 0 | isempty(ntot)),
 	      flg = 0; fprintf(2,'Error: the input is not a number. Please try again.\n');
 	   else flg = 1;
      end
	
	switch NF
	case 1,
		for (ii = 1:1:FL(1).N_level),		
				   fprintf (2,'\nFor factor %c (%s) at level %i (%s),', 64+1, FL(1).expr, ii, FL(1).level(ii).expr);
					sz = input('\nsample size is: ');
					fprintf (2,'\nProvide those %i input files:\n', sz);					
					for (rr = 1:1:sz),
						FI = FI + 1;
						GP(1, FI) = {FL(1).level(ii).expr};
						
						flg = 0;	
						while flg == 0,
						   fprintf (2,'No. %i file ', rr);
							file(FI).nm = input('is: ', 's');
							fid = fopen (file(FI).nm,'r');	
							if (fid == -1), flg = 0; fprintf(2,'Error: File %s does not exist. Please try it again. \n', file(FI).nm);
							else flg = 1; fclose (fid);	end
							if isempty(strfind(file(FI).nm, 'tlrc')) == 0
					         format = 'tlrc';
					      elseif isempty(strfind(file(FI).nm, 'orig')) == 0
					         format = 'orig';
					      else 	
					         if isempty(strfind(file(FI).nm, '1D')) == 0		% if 1D file (surface data)
								   format = '1D';
								else 					
								   while (1); fprintf(2,'Error: format of file %s is incorrect!\n', file(FI).nm);
						         fprintf(2,'Halted: Ctrl+c to exit'); pause; end
								end
					      end
						end
						
					end % for (rr = 1:1:sz)
		end % for (ii = 1:1:FL(1).N_level)		

	case 2,
		for (ii = 1:1:FL(1).N_level),		
   	   for (jj = 1:1:FL(2).N_level),		
				   fprintf (2,'\nFor factor %c (%s) at level %i (%s),', 64+1, FL(1).expr, ii, FL(1).level(ii).expr);
					fprintf (2,'\n    factor %c (%s) at level %i (%s),', 64+2, FL(2).expr, jj, FL(2).level(jj).expr);
					sz = input('\nsample size is: ');
					fprintf (2,'\nProvide those %i input files:\n', sz);					
					for (rr = 1:1:sz),
						FI = FI + 1;
						GP(1, FI) = {FL(1).level(ii).expr};
						GP(2, FI) = {FL(2).level(jj).expr};
						
						flg = 0;	
						while flg == 0,
						   fprintf (2,'No. %i file ', rr);
							file(FI).nm = input('is: ', 's');
							fid = fopen (file(FI).nm,'r');	
							if (fid == -1), flg = 0; fprintf(2,'Error: File %s does not exist. Please try it again. \n', file(FI).nm);
							else flg = 1; fclose (fid);	end
							if isempty(strfind(file(FI).nm, 'tlrc')) == 0
					         format = 'tlrc';
					      elseif isempty(strfind(file(FI).nm, 'orig')) == 0
					         format = 'orig';
					      else 	
					         if isempty(strfind(file(FI).nm, '1D')) == 0		% if 1D file (surface data)
								   format = '1D';
								else	
								   while (1); fprintf(2,'Error: format of file %s is incorrect!\n', file(FI).nm);
						         fprintf(2,'Halted: Ctrl+c to exit'); pause; end
								end
					      end
						end
						
					end % for (rr = 1:1:sz)
			end % for (jj = 1:1:FL(2).N_level)
		end % for (ii = 1:1:FL(1).N_level)
	
	case 3,
		for (ii = 1:1:FL(1).N_level),		
   	   for (jj = 1:1:FL(2).N_level),		
   	      for (kk = 1:1:FL(3).N_level),
				   fprintf (2,'\nFor factor %c (%s) at level %i (%s),', 64+1, FL(1).expr, ii, FL(1).level(ii).expr);
					fprintf (2,'\n    factor %c (%s) at level %i (%s),', 64+2, FL(2).expr, jj, FL(2).level(jj).expr);
					fprintf (2,'\n    factor %c (%s) at level %i (%s),', 64+3, FL(3).expr, kk, FL(3).level(kk).expr);
					sz = input('\nsample size is: ');
					fprintf (2,'\nProvide those %i input files:\n', sz);					
					for (rr = 1:1:sz),
						FI = FI + 1;
						GP(1, FI) = {FL(1).level(ii).expr};
						GP(2, FI) = {FL(2).level(jj).expr};
						GP(3, FI) = {FL(3).level(kk).expr};
						
						flg = 0;	
						while flg == 0,
						   fprintf (2,'No. %i file ', rr);
							file(FI).nm = input('is: ', 's');
							fid = fopen (file(FI).nm,'r');	
							if (fid == -1), flg = 0; fprintf(2,'Error: File %s does not exist. Please try it again. \n', file(FI).nm);
							else flg = 1; fclose (fid);	end
							if isempty(strfind(file(FI).nm, 'tlrc')) == 0
					         format = 'tlrc';
					      elseif isempty(strfind(file(FI).nm, 'orig')) == 0
					         format = 'orig';
					      else 	
					         if isempty(strfind(file(FI).nm, '1D')) == 0		% if 1D file (surface data)
								   format = '1D';
								else 					
								   while (1); fprintf(2,'Error: format of file %s is incorrect!\n', file(FI).nm);
						         fprintf(2,'Halted: Ctrl+c to exit'); pause; end
								end
					      end
						end
						
					end % for (rr = 1:1:sz)
				end % for (kk = 1:1:FL(3).N_level)
			end % for (jj = 1:1:FL(2).N_level)
		end % for (ii = 1:1:FL(1).N_level)

   case 4,		
		for (ii = 1:1:FL(1).N_level),		
   	   for (jj = 1:1:FL(2).N_level),		
   	      for (kk = 1:1:FL(3).N_level),
				for (ll = 1:1:FL(4).N_level),
				   fprintf (2,'\nFor factor %c (%s) at level %i (%s),', 64+1, FL(1).expr, ii, FL(1).level(ii).expr);
					fprintf (2,'\n    factor %c (%s) at level %i (%s),', 64+2, FL(2).expr, jj, FL(2).level(jj).expr);
					fprintf (2,'\n    factor %c (%s) at level %i (%s),', 64+3, FL(3).expr, kk, FL(3).level(kk).expr);
					fprintf (2,'\n    factor %c (%s) at level %i (%s),', 64+3, FL(3).expr, kk, FL(4).level(ll).expr);
					sz = input('\nsample size is: ');
					fprintf (2,'\nProvide those %i input files:\n', sz);					
					for (rr = 1:1:sz),
						FI = FI + 1;
						GP(1, FI) = {FL(1).level(ii).expr};
						GP(2, FI) = {FL(2).level(jj).expr};
						GP(3, FI) = {FL(3).level(kk).expr};
						GP(4, FI) = {FL(4).level(ll).expr};
						flg = 0;	
						while flg == 0,
						   fprintf (2,'No. %i file ', rr);
							file(FI).nm = input('is: ', 's');
							fid = fopen (file(FI).nm,'r');	
							if (fid == -1), flg = 0; fprintf(2,'Error: File %s does not exist. Please try it again. \n', file(FI).nm);
							else flg = 1; fclose (fid);	end
							if isempty(strfind(file(FI).nm, 'tlrc')) == 0
					         format = 'tlrc';
					      elseif isempty(strfind(file(FI).nm, 'orig')) == 0
					         format = 'orig';
					      else 	
					         if isempty(strfind(file(FI).nm, '1D')) == 0		% if 1D file (surface data)		   					
									format = '1D';
							   else
								   while (1); fprintf(2,'Error: format of file %s is incorrect!\n', file(FI).nm);
						         fprintf(2,'Halted: Ctrl+c to exit'); pause; end
								end
					      end
						end
						
					end % for (rr = 1:1:sz)
					end % for (ll = 1:1:FL(4).N_level)
				end % for (kk = 1:1:FL(3).N_level)
			end % for (jj = 1:1:FL(2).N_level)
		end % for (ii = 1:1:FL(1).N_level)
		
   end % switch NF
	if (ntot == FI), fprintf (2,'\n%i input files have been read in. \n', FI);
	else fprintf(2,'Error: Total number of files do not match up. \n');
	while (1); fprintf(2,'Halted: Ctrl+c to exit'); pause; end
	end	
	end % if (dsgn == 1)					
	
	if ((NF == 3 & dsgn == 3)),	
%      acc = 0;
		FI = 0; % File index
%    	if (unbalanced.type == 1),
   		for (i = 1:1:FL(1).N_level),
   	   for (j = 1:1:FL(2).N_level),
   	   for (k = 1:1:FL(3).UL(i).N_level), 	
 	      			
%   			FI = acc + (j-1)*FL(3).UL(i).N_level + k;   % file index
 	 		   % Create a matrix for group indices
% 		   	GP(1, FI) = {FL(1).level(i).expr};
%   		   GP(2, FI) = {FL(2).level(j).expr};
%   		   GP(3, FI) = {FL(3).UL(i).n(k).expr};   		
				
				for (r = 1:1:FL(NF+1).N_level),  % if there is any repeated observations
				FI = FI + 1;
				GP(1, FI) = {FL(1).level(i).expr};
   		   GP(2, FI) = {FL(2).level(j).expr};
   		   GP(3, FI) = {FL(3).UL(i).n(k).expr};
				
				flg = 0;		
				while flg == 0,
   	         fprintf (2,'\n(%i) factor combination:\n', FI);
               for (m=1:1:NF),
	               fprintf('\tfactor %c (%s) at level %s \n', 64+m, FL(m).expr, char(GP(m, FI)));   			
               end
	            fprintf('\tat repeat %i \n', r);
   	         file(FI).nm = input('is: ', 's');	
	            fid = fopen (file(FI).nm,'r');	
               if (fid == -1),
   	            flg = 0; fprintf(2,'Error: File %s does not exist. Please try it again. \n', file(FI).nm);
               else flg = 1; fclose (fid);	end
					if isempty(strfind(file(FI).nm, 'tlrc')) == 0
					   format = 'tlrc';
					elseif isempty(strfind(file(FI).nm, 'orig')) == 0
					   format = 'orig';
					else 	
					   if isempty(strfind(file(FI).nm, '1D')) == 0		% if 1D file (surface data)
							format = '1D';
						else 					
						   while (1); fprintf(2,'Error: format of file %s is incorrect!\n', file(FI).nm);
						   fprintf(2,'Halted: Ctrl+c to exit'); pause; end
						end
					end					
            end
				end % for (r = 1:1:FL(NF+1).N_level),	
 			
   	 	end  % k
 	  	   end  % j
%   		acc = FI;
   		end  % i
 	 		  		
%      end % if (unbalanced.type == 1),
   end % if ((NF == 3 & dsgn == 3))	
	
   if ((NF == 4 & dsgn == 3)),	
%      acc = 0;
%    	if (unbalanced.type == 1),
   		FI = 0; % File index
			for (i = 1:1:FL(1).N_level),
   	   for (j = 1:1:FL(2).N_level),
   	   for (k = 1:1:FL(3).N_level),
   	   for (l = 1:1:FL(4).UL(i).N_level),
 	      			
%   			FI = acc + (j-1)*FL(3).N_level*FL(4).UL(i).N_level+(k-1)*FL(4).UL(i).N_level + l;   % file index
 	 		   % Create a matrix for group indices
% 		   	GP(1, FI) = {FL(1).level(i).expr};
%   		   GP(2, FI) = {FL(2).level(j).expr};
%   		   GP(3, FI) = {FL(3).level(k).expr};
%   		   GP(4, FI) = {FL(4).UL(i).n(l).expr};
				
				for (r = 1:1:FL(NF+1).N_level),				
				FI = FI +1;
				GP(1, FI) = {FL(1).level(i).expr};
   		   GP(2, FI) = {FL(2).level(j).expr};
   		   GP(3, FI) = {FL(3).level(k).expr};
   		   GP(4, FI) = {FL(4).UL(i).n(l).expr};
				
				flg = 0;		
				while flg == 0,
   	         fprintf (2,'\n(%i) factor combination:\n', FI);
               for (m=1:1:NF),
	               fprintf('\tfactor %c (%s) at level %s \n', 64+m, FL(m).expr, char(GP(m, FI)));   			
               end
	            fprintf('\tat repeat %i \n', r);
   	         file(FI).nm = input('is: ', 's');	
	            fid = fopen (file(FI).nm,'r');	
               if (fid == -1),
   	            flg = 0; fprintf(2,'Error: File %s does not exist. Please try it again. \n', file(FI).nm);
               else flg = 1; fclose (fid);	end
					
					if isempty(strfind(file(FI).nm, 'tlrc')) == 0
					   format = 'tlrc';
					elseif isempty(strfind(file(FI).nm, 'orig')) == 0
					   format = 'orig';
					else 	
					         if isempty(strfind(file(FI).nm, '1D')) == 0		% if 1D file (surface data)
								   format = '1D';
								else 					
								   while (1); fprintf(2,'Error: format of file %s is incorrect!\n', file(FI).nm);
						         fprintf(2,'Halted: Ctrl+c to exit'); pause; end
								end
					end					
            end
				end % for (r = 1:1:FL(NF+1).N_level),	
 			
   	 	end
 		   end
 	  	   end
%   		acc = FI;
   		end
 	 		  		
%     end % if (unbalanced.type == 1),
     end % if ((NF == 4 & dsgn == 3))
	
   if ((NF == 4 & dsgn == 5)),	% CXD(AXB)
   		FI = 0; % File index
			for (i = 1:1:FL(1).N_level),
   	   for (j = 1:1:FL(2).N_level),
   	   for (k = 1:1:FL(3).N_level),
   	   for (l = 1:1:FL(4).UL(i, j).N_level),
 	      			
				for (r = 1:1:FL(NF+1).N_level),				
				FI = FI +1;
				GP(1, FI) = {FL(1).level(i).expr};
   		   GP(2, FI) = {FL(2).level(j).expr};
   		   GP(3, FI) = {FL(3).level(k).expr};
   		   GP(4, FI) = {FL(4).UL(i, j).n(l).expr};
				
				flg = 0;		
				while flg == 0,
   	         fprintf (2,'\n(%i) factor combination:\n', FI);
               for (m=1:1:NF),
	               fprintf('\tfactor %c (%s) at level %s \n', 64+m, FL(m).expr, char(GP(m, FI)));   			
               end
	            fprintf('\tat repeat %i \n', r);
   	         file(FI).nm = input('is: ', 's');	
	            fid = fopen (file(FI).nm,'r');	
               if (fid == -1),
   	            flg = 0; fprintf(2,'Error: File %s does not exist. Please try it again. \n', file(FI).nm);
               else flg = 1; fclose (fid);	end
					
					if isempty(strfind(file(FI).nm, 'tlrc')) == 0
					   format = 'tlrc';
					elseif isempty(strfind(file(FI).nm, 'orig')) == 0
					   format = 'orig';
					else 	
					         if isempty(strfind(file(FI).nm, '1D')) == 0		% if 1D file (surface data)
								   format = '1D';
								else 					
								   while (1); fprintf(2,'Error: format of file %s is incorrect!\n', file(FI).nm);
						         fprintf(2,'Halted: Ctrl+c to exit'); pause; end
								end
					end					
            end
				end % for (r = 1:1:FL(NF+1).N_level),	
 			
   	 	end
 		   end
 	  	   end
   		end
 	 		  		
     end % if ((NF == 4 & dsgn == 5))	  	
	
     %end  % (unbalanced.yes == 1),

   else  % Balanced designs	
	
	fprintf('\nPlease provide input files.');
   for (i=1:1:ntot),
   %Ziad Saad modified Matlab function ind2sub for the purpose here
	%Converting the index into multiple subscripts
   %I want vary the levels starting the last factor in stead of the first.
	   [err, file(i).v] = gind2sub (fliplr(sz), i);   %flip before subscripting. Doing flip because
   %I want vary the levels starting the last factor in stead of the first.
      scpt = fliplr(file(i).v);  %flip back to restore the original order	
	
	   flg = 0;
   	while flg == 0,
         fprintf (2,'\n(%i) factor combination:\n', i);
         for (k=1:1:NF),
	         fprintf('\tfactor %c (%s) at level %i (%s) \n', 64+k, FL(k).expr, scpt(k), FL(k).level(scpt(k)).expr);
   			GP(k, i) = {FL(k).level(scpt(k)).expr};
         end
	      fprintf('\tat repeat %i \n', scpt(NF+1));
   	   file(i).nm = input('is: ', 's');	
	      fid = fopen (file(i).nm,'r');	
         if (fid == -1),
   	      flg = 0; fprintf(2,'Error: File %s does not exist. Please try it again. \n', file(i).nm);
         else flg = 1; fclose (fid);	end
			
			if isempty(strfind(file(i).nm, 'tlrc')) == 0
			   format = 'tlrc';
			elseif isempty(strfind(file(i).nm, 'orig')) == 0
			   format = 'orig';
			else 	
					         if isempty(strfind(file(i).nm, '1D')) == 0		% if 1D file (surface data)
								   format = '1D';
								else 					
								   while (1); fprintf(2,'Error: format of file %s is incorrect!\n', file(i).nm);
						         fprintf(2,'Halted: Ctrl+c to exit'); pause; end
								end
			end					
      end	
   end
	end % if (unbalanced.yes == 1),
end % if (file_format == 1),

if (file_format == 0),   % unused now!

   if (unbalanced.yes == 1),    % Try 4-way design 3 first
   if ((NF == 4 & dsgn == 3)),	
      acc = 0;
    	if (unbalanced.type == 1),
   		for (i = 1:1:FL(1).N_level),
   	   for (j = 1:1:FL(2).N_level),
   	   for (k = 1:1:FL(3).N_level),
   	   for (l = 1:1:FL(4).UL(i).N_level),
 	      			
   			FI = acc + (j-1)*FL(3).N_level*FL(4).UL(i).N_level+(k-1)*FL(4).UL(i).N_level + l;   % file index
 	 		
 		   	GP(1, FI) = {FL(1).level(i).expr};
   		   GP(2, FI) = {FL(2).level(j).expr};
   		   GP(3, FI) = {FL(3).level(k).expr};
   		   GP(4, FI) = {FL(4).UL(i).n(l).expr}; 			
 			
   	 	end
 		   end
 	  	   end
   		acc = FI;
   		end
 		
   		for (m = 1:1:FL(1).N_level),
 	         counter(m) = 0;
      	end   %Just initiation
 		
   		for (i=1:1:file_num),

      	   flg = 0;
   			while flg == 0,
   	         fprintf('Input file for subject #%i ', i);
    	         file(i).nm = input('is: ', 's');
    			   fid = fopen (file(i).nm,'r');
   			   fprintf('\tThe corresponding level of factor 1 for subject %i ', i);
   			   file(i).F1L = input('is: ');   % level for factor 1
   			   if (isnumeric(file(i).F1L) == 0 | isempty(file(i).F1L)),
   			      flg = 0; fprintf(2,'Error: the input is not a number. Please try it again.\n');
   			   end

   			   for (m = 1:1:FL(1).N_level),
   			   if (file(i).F1L == m),
   				   counter(m) = counter(m) + 1;
 	   			   file(i).cntr(m) = counter(m);  % Any better way to implement this without the temporary array of counter(j)?
   				end   % set a counter of 4th factor level marker for later use
   			   end	
 			
   			   if (fid == -1),
      	         flg = 0; fprintf(2,'Error: File %s does not exist. Please try it again. \n', file(i).nm);
               else flg = 1; fclose (fid);	end
   		   end
   		end   % for (i=1:1:file_num),
 		
   		fprintf('\nIt is assumed that all files have the same subbrik order.');
         fprintf('\nFor each subbrik of an input file, specify the factor level of 2nd and 3rd factors. \n');

   	   for (j = 1:1:file_SB),
   	      fprintf('\n\tFor subbrik %i: \n', j-1);
   	      for (k = 1:1:(NF-2)),   % Here we only keep record of 2nd and 3rd factor levels
               fprintf('\t\tlevel of factor %c (%s) ', 65+k, FL(k+1).expr);
   	         SB(j).lv(k) = input('is: ');
   	      end   %for (j = 1:1:file_SB),
   	   end	 % for (j = 1:1:file_SB),			
 		
     end % if (unbalanced.type == 1),
     end % if ((NF == 4 & dsgn == 3))
     %end  % (unbalanced.yes == 1),
   else



   for (i=1:1:ntot),
		[err, file(i).v] = gind2sub (fliplr(sz), i);   %flip before subscripting. Doing flip because
         %I want vary the levels starting the last factor in stead of the first.
      scpt = fliplr(file(i).v);  %flip back to restore the original order	
      for (k=1:1:NF),
	      GP(k, i) = {FL(k).level(scpt(k)).expr};
      end
	end

	for (j = 1:1:FL(1).N_level),
	   counter(j) = 0;
	end   %Just initiation
   for (i = 1:1:file_num),
	  	flg = 0;
   	while flg == 0,
	      fprintf('Input file for subject #%i ', i);
	      file(i).nm = input('is: ', 's');
			fid = fopen (file(i).nm,'r');
			fprintf('\tThe corresponding level of factor 1 for subject %i ', i);
			file(i).F1L = input('is: ');   % level for factor 1
			if (isnumeric(file(i).F1L) == 0 | isempty(file(i).F1L)),
			   flg = 0; fprintf(2,'Error: the input is not a number. Please try it again.\n');
			end

			for (j = 1:1:FL(1).N_level),
			   if (file(i).F1L == j),
				   counter(j) = counter(j) + 1;
				   file(i).cntr(j) = counter(j);  % Any better way to implement this without the temporary array of counter(j)?
				end   % set a counter of 4th factor level marker for later use
			end	
			
			if (fid == -1),
   	      flg = 0; fprintf(2,'Error: File %s does not exist. Please try it again. \n', file(i).nm);
         else flg = 1; fclose (fid);	end
		end
	end
	fprintf('\nIt is assumed that all files have the same subbrik order.');
   fprintf('\nFor each subbrik of an input file, specify the factor level of 2nd and 3rd factors. \n');

	for (j = 1:1:file_SB),
	   fprintf('\n\tFor subbrik %i: \n', j-1);
	   for (k = 1:1:(NF-2)),   %Here we only keep record of 2nd and 3rd factor levels
         fprintf('\t\tlevel of factor %c (%s) ', 65+k, FL(k+1).expr);
	      SB(j).lv(k) = input('is: ');
	   end
	end 	
end % if (unbalanced.yes == 1)	
end	% if (file_format == 0),



%Obtain output file name	
flg = 0;
while flg == 0,
   OutFN = input('\nOutput file name (in bucket format): ', 's');
   OutFull = sprintf('%s+%s.HEAD', OutFN, format);
   fid2 = fopen(OutFull,'r');
   if (fid2 ~= -1),
      flg = 0; fprintf(2,'Error: File %s exists already. Please give another name. \n', OutFN);
   else flg = 1; end
end

% Gather contrast information: Everything stored in structure Contr

flg = 0;
while flg == 0,
   Contr.do = input('\nAny contrast test (1 - Yes, 0 - No)? ');
   if (Contr.do ~= 0 & Contr.do ~= 1),
      flg = 0; fprintf(2,'Error: invalid answer. Try it again. \n', OutFN);
   else flg = 1; end
end

if (Contr.do == 0),
   Contr.ord1.tot = 0;   % assign these for output later
	Contr.ord2.tot = 0;
	Contr.ord3.tot = 0;
	Contr.ord4.tot = 0;
else

	fprintf('\nNow coding contrasts:\n');
	fprintf('\n\t1. Each contrast should contain at least 2 terms;');
	fprintf('\n\t2. Each term in a contrast should have %i character(s), \n', NF);
	for (i = 1:1:NF), fprintf('\tNo. %i character corresponds to the level of factor %c (%s),\n', i, 'A'+i-1, FL(i).expr); end
	fprintf('\n\tUse 0 if a factor is collapsed.\n');
	fprintf('\n\tIf a factor level is smaller than 9, use its ordinal number;');
	fprintf('\n\tIf a factor level is bigger  than 9, use a, b, c, ... (no capitals) for 10, 11, 12, ... \n');	
	fprintf('\n\t3. All weights/coefficients in a contrast have to add up to zero.\n');

if (NF == 1),
   % 1st order contrasts
   flg = 0;
	
   while flg == 0,
      fprintf('\n1st order contrasts have %i factor(s) collapsed.\n', NF-1);
      Contr.ord1.tot = input('\nHow many 1st-order contrasts? (0 if none) ');
      if (isnumeric(Contr.ord1.tot) == 0 | Contr.ord1.tot < 0),
	      flg = 0; fprintf(2,'Error: inapproriate input. Please try again.\n');
   	else flg = 1;
	   end
   end
%	fprintf('\nUse factor level to code each term in a contrast.  For example, 1 means the factor is at first level.\n');
   for (i = 1:1:Contr.ord1.tot),
	   fprintf('\nLabel for 1st order contrast No. %i ', i);
		Contr.ord1.label(i).nm = input('is: ', 's');
		flg = 0;
		while flg == 0,
		Contr.ord1.cnt(i).NT = input('How many terms are involved in this contrast? ');   % NT = number of terms involved in this contrast
		if (isnumeric(Contr.ord1.cnt(i).NT) == 0 | Contr.ord1.cnt(i).NT < 2),
	      flg = 0; fprintf(2,'\nError: inapproriate input. Please try again.\n\n');
   	else flg = 1; end
      end
		
		flg0 = 0;
		while flg0 == 0,
		for (j = 1:1:Contr.ord1.cnt(i).NT),
		   flg = 0;
         while flg == 0,
			   fprintf('Factor index for No. %i term ', j);
			   Contr.ord1.cnt(i).code(j).str = input('is (e.g., 2): ', 's');
				if (length(Contr.ord1.cnt(i).code(j).str) ~= NF),
               flg = 0; fprintf(2,'Error: invalid input. Try it again. \n', OutFN);
            else flg = 1; end
			end
			Contr.ord1.cnt(i).coef(j) = input('Corresponding coefficient (e.g., 1 or -1): ');
		end
		if (abs(sum(Contr.ord1.cnt(i).coef)) > tol),
		   flg0 = 0; fprintf(2,'\nError: All the coefficients of a contrast have to sum up to 0! Try again...\n\n');
		else flg0 = 1; end
      end
	end	
	Contr.ord2.tot = 0;
	Contr.ord3.tot = 0;
	Contr.ord4.tot = 0;
   fprintf(1,'Done with 1st order contrast information.\n');
end

if (NF == 2),
   % 1st order contrasts ONLY at this point
   flg = 0;
   while flg == 0,
      fprintf('\n1st order contrasts have %i factor(s) collapsed.\n', NF-1);
      Contr.ord1.tot = input('\nHow many 1st-order contrasts? (0 if none) ');
      if (isnumeric(Contr.ord1.tot) == 0 | Contr.ord1.tot < 0),
	      flg = 0; fprintf(2,'\nError: inapproriate input. Please try again.\n');
   	else flg = 1;
	   end
   end
%	fprintf('\nUse factor level to code each term in a contrast.  For example, 1 means the factor is at first level.\n');
   for (i = 1:1:Contr.ord1.tot),
	   fprintf('\nLabel for 1st order contrast No. %i ', i);
		Contr.ord1.label(i).nm = input('is: ', 's');
		flg = 0;
		while flg == 0,
		Contr.ord1.cnt(i).NT = input('How many terms are involved in this contrast? ');   % NT = number of terms involved in this contrast
		if (isnumeric(Contr.ord1.cnt(i).NT) == 0 | Contr.ord1.cnt(i).NT < 2),
	      flg = 0; fprintf(2,'\nError: inapproriate input. Please try again.\n\n');
   	else flg = 1; end
      end

		flg0 = 0;
		while flg0 == 0,
		for (j = 1:1:Contr.ord1.cnt(i).NT),
		   flg = 0;
         while flg == 0,
			   fprintf('Factor index for No. %i term ', j);
			   Contr.ord1.cnt(i).code(j).str = input('is (e.g., 10): ', 's');
				if (length(Contr.ord1.cnt(i).code(j).str) ~= NF),
               flg = 0; fprintf(2,'\nError: invalid input. Try it again. \n', OutFN);
            else flg = 1; end
			end
			Contr.ord1.cnt(i).coef(j) = input('Corresponding coefficient (e.g., 1 or -1): ');
		end
		if (abs(sum(Contr.ord1.cnt(i).coef)) > tol),
		   flg0 = 0; fprintf(2,'\nError: All the coefficients of a contrast have to sum up to 0! Try again...\n\n');
		else flg0 = 1; end
      end
	end	
   fprintf(1,'Done with 1st order contrast information.\n');	
   % 2nd order contrasts
   flg = 0;
   while flg == 0,
      fprintf('\n2nd order contrasts have %i factor(s) collapsed.\n', NF-2);
	   fprintf('\nNotice: Contrasts for random factor are NOT feasible.\n');
      Contr.ord2.tot = input('\nHow many 2nd-order contrasts? (0 if none) ');
      if (isnumeric(Contr.ord2.tot) == 0 | Contr.ord2.tot < 0),
	      flg = 0; fprintf(2,'\nError: inapproriate input. Please try again.\n');
   	else flg = 1;
	   end
   end
%	fprintf('\nUse factor level to code each term in a contrast.  For example, 12 means ');
%	fprintf('\nthe 1st and 2nd factors are at their first and second level respectively.\n');
   for (i = 1:1:Contr.ord2.tot),
	   fprintf('\nLabel for 2nd order contrast No. %i: ', i);
		Contr.ord2.label(i).nm = input('is: ', 's');
		flg = 0;
		while flg == 0,
		Contr.ord2.cnt(i).NT = input('How many terms are involved in this contrast? ');   % NT = number of terms involved in this contrast
		if (isnumeric(Contr.ord2.cnt(i).NT) == 0 | Contr.ord2.cnt(i).NT < 2),
	      flg = 0; fprintf(2,'\nError: inapproriate input. Please try again.\n\n');
   	else flg = 1; end
      end

		flg0 = 0;
		while flg0 == 0,
		for (j = 1:1:Contr.ord2.cnt(i).NT),		
		   flg = 0;
         while flg == 0,
			   fprintf('Factor index for No. %i term ', j);
			   Contr.ord2.cnt(i).code(j).str = input('is (e.g., 12): ', 's');
				if (length(Contr.ord2.cnt(i).code(j).str) ~= NF),
               flg = 0; fprintf(2,'\nError: invalid input. Try it again. \n', OutFN);
            else flg = 1; end
			end
			Contr.ord2.cnt(i).coef(j) = input('Corresponding coefficient (e.g., 1 or -1): ');
		end
		if (abs(sum(Contr.ord2.cnt(i).coef)) > tol),
		   flg0 = 0; fprintf(2,'\nError: All the coefficients of a contrast have to sum up to 0! Try again...\n\n');
		else flg0 = 1; end
      end
	end	

	Contr.ord3.tot = 0;
	Contr.ord4.tot = 0;
   fprintf(1,'Done with 2nd order contrast information.\n');
end

if (NF == 3 | NF == 4),

   % 1st order contrasts
   flg = 0;
   while flg == 0,
      fprintf('\n1st order contrasts have %i factor(s) collapsed.\n', NF-1);
   	%fprintf('\t1. Simple mean for a factor level, such as [1 0 0];\n');
	   %fprintf('\t2. Difference between two factor levels, such as [1 -1 0];\n');
   	%fprintf('\t3. Linear combination of factor levels, such as [0.5 0.5 -1];\n');
	   fprintf('\nNotice: Contrasts for random factor are NOT feasible.\n');
      Contr.ord1.tot = input('\nHow many 1st-order contrasts? (0 if none) ');
      if (isnumeric(Contr.ord1.tot) == 0 | Contr.ord1.tot < 0),
	      flg = 0; fprintf(2,'\nError: inapproriate input. Please try again.\n');
   	else flg = 1;
	   end
   end
%	fprintf('\nUse factor level to code each term in a contrast.  For example, 0100 means the first, third ');
%	fprintf('\nand fourth factors are collapsed while 2nd factor is at first level.\n');
   for (i = 1:1:Contr.ord1.tot),
	   fprintf('\nLabel for 1st order contrast No. %i ', i);
		Contr.ord1.label(i).nm = input('is: ', 's');
		flg = 0;
		while flg == 0,
		Contr.ord1.cnt(i).NT = input('How many terms are involved in this contrast? ');   % NT = number of terms involved in this contrast
		if (isnumeric(Contr.ord1.cnt(i).NT) == 0 | Contr.ord1.cnt(i).NT < 2),
	      flg = 0; fprintf(2,'\nError: inapproriate input. Please try again.\n\n');
   	else flg = 1; end
      end

		flg0 = 0;
		while flg0 == 0,
		for (j = 1:1:Contr.ord1.cnt(i).NT),
		   flg = 0;
         while flg == 0,
			   fprintf('Factor index for No. %i term ', j);
			   Contr.ord1.cnt(i).code(j).str = input('is (e.g., 0020): ', 's');
				if (length(Contr.ord1.cnt(i).code(j).str) ~= NF),
               flg = 0; fprintf(2,'\nError: invalid input. Try it again. \n', OutFN);
            else flg = 1; end
			end
			Contr.ord1.cnt(i).coef(j) = input('Corresponding coefficient (e.g., 1 or -1): ');
		end
		if (abs(sum(Contr.ord1.cnt(i).coef)) > tol),
	      flg0 = 0; fprintf(2,'\nError: All the coefficients of a contrast have to sum up to 0! Try again...\n\n');
   	else flg0 = 1; end
     end	
	end
   fprintf(1,'Done with 1st order contrast information.\n');

   % 2nd order contrasts
   flg = 0;
   while flg == 0,
      fprintf('\n2nd order contrasts have %i factor(s) collapsed.\n', NF-2);
	   fprintf('\nNotice: Contrasts for random factor are NOT feasible.\n');
      Contr.ord2.tot = input('\nHow many 2nd-order contrasts? (0 if none) ');
      if (isnumeric(Contr.ord2.tot) == 0 | Contr.ord2.tot < 0),
	      flg = 0; fprintf(2,'\nError: inapproriate input. Please try again.\n');
   	else flg = 1;
	   end
   end
%	fprintf('\nUse factor level to code each term in a contrast.  For example, 0120 means both the first ');
%	fprintf('\nand fourth factors are collapsed while 2nd and 3rd factors are at first and second level respectively.\n');
   for (i = 1:1:Contr.ord2.tot),
	   fprintf('\nLabel for 2nd order contrast No. %i ', i);
		Contr.ord2.label(i).nm = input('is: ', 's');
		flg = 0;
		while flg == 0,
		Contr.ord2.cnt(i).NT = input('How many terms are involved in this contrast? ');   % NT = number of terms involved in this contrast
		if (isnumeric(Contr.ord2.cnt(i).NT) == 0 | Contr.ord2.cnt(i).NT < 2),
	      flg = 0; fprintf(2,'\nError: inapproriate input. Please try again.\n\n');
   	else flg = 1; end
      end

		flg0 = 0;
		while flg0 == 0,
		for (j = 1:1:Contr.ord2.cnt(i).NT),
		   flg = 0;
         while flg == 0,
			   fprintf('Factor index for No. %i term ', j);
			   Contr.ord2.cnt(i).code(j).str = input('is (e.g., 0100): ', 's');
				if (length(Contr.ord2.cnt(i).code(j).str) ~= NF),
               flg = 0; fprintf(2,'Error: invalid input. Try it again. \n', OutFN);
            else flg = 1; end
			end
			Contr.ord2.cnt(i).coef(j) = input('Corresponding coefficient (e.g., 1 or -1): ');
		end
		if (abs(sum(Contr.ord2.cnt(i).coef)) > tol),
	      flg0 = 0; fprintf(2,'\nError: All the coefficients of a contrast have to sum up to 0! Try again...\n\n');
	   else flg0 = 1; end
   end  %flg0 = 0;
	end

   fprintf(1,'Done with 2nd order contrast information.\n');
	
	
	% 3rd-order contrasts
	flg = 0;
   while flg == 0,
      fprintf('\n3rd order contrasts have %i factor(s) collapsed.\n', NF-3);
	   fprintf('\nNotice: Contrasts for random factor are NOT feasible.\n');
      Contr.ord3.tot = input('\nHow many 3rd-order contrasts? (0 if none) ');
      if (isnumeric(Contr.ord3.tot) == 0 | Contr.ord3.tot < 0),
	      flg = 0; fprintf(2,'\nError: inapproriate input. Please try again.\n');
   	else flg = 1;
	   end
   end
%	fprintf('\nUse factor level to code each term in a contrast.  For example, 1230 means the fourth');
%	fprintf('\nfactor is collapsed while all the other 3 factors are at first, second and third level respectively.\n');
   for (i = 1:1:Contr.ord3.tot),
	   fprintf('\nLabel for 3rd order contrast No. %i ', i);
		Contr.ord3.label(i).nm = input('is: ', 's');
		flg = 0;
		while flg == 0,
		Contr.ord3.cnt(i).NT = input('How many terms are involved in this contrast? ');   % NT = number of terms involved in this contrast
		if (isnumeric(Contr.ord3.cnt(i).NT) == 0 | Contr.ord3.cnt(i).NT < 2),
	      flg = 0; fprintf(2,'\nError: inapproriate input. Please try again.\n\n');
   	else flg = 1; end
      end

		flg0 = 0;
		while flg0 == 0,
		for (j = 1:1:Contr.ord3.cnt(i).NT),
		   flg = 0;
         while flg == 0,
			   fprintf('Factor index for No. %i term ', j);
			   Contr.ord3.cnt(i).code(j).str = input('is (e.g., 1230): ', 's');
				if (length(Contr.ord3.cnt(i).code(j).str) ~= NF),
               flg = 0; fprintf(2,'\nError: invalid input. Try it again. \n', OutFN);
            else flg = 1; end
			end
			Contr.ord3.cnt(i).coef(j) = input('Corresponding coefficient (e.g., 1 or -1): ');
		end
		if (abs(sum(Contr.ord3.cnt(i).coef)) > tol),
	      flg0 = 0; fprintf(2,'\nError: All the coefficients of a contrast have to sum up to 0! Try again...\n\n');
	   else flg0 = 1; end
      end  %flg0 = 0;
	end
	Contr.ord4.tot = 0;	
   fprintf(1,'Done with 3rd order contrast information.\n');
	
end  % if (NF == 4)

if (NF == 5),

   % 1st order contrasts
   flg = 0;
   while flg == 0,
      fprintf('\n1st order contrasts have %i factor(s) collapsed.\n', NF-1);
   	%fprintf('\t1. Simple mean for a factor level, such as [1 0 0];\n');
	   %fprintf('\t2. Difference between two factor levels, such as [1 -1 0];\n');
   	%fprintf('\t3. Linear combination of factor levels, such as [0.5 0.5 -1];\n');
	   fprintf('\nNotice: Contrasts for random factor are NOT feasible.\n');
      Contr.ord1.tot = input('\nHow many 1st-order contrasts? (0 if none) ');
      if (isnumeric(Contr.ord1.tot) == 0 | Contr.ord1.tot < 0),
	      flg = 0; fprintf(2,'\nError: inapproriate input. Please try again.\n');
   	else flg = 1;
	   end
   end
%	fprintf('\nUse factor level to code each term in a contrast.  For example, 0100 means the first, third ');
%	fprintf('\nand fourth factors are collapsed while 2nd factor is at first level.\n');
   for (i = 1:1:Contr.ord1.tot),
	   fprintf('\nLabel for 1st order contrast No. %i ', i);
		Contr.ord1.label(i).nm = input('is: ', 's');
		flg = 0;
		while flg == 0,
		Contr.ord1.cnt(i).NT = input('How many terms are involved in this contrast? ');   % NT = number of terms involved in this contrast
		if (isnumeric(Contr.ord1.cnt(i).NT) == 0 | Contr.ord1.cnt(i).NT < 2),
	      flg = 0; fprintf(2,'\nError: inapproriate input. Please try again.\n\n');
   	else flg = 1; end
      end

		flg0 = 0;
		while flg0 == 0,
		for (j = 1:1:Contr.ord1.cnt(i).NT),
		   flg = 0;
         while flg == 0,
			   fprintf('Factor index for No. %i term ', j);
			   Contr.ord1.cnt(i).code(j).str = input('is (e.g., 00200): ', 's');
				if (length(Contr.ord1.cnt(i).code(j).str) ~= NF),
               flg = 0; fprintf(2,'Error: invalid input. Try it again. \n', OutFN);
            else flg = 1; end
			end
			Contr.ord1.cnt(i).coef(j) = input('Corresponding coefficient (e.g., 1 or -1): ');
		end  % for (j = 1:1:Contr.ord1.cnt(i).NT)
		if (abs(sum(Contr.ord1.cnt(i).coef)) > tol),
	      flg0 = 0; fprintf(2,'\nError: All the coefficients of a contrast have to sum up to 0! Try again...\n\n');
	   else flg0 = 1; end
      end  %flg0 = 0;
	end	
   fprintf(1,'Done with 1st order contrast information.\n');

   % 2nd order contrasts
   flg = 0;
   while flg == 0,
      fprintf('\n2nd order contrasts have %i factor(s) collapsed.\n', NF-2);
	   fprintf('\nNotice: Contrasts for random factor are NOT feasible.\n');
      Contr.ord2.tot = input('\nHow many 2nd-order contrasts? (0 if none) ');
      if (isnumeric(Contr.ord2.tot) == 0 | Contr.ord2.tot < 0),
	      flg = 0; fprintf(2,'Error: inapproriate input. Please try again.\n');
   	else flg = 1;
	   end
   end  % while flg == 0
%	fprintf('\nUse factor level to code each term in a contrast.  For example, 0120 means both the first ');
%	fprintf('\nand fourth factors are collapsed while 2nd and 3rd factors are at first and second level respectively.\n');
   for (i = 1:1:Contr.ord2.tot),
	   fprintf('\nLabel for 2nd order contrast No. %i ', i);
		Contr.ord2.label(i).nm = input('is: ', 's');
		flg = 0;
		while flg == 0,
		Contr.ord2.cnt(i).NT = input('How many terms are involved in this contrast? ');   % NT = number of terms involved in this contrast
		if (isnumeric(Contr.ord2.cnt(i).NT) == 0 | Contr.ord2.cnt(i).NT < 2),
	      flg = 0; fprintf(2,'\nError: inapproriate input. Please try again.\n\n');
   	else flg = 1; end
      end

		flg0 = 0;
		while flg0 == 0,
		for (j = 1:1:Contr.ord2.cnt(i).NT),
		   flg = 0;
         while flg == 0,
			   fprintf('Factor index for No. %i term ', j);
			   Contr.ord2.cnt(i).code(j).str = input('is (e.g., 01200): ', 's');
				if (length(Contr.ord2.cnt(i).code(j).str) ~= NF),
               flg = 0; fprintf(2,'\nError: invalid input. Try it again. \n', OutFN);
            else flg = 1; end
			end  % while flg == 0
			Contr.ord2.cnt(i).coef(j) = input('Corresponding coefficient (e.g., 1 or -1): ');
		end  % for (j = 1:1:Contr.ord2.cnt(i).NT)
		if (abs(sum(Contr.ord2.cnt(i).coef)) > tol),
	      flg0 = 0; fprintf(2,'\nError: All the coefficients of a contrast have to sum up to 0! Try again...\n\n');
	   else flg0 = 1; end
      end  %flg0 = 0;
	end	 % for (i = 1:1:Contr.ord2.tot)
   fprintf(1,'Done with 2nd order contrast information.\n');
	
	
	% 3rd-order contrasts
	flg = 0;
   while flg == 0,
      fprintf('\n3rd order contrasts have %i factor(s) collapsed.\n', NF-3);
	   fprintf('\nNotice: Contrasts for random factor are NOT feasible.\n');
      Contr.ord3.tot = input('\nHow many 3rd-order contrasts? (0 if none) ');
      if (isnumeric(Contr.ord3.tot) == 0 | Contr.ord3.tot < 0),
	      flg = 0; fprintf(2,'\nError: inapproriate input. Please try again.\n');
   	else flg = 1;
	   end
   end
%	fprintf('\nUse factor level to code each term in a contrast.  For example, 1230 means the fourth');
%	fprintf('\nfactor is collapsed while all the other 3 factors are at first, second and third level respectively.\n');
   for (i = 1:1:Contr.ord3.tot),
	   fprintf('\nLabel for 3rd order contrast No. %i ', i);
		Contr.ord3.label(i).nm = input('is: ', 's');
		flg = 0;
		while flg == 0,
		Contr.ord3.cnt(i).NT = input('How many terms are involved in this contrast? ');   % NT = number of terms involved in this contrast
		if (isnumeric(Contr.ord3.cnt(i).NT) == 0 | Contr.ord3.cnt(i).NT < 2),
	      flg = 0; fprintf(2,'\nError: inapproriate input. Please try again.\n\n');
   	else flg = 1; end
      end

		flg0 = 0;
		while flg0 == 0,
		for (j = 1:1:Contr.ord3.cnt(i).NT),
		   flg = 0;
         while flg == 0,
			   fprintf('Factor index for No. %i term ', j);
			   Contr.ord3.cnt(i).code(j).str = input('is (e.g., 1230): ', 's');
				if (length(Contr.ord3.cnt(i).code(j).str) ~= NF),
               flg = 0; fprintf(2,'\nError: invalid input. Try it again. \n', OutFN);
            else flg = 1; end
			end
			Contr.ord3.cnt(i).coef(j) = input('Corresponding coefficient (e.g., 1 or -1): ');
		end
		if (abs(sum(Contr.ord3.cnt(i).coef)) > tol),
	      flg0 = 0; fprintf(2,'\nError: All the coefficients of a contrast have to sum up to 0! Try again...\n\n');
	   else flg0 = 1; end
      end  %flg0 = 0;
	end	
   fprintf(1,'Done with 3rd order contrast information.\n');
	
	% 4th-order contrasts
	flg = 0;
   while flg == 0,
      fprintf('\n4th order contrasts have %i factor(s) collapsed.\n', NF-4);
	   fprintf('\nNotice: Contrasts for random factor are NOT feasible.\n');
      Contr.ord4.tot = input('\nHow many 4th-order contrasts? (0 if none) ');
      if (isnumeric(Contr.ord4.tot) == 0 | Contr.ord4.tot < 0),
	      flg = 0; fprintf(2,'\nError: inapproriate input. Please try again.\n');
   	else flg = 1;
	   end
   end
%	fprintf('\nUse factor level to code each term in a contrast.  For example, 1230 means the fourth');
%	fprintf('\nfactor is collapsed while all the other 3 factors are at first, second and third level respectively.\n');
   for (i = 1:1:Contr.ord4.tot),
	   fprintf('\nLabel for 4th order contrast No. %i ', i);
		Contr.ord4.label(i).nm = input('is: ', 's');
		flg = 0;
		while flg == 0,
		Contr.ord4.cnt(i).NT = input('How many terms are involved in this contrast? ');   % NT = number of terms involved in this contrast
		if (isnumeric(Contr.ord4.cnt(i).NT) == 0 | Contr.ord4.cnt(i).NT < 2),
	      flg = 0; fprintf(2,'\nError: inapproriate input. Please try again.\n\n');
   	else flg = 1; end
      end

		flg0 = 0;
		while flg0 == 0,
		for (j = 1:1:Contr.ord4.cnt(i).NT),
		   flg = 0;
         while flg == 0,
			   fprintf('Factor index for No. %i term ', j);
			   Contr.ord4.cnt(i).code(j).str = input('is (e.g., 01230): ', 's');
				if (length(Contr.ord4.cnt(i).code(j).str) ~= NF),
               flg = 0; fprintf(2,'\nError: invalid input. Try it again. \n', OutFN);
            else flg = 1; end
			end
			Contr.ord4.cnt(i).coef(j) = input('Corresponding coefficient (e.g., 1 or -1): ');
		end
		if (abs(sum(Contr.ord4.cnt(i).coef)) > tol),
	      flg0 = 0; fprintf(2,'\nError: All the coefficients of a contrast have to sum up to 0! Try again...\n\n');
	   else flg0 = 1; end
      end  %flg0 = 0;
	end	
   fprintf(1,'Done with 4th order contrast information.\n');
	
end  % if (NF == 5)


end  % end of if (Contrast)

% Done for gathering info from the user.

%=========================================

cov.marker = 10000;  % temporary solution when no covariate exists in PreProc.m

if (cov.do),  % new design type after covariate is added: Seems this section not needed?
   switch NF
	   case 1,
		   NF = 2;   % 1-ay ANCOVA
			dsgn = 1;	
		case 2,
		   NF = 3;   % 2-ay ANCOVA
			if (dsgn == 1), dsgn = 1; end
			if (dsgn == 2), dsgn = 2; end
		case 3,
		   NF = 4;   % 3-ay ANCOVA
			if (dsgn == 1), dsgn = 1; end
			if (dsgn == 2), dsgn = 2; end
			if (dsgn == 3), dsgn = 3; end			
	end			
end

for (i = 1:1:NF),
   if (i == NF & cov.do),
	   group(NF) = {cov.vec'};    % Convert the cloumn into one row
      varnames(NF) = {cov.label};
   else
      group(i) = {GP(i,:)};
      varnames(i) = {FL(i).expr};
	end
end

if (cov.do),         % swap the covariate with another factor for some special cases: more elegant way to do this?

   switch NF
   case 2,    % 1-way ANCOVA
	
	cov.marker = 2;
	
	for (i = 1:1:Contr.ord1.tot),    %1st order contrasts
   for (j = 1:1:Contr.ord1.cnt(i).NT),
	   Contr.ord1.cnt(i).code(j).str = [Contr.ord1.cnt(i).code(j).str(1) '0'];
	end
	end
	
	case 3,    % 2-way ANCOVA
	
	cov.marker = 3;    %default
	
   if (dsgn == 2),   % for 2-way ANCOVA with design of A(F)XB(R), make covarite the 2nd factor so that 2-way ANCOVA would be 3-way ANOVA of A(F)XB(F)XC(R);
      temp1 = group{2};
		group{2} = group{3};
		group{3} = temp1;
		
		temp2 = varnames(2);
		varnames(2) = varnames(3);
		varnames(3) = temp2;
		
		temp3 = FL(2);
		FL(2) = FL(3);
		FL(3) = temp3;
		
		cov.marker = 2;    % now 2nd factor is covariate
		
	for (i = 1:1:Contr.ord1.tot),
     for (j = 1:1:Contr.ord1.cnt(i).NT),
	      Contr.ord1.cnt(i).code(j).str = [Contr.ord1.cnt(i).code(j).str(1) '0' Contr.ord1.cnt(i).code(j).str(2)];
		end
	end
	
	for (i = 1:1:Contr.ord2.tot),
      for (j = 1:1:Contr.ord2.cnt(i).NT),
	      Contr.ord2.cnt(i).code(j).str = [Contr.ord2.cnt(i).code(j).str(1) '0' Contr.ord2.cnt(i).code(j).str(2)];
	   end
	end
	
	else   % other two designs
		for (i = 1:1:Contr.ord1.tot),
      for (j = 1:1:Contr.ord1.cnt(i).NT),
	      Contr.ord1.cnt(i).code(j).str = [Contr.ord1.cnt(i).code(j).str(1) Contr.ord1.cnt(i).code(j).str(2) '0'];
		end
		end

   	for (i = 1:1:Contr.ord2.tot),
      for (j = 1:1:Contr.ord2.cnt(i).NT),
	      Contr.ord2.cnt(i).code(j).str = [Contr.ord2.cnt(i).code(j).str(1) Contr.ord2.cnt(i).code(j).str(2) '0'];
	   end
	   end			
		
	end   % case 3
		
	case 4,	
	
	cov.marker = 4;    %default
						
	if (dsgn == 2 | dsgn == 3), % for 3-way ANCOVA, swap the 3rd factor and the covariate
	   temp1 = group{3};
		group{3} = group{4};
		group{4} = temp1;
		
		temp2 = varnames(3);
		varnames(3) = varnames(4);
		varnames(4) = temp2;
		
		temp3 = FL(3);
		FL(3) = FL(4);
		FL(4) = temp3;
		
		cov.marker = 3;
		
		for (i = 1:1:Contr.ord1.tot),
      for (j = 1:1:Contr.ord1.cnt(i).NT),
	      Contr.ord1.cnt(i).code(j).str = [Contr.ord1.cnt(i).code(j).str(1) Contr.ord1.cnt(i).code(j).str(2) '0' Contr.ord1.cnt(i).code(j).str(3)];
		end
		end
		
		for (i = 1:1:Contr.ord2.tot),
      for (j = 1:1:Contr.ord2.cnt(i).NT),
	      Contr.ord2.cnt(i).code(j).str = [Contr.ord2.cnt(i).code(j).str(1) Contr.ord2.cnt(i).code(j).str(2) '0' Contr.ord2.cnt(i).code(j).str(3)];
		end
		end
		
		for (i = 1:1:Contr.ord3.tot),
      for (j = 1:1:Contr.ord3.cnt(i).NT),
	      Contr.ord3.cnt(i).code(j).str = [Contr.ord3.cnt(i).code(j).str(1) Contr.ord3.cnt(i).code(j).str(2) '0' Contr.ord3.cnt(i).code(j).str(3)];
		end
		end
				
	else
		for (i = 1:1:Contr.ord1.tot),
      for (j = 1:1:Contr.ord1.cnt(i).NT),
	      Contr.ord1.cnt(i).code(j).str = [Contr.ord1.cnt(i).code(j).str(1) Contr.ord1.cnt(i).code(j).str(2) Contr.ord1.cnt(i).code(j).str(3) '0'];
		end
		end
		
		for (i = 1:1:Contr.ord2.tot),
      for (j = 1:1:Contr.ord2.cnt(i).NT),
	      Contr.ord2.cnt(i).code(j).str = [Contr.ord2.cnt(i).code(j).str(1) Contr.ord2.cnt(i).code(j).str(2) Contr.ord2.cnt(i).code(j).str(3) '0'];
		end
		end
		
		for (i = 1:1:Contr.ord3.tot),
      for (j = 1:1:Contr.ord3.cnt(i).NT),
	      Contr.ord3.cnt(i).code(j).str = [Contr.ord3.cnt(i).code(j).str(1) Contr.ord3.cnt(i).code(j).str(2) Contr.ord3.cnt(i).code(j).str(3) '0'];
		end
		end	
	end  % case 4

	end % switch
end

N_Brik = 2^NF - 1;  % For 4way with crossed design, there are 15 effect terms: A, B, C, D, AB, AC, AD, BC, BD, ABC, ABD, BCD, and ABCD.
switch NF
   case 1,	
	
	case 2,
	
	case 3,
	   if (dsgn == 3 | dsgn == 4), N_Brik = 5; end
				
	case 4,
      if (dsgn == 3 | dsgn == 4), N_Brik = 11; end
		if (dsgn == 5), N_Brik = 9; end
		
	case 5,
	   if (dsgn == 3 | dsgn == 4), N_Brik = 23; end	
		
end

%Voxel independent stuff
[err, Qd, s, termname, nterms, sindices, dfbothSS, modw, modwo, tnames, dfterm, dfe, Contr] = PreProc(ntot, NF, group, varnames, FL, Contr, cov, unbalanced);
%[err, Qd, s, termname, nterms, sindices, dfbothSS, modw, modwo, tnames, dfterm, dfe, Contr] = PreProc(ntot, NF, group, varnames, FL, Contr, dsgn, cov);

%Use Ziad's function BrikLoad to load all 'ntot' of datasets in a column format
%since matlab function 'anovan' only accepts vectors.

Opt.Format = 'matrix';

t0 = clock;

[err, FileInfo] = BrikInfo(file(1).nm);
if (data_type == 0),
   slices = FileInfo.DATASET_DIMENSIONS(3);    % Get the number of slices along Z axis from file header
elseif (data_type == 1 | data_type == 2 ),
   if (Frame_N == 1), Opt.method = 1; % if being purified
	else Opt.method = 2; end
   if (data_type == 2 ),
      %finalize NI
      %I would have liked to do it when the user enters the number or the
      %filename but the function: newid() seems to hang when the users
      %use cut and paste to answer all the questions at once - don't know
      %why that is.
      node_n=str2num(node_file);
      NI = [];
      if (isempty(node_n)),
         ropt.method = 3;
         ropt.verb = 2;
         NI = Read_1D(node_file,ropt);
         node_n = length(NI);
      else
         NI = [0:1:node_n-1];
      end
      if (isempty(NI)),
	      flg = 0; fprintf(2,'Error: the input is not a number, or 1D file reading failed. \n');
         return;
	   end
   end

   Opt.SliceSize_1D = 50000;  % each time run 50000 nodes due to memory limit
	slices =  ceil(node_n/Opt.SliceSize_1D);
end

fprintf(1, '\nTotal slices along Z axis: %i - You can estimate the total runtime\n', slices);
fprintf(1, '\tby multiplying the runtime for each slice from down below.\n');
fprintf(1, '\nRunning analysis on slice:\n');

for (sn = 1:1:slices),
	tic,

	fprintf('\t#%d... ', sn);

	if (file_format == 1),   % Each file contains only single subbrik
%   fprintf(1,'\nReading %d input files...', ntot);
		Opt.Slices = sn;   % Right now just try once slice
%		Opt.Frames = 1;
	   for (i=1:1:ntot),      % Read in each file
		 %V(i) is supposed to be a N X M X K matrix? The information regarding the dimensions
		 %should be in Info: Info.DATASET_DIMENSIONS(1), Info.DATASET_DIMENSIONS(2),
		 %Info.DATASET_DIMENSIONS(3), Info.DATASET_RANK(2)?
%      [err, X(:, :, :, i), Info, ErrMessage] = BrikLoad(file(i).nm, Opt);
			
			if (data_type == 0),  % Volume data
			   [err, X(:, :, i), Info, ErrMessage] = BrikLoad(file(i).nm, Opt);
			elseif (data_type == 1 | data_type == 2),	% Surface data
			   Opt.Frames = Frame_N;
				[err, Z, Info, ErrMessage] = BrikLoad(file(i).nm, Opt);
				if (sn ~= 1), X(:,1,i) = zeros(size(X(:,1,i))); end
				% For the last chunk of nodes, which are not a whole set of 50,000.
				% Not an elegant solution: This would create some dangling 0's in
            % the output file!
				X(1:size(Z,1), 1, i) = Z;
				clear Z;
			end	
				
	   	if (err == 1),
         	fprintf(2,'Error: failure on loading file %s -- %s. \n', file(i).nm, ErrMessage);
			while (1); fprintf(2,'Halted: Ctrl+c to exit'); pause; end
      	end	
   	end
	end

	if (file_format == 0),   % Each file contains mutliple subbriks from one subject
%   fprintf(1,'\nReading slice #%d from %d input files... ', sn, file_num);
		Opt.Slices = sn;   % Right now just try once slice
	% Create a subbrik array. Here we assume the first file_SB subbriks are for ANOVA. Need to change for general situation?
		for (i = 1:1:file_SB), Opt.Frames = [Opt.Frames, i]; end	
	
	   for (i=1:1:file_num),      % Read in each file. For each subject (file), there should have FL(1).N_level*FL(4).N_level files
%      [err, tmp(:, :, :, (i-1)*file_SB+1:i*file_SB), Info, ErrMessage] = BrikLoad(file(i).nm,Opt);
   	   [err, tmp(:, :, (i-1)*file_SB+1:i*file_SB), Info, ErrMessage] = BrikLoad(file(i).nm,Opt);
		   if (err == 1),
	   	   fprintf(2,'Error: failure on loading file %s -- %s. \n', file(i).nm, ErrMessage);
			   while (1); fprintf(2,'Halted: Ctrl+c to exit'); pause; end
	      end
%	   for (j = 1:1:file_SB),       % For each subbrik, there should have FL(2).N_level*FL(3).N_level subbriks
%	      X(:, :, :, (file(i).F1L-1)*file_SB*FL(4).N_level + (SB(j).lv(1)-1)*FL(3).N_level*FL(4).N_level + ...
%			   (SB(j).lv(2)-1)*FL(4).N_level + file(i).cntr(file(i).F1L)) = tmp(:, :, :, (i-1)*file_SB+j);
		
		   if (NF == 4),				
			for (j = 1:1:file_SB),       % For each subbrik, there should have FL(2).N_level*FL(3).N_level subbriks
	      	X(:, :, (file(i).F1L-1)*file_SB*FL(4).N_level + (SB(j).lv(1)-1)*FL(3).N_level*FL(4).N_level + ...
			   	(SB(j).lv(2)-1)*FL(4).N_level + file(i).cntr(file(i).F1L)) = tmp(:, :, (i-1)*file_SB+j);		
		   end
			else fprintf(2,'Sorry, I cannot handle mutliple subbriks for this case. \n');
			   while (1); fprintf(2,'Halted: Ctrl+c to exit'); pause; end
	      end
					
  	    end

	% Here we run 4-way ANOVA one slice a time due to swap memory issue. In matrix tmp, the first two
	% dimensions are voxels along X and Y axes while the 4th dimension is subbrik tacked together for
	% each subject. We need to shuffle the 4th dimension in tmp so that the 4th dimention is arranged
	% in the following order: varying 4th factor first, then 3rd, then 2nd and 1st factors.
	
	% 1st factor: (file(i).F1L-1)*file_SB*FL(4).N_level
	% 2nd factor: (SB(j).lv(1)-1)*FL(3).N_level*FL(4).N_level
	% 3rd factor: (SB(j).lv(2)-1)*FL(4).N_level
	% 4th factor: file(i).cntr(file(i).F1L)
	
%			X(:, :, :, (i-1)*file_SB*FL(1).N_level + (j-1)*FL(3).N_level*FL(4).N_level + (k-1)*FL(4).N_level + l)	
% Here i, j, k, and l are indeces for 1st, 2nd, 3rd and 4th factor levels.

% file_SB = FL(3).N_level*FL(2).N_level   Maybe I should set a checking condition to make sure this equality does hold?

% Is there a better way to load those subbriks directly into matrix X instead of relaying them into
% matrix tmp???
			
	end

	D1 = Info.DATASET_DIMENSIONS(1);      % X
	D2 = Info.DATASET_DIMENSIONS(2);		% Y
	D3 = Info.DATASET_DIMENSIONS(3);		% Z
	dim = size(X);
	
%fstat = zeros(D1, D2, D3, N_Brik);
	fstat = zeros(D1, D2, N_Brik);
%intensity = zeros(D1, D2, D3, N_Brik);
	intensity = zeros(D1, D2, N_Brik);

%Y = reshape(X, ntot, D1, D2);  % Have not figured out how to change the shape of the matrix yet!
	for (i = 1:1:D1),
	   for (j = 1:1:D2),
%for (k = 1:1:D3),
%   [err,fstat(i, j, k, :), intensity(i, j, k, :), dfterm_new, dfdenom, tnames_new] = SumsOfSquares(reshape(X(i, j, k, :), ntot, 1), ...
%	   nterms, Qd, s, sindices, dfbothSS, modw, modwo, tnames, dfterm, dfe, Nest);
		   [err, fstat(i, j, :), intensity(i, j, :), dfterm_new, dfdenom, tnames_new, LC(i, j)] = SumsOfSquares(reshape(X(i, j, :), ntot, 1), ...
	   	   NF, FL, ntot, nterms, Qd, s, sindices, dfbothSS, modw, modwo, tnames, dfterm, dfe, dsgn, N_Brik, Contr);	
%end
   	end
	end


%====================================

%M = zeros(D1, D2, D3, 2*N_Brik);	%initialization

%   if (cov.do == 1),
%      M = zeros(D1, D2, 2*(N_Brik+Contr.ord1.tot+Contr.ord2.tot+Contr.ord3.tot+1));	%initialization
%	else M = zeros(D1, D2, 2*(N_Brik+Contr.ord1.tot+Contr.ord2.tot+Contr.ord3.tot));	%initialization
%	end

   M = zeros(D1, D2, 2*(N_Brik+Contr.ord1.tot+Contr.ord2.tot+Contr.ord3.tot+Contr.ord4.tot));	%initialization

% Assemble a matrix with intensity and F value sandwiched with each other
	for (i = 1:1:D1),
		for (j = 1:1:D2),
%for (k = 1:1:D3),
		   for (l = 1:1:N_Brik),
%   M(i, j, k, 2*l-1) = intensity(i, j, k, l)*isfinite(intensity(i, j, k, l));
			   M(i, j, 2*l-1) = intensity(i, j, l)*isfinite(intensity(i, j, l));
%	   M(i, j, k, 2*l) = fstat(i, j, k, l)*isfinite(fstat(i, j, k, l));
	   		M(i, j, 2*l) = fstat(i, j, l)*isfinite(fstat(i, j, l));
		   end
			if (Contr.ord1.tot>0),
   		for (l = 1:1:Contr.ord1.tot),
			   M(i, j, 2*(l+N_Brik)-1) = LC(i, j).t1(l).value*isfinite(LC(i, j).t1(l).value);
	   		M(i, j, 2*(l+N_Brik)) = LC(i, j).t1(l).t*isfinite(LC(i, j).t1(l).t);
   		end
			end
			if (Contr.ord2.tot>0),
			for (l = 1:1:Contr.ord2.tot),
			   M(i, j, 2*(l+N_Brik+Contr.ord1.tot)-1) = LC(i, j).t2(l).value*isfinite(LC(i, j).t2(l).value);
	   		M(i, j, 2*(l+N_Brik+Contr.ord1.tot)) = LC(i, j).t2(l).t*isfinite(LC(i, j).t2(l).t);
   		end
			end
			if (Contr.ord3.tot>0),
			for (l = 1:1:Contr.ord3.tot),
			   M(i, j, 2*(l+N_Brik+Contr.ord1.tot+Contr.ord2.tot)-1) = LC(i, j).t3(l).value*isfinite(LC(i, j).t3(l).value);
	   		M(i, j, 2*(l+N_Brik+Contr.ord1.tot+Contr.ord2.tot)) = LC(i, j).t3(l).t*isfinite(LC(i, j).t3(l).t);
   		end
			end
			if (Contr.ord4.tot>0),
			for (l = 1:1:Contr.ord4.tot),
			   M(i, j, 2*(l+N_Brik+Contr.ord1.tot+Contr.ord2.tot+Contr.ord3.tot)-1) = LC(i, j).t4(l).value*isfinite(LC(i, j).t4(l).value);
	   		M(i, j, 2*(l+N_Brik+Contr.ord1.tot+Contr.ord2.tot+Contr.ord3.tot)) = LC(i, j).t4(l).t*isfinite(LC(i, j).t4(l).t);
   		end
			end	
		end
	end

% Get df for 1st order  contrasts
%	if ((NF == 3 | NF == 4) & Contr.ord1.tot > 0),
	if (Contr.ord1.tot > 0),		
		for (i = 1:1:Contr.ord1.tot),
		   Contr.ord1.df(i) = dfdenom(Contr.ord1.cnt(i).idx1);
      end
	end

% Get df for 2nd order contrasts	

   if ((NF == 2) & Contr.ord2.tot > 0),
   for (i = 1:1:Contr.ord2.tot),
		   switch Contr.ord2.cnt(i).idx1   % 1st factor whose level is fixed
		   case 1,
			   switch Contr.ord2.cnt(i).idx2    % 2nd factor whose level is fixed
				   case 2, Contr.ord2.df(i) = dfdenom(3);  % MSAB
				end	
			case 2,
			   fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;
		   end
   end
   end	

   if ((NF == 3) & Contr.ord2.tot > 0),
   for (i = 1:1:Contr.ord2.tot),
		   switch Contr.ord2.cnt(i).idx1   % 1st factor whose level is fixed
		   case 1,
			   switch Contr.ord2.cnt(i).idx2    % 2nd factor whose level is fixed
				   case 2, Contr.ord2.df(i) = dfdenom(4);  % MSAB
					case 3, Contr.ord2.df(i) = dfdenom(5) * (dsgn == 1 | dsgn == 2) + dfdenom(3) * (dsgn == 4);  % MSAC
				end	
			case 2,
			   if (Contr.ord2.cnt(i).idx2 == 3), Contr.ord2.df(i) = dfdenom(6) * (dsgn == 1 | dsgn == 2) + dfdenom(5) * (dsgn == 3 | dsgn == 4);  % MSBC
				else fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;	
				end		
			case 3,
			   fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;
		   end
   end
   end	
	
	if ((NF == 4) & Contr.ord2.tot > 0),
   for (i = 1:1:Contr.ord2.tot),
%		if (dsgn == 3),
		   switch Contr.ord2.cnt(i).idx1
		   case 1,
			   switch Contr.ord2.cnt(i).idx2
				   case 2, Contr.ord2.df(i) = dfdenom(5);  % MSAB
					case 3, Contr.ord2.df(i) = dfdenom(6);  % MSAC
					case 4, Contr.ord2.df(i) = dfdenom(7)*(dsgn == 1 | dsgn == 2) + dfdenom(4) * (dsgn == 4);  % MSAD
				end	
			case 2,
			   switch Contr.ord2.cnt(i).idx2
				   case 3, Contr.ord2.df(i) = dfdenom(8) * (dsgn == 1 | dsgn == 2) + dfdenom(7) * (dsgn == 3 | dsgn == 4 | dsgn == 5);  % MSBC
					case 4, Contr.ord2.df(i) = dfdenom(9) * (dsgn == 1 | dsgn == 2) + dfdenom(8) * (dsgn == 3 | dsgn == 4 | dsgn == 5);  % Less likely occur: MSBD	
				end		
			case 3,   % Less likely occur
			   switch Contr.ord2.cnt(i).idx2MSE
					case 4, Contr.ord2.df(i) = dfdenom(10) * (dsgn == 1 | dsgn == 2) + dfdenom(9)*(dsgn == 3 | dsgn == 4 | dsgn == 5); % Less likely occur: MSCD		
		      end
		   end
%		end
   end
   end	
	
	if ((NF == 5) & Contr.ord2.tot > 0),   % Refer SumsOfSquares.m for details
   for (i = 1:1:Contr.ord2.tot),
		   switch Contr.ord2.cnt(i).idx1
		   case 1,
			   switch Contr.ord2.cnt(i).idx2
				   case 2, Contr.ord2.df(i) = dfdenom(6);  % MSAB
					case 3, Contr.ord2.df(i) = dfdenom(7);  % MSAC
					case 4, Contr.ord2.df(i) = dfdenom(8);  % MSAD
					case 5, Contr.ord2.df(i) = dfdenom(9)*(dsgn == 1 | dsgn == 2);  % MSAD
				end	
			case 2,
			   switch Contr.ord2.cnt(i).idx2
				   case 3, Contr.ord2.df(i) = dfdenom(10) * (dsgn == 1 | dsgn == 2) + dfdenom(9) * (dsgn == 3 | dsgn == 4);  % MSBC
					case 4, Contr.ord2.df(i) = dfdenom(11) * (dsgn == 1 | dsgn == 2) + dfdenom(10) * (dsgn == 3 | dsgn == 4);  % Less likely occur: MSBD
					case 5, Contr.ord2.df(i) = dfdenom(12) * (dsgn == 1 | dsgn == 2) + dfdenom(11) * (dsgn == 3 | dsgn == 4);  % Less likely occur: MSBE
				end		
			case 3,   % Less likely occur
			   switch Contr.ord2.cnt(i).idx2
					case 4, Contr.ord2.df(i) = dfdenom(13) * (dsgn == 1 | dsgn == 2) + dfdenom(12)*(dsgn == 3 | dsgn == 4); % Less likely occur: MSCD	
					case 5, Contr.ord2.df(i) = dfdenom(14) * (dsgn == 1 | dsgn == 2) + dfdenom(13)*(dsgn == 3 | dsgn == 4); % Less likely occur: MSCE		
		      end
			case 4,   % Less likely occur
			   switch Contr.ord2.cnt(i).idx2
					case 5, Contr.ord2.df(i) = dfdenom(15) * (dsgn == 1 | dsgn == 2) + dfdenom(14)*(dsgn == 3 | dsgn == 4); % Less likely occur: MSDE		
		      end	
		   end  % switch Contr.ord2.cnt(i).idx1
%		end
   end
   end	% if ((NF == 5) & Contr.ord2.tot > 0)
	
	

% Get df for 3rd order  contrasts

	if ((NF == 3) & Contr.ord3.tot > 0),
   for (i = 1:1:Contr.ord3.tot),
		   switch Contr.ord3.cnt(i).idx1
		   case 1,
			   switch Contr.ord3.cnt(i).idx2
				   case 2,
					   if (Contr.ord3.cnt(i).idx3 == 3), Contr.ord3.df(i) = dfdenom(7)*(dsgn == 1 | dsgn == 2);  % MSABC
							else fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;
						end	
					case 3, fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;
				end	
			case 2, fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;	
			case 3, fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;	
		   end
   end
   end	

	if ((NF == 4) & Contr.ord3.tot > 0),
   for (i = 1:1:Contr.ord3.tot),
%		if (dsgn == 3),
		   switch Contr.ord3.cnt(i).idx1
		   case 1,
			   switch Contr.ord3.cnt(i).idx2
				   case 2,
					   switch Contr.ord3.cnt(i).idx3
						   %case 3, Contr.ord3.df(i) = dfdenom(11)*(dsgn == 1 | dsgn == 2) + dfdenom(10)* (dsgn == 3 | dsgn == 4 | dsgn == 5);  % MSABC
							case 3, Contr.ord3.df(i) = dfdenom(11)*(dsgn == 1 | dsgn == 2) + dfdenom(10)* (dsgn == 3 | dsgn == 4) + dfdenom(9)* (dsgn == 5);  % MSABC
							case 4, Contr.ord3.df(i) = dfdenom(12)*(dsgn == 1 | dsgn == 2);  % MSABD not exist for (dsgn == 3 | dsgn == 4 | dsgn == 5)			
						end	
					case 3,
					   if (Contr.ord3.cnt(i).idx3 == 4), Contr.ord3.df(i) = dfdenom(13)*(dsgn == 1 | dsgn == 2);  % MSACD
						else fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;
					   end
				end	
			case 2,
			   switch Contr.ord3.cnt(i).idx2
				   case 3,
					   if (Contr.ord3.cnt(i).idx3 == 4), Contr.ord3.df(i) = dfdenom(14)*(dsgn == 1 | dsgn == 2) + dfdenom(11)* (dsgn == 3 | dsgn == 4 | dsgn == 5);   % MSBCD
						else fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;end
				   case 4,
					   fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;	
				end
			case 3,
			   fprintf('\nSomething is wrong in the contrast coding!\n');
		      fprintf(2,'Halted: Ctrl+c to exit'); pause;	
			case 4,
			   fprintf('\nSomething is wrong in the contrast coding!\n');
		      fprintf(2,'Halted: Ctrl+c to exit'); pause;			
		   end
%		end
   end
   end
	
	if ((NF == 5) & Contr.ord3.tot > 0),
   for (i = 1:1:Contr.ord3.tot),
		   switch Contr.ord3.cnt(i).idx1
		   case 1,
			   switch Contr.ord3.cnt(i).idx2
				   case 2,
					   switch Contr.ord3.cnt(i).idx3
						   case 3, Contr.ord3.df(i) = dfdenom(16)*(dsgn == 1 | dsgn == 2) + dfdenom(15)* (dsgn == 3 | dsgn == 4);  % MSABC
							case 4, Contr.ord3.df(i) = dfdenom(17)*(dsgn == 1 | dsgn == 2) + dfdenom(16)* (dsgn == 3 | dsgn == 4);  % MSABD
							case 5, Contr.ord3.df(i) = dfdenom(18)*(dsgn == 1 | dsgn == 2);  % MSABE			
						end	
					case 3,
					   switch Contr.ord3.cnt(i).idx3
						   case 4, Contr.ord3.df(i) = dfdenom(19)*(dsgn == 1 | dsgn == 2) + dfdenom(17)* (dsgn == 3 | dsgn == 4);  % MSACD
							case 5, Contr.ord3.df(i) = dfdenom(20)*(dsgn == 1 | dsgn == 2);  % MSACE
						end
					case 4,	
						if (Contr.ord3.cnt(i).idx3 == 5), Contr.ord3.df(i) = dfdenom(21)*(dsgn == 1 | dsgn == 2);  % MSADE
						else fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;
					   end
				end	% switch Contr.ord3.cnt(i).idx2
			case 2,
			   switch Contr.ord3.cnt(i).idx2
				   case 3,
					   switch Contr.ord3.cnt(i).idx3
						   case 4, Contr.ord3.df(i) = dfdenom(22)*(dsgn == 1 | dsgn == 2) + dfdenom(18)* (dsgn == 3 | dsgn == 4);   % MSBCD
						   case 5, Contr.ord3.df(i) = dfdenom(23)*(dsgn == 1 | dsgn == 2) + dfdenom(19)* (dsgn == 3 | dsgn == 4);   % MSBCE
						end
					case 4,	
						if (Contr.ord3.cnt(i).idx3 == 5), Contr.ord3.df(i) = dfdenom(24)*(dsgn == 1 | dsgn == 2) + dfdenom(20)* (dsgn == 3 | dsgn == 4);   % MSBDE
						else fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;end
				   case 4,
					   fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;	
				end   % switch Contr.ord3.cnt(i).idx2
			case 3,
			   if (Contr.ord3.cnt(i).idx2 == 4 & Contr.ord3.cnt(i).idx3 == 5),
				if (dsgn == 1 | dsgn == 2),
				   Contr.ord3.df(i) = dfdenom(25);
				elseif (dsgn == 3),
				 	dfdenom(21);   % MSCDE
				end	
				else fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;end	
					
			case 4,
			   fprintf('\nSomething is wrong in the contrast coding!\n');
		      fprintf(2,'Halted: Ctrl+c to exit'); pause;			
		end % switch Contr.ord3.cnt(i).idx1
   end
   end
	
	if ((NF == 5) & Contr.ord4.tot > 0),		
	   for (i = 1:1:Contr.ord4.tot),
		   switch Contr.ord4.cnt(i).idx1
			   case 1,
			   switch Contr.ord4.cnt(i).idx2		
					case 2,
		         switch Contr.ord4.cnt(i).idx3
			         case 3,
						switch Contr.ord4.cnt(i).idx4
						   case 4,
							if (dsgn == 1 | dsgn == 2),
							   Contr.ord4.df(i) = dfdenom(26);
							elseif (dsgn == 3  | dsgn == 4),
							   Contr.ord4.df(i) = dfdenom(22);  % MSABCD
							end	
							case 5,
							if (dsgn == 1 | dsgn == 2),
							   Contr.ord4.df(i) = dfdenom(27);
							end   % MSABCE
						end  % switch Contr.ord4.cnt(i).idx4
						case 4,
						if (Contr.ord4.cnt(i).idx4 == 5),
						if (dsgn == 1 | dsgn == 2),
						   Contr.ord4.df(i) = dfdenom(28);   % MSABDE
						end	
						else fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;	
	               end % if (Contr.ord4.cnt(i).idx4 == 5)
						case 5, fprintf('\nSomething is wrong in the contrast coding!\n');
	                  fprintf(2,'Halted: Ctrl+c to exit'); pause;
					end  % switch Contr.ord4.cnt(i).idx3
					case 3,
					switch Contr.ord4.cnt(i).idx3
					   case 4,
						if (Contr.ord4.cnt(i).idx4 == 5),
						if (dsgn == 1 | dsgn == 2),
						   Contr.ord4.df(i) = dfdenom(29);   % MSACDE
						end
						else fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;
						end
						case 5, fprintf('\nSomething is wrong in the contrast coding!\n');
	                  fprintf(2,'Halted: Ctrl+c to exit'); pause;
					end  % switch Contr.ord4.cnt(i).idx3
					case 4, fprintf('\nSomething is wrong in the contrast coding!\n');
	               fprintf(2,'Halted: Ctrl+c to exit'); pause;
	         end % switch Contr.ord4.cnt(i).idx2
				MSE
				case 2,
				switch Contr.ord4.cnt(i).idx2
				   case 3,
					switch Contr.ord4.cnt(i).idx3
					   case 4,
						if (Contr.ord4.cnt(i).idx4 == 5),
						if (dsgn == 1 | dsgn == 2),
						   Contr.ord4.df(i) = dfdenom(30);   % MSBCDE
						end	
						else fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;	
	               end
						case 5, fprintf('\nSomething is wrong in the contrast coding!\n');
	                  fprintf(2,'Halted: Ctrl+c to exit'); pause;
					end % switch Contr.ord4.cnt(i).idx3
					case 4,
					fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;
				end % switch Contr.ord4.cnt(i).idx2
				case {3,4,5},
				   fprintf('\nSomething is wrong in the contrast coding!\n');
	            fprintf(2,'Halted: Ctrl+c to exit'); pause;
				
			end	% switch Contr.ord4.cnt(i).idx1			
		end  % for (i = 1:1:Contr.ord4.tot)
	end	% 	if ((NF == 5) & Contr.ord4.tot > 0)	
	
%	if (cov.do == 1), covoutdf = dfe;  end % df of MSE

%Set up options for writing out the results
	
   Opt.Prefix = sprintf('%s', OutFN);
   Opt.NoCheck = 0;
   Opt.Scale = 1;
   if (data_type == 0), Opt.View = sprintf('+%s', format);
	elseif (data_type == 1) Opt.View = '.1D.dset';
   elseif (data_type == 2) Opt.View = '.niml.dset'; end

   Info.TYPESTRING = '3DIM_HEAD_FUNC';
   %Info.DATASET_RANK(2) = 2*N_Brik;    % Number of subbriks
%	if (cov.do == 1),
%	   Info.DATASET_RANK(2) = 2*(N_Brik + Contr.ord1.tot + Contr.ord2.tot + Contr.ord3.tot +1);
%   else Info.DATASET_RANK(2) = 2*(N_Brik + Contr.ord1.tot + Contr.ord2.tot + Contr.ord3.tot);    % Number of subbriks
%	end
   Info.DATASET_RANK(2) = 2*(N_Brik + Contr.ord1.tot + Contr.ord2.tot + Contr.ord3.tot + Contr.ord4.tot);    % Number of subbriks

   Info.DATASET_DIMENSIONS = [D1 D2 D3 0 0];
   Info.BRICK_STATS = [];
%   Info.BRICK_TYPES = [];
   Info.BRICK_FLOAT_FACS = [];
   Info.BRICK_STATAUX = [];   %initialization
   Info.BRICK_TYPES = [];   %initialization: generate results in float with 3
   Info.TypeBytes = 4; % Force the output to be float. Why?

% For all terms
   for (i = 1:1:N_Brik),
      Info.BRICK_STATAUX = [Info.BRICK_STATAUX 2*i-1 4 2 dfterm_new(i) dfdenom(i)];   % 4 is for F stat defined in 3ddata.h
	   Info.BRICK_TYPES = [Info.BRICK_TYPES 3 3];  % 3 for float; 1 for short
   end

% For 1st oder contrasts
   if (Contr.ord1.tot>0),
   for (i = 1:1:Contr.ord1.tot),
      Info.BRICK_STATAUX = [Info.BRICK_STATAUX 2*(i+N_Brik)-1 3 1 Contr.ord1.df(i)];   % 3 is for t stat defined in 3ddata.h
   	Info.BRICK_TYPES = [Info.BRICK_TYPES 3 3];  % 3 for float; 1 for short
   end
	end
	
% For 2nd oder contrasts
   if (Contr.ord2.tot>0),
   for (i = 1:1:Contr.ord2.tot),
      Info.BRICK_STATAUX = [Info.BRICK_STATAUX 2*(i+N_Brik+Contr.ord1.tot)-1 3 1 Contr.ord2.df(i)];   % 3 is for t stat defined in 3ddata.h
   	Info.BRICK_TYPES = [Info.BRICK_TYPES 3 3];  % 3 for float; 1 for short
   end
	end	
	
% For 3rd oder contrasts
   if (Contr.ord3.tot>0),
   for (i = 1:1:Contr.ord3.tot),
      Info.BRICK_STATAUX = [Info.BRICK_STATAUX 2*(i+N_Brik+Contr.ord1.tot+Contr.ord2.tot)-1 3 1 Contr.ord3.df(i)];   % 3 is for t stat defined in 3ddata.h
   	Info.BRICK_TYPES = [Info.BRICK_TYPES 3 3];  % 3 for float; 1 for short
   end
	end
	
% For 3rd oder contrasts
   if (Contr.ord4.tot>0),
   for (i = 1:1:Contr.ord4.tot),
      Info.BRICK_STATAUX = [Info.BRICK_STATAUX 2*(i+N_Brik+Contr.ord1.tot+Contr.ord2.tot+Contr.ord3.tot)-1 3 1 Contr.ord4.df(i)];   % 3 is for t stat defined in 3ddata.h
   	Info.BRICK_TYPES = [Info.BRICK_TYPES 3 3];  % 3 for float; 1 for short
   end
	end

% ANOCVA	
%	if (cov.do == 1),
%	   Info.BRICK_STATAUX = [Info.BRICK_STATAUX 2*(i+N_Brik+Contr.ord1.tot+Contr.ord2.tot+Contr.ord3.tot)-1 3 1 covoutdf];   % 3 is for t stat defined in 3ddata.h
%   	Info.BRICK_TYPES = [Info.BRICK_TYPES 3 3];  % 3 for float; 1 for short
%	end
	
% For all terms
   Info.BRICK_LABS =  [];
   for (i = 1:1:N_Brik),  %create labels for all terms
      Info.BRICK_LABS =  [Info.BRICK_LABS cell2mat(tnames_new(i)) '~'];   % Label for intensity
	   Info.BRICK_LABS =  [Info.BRICK_LABS cell2mat(tnames_new(i)) ' ' 'F' '~'];   % Label for F
   end

% For 1st order contrasts
   if (Contr.ord1.tot>0),
   for (i = 1:1:Contr.ord1.tot),  %create labels for contrasts
      Info.BRICK_LABS =  [Info.BRICK_LABS Contr.ord1.label(i).nm '~'];   % Label for intensity
	   Info.BRICK_LABS =  [Info.BRICK_LABS Contr.ord1.label(i).nm ' ' 't' '~'];   % Label for t
   end
	end

% For 2nd order contrasts
   if (Contr.ord2.tot>0),
   for (i = 1:1:Contr.ord2.tot),  %create labels for contrasts
      Info.BRICK_LABS =  [Info.BRICK_LABS Contr.ord2.label(i).nm '~'];   % Label for intensity
	   Info.BRICK_LABS =  [Info.BRICK_LABS Contr.ord2.label(i).nm ' ' 't' '~'];   % Label for t
   end
	end

% For 3rd order contrasts
   if (Contr.ord3.tot>0),
   for (i = 1:1:Contr.ord3.tot),  %create labels for contrasts
      Info.BRICK_LABS =  [Info.BRICK_LABS Contr.ord3.label(i).nm '~'];   % Label for intensity
	   Info.BRICK_LABS =  [Info.BRICK_LABS Contr.ord3.label(i).nm ' ' 't' '~'];   % Label for t
   end
	end
	
% For 4th order contrasts
   if (Contr.ord4.tot>0),
   for (i = 1:1:Contr.ord4.tot),  %create labels for contrasts
      Info.BRICK_LABS =  [Info.BRICK_LABS Contr.ord4.label(i).nm '~'];   % Label for intensity
	   Info.BRICK_LABS =  [Info.BRICK_LABS Contr.ord4.label(i).nm ' ' 't' '~'];   % Label for t
   end
	end	
	
	% ANOCVA	
%	if (cov.do == 1),
%	   Info.BRICK_LABS =  [Info.BRICK_LABS cov.label '~'];   % Label for beta
%	   Info.BRICK_LABS =  [Info.BRICK_LABS cov.label ' ' 't' '~'];   % Label for t
%	end

   Info.SCENE_DATA(2) = 11;   % 11 means Bucket Func type, defined in AFNI doc README.attributes
   Info.SCENE_DATA(3) = 1;    % 1 means 3DIM_HEAD_FUNC, defined in AFNI doc README.attributes
	
   if (isfield(Info,{'TAXIS_NUMS'})),
      Info = rmfield(Info, {'TAXIS_NUMS'});
   end
   if (isfield(Info,{'TAXIS_FLOATS'})),
      Info = rmfield(Info, {'TAXIS_FLOATS'});
   end
   if (isfield(Info,{'TAXIS_OFFSETS'})),
      Info = rmfield(Info, {'TAXIS_OFFSETS'});
   end
   if (isfield(Info,{'BRICK_FLOAT_FACS'})),
      Info = rmfield(Info, {'BRICK_FLOAT_FACS'});
   end


   Opt.Frames = [];  %Because it might have been set as the frame list in the case of input files with multiple subbriks during loading

% Write output: reshape because the 3rd dimension is supposed to be Z in the normal situation.
   OptW = Opt;
   OptW.Scale = 0;
   if (data_type == 0),
      [err2, ErrMessage, NewInfo] = WriteBrik(reshape(M, D1, D2, 1, Info.DATASET_RANK(2)), Info, OptW);
	elseif (data_type == 1), % Collapse the 2nd dimsion, which is 1 for surface data.
      [err2, ErrMessage, NewInfo] = WriteBrik(squeeze(reshape(M, D1, D2, 1, Info.DATASET_RANK(2))), Info, OptW);
	elseif (data_type == 2),
      [sstt, uid] = system('3dnewid -fun');
      ftmp = sprintf('__GroupAna_tmp%d_%s.mat', sn, uid);
      Mo = squeeze(reshape(M, D1, D2, 1, Info.DATASET_RANK(2)));
      NIo = NI((1+(sn-1)*Opt.SliceSize_1D):min([length(NI) sn*Opt.SliceSize_1D]));
      save(ftmp, 'Mo', 'NIo'); clear('Mo'); clear('NIo');
      flist(sn) = cellstr(ftmp);
   end	
	fprintf(1, 'done in %f seconds\n', toc);	
end  % end of the big loop of running ANOVA and writing up: One slice a time
if (data_type == 2),
   clear('LC'); clear ('M');
   for (it=1:1:length(flist)),
      F = load(char(flist(it)));
      if (it==1),
         M = F.Mo;
      else
         M = [M ; F.Mo];
      end
      clear('F');
      delete(char(flist(it)));
   end
   Info.FileFormat='NIML';
   Info.NodeIndices = NI;
   OptW.Slices = [];
   [err2, ErrMessage, NewInfo] = WriteBrik(M, Info, OptW);
end
fprintf(1, '\nCongratulations, job is successfully done!!! Total runtime: %f minutes...', etime(clock,t0)/60);	
fprintf(1, '\nOutput files are %s%s.*\n\n', Opt.Prefix, Opt.View);	
err = 0;
return;
