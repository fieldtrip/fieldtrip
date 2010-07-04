function [err] = m3dReorder (Input, Prefix, Mapfile, Opt)
%
%   [err] = m3dReorder (Input, Prefix, Mapfile, Opt)
%
%Purpose:
%   Reorders a time series data set (3D+time) a la AFNI plugin Reorder   
%   
%   
%Input Parameters:
%
%   Input: Input of 3d+time brick
%   Prefix: prefix of output data set
%   Mapfile: Name of map file
%   Opt is the options structure with the following fields
%   	.Dup : [Col]/Ave Collate(default) or Average the multiple
%             instances in the map file. 
%     .Detrend : 0/1/[2] Linear trend removal. 0 for none, 
%                1 for mean only, 2 for linear trend (default)  
%     .Verbose : [0]/1 verbosity of function ...
%
%Output Parameters:
%   err : 0 No Problem
%       : 1 Mucho Problems
%   
%   
%      
%Key Terms:
%   
%More Info :
%   
%   help button in Reorder plugin
%   
%
%     Author : Ziad Saad
%     Date : Tue Sep 18 10:39:36 EDT 2001
%     LBC/NIMH/ National Institutes of Health, Bethesda Maryland

%Define the function name for easy referencing
FuncName = 'm3dReorder';

%Debug Flag
DBG = 1;

%initailize return variables
err = 1;

%check for options	
	if (nargin < 4), Opt.Dup = 'Col'; Opt.Detrend = 2; Opt.Verbose = 0; end
	
	if (~isfield(Opt,'Dup') | isempty(Opt.Dup)), Opt.Dup = 'Col'; end
	if (~isfield(Opt,'Detrend') | isempty(Opt.Detrend)), Opt.Detrend = 2; end
	if (~isfield(Opt,'Verbose') | isempty(Opt.Verbose)), Opt.Verbose = 0; end
	
%check for valid parameters
	if (~eq_str(Opt.Dup,'Col') & ~eq_str(Opt.Dup,'Ave')),
		err = ErrEval(FuncName,'Err_Bad value for Opt.Dup');
		return;
	end

%check the prefix
	[status, Prefix] = PrefixStatus(Prefix); 
	if (status == 0), err = ErrEval(FuncName,'Err_Bad Prefix'); return; end
	
%check and read the input file
	if (Opt.Verbose), fprintf(1,'\nLoading & Sorting Mapfile ...'); end
	fid = fopen (Mapfile,'ro');
	sall = fscanf(fid,'%c');
	fclose (fid);
	
%load the file into a cell string
	%search for Newlines
	inl = find (sall == 10);
	N_inl = length(inl);
	%fill in sall to call, ignoring #
	icell = 1;
	ipos = 1;
	vsall = '';
	for (i = 1:1:N_inl),
		iend = inl(i) -1;
		sword = sall(ipos:iend);
		if (sall(ipos) ~= '#' & ~isempty(sword)), %Not a comment 
				%vertically concatenate the strings
				vsall = strvcat(vsall, sword);
		end
		ipos = inl(i)+1;
	end

% sort the results
	[vsall_sort,i_sort] = sortrows (vsall);

%find the ignore points
	ikeep = find (vsall_sort(:,1) ~= '-');
	imap = i_sort(ikeep);
	%vsall_sort(1:20,:)
	%ikeep(1:5)
	%length(imap) 
	
%load the data set
	%make sure that this is a 3D+time
	[err, Info] = BrikInfo(Input);
	if (err), 
		err = ErrEval(FuncName,'Err_error in BrikInfo. Check Input Brick Name'); return; 
	end
	if (Info.SCENE_DATA(2) ~= 2), 
		err = ErrEval(FuncName,sprintf('Err_%s is not an EPI type', Input)); return;
	end
	%make sure length of time series is OK
	if (Info.DATASET_RANK(2) ~= size(vsall,1)),
		err = ErrEval(FuncName,...
		sprintf('Err_Length mismatch between time series (%g) and Mapfile(%g)',...
		 Info.DATASET_RANK(2),size(vsall,1))); return;  
	end
	
	
	%OK, load data
	if (Opt.Verbose), fprintf(1,'\nLoading Brick ...'); end
	OptR.Format = 'vector';
	[err, V, Info, ErrMessage] = BrikLoad (Input, OptR);
	if (err), err = ErrEval(FuncName,'Err_error in BrikLoad'); return; end
	
	%figure(1);clf; subplot (211);  plot (V(1,:), 'b'); 
	
	%detrend ?
	if (Opt.Detrend == 1),
		if (Opt.Verbose), fprintf(1,'\nDetrending ...'); end
		V = detrend(V','constant'); V = V';
	elseif (Opt.Detrend == 2),
		if (Opt.Verbose), fprintf(1,'\nDetrending ...'); end
		V = detrend(V','linear'); V = V';
	elseif (Opt.Detrend ~= 0),
		err = ErrEval(FuncName, 'Err_Bad value for Opt.Detrend'); return;
	end	
	%plot (V(1,:), 'r'); hold on; drawnow
	%pause

%Now the output brick is:
	Vo = V(:,imap);
	%subplot (212); plot (Vo(1,:)); drawnow
	%pause
%update the header
	Info_o = Info;
	Info_o.RootName = '';
	Info_o.IDCODE_STRING = '';
	Info_o.IDCODE_DATE = '';
	Info_o.TAXIS_OFFSETS = [];
	Info_o.BRICK_FLOAT_FACS = [];
	Info_o.BRICK_STATS = [];
	Info_o.BRICK_TYPES = Info.BRICK_TYPES(imap);
	Info_o.DATASET_RANK(2) = length(imap);
	Info_o.TAXIS_NUMS(1:2) = [length(imap) 0]; %remove time offset because it's meaningless when you scramble the data

%Write out results
	if (Opt.Verbose), fprintf(1,'\nWriting Brick ...'); end
	OptW.Scale = 1;
	OptW.Prefix = Prefix;
	OptW.AppendHistory = 1;
	OptW.NoCheck = 1;
	[err, ErrMessage, Info] = WriteBrik (Vo, Info_o, OptW);
	if (err), err = ErrEval(FuncName,'Err_Error writing brick to disk'); return; end

err = 0;
return;

