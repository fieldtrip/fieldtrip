function [err,S] = HistoryTrace ( Opt)
%
%   [err,S] = HistoryTrace (Opt)
%
%Purpose:
%   returns a signature/history sting, useful for keeping a log
%
%   The trace string is formed as such
%      For example if script MainScript calls at line 12 
%          function F1 which at line 6 calls F2 wich at line 192 calls 
%          F3 which at line 56 calls HistoryTrace then the trace string is:
%         MainScript(ln12)->F1(ln6)->F2(ln192)->F3(ln56)
%   
%Input Parameters:
%   Opt is an optional options structure with the following optional fields
%      .NoPath : (0/[1]) strips the paths of the function names in the call trace
%      .PerSig : a string to be used instead of ['Ziad S. Saad LBC/NIMH/NIH']
%      .AFNI_Format : ([0]/1) formats S to fit AFNI's History field
%      .MaxTraceLevel : ([3]) Maximum number of functions to display
%          from the top of the stack. In the example for the trace string, 
%          F3 will not be shown. If .MaxTraceLevel is negative then counting
%          is done from the bottom of the stack. In the example, MainScript 
%          would not show up.
%
%Output Parameters:
%   err : 0 No Problem
%       : 1 Mucho Problems
%   S the string wth PerSig, function call trace with the line numbers in 
%     parenthesis, machine name and current working directory and date
%   
%   
%      
%Key Terms:
%   
%More Info :
%   plotsign2
%   tross_Encode_String
%   
%
%     Author : Ziad Saad
%     Date : Tue Apr 10 17:29:11 PDT 2001
%     LBC/NIMH/ National Institutes of Health, Bethesda Maryland


%Define the function name for easy referencing
FuncName = 'HistoryTrace';

%Debug Flag
DBG = 1;

%initailize return variables
err = 1;

if (nargin == 0), Opt = []; end

if (~isfield(Opt,'NoPath') | isempty(Opt.NoPath)), Opt.NoPath = 1; end
if (~isfield(Opt,'AFNI_Format') | isempty(Opt.AFNI_Format)), Opt.AFNI_Format = 0; end
if (~isfield(Opt,'MaxTraceLevel') | isempty(Opt.MaxTraceLevel)), Opt.MaxTraceLevel = 3; end

if (~isfield(Opt,'PerSig') | isempty(Opt.PerSig)), 
	Opt.PerSig = sprintf ('Ziad Saad LBC/NIMH/NIH'); 
end

	[ST,I] = dbstack; 
	N_ST = length(ST);
	if (N_ST > 1),
		if (Opt.AFNI_Format),	stk = sprintf('\n\t'); else stk = sprintf('\n');end
		if (Opt.MaxTraceLevel > 0),
			strt = N_ST; stp = max([(N_ST - Opt.MaxTraceLevel + 1) 2]);
		else
			strt = min([(1-Opt.MaxTraceLevel) N_ST]); stp = 2;
		end
		for (i=strt:-1:stp), 
			if (Opt.NoPath), [err,PathString,ST(i).name] = GetPath (ST(i).name); end
			stk = sprintf('%s%s(ln%d)->', stk, ST(i).name, ST(i).line);
		end
		stk = stk(1:1:length(stk)-2);
	else
		stk = '';
	end		
	[tmp, smach] = unix('hostname');
	%remove this annoying tset message (some bug ....)
	[err, snl, Nlines] = GetNextLine(smach, 2);
	if (Nlines >= 2),
		[err, smach] = GetNextLine(smach,Nlines);
	end 
	if (tmp), 
		smach = sprintf('Dunno');	
		else
		smach = zdeblank(smach); 
	end
	c=datevec(now);

	if (Opt.AFNI_Format),
		S = sprintf ('[%s@%s: %s %s:%s:%s]\n\t%s%s',...
				 Opt.PerSig,smach,...
				 date,...
				 pad_strn(num2str(c(4)), '0', 2, 1),...
				 pad_strn(num2str(c(5)), '0', 2, 1),...
				 pad_strn(num2str(round(c(6))), '0', 2, 1),...
				 pwd, ...
				 stk);

		[err, S] = tross_Encode_String(S);
	else
		S = sprintf ('%s: %s %s:%s:%s\n%s (%s)%s',...
					 Opt.PerSig,...
					 date,...
					 pad_strn(num2str(c(4)), '0', 2, 1),...
					 pad_strn(num2str(c(5)), '0', 2, 1),...
					 pad_strn(num2str(round(c(6))), '0', 2, 1),...
					 smach, pwd, ...
					 stk);
	end
	
err = 0;
return;

