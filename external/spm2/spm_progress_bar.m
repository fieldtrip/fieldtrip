function spm_progress_bar(action,arg1,arg2,arg3,arg4)
% Display a 'Progress Bar'
% FORMAT spm_progress_bar('Init',height,xlabel,ylabel,flgs)
% Initialises the bar in the 'Interactive' window.
% If flgs contains a 't', then use tex interpreter for labels.
%
% FORMAT spm_progress_bar('Set',value)
% Sets the height of the bar itself.
%
% FORMAT spm_progress_bar('Clear')
% Clears the 'Interactive' window.
%
%-----------------------------------------------------------------------
% @(#)spm_progress_bar.m	2.2 John Ashburner 02/08/01

global pb_pointer pb_name ax
%-----------------------------------------------------------------------
if nargin == 0,
	spm_progress_bar('Init');
else,
	action = lower(action);
	% initialize
	%---------------------------------------------------------------
	if strcmp(action,'init'),
		if nargin<5,
			arg4 = ' ';
			if nargin<4,
				arg3 = '';
				if nargin<3,
					arg2 = 'Computing';
					if nargin<2,
						arg1 = 1;
					end;
				end;
			end;
		end;
		if any(arg4=='t'), interp = 'tex'; else, interp = 'none'; end;
		fg = spm_figure('FindWin','Interactive');
		if ~isempty(fg),
			pb_pointer = get(fg,'Pointer');
			pb_name    = get(fg,'Name');
			spm_progress_bar('Clear');
			set(fg,'Pointer','watch');
			set(fg,'Name',pb_name);
			ax = axes('Position', [0.45 0.2 0.05 0.6],...
				'XTick',[],...
				'Xlim', [0 1],...
				'Ylim', [0 max([arg1 eps])],...
				'Box', 'on',...
				'Parent',fg);
			lab = get(ax,'Xlabel');
			set(lab,'string',arg2,'FontSize',10,'Interpreter',interp);
			lab = get(ax,'Ylabel');
			set(lab,'string',arg3,'FontSize',10,'Interpreter',interp);
			lab = get(ax,'Title');
			set(lab,'string','0% Complete','Interpreter',interp);
			tim = clock;
			tim = tim(4:6);
			str = sprintf('Began %2.0f:%2.0f:%2.0f',tim(1),tim(2),tim(3));
			t1=text(2,arg1/2,0,str,'FontSize',10,'Parent',ax);
			set(t1,'Tag','StartTime');
			line('Xdata',[0.5 0.5], 'Ydata',[0 0],...
				'LineWidth',8, 'Color', [1 0 0],'Tag','ProgressBar',...
				'Parent',ax);
			drawnow;
		end;
	% reset
	%---------------------------------------------------------------
	elseif strcmp(action,'set'),
		if nargin<2,
			arg1 = 0;
		end;
		F = spm_figure('FindWin','Interactive');
		if ~isempty(F),
			br = findobj(F,'Tag','ProgressBar');
			if (~isempty(br))
				set(br,'Ydata',[0 arg1]);
				lim = get(get(br,'Parent'),'Ylim');lim=lim(2);
				lab = get(ax,'Title'); 
				set(lab,'string',sprintf('%.0f%% Complete',100*arg1/lim)); 
				drawnow;
			end;
		end;
	% clear
	%---------------------------------------------------------------
	elseif strcmp(action,'clear'),
		fg = spm_figure('FindWin','Interactive');
		if ~isempty(fg),
			spm_figure('Clear',fg);
			if ~isempty(pb_pointer),
				set(fg,'Pointer',pb_pointer);
				set(fg,'Name',pb_name);
			end;
			drawnow;
		end;
	end;
end;
