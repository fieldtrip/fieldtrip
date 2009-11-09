function data = spm_chi2_plot(action,arg1,arg2,arg3,arg4)
% Display a plot showing convergence of an optimization routine.
% FORMAT spm_chi2_plot('Init',title,xlabel,ylabel)
% Initialises the plot in the 'Interactive' window.
%
% FORMAT spm_chi2_plot('Set',value)
% Updates the plot.
%
% FORMAT spm_chi2_plot('Clear')
% Clears the 'Interactive' window.
%
%-----------------------------------------------------------------------
% @(#)spm_chi2_plot.m	1.2 John Ashburner 97/09/11

global pb_pointer pb_name ax
%-----------------------------------------------------------------------
if (nargin == 0)
	spm_chi2_plot('Init');
else
	% initialize
	%---------------------------------------------------------------
	if (strcmp(lower(action),'init'))
		if (nargin<4)
			arg3 = 'Iteration #';
			if (nargin<3)
				arg2 = 'Chi-squared';
				if (nargin<2)
					arg1 = 'Optimizing';
				end
			end
		end
		fg = spm_figure('FindWin','Interactive');
		if ~isempty(fg)
			pb_pointer = get(fg,'Pointer');
			pb_name    = get(fg,'Name');
			spm_chi2_plot('Clear');
			set(fg,'Pointer','watch');
			set(fg,'Name',pb_name);
			ax = axes('Position', [0.15 0.1 0.8 0.75],...
				'Box', 'on','Parent',fg);
			lab = get(ax,'Xlabel');
			set(lab,'string',arg3,'FontSize',10);
			lab = get(ax,'Ylabel');
			set(lab,'string',arg2,'FontSize',10);
			lab = get(ax,'Title');
			set(lab,'string',arg1);
			line('Xdata',[], 'Ydata',[],...
				'LineWidth',2,'Tag','Chi2Plot','Parent',ax);
			drawnow;
		end

	% reset
	%---------------------------------------------------------------
	elseif (strcmp(lower(action),'set'))
		if (nargin<2)
			arg1 = 0;
		end
		F = spm_figure('FindWin','Interactive');
		br = findobj(F,'Tag','Chi2Plot');
		if (~isempty(br))
			xd = get(br,'Xdata');
			yd = [get(br,'Ydata') arg1];
			xd = [xd (length(xd)+1)];
			set(br,'Ydata',yd,'Xdata',xd);
			drawnow;
		end

	% clear
	%---------------------------------------------------------------
	elseif (strcmp(lower(action),'clear'))
		fg = spm_figure('FindWin','Interactive');
		spm_figure('Clear',fg);
		set(fg,'Pointer',pb_pointer);
		set(fg,'Name',pb_name);
		drawnow;
	elseif (strcmp(lower(action),'data'))
		F = spm_figure('FindWin','Interactive');
		br = findobj(F,'Tag','Chi2Plot');
		if (~isempty(br))
			data = get(br,'Ydata');
		else
			data = [];
		end
	end
end
