function spm_plot_convergence(action,varargin)
% Display a plot showing convergence of an optimisation routine.
% FORMAT spm_plot_convergence('Init',title,ylabel,xlabel)
% Initialise the plot in the 'Interactive' window.
%
% FORMAT spm_plot_convergence('Set',value)
% Update the plot.
%
% FORMAT spm_plot_convergence('Clear')
% Clear the 'Interactive' window.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_plot_convergence.m 4147 2010-12-24 13:50:06Z guillaume $

if ~nargin, action = 'Init'; end

% Find the interactive window and exit if not
%--------------------------------------------------------------------------
Finter = spm_figure('FindWin','Interactive');
if isempty(Finter), return; end

switch lower(action)
    % Initialise
    %----------------------------------------------------------------------
    case 'init'
        if nargin > 1, arg1 = varargin{1}; else arg1 = 'Optimising';  end
        if nargin > 2, arg2 = varargin{2}; else arg2 = 'Chi-squared'; end
        if nargin > 3, arg3 = varargin{3}; else arg3 = 'Iteration';   end
        pb = struct('pointer',get(Finter,'Pointer'),...
                    'name',   get(Finter,'Name'),...
                    'buffer', get(Finter,'DoubleBuffer'));
        spm_plot_convergence('Clear');
        set(Finter,'Pointer','watch');
        set(Finter,'Name',pb.name);
        set(Finter,'DoubleBuffer','on');
        pb.ax = axes('Position',[0.15 0.1 0.8 0.75],...
                     'Box',     'on',...
                     'Parent',  Finter);
        set(get(pb.ax,'Xlabel'), 'string',arg3, 'FontSize',10);
        set(get(pb.ax,'Ylabel'), 'string',arg2, 'FontSize',10);
        set(get(pb.ax,'Title'),  'string',arg1);
        l = line('Xdata',     [],...
                 'Ydata',     [],...
                 'LineWidth', 2,...
                 'Tag',       'PlotConvOptim',...
                 'Parent',    pb.ax);
        set(l,'UserData',pb);
        
    % Set
    %----------------------------------------------------------------------
    case 'set'
        if nargin > 1, arg1 = varargin{1}; else arg1 = 0; end
        br = findobj(Finter,'Tag','PlotConvOptim');
        if ~isempty(br)
            xd = get(br,'Xdata');
            yd = [get(br,'Ydata') arg1];
            xd = [xd (length(xd)+1)];
            set(br,'Ydata',yd, 'Xdata',xd);
        end
        
    % Clear
    %----------------------------------------------------------------------
    case 'clear'
        pb = get(findobj(Finter,'Tag','PlotConvOptim'),'UserData');
        spm_figure('Clear',Finter);
        if isstruct(pb)
            set(Finter,'Pointer',     pb.pointer);
            set(Finter,'Name',        pb.name);
            set(Finter,'DoubleBuffer',pb.buffer);
        end
        
    % Error
    %----------------------------------------------------------------------
    otherwise
        error('Unknown action string');
end

drawnow;
