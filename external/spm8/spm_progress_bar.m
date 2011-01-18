function spm_progress_bar(action,varargin)
% Display a 'Progress Bar' in the 'Interactive' window
% FORMAT spm_progress_bar('Init',height,xlabel,ylabel,flgs)
% Initialises the bar in the 'Interactive' window.
% If flgs contains a 't', then use tex interpreter for labels.
%
% FORMAT spm_progress_bar('Set',value)
% Sets the height of the bar itself.
%
% FORMAT spm_progress_bar('Clear')
% Clears the 'Interactive' window.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_progress_bar.m 2914 2009-03-20 18:30:31Z guillaume $

if ~nargin, action = 'Init'; end

% Find the interactive window and exit if not
%-----------------------------------------------------------------------
Finter = spm_figure('FindWin','Interactive');
if isempty(Finter), return; end

switch lower(action)
    % Initialise
    %-------------------------------------------------------------------
    case 'init'
        error(nargchk(0,5,nargin));
        if nargin > 1, arg1 = varargin{1}; else arg1 = 1;           end
        if nargin > 2, arg2 = varargin{2}; else arg2 = 'Computing'; end
        if nargin > 3, arg3 = varargin{3}; else arg3 = '';          end
        if nargin > 4, arg4 = varargin{4}; else arg4 = ' ';         end
        if any(arg4 == 't'), interp = 'tex'; else interp = 'none';  end
        pb = struct('pointer',get(Finter,'Pointer'),...
                    'name'   ,get(Finter,'Name'),...
                    'buffer', get(Finter,'DoubleBuffer'));
        spm_progress_bar('Clear');
        set(Finter,'Pointer','watch');
        set(Finter,'Name',pb.name);
        set(Finter,'DoubleBuffer','on');
        pb.ax = axes('Position', [0.45 0.2 0.05 0.6],...
                     'XTick',    [],...
                     'Xlim',     [0 1],...
                     'Ylim',     [0 max([arg1 eps])],...
                     'Box',      'on',...
                     'Parent',   Finter);
        lab = get(pb.ax,'Xlabel');
        set(lab,'string',arg2,'FontSize',10,'Interpreter',interp);
        lab = get(pb.ax,'Ylabel');
        set(lab,'string',arg3,'FontSize',10,'Interpreter',interp);
        lab = get(pb.ax,'Title');
        set(lab,'string','0% Complete','Interpreter',interp);
        t = clock;
        str = sprintf('Began %2.0f:%02.0f:%02.0f',t(4),t(5),t(6));
        text(2,arg1/2,0,str,'FontSize',10,'Parent',pb.ax);
        l = line('Xdata',     [0.5 0.5],...
                 'Ydata',     [0 0],...
                 'LineWidth', 8,...
                 'Color',     [1 0 0],...
                 'Tag',       'ProgressBar',...
                 'Parent',    pb.ax);
        set(l,'UserData',pb);
        drawnow;
        
    % Set
    %-------------------------------------------------------------------
    case 'set'
        error(nargchk(1,2,nargin));
        if nargin == 1, value = 0; else  value = varargin{1}; end
        br = findobj(Finter,'Tag','ProgressBar');
        if ~isempty(br)
            pb = get(br,'UserData');
            set(br,'Ydata',[0 value]);
            lim = get(get(br,'Parent'),'Ylim');lim=lim(2);
            lab = get(pb.ax,'Title'); 
            set(lab,'string',sprintf('%.0f%% Complete',100*value/lim)); 
            drawnow;
        end 
    
    % Clear
    %-------------------------------------------------------------------
    case 'clear'
        error(nargchk(1,1,nargin));
        pb = get(findobj(Finter,'Tag','ProgressBar'),'UserData');
        spm_figure('Clear',Finter);
        if isstruct(pb)
            set(Finter,'Pointer',     pb.pointer);
            set(Finter,'Name',        pb.name);
            set(Finter,'DoubleBuffer',pb.buffer);
        end
        drawnow;
    
    % Error
    %-------------------------------------------------------------------
    otherwise
        error('Unknown action string');
end
