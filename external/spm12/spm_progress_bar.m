function spm_progress_bar(action,varargin)
% Display a 'Progress Bar' in the 'Interactive' window
% FORMAT spm_progress_bar('Init',height,xlabel,ylabel,flgs)
% Initialise the bar in the 'Interactive' window.
% If flgs contains a 't', then use tex interpreter for labels.
%
% FORMAT spm_progress_bar('Set',value)
% Set the height of the bar itself.
%
% FORMAT spm_progress_bar('Set','xlabel',xlabel)
% FORMAT spm_progress_bar('Set','ylabel',ylabel)
% Set the progress bar labels.
%
% FORMAT spm_progress_bar('Set','height',height)
% Set the height of the progress bar.
%
% FORMAT spm_progress_bar('Clear')
% Clear the 'Interactive' window.
%__________________________________________________________________________
% Copyright (C) 1996-2015 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_progress_bar.m 6383 2015-03-19 17:20:41Z guillaume $


persistent pbar;

if ~nargin, action = 'Init'; end

switch lower(action)
    
    % Initialise
    %======================================================================
    case 'init'
        Finter = spm_figure('FindWin','Interactive');
        if isempty(Finter), pbar = []; return; end
        if nargin > 1, arg1 = varargin{1}; else arg1 = 1;           end
        if nargin > 2, arg2 = varargin{2}; else arg2 = 'Computing'; end
        if nargin > 3, arg3 = varargin{3}; else arg3 = '';          end
        if nargin > 4, arg4 = varargin{4}; else arg4 = ' ';         end
        if any(arg4 == 't'), interp = 'tex'; else interp = 'none';  end
        pb = struct('pointer',get(Finter,'Pointer'),...
                    'name',   get(Finter,'Name'),...
                    'buffer', get(Finter,'DoubleBuffer'));
        spm_progress_bar('Clear');
        set(Finter,'Pointer','watch');
        set(Finter,'Name',pb.name);
        set(Finter,'DoubleBuffer','on'); % no effect since R2013a
        if ischar(arg2), arg2 = repmat({arg2},1,numel(arg1)); end
        if ischar(arg3), arg3 = repmat({arg3},1,numel(arg1)); end
        
        for i=1:numel(arg1)
            %-
            pb.ax(i) = axes(...
                'Position', [((i/(numel(arg1)+1))-0.05) 0.2 0.05 0.6],...
                'XTick',    [],...
                'Xlim',     [0 1],...
                'Ylim',     [0 max([arg1(i) eps])],...
                'Box',      'on',...
                'Parent',   Finter);
            try, set(pb.ax(i),'ClippingStyle','rectangle'); end
            
            %-XLabel
            lab = get(pb.ax(i),'Xlabel');
            if numel(arg2) < i, arg2{i} = ''; end
            set(lab,'string',arg2{i},'FontSize',10,'Interpreter',interp);
            
            %-YLabel
            lab = get(pb.ax(i),'Ylabel');
            if numel(arg3) < i, arg3{i} = ''; end
            set(lab,'string',arg3{i},'FontSize',10,'Interpreter',interp);
            
            %-Title
            pb.t(i) = get(pb.ax(i),'Title');
            set(pb.t(i),'string','0% Complete','Interpreter',interp);
            
            %-Began...
            t = clock;
            if numel(arg1) == 1, opts = {};
            else opts = {'Rotation',90,'HorizontalAlignment','center'}; end
            str = sprintf('Began %2.0f:%02.0f:%02.0f',t(4),t(5),t(6));
            pb.b(i) = text(2,arg1(i)/2,0,str,...
                'FontSize',10,...
                'Parent',pb.ax(i),...
                opts{:});
            
            %-Progress bar
            pb.l(i) = line(...
                'Xdata',     [0.5 0.5],...
                'Ydata',     [0 0],...
                'LineWidth', 8,...
                'Color',     [1 0 0],...
                'Parent',    pb.ax(i));
        end
        
        pbar = pb;
        drawnow;
        
    % Set
    %======================================================================
    case 'set'
        if isempty(pbar) || ~all(ishandle(pbar.l)), pbar = []; return; end
        if nargin == 1, value = 0; else value = varargin{1}; end
        
        if ischar(value)
            if nargin == 2, str = ''; else str = varargin{2}; end
            if nargin == 3, p   = 1;  else p   = varargin{3}; end
            switch lower(value)
                case {'xlabel','ylabel'}
                    set(get(pbar.ax(p),value),'String',str);
                case 'height'
                    t = clock;
                    bstr = sprintf('Began %2.0f:%02.0f:%02.0f',t(4),t(5),t(6));
                    set(pbar.b(p),'String',bstr);
                    set(pbar.ax(p),'YLim',[0 max([str eps])]);
                otherwise
                    error('Unknown action.');
            end
        else
            if nargin == 2, p = 1; else p = varargin{2}; end
            set(pbar.l(p),'Ydata',[0 value]);
            lim = get(pbar.ax(p),'Ylim');lim=lim(2);
            set(pbar.t(p),'string',sprintf('%.0f%% Complete',100*value/lim));
        end
        try, drawnow limitrate; catch, drawnow; end
        
    % Clear
    %======================================================================
    case 'clear'
        Finter = spm_figure('FindWin','Interactive');
        if isempty(Finter), pbar = []; return; end
        spm_figure('Clear',Finter);
        if isstruct(pbar)
            set(Finter,'Pointer',     pbar.pointer);
            set(Finter,'Name',        pbar.name);
            set(Finter,'DoubleBuffer',pbar.buffer);
        end
        pbar = [];
        drawnow;
    
    % Error
    %======================================================================
    otherwise
        error('Unknown action string');
end
