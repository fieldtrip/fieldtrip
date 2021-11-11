function pBar = bst_progress(varargin)
% bst_progress: Manage the Brainstorm progress bar
%
% USAGE : pBar = bst_progress('start', title, msg, valStart, valStop) : Create a progress bar with start and stop bounds
%         pBar = bst_progress('start', title, msg)                    : Create a progress bar (unlimited)
%         pBar = bst_progress('stop')        : stop and hide progress bar
%         pBar = bst_progress('inc', valInc) : increment of 'valInc' the position of the progress bar
%         pBar = bst_progress('set', pos)    : set the position
%          pos = bst_progress('get')         : get the position
%         pBar = bst_progress('text', txt)   : set the text
%    isVisible = bst_progress('isvisible')   : return 1 if progress bar is visible, 0 else
%         pBar = bst_progress('show')        : display previously defined progress bar
%         pBar = bst_progress('hide')        : hide progress bar
%         pBar = bst_progress('setimage', imagefile) : display an image in the wait bar
%         pBar = bst_progress('setlink', url)        : clicking on the image opens a browser to display the url
%         pBar = bst_progress('removeimage')         : Remove the image from the wait bar

% NOTES : The window is created once, and then never deleted, just hidden.
%         Progress bar is represented by a structure: 
%            |- jWindow
%            |- jLabel
%            |- jProgressBar
%
% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2020 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Francois Tadel, 2008-2013

% JAVA imports
import org.brainstorm.icon.*;
import java.awt.Dimension;

global GlobalData;


%% ===== PARSE INPUTS =====
if ((nargin >= 1) && ischar(varargin{1}))
    commandName = varargin{1};
else
    error('Usage : bst_progress(commandName, parameters)');
end


%% ===== GET OR CREATE PROGRESS BAR =====
% Do nothing in case of server mode
if ~isempty(GlobalData) && ~isempty(GlobalData.Program) && isfield(GlobalData.Program, 'GuiLevel') && (GlobalData.Program.GuiLevel == -1)
    if ismember(lower(commandName), {'pos','isvisible'})
        pBar = 0;
    else
        pBar = [];
    end
    return;
end
% If running in NOGUI mode: just display the message in the command window
if ~bst_get('isGUI')
    switch lower(commandName)
        case 'start',     disp(['PROGRESS> ' varargin{2} ': ' varargin{3}]); pBar = [];
        case 'text',      disp(['PROGRESS> ' varargin{2}]); pBar = [];
        case 'get',       pBar = 1;
        case 'isvisible', pBar = 0;
    end
    return;
end
% Get progress bar
if ~isempty(GlobalData) && ~isempty(GlobalData.Program) && isfield(GlobalData.Program, 'ProgressBar') && ~isempty(GlobalData.Program.ProgressBar)
    pBar = GlobalData.Program.ProgressBar;
else
    pBar = [];
end
% Get Brainstorm GUI context
jBstFrame = bst_get('BstFrame');
if isempty(jBstFrame)
    return
end

% Default window size
if isempty(pBar) || strcmpi(commandName, 'removeimage') || (strcmpi(commandName, 'stop') && pBar.isImage)
    DefaultSize = java_scaled('dimension', 350, 130);
end

% If progress bar was not created yet : create it
if isempty(pBar)
    % If action=isvisible: no need to create the progress bar
    if strcmpi(commandName, 'isvisible')
        pBar = 0;
        return
    end
    % Create a JDialog, if possible dependent of the main Brainstorm JFrame
    pBar.jWindow = java_create('javax.swing.JDialog');
    % Set icon
    try
        pBar.jWindow.setIconImage(IconLoader.ICON_APP.getImage());
    catch
        % Old matlab... just ignore...
    end
    % Set as always-on-top / non-focusable
    pBar.jWindow.setAlwaysOnTop(1);
    pBar.jWindow.setFocusable(0);
    pBar.jWindow.setFocusableWindowState(0);
    % Non-modal
    pBar.jWindow.setModal(0);
    
    % Closing callback
%     if bst_iscompiled()
        pBar.jWindow.setDefaultCloseOperation(pBar.jWindow.HIDE_ON_CLOSE);
%     else
%         pBar.jWindow.setDefaultCloseOperation(pBar.jWindow.DO_NOTHING_ON_CLOSE);
%         java_setcb(pBar.jWindow, 'WindowClosingCallback', @(h,ev)CloseCallback);
%     end

    % Configure window
    pBar.jWindow.setPreferredSize(DefaultSize);
    % Main panel
    pBar.jPanel = java_create('javax.swing.JPanel');
    pBar.jPanel.setLayout(java_create('java.awt.GridBagLayout'));

    % Create objects
    pBar.isImage = 0;
    pBar.jImage = java_create('javax.swing.JLabel');
    pBar.jLabel = java_create('javax.swing.JLabel', 'Ljava.lang.String;', '...');
    pBar.jLabel.setFont(bst_get('Font'));
    pBar.jProgressBar = java_create('javax.swing.JProgressBar', 'II', 0, 99);
    % Update constraints
    UpdateConstraints(0);
    
    % Add the main Panel
    pBar.jWindow.getContentPane.add(pBar.jPanel);
    pBar.jWindow.pack();
    % Set window size and location
    %pBar.jWindow.setLocationRelativeTo(pBar.jWindow.getParent());
    jLoc = jBstFrame.getLocation();
    jSize = jBstFrame.getSize();
    pos = [jLoc.getX() + ((jSize.getWidth() - DefaultSize.getWidth()) / 2), ...
           jLoc.getY() + ((jSize.getHeight() - DefaultSize.getHeight()) / 2)];
    pBar.jWindow.setLocation(pos(1), pos(2));
    % Save progress bar
    GlobalData.Program.ProgressBar = pBar;
end

% Linux: need to print something on the command window (don't know why...)
if strcmpi(commandName, 'stop') && ismember(computer('arch'), {'glnx86', 'glnxa64'})
    drawnow();
    fprintf(' ');
    fprintf('\b');
    drawnow();
    % The dialog needs to be displayed for a short period before being hidden
    pause(0.05);
end


%% ===== SWITCH BETWEEN COMMANDS =====
switch (lower(commandName))
    % ==== START ====
    case 'start'
        % Set as "always on top"
        java_call(pBar.jWindow, 'setAlwaysOnTop', 'Z', 1);
        java_call(pBar.jWindow, 'setFocusable',   'Z', 0);
        java_call(pBar.jWindow, 'setFocusableWindowState', 'Z', 0);
        % Call: bst_progress(''start'', title, msg)
        if (nargin == 3) && ischar(varargin{3}) && ischar(varargin{2})
            % Get window title and comment
            wndTitle    = varargin{2};
            msg = varargin{3};
            % Set Progress bar in inderminate mode
            pBar.jProgressBar.setIndeterminate(1);
            pBar.jProgressBar.setStringPainted(0);
            % Set progress bar bounds
            pBar.jProgressBar.setMinimum(0);
            pBar.jProgressBar.setMaximum(100);
            % Set initial value to start
            pBar.jProgressBar.setValue(0);
            
        % Call: bst_progress(''start'', title, msg, start, stop)
        elseif ((nargin == 5) && ischar(varargin{2}) && ischar(varargin{3}) && isnumeric(varargin{4}) && isnumeric(varargin{5}))
            % Get window title and comment
            wndTitle = varargin{2};
            msg = varargin{3};
            % Set Progress bar in derminate mode
            pBar.jProgressBar.setIndeterminate(0);
            pBar.jProgressBar.setStringPainted(1);
            % Get progress bar bounds
            valStart = varargin{4};
            valStop  = varargin{5};
            % Test bounds
            if ( (valStart >= valStop) || (valStop <= 0) )
                % Set indeterminate bounds
                pBar.jProgressBar.setIndeterminate(1);
                pBar.jProgressBar.setStringPainted(0);
                valStart = 0;
                valStop  = 100;
            end
            % Set progress bar bounds
            pBar.jProgressBar.setMinimum(valStart);
            pBar.jProgressBar.setMaximum(valStop);
            % Set initial value to start
            pBar.jProgressBar.setValue(valStart);
            
        else
            error(['Usage : bst_progress(''start'', title, comment) ' 10 '        bst_progress(''start'', title, comment, valStart, valStop)']);
        end
        % Set window title
        pBar.jWindow.setTitle(wndTitle);
        % Set window comment (central label)
        pBar.jLabel.setText(msg);
        % Show window
        java_call(pBar.jWindow, 'setVisible', 'Z', 1);
        % Set watch cursor
        jBstFrame.setCursor(java_create('java.awt.Cursor', 'I', java.awt.Cursor.WAIT_CURSOR));
        
    % ==== STOP ====
    case 'stop'
        % Remove the "always on top" status
        java_call(pBar.jWindow, 'setAlwaysOnTop', 'Z', 0);
        java_call(pBar.jWindow, 'setFocusable',   'Z', 1);
        java_call(pBar.jWindow, 'setFocusableWindowState', 'Z', 1);
        % Hide window
        java_call(pBar.jWindow, 'setVisible', 'Z', 0);
        % Restore cursor
        jBstFrame.setCursor([]);
        % Remove image
        if pBar.isImage
            GlobalData.Program.ProgressBar.isImage = 0;
            pBar.jImage.setIcon([]);
            java_setcb(pBar.jImage, 'MouseClickedCallback', []);
            UpdateConstraints(0);
            pBar.jWindow.setPreferredSize(DefaultSize);
            pBar.jWindow.pack();
        end
        
    % ==== INCREMENT ====
    case 'inc'
        % Parse arguments
        if ((nargin == 2) && isnumeric(varargin{2}))
            valInc = varargin{2};
        else
            error('Usage : bst_progress(''inc'', valInc)');
        end
        % Get current value
        curValue = pBar.jProgressBar.getValue();
        % Get the incremented progress bar position
        newVal = min(curValue + valInc, pBar.jProgressBar.getMaximum());
        pBar.jProgressBar.setValue(newVal);
        
    % ==== SET POSITION ====
    case 'set'
        % Parse arguments
        if ((nargin == 2) && isnumeric(varargin{2}))
            newVal = varargin{2};
        else
            error('Usage : bst_progress(''set'', pos)');
        end
        % Get current value
        curValue = pBar.jProgressBar.getValue();
        % Get the incremented progress bar position
        if (curValue ~= newVal)
            newVal = min(newVal, pBar.jProgressBar.getMaximum());
            pBar.jProgressBar.setValue(newVal);
        end
        
    % ==== GET POSITION ====
    case 'get'
        % Get the incremented progress bar position
        pBar = pBar.jProgressBar.getValue();
        
    % ==== SET TEXT ====
    case 'text'
        % Parse arguments
        if ((nargin == 2) && ischar(varargin{2}))
            % Set new label
            pBar.jLabel.setText(varargin{2});
        else
            % Get label
            pBar.jLabel.getText();
        end
        
    % ==== IS VISIBLE ====
    case 'isvisible'
        pBar = pBar.jWindow.isVisible();
        
    % ==== SHOW ====
    case 'show'
        % Set as "always on top"
        java_call(pBar.jWindow, 'setAlwaysOnTop', 'Z', 1);
        java_call(pBar.jWindow, 'setFocusable',   'Z', 0);
        java_call(pBar.jWindow, 'setFocusableWindowState', 'Z', 0);
        % Show window
        java_call(pBar.jWindow, 'setVisible', 'Z', 1);
        % Set watch cursor
        jBstFrame.setCursor(java_create('java.awt.Cursor', 'I', java.awt.Cursor.WAIT_CURSOR));
    
    % ==== HIDE ====
    case 'hide'
        % Remove the "always on top" status
        java_call(pBar.jWindow, 'setAlwaysOnTop', 'Z', 0);
        java_call(pBar.jWindow, 'setFocusable',   'Z', 1);
        java_call(pBar.jWindow, 'setFocusableWindowState', 'Z', 1);
        % Hide window
        java_call(pBar.jWindow, 'setVisible', 'Z', 0);
        % Restore cursor
        jBstFrame.setCursor([]);
        
    % ==== SET IMAGE ====
    case 'setimage'
        % Get image path
        imagefile = varargin{2};
        if ~file_exist(imagefile)
            imagefile = bst_fullfile(bst_get('BrainstormDocDir'), imagefile);
        end
        if ~file_exist(imagefile)
            warning(['Image not found: ' imagefile]);
            return
        end
        % Image in label
        pBar.jImage.setIcon(javax.swing.ImageIcon(imagefile));
        GlobalData.Program.ProgressBar.isImage = 1;
        % Extend size of the frame
        UpdateConstraints(1);
        pBar.jWindow.setPreferredSize([]);
        pBar.jWindow.pack();
        
    % ==== SET LINK ====
    case 'setlink'
        url = varargin{2};
        java_setcb(pBar.jImage, 'MouseClickedCallback', @(h,ev)web(url, '-browser'));
        
    % ==== REMOVE IMAGE ====
    case 'removeimage'
        % Remove image
        GlobalData.Program.ProgressBar.isImage = 0;
        pBar.jImage.setIcon([]);
        java_setcb(pBar.jImage, 'MouseClickedCallback', []);
        UpdateConstraints(0);
        pBar.jWindow.setPreferredSize(DefaultSize);
        pBar.jWindow.pack();
end

%% ===== ADD COMPONENTS =====
    function UpdateConstraints(isImage)
        import java.awt.GridBagConstraints;
        import java.awt.Insets;
        % Remove all components
        pBar.jPanel.removeAll();
        % Generic constraints
        c = GridBagConstraints();
        c.fill = GridBagConstraints.BOTH;
        c.gridx = 1;
        c.weightx = 1;
        % IMAGE
        c.gridy = 1;
        c.weighty = isImage;
        c.insets = Insets(0,0,0,0);
        pBar.jPanel.add(pBar.jImage, c);
        % TEXT
        c.gridy = 2;
        c.weighty = ~isImage;
        c.insets = Insets(5,12,5,12);
        pBar.jPanel.add(pBar.jLabel, c);
        % PROGRESS BAR
        c.gridy = 3;
        c.weighty = 0;
        c.insets = Insets(0,12,9,12);
        pBar.jPanel.add(pBar.jProgressBar, c);
%         % CANCEL BUTTON
%         c.gridy = 4;
%         c.weighty = 0;
%         c.insets = Insets(0,12,9,12);
%         c.weightx = 0;
%         pBar.jPanel.add(pBar.jButtonCancel, c);
    end



%     %% ===== CLOSE CALLBACK =====
%     function CloseCallback()
%         % Hide progress bar
%         %java_call(pBar.jWindow, 'setVisible', 'Z', 0);
%         bst_progress('stop');
%         
%         if bst_iscompiled()
%             try 
%                 % Get command window
%                 cmdWindow = com.mathworks.mde.cmdwin.CmdWin.getInstance();
%                 cmdWindow.grabFocus();
%                 %2) Wait for focus transfer to complete (up to 2 seconds)
%                 focustransferTimer = tic;
%                 while ~cmdWindow.isFocusOwner
%                     pause(0.1);  %Pause some small interval
%                     if (toc(focustransferTimer) > 2)
%                         error('Error transferring focus for CTRL+C press.')
%                     end
%                 end
% 
%                 %3) Use Java robot to execute a CTRL+C in the (now focused) command window.
% 
%                 %3.1)  Setup a timer to relase CTRL + C in 0.3 second
%                 %  Try to reuse an existing timer if possible (this would be a holdover
%                 %  from a previous execution)
%                 t_all = timerfindall;
%                 releaseTimer = [];
%                 ix_timer = 1;
%                 while isempty(releaseTimer) && (ix_timer<= length(t_all))
%                     if isequal(t_all(ix_timer).TimerFcn, @releaseCtrl_C)
%                         releaseTimer = t_all(ix_timer);
%                     end
%                     ix_timer = ix_timer+1;
%                 end
%                 if isempty(releaseTimer)
%                     releaseTimer = timer;
%                     releaseTimer.TimerFcn = @(h,ev)releaseCtrl_C;
%                 end
%                 releaseTimer.StartDelay = 0.3;
%                 start(releaseTimer);
% 
%                 %3.2)  Press press CTRL+C
%                 pressCtrl_C();
%             catch
%                 disp('BST> Could not post a CTRL+C signal in the command window.');
%             end
%         end
%     end
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     function pressCtrl_C()
%         SimKey = java.awt.Robot;
%         SimKey.keyPress(java.awt.event.KeyEvent.VK_CONTROL);
%         SimKey.keyPress(java.awt.event.KeyEvent.VK_C);
%     end
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     function releaseCtrl_C()
%         SimKey = java.awt.Robot;
%         SimKey.keyRelease(java.awt.event.KeyEvent.VK_CONTROL);
%         SimKey.keyRelease(java.awt.event.KeyEvent.VK_C);
%         jBstFrame.setVisible(1);
%     end

end



