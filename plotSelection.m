function plotSelection(userData, buttonState)

% PLOTSELECTION is the callback function for interactive ER or TFR plots.
% The function can be called by one of the following functions:
% multiplotER, multiplotTFR, singleplotER, singleplotTFR, and topoplot.
%
% This function is only to be used as a callback function, i.e. you should
% not call it directly.
%
% This function updates the selected areas in the plot, and produces a new
% interactive plot when a selected area is clicked. Multiple areas can be
% selected by holding down the SHIFT key. The newly produced plot is the
% outcome of a statistical function applied to the selected region(s),
% based on the calling plot type (multi, single, or topo). The statistical
% functions used by default are:
%	    for multiplots   : mean of selected channels
%	    for singleplots  : none
%	    for topoplots    : none
%
% userData is a structure array linked to the 'UserData' property of the
% figure that contains the interactive ER/TFR plot. It contains the
% following fields:
%     hFigure          = figure handle
%     hAxes            = axes handle of the ER/TFR plot
%     hSelection{:}    = cell array of line handles of selection rectangles
%     iSelection       = cell array element index of current selection
%                        rectangle
%     selecting        = value specifying the current selecting status:
%                          0: not selecting
%                          1: starting a new selection (mouse button just
%                             pressed)
%                          2: continuing current selection (dragging
%                             rectangle)
%     selectionType    = string specifying whether to extend current
%                        selection (when SHIFT key is held down), or to start
%                        a new selection (otherwise).
%     selectAxes       = string specifying which axes to select; any
%                        combination of x,y,z. For example, 'xy' means the
%                        user can select a range on the x- and y-axes.
%     lastClick        = mouse pointer coordinates at the latest button down
%                        event
%     cfg              = copy of configuration used for the plotting
%     data             = copy of data structure used for the plotting
%     chanX,chanY      = [only used by multiplot routines] x,y-coordinates of
%                        the plotted channels (centre coordinates)
%     chanLabels       = [only used by multiplot routines] label names of the
%                        plotted channels
%
% buttonState can be one of the following values:
%	    0: button is not going up or down (it may be HELD down, though)
%	    1: button is going down
%	    2: button is going up
%

% Undocumented local options:
% cfg.channelname
% cfg.interactive
% cfg.layout
% cfg.maplimits
% cfg.xlim
% cfg.xparam
% cfg.ylim
% cfg.yparam
% cfg.zlim

% Copyright (C) 2006, Dennis Pasveer
%
% $Log: plotSelection.m,v $
% Revision 1.20  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.19  2008/09/22 14:23:25  roboos
% fixed some lint warnings, i.e. no functional change, only efficiency
%
% Revision 1.18  2007/10/01 07:37:54  roboos
% copy cfg.rotate along in the new_cfg, for subsequent topoplots (thanks to Vladimir)
%
% Revision 1.17  2007/09/25 11:30:10  ingnie
% set cfg.baseline = 'no' to prevent multiple calls to baseline correction in interactive plotting
%
% Revision 1.16  2007/07/02 13:42:35  roboos
% only pass cfg.maskparameter on to next function if present
%
% Revision 1.15  2007/06/19 13:59:00  ingnie
% give cfg.maskparameter along in new_cfg
%
% Revision 1.14  2007/06/05 15:05:12  ingnie
% removed dimord repl_chan_time, added dimord rpt_chan_time
%
% Revision 1.13  2007/04/19 10:56:16  roboos
% pass the cfg.interpolation option along to subsequent interactive plots
%
% Revision 1.12  2007/03/14 11:21:29  chrhes
% fixed same bug as in previous revision but on line 252
%
% Revision 1.11  2007/03/14 11:09:25  chrhes
% fixed a small bug to do with copying a possibly nonexistent field 'layout'
% from the userData.cfg structure to new_cfg (on line 222)
%
% Revision 1.10  2007/01/09 11:12:38  roboos
% pass cfg.renderer on to new_cfg for next plot function
%
% Revision 1.9  2007/01/09 10:41:43  roboos
% Added the option cfg.renderer, default is opengl. Interactive plotting
% on linux/VNC sometimes does not work, using cfg.renderer='painters'
% seems to fix it.
%
% Revision 1.8  2006/07/27 16:11:19  roboos
% also recognise timelocked data with keeptrials
% pass the xparam,zparam and layout along to the next plot in case of interactive
%
% Revision 1.7  2006/06/06 16:24:06  ingnie
% replaced cfg.channelindex and cfg.channelname with cfg.channel for consistency
% with other functions
%
% Revision 1.6  2006/05/30 20:57:15  roboos
% updated documentation
%
% Revision 1.5  2006/04/27 16:01:30  ingnie
% added dimord subj_chan_time
%
% Revision 1.4  2006/04/20 09:58:34  roboos
% updated documentation
%
% Revision 1.3  2006/03/22 16:22:34  jansch
% implemented the possibility to take cohspctrm from a multiplot to a singleplot
%
% Revision 1.2  2006/03/13 14:09:28  denpas
% Added semicolon to prevent display of new_cfg assignment.
%
% Revision 1.1  2006/03/03 10:03:43  roboos
% moved plotSelection from private into FT main, since a figure callback
% to a private does not work
%
% Revision 1.3  2006/03/02 12:39:55  jansch
% fixed typo
%
% Revision 1.2  2006/03/02 10:51:39  jansch
% corrected updating of xlim and ylim when creating topoplotTFR from singleplotTFR
%
% Revision 1.1  2006/02/27 15:03:47  denpas
% new function, is used as callback function for interactive plotting in multi-, single- and topoplot
%

fieldtripdefs

% Get cursor coordinates [x y]:
p = get(userData.hAxes, 'CurrentPoint');
p = p(1,1:2);

% Get physical axes limits:
xLim = get(userData.hAxes, 'XLim');
yLim = get(userData.hAxes, 'YLim');

% Update pointer:
if ~userData.selecting
   if buttonState==0
      % Mouse pointer becomes normal when outside plot boundaries, or when in selected area (=clickable):
      if (p(1) < xLim(1)) || (p(1) > xLim(2)) || (p(2) < yLim(1)) || (p(2) > yLim(2)) || inSelection(p, userData.range)
         set(userData.hFigure, 'Pointer', 'default');
      else
         set(userData.hFigure, 'Pointer', 'crosshair');
      end
   end
else
   % While selecting, mouse pointer is always crosshair:
   set(userData.hFigure, 'Pointer', 'crosshair');
end

% Button down starts a selection:
if buttonState==1
   % Check cursor coordinates:
   if (p(1) < xLim(1)) || (p(1) > xLim(2)) || (p(2) < yLim(1)) || (p(2) > yLim(2))
      return;
   end

   if strcmp(userData.selectAxes, 'x')
      p(2) = yLim(1);
   end

   % Starting new selection rectangle:
   userData.selecting = 1;

   % Keep button down coordinates and selection type:
   userData.lastClick = p;
   userData.selectionType = get(userData.hFigure, 'SelectionType');

   % Update figure's UserData:
   set(userData.hFigure, 'UserData', userData);

elseif userData.selecting
   % Distinguish between clicking and dragging:
   clicked = sum(p == userData.lastClick) == 2;
   if strcmp(userData.selectAxes, 'x')
      clicked = p(1) == userData.lastClick(1);
   end

   % Limit cursor coordinates:
   if p(1)<xLim(1), p(1)=xLim(1); end;
   if p(1)>xLim(2), p(1)=xLim(2); end;
   if p(2)<yLim(1), p(2)=yLim(1); end;
   if p(2)>yLim(2), p(2)=yLim(2); end;

   if strcmp(userData.selectAxes, 'xy')
      % Determine selected x,y-range:
      if size(userData.range{1}, 2) == 4
         set(userData.hFigure, 'Name', sprintf('%0s=[%.3g %.3g]/%0s=[%.3g %.3g]', userData.cfg.xparam, min(userData.range{1}([1 3])), max(userData.range{1}([1 3])), userData.cfg.yparam, min(userData.range{1}([2 4])), max(userData.range{1}([2 4]))));
      end
   elseif strcmp(userData.selectAxes, 'z')
      % Determine selected channels:
      selChannels = [];
      for k=1:length(userData.chanLabels)
         if ~isempty(userData.chanLabels{k})
            if inSelection([userData.chanX(k) userData.chanY(k)], userData.range)
               selChannels = [selChannels k];
            end
         end
      end

      % Update window caption:
      set(userData.hFigure, 'Name', sprintf('%.0f channels selected', length(selChannels)));
   elseif strcmp(userData.selectAxes, 'x')
      % Determine selected x-range:
      if size(userData.range{1}, 2) == 4
         set(userData.hFigure, 'Name', sprintf('%0s=[%.3g %.3g]', userData.cfg.xparam, min(userData.range{1}([1 3])), max(userData.range{1}([1 3]))));
      end

      % Set y to the vertical limit:
      p(2) = yLim(2);
   end

   % Launch a new figure if a selected area has been clicked:
   if (buttonState == 2) ...
         && (userData.selecting == 1) ...
         && (clicked) ...
         && (inSelection(p, userData.range))

      % Button has been released, so stop selecting:
      userData.selecting = 0;

      % Update figure's UserData:
      set(userData.hFigure, 'UserData', userData);

      % Check data dimord for ERP data:
      dimord = getDimord(userData);
      if strcmp(dimord, 'chan_time') || strcmp(dimord, 'chan_freq') || strcmp(dimord, 'subj_chan_time') || strcmp(dimord, 'rpt_chan_time')
         if strcmp(userData.selectAxes, 'z')
            if ~isempty(selChannels)
               % Launch ER singleplot figure:
               new_cfg = [];
               if isfield(userData.cfg, 'layout')
                  new_cfg.layout = userData.cfg.layout;
               end
               if isfield(userData.cfg, 'rotate')
                  new_cfg.rotate = userData.cfg.rotate;
               end
               if isfield(userData.cfg, 'interpolation')
                  new_cfg.interpolation = userData.cfg.interpolation;
               end
               if isfield(userData.cfg, 'maskparameter')
                 new_cfg.maskparameter = userData.cfg.maskparameter;
               end
               new_cfg.xparam = userData.cfg.xparam;
               new_cfg.zparam = userData.cfg.zparam;
               new_cfg.interactive = 'yes';
               new_cfg.xlim = userData.cfg.xlim;
               new_cfg.ylim = userData.cfg.ylim;
               new_cfg.renderer = userData.cfg.renderer;
               new_cfg.baseline = 'no'; %to prevent baseline correction to take place multiple times

               % Set xlim to 'maxmin' when going from topoplot to singleplot:
               if strcmp(userData.plotType, 'topoplot')
                  new_cfg.xlim = 'maxmin';
               end

               new_cfg.channel = userData.chanLabels(selChannels);

               figure;
               copyPosition(userData.hFigure);

               if iscell(userData.data)
                  singleplotER(new_cfg, userData.data{:});
               else
                  singleplotER(new_cfg, userData.data);
               end
            end
         else
            % Configure ER topoplot figure:
            new_cfg = [];
            if isfield(userData.cfg, 'layout')
               new_cfg.layout = userData.cfg.layout;
            end
             if isfield(userData.cfg, 'rotate')
                new_cfg.rotate = userData.cfg.rotate;
             end
            if isfield(userData.cfg, 'interpolation')
               new_cfg.interpolation = userData.cfg.interpolation;
            end
            new_cfg.xparam = userData.cfg.xparam;
            new_cfg.zparam = userData.cfg.zparam;
            new_cfg.interactive = 'yes';
            new_cfg.maplimits = 'maxmin';
            new_cfg.ylim = userData.cfg.ylim;
            new_cfg.renderer = userData.cfg.renderer;
            new_cfg.baseline = 'no'; %to prevent baseline correction to take place multiple times

            % Calculate xlim:
            new_cfg.xlim = sort(userData.range{1}([1 3]));

            % Produce topoplot:
            figure;
            copyPosition(userData.hFigure);

            if iscell(userData.data)
               topoplotER(new_cfg, userData.data{:});
            else
               topoplotER(new_cfg, userData.data);
            end
         end

         % Check data dimord for TFR data:
      elseif strcmp(dimord, 'chan_freq_time')
         if strcmp(userData.selectAxes, 'z')
            if ~isempty(selChannels)
               % Launch TFR singleplot figure:
               new_cfg      = userData.cfg;
               new_cfg.xlim = userData.cfg.xlim;
               new_cfg.ylim = userData.cfg.ylim;
               new_cfg.zlim = userData.cfg.zlim;
               new_cfg.baseline = 'no'; %to prevent baseline correction to take place multiple times

               % Set xlim and ylim to 'maxmin' when going from topoplot to singleplot:
               if strcmp(userData.plotType, 'topoplot')
                  new_cfg.xlim = 'maxmin';
                  new_cfg.ylim = 'maxmin';
               end

               % Get channel names of selected channels:
               new_cfg.channel = userData.chanLabels(selChannels);

               figure;
               copyPosition(userData.hFigure);
               singleplotTFR(new_cfg, userData.data);
            end
         else
            % Launch TFR topoplot figure:
            new_cfg           = userData.cfg;
            new_cfg.maplimits = 'maxmin';
            new_cfg.xlim = userData.data.time([ ...
               nearest(userData.data.time, min((userData.range{1}([1 3])))) ...
               nearest(userData.data.time, max((userData.range{1}([1 3])))) ...
               ]);
            new_cfg.ylim = userData.data.freq([ ...
               nearest(userData.data.freq, min((userData.range{1}([2 4])))) ...
               nearest(userData.data.freq, max((userData.range{1}([2 4])))) ...
               ]);
            new_cfg.baseline = 'no'; %to prevent baseline correction to take place multiple times

            % Produce topoplot:
            figure;
            copyPosition(userData.hFigure);
            topoplotTFR(new_cfg, userData.data);
         end

         % Check data dimord for FREQ data:
      elseif (isfield(userData.cfg, 'xparam') && strcmp(userData.cfg.yparam, 'freq'))


      else
         userData.cfg
         error('Unable to produce figure: unidentified data.dimord value.');
      end

      return;
   end

   % If starting new selection rectangle, determine rectangle index:
   if userData.selecting ~= 2
      if strcmp(userData.selectionType, 'normal')
         % If previous selection was finished, clear selection:
         for i=1:length(userData.hSelection)
            set(userData.hSelection{i}, 'Visible', 'off');
            userData.range{i} = [];
         end

         % Update window caption for empty selection:
         if strcmp(userData.selectAxes, 'xy')
            set(userData.hFigure, 'Name', 'No selection');
         elseif strcmp(userData.selectAxes, 'z')
            set(userData.hFigure, 'Name', '0 channels selected');
         elseif strcmp(userData.selectAxes, 'x')
            set(userData.hFigure, 'Name', 'No selection');
         end

         userData.iSelection = 1;
      else
         % Proceed to next selection rectangle:
         userData.iSelection = userData.iSelection + 1;
         if userData.iSelection > length(userData.hSelection)
            userData.iSelection = 1;
         end
      end

      % Set first corner coordinates:
      userData.range{userData.iSelection} = userData.lastClick;

      % We're now selecting:
      userData.selecting = 2;
   end

   % Set rectangle range:
   userData.range{userData.iSelection} = [userData.range{userData.iSelection}([1 2]) p];

   % Keep rectangle coordinates:
   xData = userData.range{userData.iSelection}([1 3 3 1 1]);
   yData = userData.range{userData.iSelection}([2 2 4 4 2]);

   % Plot selection rectangle:
   set(userData.hSelection{userData.iSelection}, 'XData', xData);
   set(userData.hSelection{userData.iSelection}, 'YData', yData);
   set(userData.hSelection{userData.iSelection}, 'Color', [0 0 0]);
   set(userData.hSelection{userData.iSelection}, 'EraseMode', 'xor');
   set(userData.hSelection{userData.iSelection}, 'LineStyle', '--');
   set(userData.hSelection{userData.iSelection}, 'LineWidth', 1.5);
   set(userData.hSelection{userData.iSelection}, 'Visible', 'on');

   % On buttonUp, the selection rectangle is done:
   if buttonState == 2
      userData.selecting = 0;
   end

   % Update figure's UserData:
   set(userData.hFigure, 'UserData', userData);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = inSelection(aPoint, selectionRanges)
% Function returns selection rectangle index number of the rectangle of a
% certain coordinate.

a = 0;
for i=1:length(selectionRanges)
   % Look only in complete selections (containing 4 elemens [left top right bottom]):
   if (length(selectionRanges{i}) == 4) ...
         && (aPoint(:,1) >= min(selectionRanges{i}([1 3]))) ...
         && (aPoint(:,1) <= max(selectionRanges{i}([1 3]))) ...
         && (aPoint(:,2) >= min(selectionRanges{i}([2 4]))) ...
         && (aPoint(:,2) <= max(selectionRanges{i}([2 4])))
      a = i;
      return;
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function copyPosition(hFigure)
% This function copies the window position from a given figure to the
% current figure.

set(gcf, 'position', get(hFigure, 'position'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dimord = getDimord(userData)
% This function reads the dimord value from the userData.data structure.
% Since userData.data may be a cell array of multiple data sets, this
% function first detects the data type.

if iscell(userData.data)
   dimord = userData.data{1}.dimord;
else
   dimord = userData.data.dimord;
end
