function [cs] = NewCs(Name, win, p1, p2, p3, p4, p5)
%
%   [cs] = NewCs (Name, Con, P1, P2)
%
%Purpose:
%   Create a command structure for passing to TellAfni
%   
%   
%Input Parameters:
%   Name: Name of command, see AFNI's README.driver for all possibilities
%         Example: SET_THRESHOLD 
%         To see a list of commands, try the function TellAfni_Commands
%         To start AFNI, use the special Name: 'start_afni'
%   Con : Index of controller where command is applied.
%         Choose from 'A' to 'J'. Default if win = '' is controller 'A'
%   P1  : First parameter passed to command. This can include the
%         controller index in the form A.  See AFNI's README.driver
%         for information and Test_TellAfni for examples.
%   P2  : All other parameters to pass to command. This would include
%         any second parameters plus options.
%Output Parameters:
%   cs  : Returned structure to feed to AFNI with TellAfni
%       .err: Error flag
%       .w  : AFNI controller
%       .c  : AFNI command (Name)
%       .v  : AFNI command value
%   
%      
%More Info :
%    TellAfni
%    TellAfni_Commands
%    Test_TellAfni
%    AFNI's README.driver
%
%     Author : Ziad Saad
%     Date : Tue Dec 6 12:05:09 EST 2005
%     SSCC/NIMH/ National Institutes of Health, Bethesda Maryland


%Define the function name for easy referencing
FuncName = 'NewCs';
%Debug Flag
DBG = 1;

%initailize return variables
cs.err = 1;
cs.c = '';
cs.v = '';
cs.w = '';

np = nargin - 2; %number of p parameters

if (np == -1), win = ''; end
if (np == 0), p1 = ''; p2 = ''; p3 = ''; p4 = ''; p5 = ''; end
if (np == 1),          p2 = ''; p3 = ''; p4 = ''; p5 = ''; end
if (np == 2),                   p3 = ''; p4 = ''; p5 = ''; end
if (np == 3),                            p4 = ''; p5 = ''; end
if (np == 4),                                     p5 = ''; end

if (~isempty(win)), 
   win = upper(win);
   uwin = win; 
else 
   win = 'A';
   uwin = '';
end

if (length(win) ~= 1 | win < 'A' | win > 'J'),
   fprintf(2,'Bad window specifer %s\n', win); 
   return;           
end
  
Name = upper(Name);
switch (Name),
   case 'START_AFNI',
      cs.c = sprintf('%s', Name);
      if (~isempty(p1)),
         if (NewCs_CheckFnames(p1) > 0),
            cs.v = p1;
         else
            fprintf(2,'Command <%s> requires an existing directory or files for the first parameter\nSomething in %s not found\n', Name, p1)
         end
      end
      cs.v = sprintf('%s %s', cs.v, p2);;
      
   case 'ADD_OVERLAY_COLOR',
      if (np ~= 2 | isnumeric(p1) | isnumeric(p2)),
         fprintf(2,'Command <%s> requires two string parameters.\nHave %d in %s\n',...
             Name, np); 
         return;
      end
      cs.c = Name;
      cs.w = win;
      cs.v = sprintf('%s %s', p1, p2);
   case 'SET_THRESHOLD',
      if (np < 1 | np > 2),
         fprintf(2,'Command <%s> needs 1 or 2 parameters.\nHave %d\n', Name, np);
         return;
      end
      
      if (~isnumeric(p1)),
         if (upper(p1(1)) >= 'A' & upper(p1(1)) <= 'J'),
            win = p1(1);
            p1 = p1(2:length(p1));
         end
         p1 = str2double(p1);
      end
      if (p1 < 0 | p1 >= 1.0),
         fprintf(2,'Command <%s> has a bad value for threshold (%s).\nThreshold must be between 0 and 1\n',...
                     Name, p1);
         return;
      end
      if (np == 2),
         if (p2 < 0 | p2 > 4),
            fprintf(2,'Command <%s> has a bad value for decimal (%f).\nDecimal must be between 0 and 4\n',...
                     cs.c, p2);
         end 
         return;
      end
      cs.v = sprintf('%c.%d %d', win, p1*10000, p2);
      cs.w = win;
      cs.c = Name;
   case 'SET_THRESHNEW',
      val = p1;
      isp = '';
      if (np == 1), dec = ''; 
      elseif (np==2), 
         if (~isempty(find(p2 == '*'))), dec = '*';
         else dec = '';
         end
         if (~isempty(find(p2 == 'p'))), 
            isp = 'p';
            if (val < 0.0 | val > 1.0),
               fprintf(2,'Command <%s> requires pvalues between 0 and 1.0\n',...
                  Name); 
            end
         else 
            isp = '';
         end
      end
      cs.w = win;
      cs.c = Name;
      cs.v = sprintf('%c %.14f %c%c', cs.w, val, dec, isp); 
   case 'SET_PBAR_NUMBER',
      if (np ~= 1),
         fprintf(2,'Command <%s> requires 1 parameter\n', Name);
         return;
      end
      if (~isnumeric(p1)),
         if (upper(p1(1)) >= 'A' & upper(p1(1)) <= 'J'),
            win = p1(1);
            p1 = p1(2:length(p1));
         end
         p1 = str2double(p1);
      end
      if (p1 < 2 | p1 > 20),
         fprintf(2,'Command <%s> requires a number between 2 and 20\nHave %f', Name, p1);
         return;
      end
      cs.w = win;
      cs.c = Name;
      cs.v = sprintf('%c.%d', cs.w, p1);
   case 'SET_PBAR_SIGN',
      if (np ~= 1),
         fprintf(2,'Command <%s> requires 1 parameter\n', Name);
         return;
      end      
      [win, p1, p1n] = NewCs_GetWin(p1, uwin);
      if (p1(1) ~= '+' & p1(1) ~= '-'),
         fprintf(2,'parameter must be either + or -\n', Name);
         return;
      end 
      cs.w = win; cs.c = Name; cs.v = sprintf('%c.%s', win, p1);
   case 'SET_PBAR_ALL',
      if (np ~= 2),
         fprintf(2,'Command <%s> requires 2 parameters\n', Name);
         return;
      end
      [win, p1, p1n] = NewCs_GetWin(p1, uwin);
      ps = p1(1); p1 = p1(2:length(p1));
      if (ps ~= '+' & ps ~= '-'),
         fprintf(2,'Missing sign for pbar\n', Name);
         return;
      end
      ncol = str2num(p1);
      if (ncol < 1 | ncol > 99),
         fprintf(2,'Command <%s> needs num between 1 and 99\n', Name);
         return;
      end
      if (ncol < 99),
         %find out how many colors are specified in p2   
         if (WordCount(p2) ~= ncol),
            fprintf(2,'Command <%s>: Mismatch between num (%d) and number of val=color strings (%d)\n', Name, ncol, WordCount(p2));
            return;
         end
         cs.w = win; cs.c = Name; cs.v = sprintf('%c.%c%d %s', win, ps, ncol, p2);
      else
         %check on the topval 
         np2 = WordCount(p2);
         if (np2 < 2 ),
            fprintf(2,'Command <%s> needs 2 parameters in second string (topval colorscale_name)\n', Name);
            return;
         end
         [err, topval] = GetWord(p2,1,' ');
         topval = str2num(topval);
         if (topval <= 0),
            fprintf(2,'Command <%s> needs a positive topval\n', Name);
            return;
         end
         cs.w = win; cs.c = Name; cs.v = sprintf('%c.%c99 %s', win, ps, p2);
      end
   case 'PBAR_ROTATE',
      if (np ~= 1),
         fprintf(2,'Command <%s> requires 1 parameter\n', Name);
         return;
      end
      [win, p1, p1n] = NewCs_GetWin(p1, uwin);
      if (p1(1) ~= '+' & p1(1) ~= '-'),
         fprintf(2,'parameter must be either + or -\n', Name);
         return;
      end 
      cs.w = win; cs.c = Name; cs.v = sprintf('%c.%s', win, p1);
   case 'DEFINE_COLORSCALE',
      if (np ~= 2),
         fprintf(2,'Command <%s> requires 2 parametera\n', Name);
         return;
      end
      cs.w = win; cs.c = Name; cs.v = sprintf('%s %s', p1, p2);
   case 'SET_FUNC_AUTORANGE',
      if (np ~= 1),
         fprintf(2,'Command <%s> requires 1 parameter\n', Name);
         return;
      end
      [win, p1, p1n] = NewCs_GetWin(p1, uwin);
      if (p1(1) ~= '+' & p1(1) ~= '-'),
         fprintf(2,'parameter must be either + or -\n', Name);
         return;
      end 
      cs.w = win; cs.c = Name; cs.v = sprintf('%c.%s', win, p1);
   case 'SET_FUNC_RANGE',
      if (np ~= 1),
         fprintf(2,'Command <%s> requires 1 parameter\n', Name);
         return;
      end
      [win, p1, p1n] = NewCs_GetWin(p1, uwin);
      if (p1n < 0.0),
         fprintf(2,'Command <%s> requires a positive parameter\n', Name);
         return;
      end
      cs.w = win; cs.c = Name; cs.v = sprintf('%c.%f', win, p1n);
   case { 'SET_FUNC_VISIBLE', 'SEE_OVERLAY', 'SEE_FUNCTION'},
      if (np ~= 1),
         fprintf(2,'Command <%s> requires 1 parameter\n', Name);
         return;
      end
      [win, p1, p1n] = NewCs_GetWin(p1, uwin);
      if (p1(1) ~= '+' & p1(1) ~= '-'),
         fprintf(2,'parameter must be either + or -\n', Name);
         return;
      end 
      cs.w = win; cs.c = Name; cs.v = sprintf('%c.%s', win, p1);
   case 'SET_FUNC_RESAM',
      if (np ~= 1),
         fprintf(2,'Command <%s> requires 1 parameter\n', Name);
         return;
      end
      [win, p1, p1n] = NewCs_GetWin(p1, uwin);
      [err,l1] = GetWord(p1,1,'.');
      if (isempty(l1)),
         fprintf(2,'Command <%s> has empty parameter\n', Name);
         return;
      end
      err = NewCs_OKresam(l1);
      if (err),
         fprintf(2,'Function resampling mode %s not recognized.\nChoose from:\n%s\n',...
                      l1, NewCs_OKresam());
         return;
      end
      [err, l2] = GetWord(p1,2,'.');
      if (~isempty(l2)),
         err = NewCs_OKresam(l2);
         if (err),
            fprintf(2,'Resampling mode %s not recognized.\nChoose from:\n%s\n',...
                      l2, NewCs_OKresam());
            return;
         end
      end
      cs.w = win; cs.c = Name; cs.v = sprintf('%c.%s.%s', win, l1, l2);
   case 'OPEN_PANEL',
      if (np ~= 1),
         fprintf(2,'Command <%s> requires 1 parameter\n', Name);
         return;
      end
      [win, p1, p1n] = NewCs_GetWin(p1, uwin);
      err = NewCs_OKpanel(p1);
      if (err),
         fprintf(2,'Panel %s not recognized.\nChoose from:\n%s\n',...
                      p1, NewCs_OKpanel());
         return;
      end
      cs.w = win; cs.c = Name; cs.v = sprintf('%c.%s',win,p1);
   case 'SYSTEM',
      if (np < 1),
         fprintf(2,'Command <%s> requires at least 1 parameter\n', Name);
         return;
      end
      cs.w = win; cs.c = Name; cs.v = sprintf('%s %s', p1, p2);
   case 'CHDIR',
      if (np ~= 1),
         fprintf(2,'Command <%s> requires 1 parameter\n', Name);
         return;
      end
      cs.w = win; cs.c = Name; cs.v = sprintf('%s', p1);
   case 'RESCAN_THIS',
      if (np > 1),
         fprintf(2,'Command <%s> requires 1 parameter\n', Name);
         return;
      end
      if (np),
         if (p1(1) < 'A' | p1(1) > 'J'),
            fprintf(2,'Command <%s> requires 1 parameter between A and J\n', Name);
            return;
         end
         win = p1(1);      
      end
      cs.w = win; cs.c = Name;  cs.v = sprintf('%c', win);
   case {'SET_SESSION', 'SWITCH_SESSION', 'SWITCH_DIRECTORY'},
      if (np ~= 1),
         fprintf(2,'Command <%s> requires 1 parameter\n', Name);
         return;
      end
      [win, p1, p1n] = NewCs_GetWin(p1, uwin);
      cs.w = win; cs.c = Name; cs.v = sprintf('%c.%s',win,p1);
   case {'SET_SUBBRICKS','SET_SUB_BRICKS'},
      if (np < 1),
         fprintf(2,'Command <%s> requires at least 1 parameter\n', Name);
         return;
      end
      if (np == 2), win = p1; p1 = p2; end
      val = str2num(p1);
      if (length(val) ~= 3),
         fprintf(2,'Second string must contain 3 values.\n', Name);
         return;
      end
      cs.w = win; cs.c = Name; cs.v = sprintf('%c %d %d %d', win, val(1), val(2), val(3));
      
   case {'SET_ANATOMY', 'SWITCH_ANATOMY', 'SWITCH_UNDERLAY'},
      if (np < 1),
         fprintf(2,'Command <%s> requires at least 1 parameter\n', Name);
         return;
      end
      [win, p1, p1n] = NewCs_GetWin(p1, uwin);
      if (np == 2),
         sb = str2num(p2);
         if (length(sb) ~= 1),
            fprintf(2,'Command <%s> requires one integer in the second parameter');
            return;
         end
      end
      cs.w = win; cs.c = Name; cs.v = sprintf('%c.%s %s', win, p1, p2);
   case {'SET_FUNCTION', 'SWITCH_FUNCTION', 'SWITCH_OVERLAY'},
      if (np < 1),
         fprintf(2,'Command <%s> requires at least 1 parameter\n', Name);
         return;
      end
      if (np == 2),
         sb = str2num(p2);
         if (length(sb) ~= 1 & length(sb) ~= 2),
            fprintf(2,'Command <%s> requires one or two integers in the second parameter');
            return;
         end
      end
      [win, p1, p1n] = NewCs_GetWin(p1, uwin);
      cs.w = win; cs.c = Name; cs.v = sprintf('%c.%s %s', win, p1, p2);
   case 'OPEN_WINDOW',
      if (np < 1),
         fprintf(2,'Command <%s> requires 1 parameter at least\n', Name);
         return;
      end
      [win, p1, p1n] = NewCs_GetWin(p1, uwin);
      err = NewCs_OKwindowname(p1);
      if (err),
         fprintf(2,'Window name %s not recognized.\nChoose from:\n%s\n',...
                      p1, NewCs_OKwindowname());
         return;
      end
      cs.w = win; cs.c = Name; cs.v = sprintf('%c.%s %s', win, p1, p2);   
   case 'CLOSE_WINDOW',
      if (np ~= 1),
         fprintf(2,'Command <%s> requires 1 parameter \n', Name);
         return;
      end
      [win, p1, p1n] = NewCs_GetWin(p1, uwin);
      err = NewCs_OKwindowname(p1);
      if (err),
         fprintf(2,'Window name %s not recognized.\nChoose from:\n%s\n',...
                      p1, NewCs_OKwindowname());
         return;
      end
      cs.w = win; cs.c = Name; cs.v = sprintf('%c.%s', win, p1);  
   case 'SAVE_JPEG',
      if (np ~= 2),
         fprintf(2,'Command <%s> requires 2 parameters\n', Name);
         return;
      end
      [win, p1, p1n] = NewCs_GetWin(p1, uwin);
      err = NewCs_OKwindowname(p1);
      if (err),
         fprintf(2,'Window name %s not recognized.\nChoose from:\n%s\n',...
                      p1, NewCs_OKwindowname());
         return;
      end
      cs.w = win; cs.c = Name; cs.v = sprintf('%c.%s %s', win, p1, p2);
   case { 'SET_DICOM_XYZ', 'SET_SPM_XYZ', 'SET_IJK'}, 
      if (np < 1),
         fprintf(2,'Command <%s> requires at least 1 parameter\n', Name);
         return;
      end
      if (np == 2), win = p1; p1 = p2; end
      val = str2num(p1);
      if (length(val) ~= 3),
         fprintf(2,'Second string must contain 3 coordinates.\n', Name);
         return;
      end
      cs.w = win; cs.c = Name; cs.v = sprintf('%c %f %f %f', win, val(1), val(2), val(3));
   case 'SET_XHAIRS',
      if (np ~= 1),
         fprintf(2,'Command <%s> requires 1 parameter \n', Name);
         return;
      end
      [win, p1, p1n] = NewCs_GetWin(p1, uwin);
      err = NewCs_OKxhaircode(p1);
      if (err),
         fprintf(2,'Xhair code %s not recognized.\nChoose from:\n%s\n',...
                      p1, NewCs_OKxhaircode());
         return;
      end
      cs.w = win; cs.c = Name; cs.v = sprintf('%c.%s', win, p1);  
   case 'PURGE_MEMORY',
      cs.w = win; cs.c = Name; cs.v = sprintf('%s', p1);  
   case 'QUIT',
      cs.w = win; cs.c = Name; cs.v = '';
   case 'SETENV',
      cs.w = win; cs.c = Name; cs.v = sprintf('%s %s', p1, p2);  
   case 'REDISPLAY',
      cs.w = win; cs.c = Name; cs.v = '';
   case 'SLEEP',
      if (np ~= 1),
         fprintf(2,'Command <%s> requires 1 parameter \n', Name);
         return;
      end
      if (~isnumeric(p1)),
         fprintf(2,'Command <%s> requires at 1 numeric value\n', Name);
         return;
      end
      cs.w = win; cs.c = Name; cs.v = sprintf('%f', p1);
   case 'OPEN_GRAPH_XY',
      if (np ~= 2),
         fprintf(2,'Command <%s> requires 2 parameters \n', Name);
         return;
      end
      cs.w = win; cs.c = Name; cs.v = sprintf('%s %s', p1, p2);
   case 'CLOSE_GRAPH_XY',
      if (np ~= 1),
         fprintf(2,'Command <%s> requires 1 parameter \n', Name);
         return;
      end
      cs.w = win; cs.c = Name; cs.v = sprintf('%s', p1);
   case 'CLEAR_GRAPH_XY',
      if (np ~= 1),
         fprintf(2,'Command <%s> requires 1 parameter \n', Name);
         return;
      end
      cs.w = win; cs.c = Name; cs.v = sprintf('%s', p1);
   case 'ADDTO_GRAPH_XY',
      if (np ~= 2),
         fprintf(2,'Command <%s> requires 2 parameters \n', Name);
         return;
      end
      cs.w = win; cs.c = Name; cs.v = sprintf('%s %s', p1, p2);
   case 'OPEN_GRAPH_1D',
      if (np ~= 2),
         fprintf(2,'Command <%s> requires 2 parameters \n', Name);
         return;
      end
      cs.w = win; cs.c = Name; cs.v = sprintf('%s %s', p1, p2);
   case 'CLOSE_GRAPH_1D',
      if (np ~= 1),
         fprintf(2,'Command <%s> requires 1 parameter \n', Name);
         return;
      end
      cs.w = win; cs.c = Name; cs.v = sprintf('%s', p1);
   case 'CLEAR_GRAPH_1D',
      if (np ~= 1),
         fprintf(2,'Command <%s> requires 1 parameter \n', Name);
         return;
      end
      cs.w = win; cs.c = Name; cs.v = sprintf('%s', p1);
   case 'ADDTO_GRAPH_1D',
      if (np ~= 2),
         fprintf(2,'Command <%s> requires 2 parameters \n', Name);
         return;
      end
      cs.w = win; cs.c = Name; cs.v = sprintf('%s %s', p1, p2);
   case 'SET_GRAPH_GEOM',
      if (np ~= 2),
         fprintf(2,'Command <%s> requires 2 parameters \n', Name);
         return;
      end
      cs.w = win; cs.c = Name; cs.v = sprintf('%s %s', p1, p2);
   case 'MADMAG',
      fprintf(2,'Command <%s> not implemented yet.\n', Name);
      return;
   otherwise,
      fprintf(2,'Command <%s> not understood.\nTry function TellAfni_Commands for available commands.\n', Name);
      return;
end

cs.err = 0;
return;

%%%%% Supporting functions below %%%%%%%%%%%%%%%%%
%figure out the window deal
function [win,p1,p1n] = NewCs_GetWin(p1, uwin)

   win = uwin;
   p1n = 0.0;

   if (~isnumeric(p1)),
      if (length(p1) > 1),
         if ((p1(1)) >= 'A' & (p1(1)) <= 'J' & p1(2) == '.'),
            win = p1(1);
            if (~isempty(uwin) & uwin ~= win),
               fprintf(2,...
                  'Warning: Conflict in window specification (%c vs %c).\nChoosing one in parameter (%c).\n',...
                  uwin, win, win);
            end
            p1 = p1(3:length(p1));
         end
      end
      p1n = str2double(p1);
   else 
      p1n = p1;
   end
   if (isempty(win)) win = 'A'; end
return

%check on windowname
function [err] = NewCs_OKwindowname(p1)
   
   lcool = {'axialimage', 'sagittalimage', 'coronalimage',...
            'axialgraph', 'sagittalgraph', 'coronalgraph'};
   if (nargin == 1),   
      switch (p1),
         case lcool,
            err = 0;
            return;
         otherwise,
            err = 1;
            return;
      end
   else
      err = '';
      for (i=1:1:length(lcool)),
         err = sprintf('%s %s\n', err, char(lcool(i)));
      end
      return;
   end

return

%Check on xhaircode
function [err] = NewCs_OKxhaircode(p1)
   
   lcool = {'OFF', 'SINGLE', 'MULTI',...
            'LR_AP', 'LR_IS', 'AP_IS', 'LR', 'AP', 'IS'};
   if (nargin == 1),   
      switch (p1),
         case lcool,
            err = 0;
            return;
         otherwise,
            err = 1;
            return;
      end
   else
      err = '';
      for (i=1:1:length(lcool)),
         err = sprintf('%s %s\n', err, char(lcool(i)));
      end
      return;
   end

return

%Check on resample
function [err] = NewCs_OKresam(p1)
   
   lcool = {'NN', 'Li', 'Cu',...
            'Bk'};
   if (nargin == 1),   
      switch (p1),
         case lcool,
            err = 0;
            return;
         otherwise,
            err = 1;
            return;
      end
   else
      err = '';
      for (i=1:1:length(lcool)),
         err = sprintf('%s %s\n', err, char(lcool(i)));
      end
      return;
   end

return

%Check on panel
function [err] = NewCs_OKpanel(p1)
   
   lcool = {'Define_Overlay', 'Define_Datamode', 'Define_Markers'};
   if (nargin == 1),   
      switch (p1),
         case lcool,
            err = 0;
            return;
         otherwise,
            err = 1;
            return;
      end
   else
      err = '';
      for (i=1:1:length(lcool)),
         err = sprintf('%s %s\n', err, char(lcool(i)));
      end
      return;
   end

return

%check on input file names
function [k] = NewCs_CheckFnames(s)
   [si, s] = strtok(s,' ');
   k = 0;
   if (exist(si,'dir') == 7 || filexist(si)), 
      k = k + 1;
   else
      k = -1;
   end
return
