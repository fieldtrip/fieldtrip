function ht = plotsign2 (h,s,Opt);
%	ht = plotsign2 (h,s,[Opt]);
%	h is the figure handle, if you need to specify a subplot within the figure
%    pass the figure handle and the subplot number in h, like: [1 213] 
%    or [1 2 1 3]
%	s is the text string, if s is not supplied, a default string is used
%    the default is : machine : date \n path:filename of calling matlab script
%
% Opt is an optional structure with the following optional fields
%  .Place is an approximate position for placing the text.
%     .Place is a three character string where the first character 
%      specifies whether the placement is relative to the Figure 'F'
%      or the plot 'P'. The second character 
%      specifies the Y position, and the second character the 
%      X position. For example 'FTC' means Figure Top Center.
%      Your options for the relative positioning are :
%              F for figure and P for plot
%      Your options for the Y position (first letter) are:
%              T or C or B (Top or Center or Bottom) 
%      Your options for the X position (second letter) are:
%       			L or C or R (Left or Center or Right)
%      text alignment is done automatically, unless you specify a .Align option.
%     The default is 'FBR'  
%  .Pos (2x1 or 3x1 vector) pins the X Y postion (relative the X Y axis) of the plot
%    .Pos overrides .Place.
%   NOTE : The units of .Pos are relative to the plot 
%  .Font is the fontsize to use (a number or string). (like 8, or 10 or 12 etc).
%     The default is 2 numbers less the current axis on the figure (usually 8)
%     You can also specify 'l' or 's' for large or small (2 less or 2 more than
%     current axis). 'n' for a font equal to that of current axis.
%  .FontName is an optional string specifying a font name
%  .Color : choose colors that you normally specify in text function, like 'r'
%       the default is 'k'. You can also specify an rgb vector 
%       (values between 1 and 0)
%  .ReverseVideo: 'y' puts a background behind text to ensure it is readable.
%                 default is to do nothing.
%  .HAlign : string specifying the horizontal alignment of text, options are
%       'left', 'center', 'right', '' for default
%  .VAlign : string specifying the vertical alignment of text, options are 
%       'top', 'middle', 'bottom', '' for default
%
%  Once you see the object, you can move it with the mouse
%  ht is the handle to the text object
%
%  example:
%  t=0:0.1:50; Y = sin(0.5.*t); cf = figure(1); plot (t,Y); 
%  Ssign = sprintf ('Z.S.Saad %s\n image:%s/%s', date, pwd, mfilename );
%  plotsign2(cf,Ssign);
%  or 
%  Opt.Place = 'FTC'; Opt.Color = 'r'; Opt.Font = 16;
%  plotsign2(cf,Ssign,Opt);
%  or
%  Opt.Pos = [t(400) Y(400)];Opt.Color = 'b'; Opt.Font = 12;
%  plotsign2(cf,Ssign,Opt);
%
% See also
%   MoveFigObject
%   HistoryTrace
%		Ziad Saad Dec 1 97, latest fix Wed May 05 21:41:30 CDT 1999
%

ht = 0; %initialize in case you return early

figure (h(1));
	if (length(h) > 1),
		if (length(h) == 2), subplot (h(2));  %in case you want to specify a subplot
			elseif (length(h) == 4), subplot (h(2), h(3), h(4));
			else err = ErrEval(FuncName,'Err_Bad size for input parameter h');
				return;
		end
	end
	

if (nargin == 1),
	[ST,I] = dbstack; 
	if (length(ST) > 1), %use the calling function name as a signature
		[jnk1, CallPath, CallFname] = GetPath(ST(2).name);
      if (CallPath(1)=='.') CallPath = pwd; end
	else %this function is called directly from command line
		CallPath = pwd;
		CallFname = mfilename;
      if (strcmp(CallFname,'plotsign2')), CallFname='ShellPrompt'; end;
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
			smach = deblank(smach); smach = smach(1:length(smach)-1);
	end
	c=datevec(now);
	s = sprintf ('Ziad Saad SSCC/NIMH/NIH\n%s : %s %s:%s:%s\n%s/%s', smach, date,...
				 pad_strn(num2str(c(4)), '0', 2, 1),...
				 pad_strn(num2str(c(5)), '0', 2, 1),...
				 pad_strn(num2str(round(c(6))), '0', 2, 1),...
				 CallPath, CallFname);
	Opt.Place = 'FBR';
elseif (nargin == 2),
	Opt.Place = 'FBR';
end

if (~isfield(Opt,'Place') | isempty(Opt.Place)),
	Opt.Place = 'FBR';
end	
if (~isfield(Opt,'HAlign')),
	Opt.HAlign = '';
end	
if (~isfield(Opt,'VAlign')),
	Opt.VAlign = '';
end	
if (~isfield(Opt,'Position')),
	Opt.Position = [];
end	

if (~isfield(Opt,'Pos') | isempty(Opt.Pos)),
	spec = 0;
else
	if (length(Opt.Pos) ~= 2 & length(Opt.Pos) ~= 3),
		ErrEval(mfilename,'Err_Bad size for Opt.Pos parameter');
		return;
	end
	spec = 1;
end

if (~isfield(Opt,'Color') | isempty(Opt.Color)),
	Opt.Color = 'k';
end

if (~isfield(Opt,'ReverseVideo') | isempty(Opt.ReverseVideo)),
   Opt.ReverseVideo = '';
end

if (Opt.ReverseVideo == 'y'), 
   if (isnumeric(Opt.Color)) Opt.BackgroundColor = 1 - Opt.Color;
   else Opt.BackgroundColor = [0 0 0];
   end
else
   Opt.BackgroundColor = [];
end

% Fontsize for text
	if (~isfield(Opt,'Font') | isempty(Opt.Font)),
		Opt.Font = 's';
	end
	if (ischar(Opt.Font)),
		if (eq_str(Opt.Font,'l')),
			fs = get(gcf,'defaultaxesfontsize')+2;
		elseif (eq_str(Opt.Font,'s')),
			fs = get(gcf,'defaultaxesfontsize')-2;
		elseif (eq_str(Opt.Font,'n')),
			fs = get(gcf,'defaultaxesfontsize');	
		end
	else
		fs = Opt.Font;
	end
	


%place text anywhere	
	ht=text(0,0,s, 'color',Opt.Color);
   if (~isempty(Opt.ReverseVideo)), 
      col = get(ht,'color');
      if (0),
         col = 1-col; %not nice
      else,
         [jjj,iii1] = max(col);
         col(iii1) = 1- col(iii1);
         [jjj,iii2] = min(col);
         col(iii2) = min([col(iii2)+0.5,1]);
         iii3=setdiff([1 2 3],[iii1 iii2]);
         col(iii3) = 1 - col(iii3);
      end
      set(ht,'BackgroundColor',col);
   end
%set fonts
	if (isfield(Opt,'FontName') & ~isempty(Opt.FontName)),
		set(ht,'FontName', Opt.FontName);
	end
%set font size
	set(ht,'fontsize',fs); 
%set interpreter
   if (~isfield(Opt,'Interpreter') | isempty(Opt.Interpreter)),
	   Opt.Interpreter = 'none';
   end
   
   set(ht,'Interpreter',Opt.Interpreter);

%store current units
	tmpunt = get(ht,'Units');

	
%set position and alignment as fit	
	if (spec), %user specified location
		set(ht,'Units','data');
		set(ht,'Position',Opt.Pos);
		
		%override auto-alignment if specified
			if (~isempty(Opt.VAlign)),	AlgnV = Opt.VAlign;	
				else	AlgnV = 'middle'; end
			if (~isempty(Opt.HAlign)),	AlgnH = Opt.HAlign;	
				else	AlgnH = 'center'; end
		
		set(ht,'horizontalalignment',AlgnH);
		set(ht,'verticalalignment',AlgnV);
	else %auto placement
		%get the size in normalized units of the axes
			UntAxes = get(gca,'Units');
			if (~eq_str(UntAxes,'normalized')),
				ErrEval(mfilename,'Err_Sorry, must have normalized units');
				return;
			end
			
			AxSize = get(gca,'Position');
			%if the axes does not fill the figure window, find out what would fill it
			Mult = [1./AxSize(3) 1./AxSize(4)];
		
		Vpos = [0 0];
		Unt = 'normalized';
		switch Opt.Place(1), %relative to figure or plot
			case 'F', %figure
				DoF = 1;
			case 'P', %plot
				DoF = 0;
			otherwise,
				ErrEval(mfilename,'Err_Cannot interpret Opt.Place(1)');
				return;
		end
		switch Opt.Place(2), %the Y axis
			case 'B', %bottom
				if (DoF), 
					Vpos(2) = (-AxSize(2) + 0.02) .* Mult(2);
					AlgnV = 'bottom';
				else	Vpos(2) = (-0.05); 
					AlgnV = 'top';
				end
				
			case 'C', %center
				if (DoF), Vpos(2) = (-AxSize(2) + 0.5 + 0.02) .* Mult(2);
				else	Vpos(2) = 0.5;	end
				AlgnV = 'middle';
				
			case 'T',
				if (DoF), Vpos(2) = (-AxSize(2) - 0.02) .* Mult(2) + Mult(2)  ;
				else	Vpos(2) = 0.95; end
				AlgnV = 'top';
			otherwise,
				ErrEval(mfilename,'Err_Cannot interpret Opt.Place(2)');
				return;
		end
		switch Opt.Place(3), %the X axis
			case 'L', %bottom
				if (DoF), Vpos(1) = (-AxSize(1) + 0.02) .* Mult(1);
				else	Vpos(1) = (-0.05); end
				AlgnH = 'left';
			case 'C', %center
				if (DoF), Vpos(1) = (-AxSize(1) + 0.5 + 0.02).* Mult(1);
				else	Vpos(1) = 0.5; end
				AlgnH = 'center';
			case 'R',
				if (DoF), Vpos(1) = (-AxSize(1) - 0.02) .* Mult(1) + Mult(1);
				else Vpos(1) = 0.95; end
				AlgnH = 'right';
			otherwise,
				ErrEval(mfilename,'Err_Cannot interpret Opt.Place(3)');
				return;
		end
		
		%override auto-alignment if specified
			if (~isempty(Opt.VAlign)),	AlgnV = Opt.VAlign;	end
			if (~isempty(Opt.HAlign)),	AlgnH = Opt.HAlign;	end
		
		set(ht,'horizontalalignment',AlgnH);
		set(ht,'verticalalignment',AlgnV);
		set(ht,'Units',Unt);
		set(ht,'Position',Vpos);
	end

set(ht,'Units',tmpunt);

%Now you need to store the handle to the object and give it a number
%By giving each object a number, you can store multiple objects on the same
%figure and manipulate them all individually

	%get the UserData of the figure
		ud = get (gcf,'UserData');
	%make sure there has been handle stored previously, otherwise add some
		if (~isfield(ud,'MaxHandle') | isempty(ud.MaxHandle)),
			ud.MaxHandle = 0;
		end
	%Increment the maximum number of handles
		ud.MaxHandle = ud.MaxHandle + 1;
	%Now add the handle in the MoveHandle vector
		ud.MoveHandle(ud.MaxHandle) = ht;
	%store the userdata back in the figure
		set(gcf,'UserData',ud);
	%create the button callback string
		stmp = sprintf('MoveFigObject(%g)',ud.MaxHandle);
	%set the button down function	
		set(ht,'ButtonDownFcn',stmp);

return;
