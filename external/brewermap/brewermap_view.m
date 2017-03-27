function [map,scheme] = brewermap_view(N,scheme)
% An interactive figure for ColorBrewer colormap selection. With demo!
%
% (c) 2017 Stephen Cobeldick
%
% View Cynthia Brewer's ColorBrewer color schemes in a figure.
%
% * Two colorbars give the color scheme in color and grayscale.
% * A button toggles between 3D-cube and 2D-lineplot of the RGB values.
% * A button toggles an endless cycle through the color schemes.
% * A button reverses the colormap.
% * 35 buttons select any ColorBrewer color scheme.
% * Text with the color scheme's type (Diverging/Qualitative/Sequential)
% * Text with the color scheme's number of nodes (defining colors).
%
%%% Syntax:
%  brewermap_view
%  brewermap_view(N)
%  brewermap_view(N,scheme)
%  brewermap_view([],...)
%  brewermap_view({axes/figure handles},...) % see "Adjust External Colormaps"
%  [map,scheme] = brewermap_view(...)
%
% Calling the function with an output argument blocks MATLAB execution until
% the figure is deleted: the final colormap and scheme are then returned.
%
% See also BREWERMAP CUBEHELIX RGBPLOT COLORMAP COLORMAPEDITOR COLORBAR UICONTROL ADDLISTENER
%
%% Adjust Colormaps of Other Figures or Axes %%
%
%%% Example:
%
% load spine
% image(X)
% brewermap_view({gca})
%
% Very useful! Simply provide a cell array of axes or figure handles when
% calling this function, and their colormaps will be updated in real-time:
% note that MATLAB versions <=2010 only support axes handles for this!
%
%% Input and Output Arguments %%
%
%%% Inputs (*=default):
%  N  = NumericScalar, an integer to define the colormap length.
%     = *[], colormap length of one hundred and twenty-eight (128).
%     = {axes/figure handles}, their colormaps will be updated by BREWERMAP_VIEW.
%  scheme = String, a ColorBrewer color scheme name.
%
%%% Outputs (these block execution until the figure is deleted!):
%  map    = NumericMatrix, the colormap defined when the figure is closed.
%  scheme = StringToken, the name of the color scheme given in <map>.
%
% [map,scheme] = brewermap_view(N,scheme)

%% Input Wrangling %%
%
persistent H
%
xtH = {};
% Parse colormap size:
if nargin<1 || isnumeric(N)&&isempty(N)
	N = 128;
elseif iscell(N)&&numel(N)
	ish = all(1==cellfun('prodofsize',N)&cellfun(@ishghandle,N));
	assert(ish,'Input <N> may be a cell array of scalar axes or figure handles.')
	xtH = N;
	N = size(colormap(xtH{1}),1);
else
	assert(isnumeric(N)&&isscalar(N),'Input <N> must be a scalar numeric.')
	assert(isreal(N)&&fix(N)==N&&N>0,'Input <N> must be positive real integer: %g+%gi',N,imag(N))
	N = double(N);
end
%
[mcs,mun,pyt] = brewermap('list');
%
% Parse scheme name:
if nargin<2
	scheme = mcs{1+rem(round(now*1e7),numel(mcs))};
else
	assert(ischar(scheme)&&isrow(scheme),'Second input <scheme> must be a string.')
end
% Check if a reversed colormap was requested:
isR = strncmp(scheme,'*',1);
scheme = scheme(1+isR:end);
%
%% Create Figure %%
%
% LHS and RHS slider bounds/limits, and slider step sizes:
lbd = 1;
rbd = 128;
%
% Define the 3D cube axis order:
xyz = 'RGB';
[~,xyz] = ismember(xyz,'RGB');
%
if isempty(H) || ~ishghandle(H.fig)
	% Check brewermap version:
	ers = 'The function BREWERMAP returned an unexpected %s.';
	assert(all(35==[numel(mcs),numel(mun),numel(pyt)]),ers,'array size')
	tmp = find(any(diff(+char(pyt)),2));
	assert(numel(tmp)==2&&all(tmp==[9;17]),ers,'scheme name sequence')
	%
	% Create a new figure:
	ClBk = struct('bmvChgS',@bmvChgS, 'bmvRevM',@bmvRevM,...
		'bmv2D3D',@bmv2D3D, 'bmvDemo',@bmvDemo, 'bmvSldr',@bmvSldr);
	H = bmvPlot(N,scheme, mcs, lbd, rbd, xyz, ClBk);
end
%
bmvUpDt()
%
if nargout
	waitfor(H.fig);
else
	clear map
end
%
%% Nested Functions %%
%
	function bmvUpDt()
		% Update all graphics objects in the figure.
		%
		% Get ColorBrewer colormap and grayscale equivalent:
		[map,num,typ] = brewermap(N,[char(42*ones(1,isR)),scheme]);
		mag = sum(map*[0.298936;0.587043;0.114021],2);
		%
		% Update colorbar values:
		set(H.cbAx, 'YLim', [0,abs(N)+(N==0)]+0.5);
		set(H.cbIm(1), 'CData',reshape(map,[],1,3))
		set(H.cbIm(2), 'CData',repmat(mag,[1,1,3]))
		%
		% Update 2D line / 3D patch values:
		if  get(H.D2D3, 'Value') % 2D
			set(H.ln2D, 'XData',linspace(0,1,abs(N)));
			set(H.ln2D, {'YData'},num2cell([map,mag],1).');
		else % 3D
			set(H.pt3D,...
				'XData',map(:,xyz(1)),...
				'YData',map(:,xyz(2)),...
				'ZData',map(:,xyz(3)), 'FaceVertexCData',map)
		end
		%
		% Update reverse button:
		set(H.bRev, 'Value',isR)
		%
		% Update warning text:
		str = {typ;sprintf('%d Nodes',num)};
		set(H.warn,'String',str);
		%
		% Update parameter value text:
		set(H.vTxt(1), 'String',sprintf('N = %.0f',N));
		%
		% Update external axes/figure:
		for k = find(cellfun(@ishghandle,xtH))
			colormap(xtH{k},map);
		end
	end
%
	function bmv2D3D(h,~)
		% Switch between 2D-line and 3D-cube representation.
		%
		if get(h,'Value') % 2D
			set(H.ax3D, 'HitTest','off', 'Visible','off')
			set(H.ax2D, 'HitTest','on')
			set(H.pt3D, 'Visible','off')
			set(H.ln2D, 'Visible','on')
		else % 3D
			set(H.ax2D, 'HitTest','off')
			set(H.ax3D, 'HitTest','on', 'Visible','on')
			set(H.ln2D, 'Visible','off')
			set(H.pt3D, 'Visible','on')
		end
		%
		bmvUpDt();
	end
%
	function bmvChgS(~,e)
		% Change the color scheme.
		%
		scheme = get(e.NewValue,'String');
		%
		bmvUpDt()
	end
%
	function bmvRevM(h,~)
		% Reverse the colormap.
		%
		isR = get(h,'Value');
		%
		bmvUpDt()
	end
%
	function bmvSldr(~,~)
		% Update the slider position.
		%
		N = round(get(H.vSld,'Value'));
		%
		bmvUpDt()
	end
%
	function bmvDemo(h,~)
		% Display all ColorBrewer schemes sequentially.
		%
		while ishghandle(h)&&get(h,'Value')
			%
			ids = 1+mod(find(strcmpi(scheme,mcs)),numel(mcs));
			set(H.bGrp,'SelectedObject',H.bEig(ids));
			scheme = mcs{ids};
			%
			bmvUpDt();
			%
			% Faster/slower:
			pause(1.2);
		end
		%
	end
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%brewermap_view
function H = bmvPlot(N,scheme, mcs, lbd, rbd, xyz, ClBk)
% Draw a new figure with RGBplot axes, ColorBar axes, and uicontrol sliders.
%
M = 9; % buttons per column
gap = 0.01; % gaps
bth = 0.04; % demo height
btw = 0.09; % demo width
uih = 0.40; % height of UI control group
cbw = 0.21; % width of both colorbars
axh = 1-uih-2*gap; % axes height
wdt = 1-cbw-2*gap; % axes width
stp = [1,10]; % slider step
%
H.fig = figure('HandleVisibility','callback', 'Color','white',...
	'IntegerHandle','off', 'NumberTitle','off',...
	'Name','ColorBrewer Interactive Scheme Selector');
%
% Add 2D lineplot:
H.ax2D = axes('Parent',H.fig, 'Position',[gap, uih+gap, wdt, axh],...
	'ColorOrder',[1,0,0; 0,1,0; 0,0,1; 0.6,0.6,0.6], 'HitTest','off',...
	'Visible','off', 'XLim',[0,1], 'YLim',[0,1], 'XTick',[], 'YTick',[]);
H.ln2D = line([0,0,0,0;1,1,1,1],[0,0,0,0;1,1,1,1], 'Parent',H.ax2D, 'Visible','off');
%
% Add 3D scatterplot:
H.ax3D = axes('Parent',H.fig, 'OuterPosition',[0, uih, wdt+2*gap, 1-uih],...
	'Visible','on', 'XLim',[0,1], 'YLim',[0,1], 'ZLim',[0,1], 'HitTest','on');
H.pt3D = patch('Parent',H.ax3D, 'XData',[0;1], 'YData',[0;1], 'ZData',[0;1],...
	'Visible','on', 'LineStyle','none', 'FaceColor','none', 'MarkerEdgeColor','none',...
	'Marker','o', 'MarkerFaceColor','flat', 'MarkerSize',10, 'FaceVertexCData',[1,1,0;1,0,1]);
view(H.ax3D,3);
grid(H.ax3D,'on')
lbl = {'Red','Green','Blue'};
xlabel(H.ax3D,lbl{xyz(1)})
ylabel(H.ax3D,lbl{xyz(2)})
zlabel(H.ax3D,lbl{xyz(3)})
%
% Add warning text:
H.warn = text('Parent',H.ax2D, 'Units','normalized', 'Position',[1,1],...
	'HorizontalAlignment','right', 'VerticalAlignment','top', 'Color','k');
%
% Add demo button:
H.demo = uicontrol(H.fig, 'Style','togglebutton', 'Units','normalized',...
	'Position',[gap,uih+gap+0*bth,btw,bth], 'String','Demo',...
	'Max',1, 'Min',0, 'Callback',ClBk.bmvDemo);
% Add 2D/3D button:
H.D2D3 = uicontrol(H.fig, 'Style','togglebutton', 'Units','normalized',...
	'Position',[gap,uih+gap+1*bth,btw,bth], 'String','2D / 3D',...
	'Max',1, 'Min',0, 'Callback',ClBk.bmv2D3D);
% Add reverse button:
H.bRev = uicontrol(H.fig, 'Style','togglebutton', 'Units','normalized',...
	'Position',[gap,uih+gap+2*bth,btw,bth], 'String','Reverse',...
	'Max',1, 'Min',0, 'Callback',ClBk.bmvRevM);
%
% Add colorbars:
C = reshape([1,1,1],1,[],3);
H.cbAx(1) = axes('Parent',H.fig, 'Visible','off', 'Units','normalized',...
	'Position',[1-cbw/1,gap,cbw/2-gap,1-2*gap], 'YLim',[0.5,1.5],...
	'YDir','reverse', 'HitTest','off');
H.cbAx(2) = axes('Parent',H.fig, 'Visible','off', 'Units','normalized',...
	'Position',[1-cbw/2,gap,cbw/2-gap,1-2*gap], 'YLim',[0.5,1.5],...
	'YDir','reverse', 'HitTest','off');
H.cbIm(1) = image('Parent',H.cbAx(1), 'CData',C);
H.cbIm(2) = image('Parent',H.cbAx(2), 'CData',C);
%
% Add parameter slider, listener, and corresponding text:
sv = max(lbd,min(rbd,N));
H.vTxt = uicontrol(H.fig,'Style','text', 'Units','normalized',...
	'Position',[gap,uih-bth,btw,bth], 'String','X');
H.vSld = uicontrol(H.fig,'Style','slider', 'Units','normalized',...
	'Position',[gap,gap,btw,uih-bth], 'Min',lbd(1), 'Max',rbd(1),...
	'SliderStep',stp(1,:)/(rbd(1)-lbd(1)), 'Value',sv(1));
addlistener(H.vSld, 'Value', 'PostSet',ClBk.bmvSldr);
%
% Add scheme button group:
H.bGrp = uibuttongroup('Parent',H.fig, 'BorderType','none', 'Units','normalized',...
	'BackgroundColor','white', 'Position',[2*gap+btw,gap,wdt-btw-gap,uih-gap]);
% Determine button locations:
Z = 1:numel(mcs);
Z = Z+(Z>17);
C = (ceil(Z/M)-1)/4;
R = (M-1-mod(Z-1,M))/M;
% Add scheme buttons to group:
for k = numel(mcs):-1:1
	H.bEig(k) = uicontrol('Parent',H.bGrp, 'Style','Toggle', 'String',mcs{k},...
		'Unit','normalized', 'Position',[C(k),R(k),1/4,1/M]);
end
set(H.bGrp,'SelectedObject',H.bEig(strcmpi(scheme,mcs)));
set(H.bGrp,'SelectionChangeFcn',ClBk.bmvChgS);
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%bmvPlot
%
% Copyright (c) 2017 Stephen Cobeldick
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
% http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and limitations under the License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%license