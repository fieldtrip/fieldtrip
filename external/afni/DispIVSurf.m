function [err,FigHndls] = DispIVSurf (NodeList,FaceSetList,UsedNodeSet,ITvect,FigHandle,DispOpt)
%
% [err,FigHndls] = DispIVSurf (NodeList,FaceSetList,UsedNodeSet,ITvect,FigHandle,[DispOpt])
%
% NodeList is a Nx3 matrix containing the node index and the
%  XYZ coordinates of each of the N nodes that form the surface model.
%  You can also send in a row major vector of Nx3 elements.
%  IT USED TO BE (Pre Oct 21 04):
%     NodeList is a Nx4 matrix containing the node index and the
%        XYZ coordinates of each of the N nodes that form the surface model.
%     but that is no longer useful although it is still supported
%
% FaceSetList is an Mx3 matrix, where M is the total number of
%  facesets forming the surface. Each row holds a triplet of
%  node indices that form the triangular faceset. You can also send
%  in a row major vector
%
% ITvect is the intensity (.it file) vector that assigns to each node
%  in NodeList an intensity. If you have no such info, just pass one value
%  which will be assigned to all nodes.
%  For unknown reasons, we end up with some meshes having an extra node,
%  Normally, the function should refuse to display this, but for now
%  it will add a 0 to the end of ITvect (or remove the last point
%   if ITvect is too big) to make ITvect match in size with NodeList.
%   That's what the Inventor tools that Seth wrote did anyway. The function
%   will issue a warning when that happens, and hopefully we'll get ridd
%   of this mess.
%
% UsedNodeSet is a vector that holds the indices of the nodes
%  that take part of the FaceSetList. It's only required with the
%  stupid 3D plot display. You can send in [] if you don't want to use plot3
%
% FigHandle is the number of the figure to display in, if you set it to zero
%  the function chooses the next available one.
%  to specify a subplot, use two values like [1 212] meaning figure 1, subplot 212
%
% DispOpt is an optional structure specifying some display parameters and options
% DispOpt has the following fields
%  .GraphType  GraphType is a string that specifies the display type required
%             choose from the following (match the case for the strings):
%             Mesh : Displays the mesh with colors specified by ITvect
%             Surf : Displays the patched mesh with colors specified by ITvect
%             The default is Surf.
%  .Shade  choose from 'flat' or 'interp' when not using .OpenGL
%          .OpenGL forces the use of 'interp'. This field is meaningless
%          for Mesh type display. default is 'interp'
%
%  .AxesSize  (specify [left bottom width height])
%                 (default is [])
%  .TrueColor flag (0/1) set to 1 if you are using true color.
%    (ie ITvect is an Nx3 matrix with each row specifying a true color (RGB, 0-1).
%    (actually, you should fix this to have the option of turning ITvect
%     into a true color vector even if it is sent as an Nx1. It's easy, using
%     ScaleToMap function)
%    default is 0.
%  .ColBarSize (specify [left bottom width height])
%                 (default is [])
%  .ColMap : An Nx3 matrix specifying N different colors. (just like what
%            you get with cmap = colormap; set to [] for using the current colour map
%                 (default is [])
%  .ColMapRange : Index determining the range of colors to use in ColMap like [1 64]
%            to go from the first colour to 64. Set to [] if you want to use the whole colourmap,
%                 (default is [])
%  .MaskValue : 1x1 If you send an ITvect vector (not just one value), you can specify
%      a value to be a mask.(default is [])
%  .MaskColor : 1x3 vector. That is a color assigned to values in ITvect that should be
%              masked (.MaskValue) . This color is added to the bottom of the color map.
%              This value has no effect unless .MaskValue is used. (default is [])
%  .DataRange : 1x2 vector. If ITvect contains data to be pseudocolored,
%               you might not want to have the highest and lowest values
%               to be mapped to the last and first colours in the colormap
%               If you do so, then for different ITvect with different values,
%               The colors will look different. So, to keep the color significance
%               constant across data sets you can specify [min max] in
%               .DataRange so that the first color is always min and the
%               last color is always max.
%               NOTE: The limits in DataRange are overrun by values in the
%               Data set that are lower than min or higher than max.
%               (default is [])
%  .AxesUnits (choose from : pixels normalized inches centimeters points)(default is '')
%  .Xlim : [min max] limits for x axis
%  .Ylim : [min max]
%  .Zlim : [min max]
%     if all the above fields are empty or non specified, then matlab's automatic settings are used
%  .Hold      (choose from: on off ' ' <- a space character for toggling )(default is '')
%          This is applied after the figure is redrawn
%  .Zoom      (choose from: (factor) on off reset out xon yon ' ' <- a space character
%            see help zoom for more info) (default is '')
%  .Dim       Choose from 3D or 2D, in 3D you can't have zoom on  (default is '2D')
%  .View  if specified, this option supersedes .Dim
%      pass along any parameters that are normally passed to matlab's view function
%
%  .DispMode (choose from 'New' or 'Add').
%          If you choose New, then old figure is cleared and colorbars are updated.
%          If you choose Add, then new figure is placed on top of old one. Colorbars are not modified
%          default is 'New'
%  .OpenGL flag 0/1, uses OpenGL for rendering. Default is 0
%  .ShowEdge flag 0/1, shows the edges of the facesets. Default is 1
%  .Grid (flag vector  [0/1 0/1 0/1]), puts an X Y Z grid on the map, default is [1 1 1]
%  .verbose (0/1) default is 1
%
% If you do not wish to set some of the previous properties, you set their values either to
% '' or [] depending on their types.
%
%  see ExtractIVdata.m for more info
%
%  err = 0 no problem
%  err = 1 Error somewhere
%
% FigHndls is a structure containing handles to the various relevant obejects in the figure window
%  The fields of FigHndls are :
%     Fig (The whole window)
%     Patch (The patches handle)
%		Axes (The axes handle)
%		ColBar (The colorbar handle)
%
%		Ziad Saad March 1 98
%

%initializations
FuncName = 'DispIVSurf';
err = 1;
FigHndls = [];

%setup defualts
	if (nargin < 6),
		DispOpt = [];
	end
	
	if (~isfield(DispOpt,'verbose') | isempty(DispOpt.verbose)),
		DispOpt.verbose = 1;
	end

	if (~isfield(DispOpt,'TrueColor') | isempty(DispOpt.TrueColor)),
		DispOpt.TrueColor = 0;
	end

	if (~isfield(DispOpt,'OpenGL') | isempty(DispOpt.OpenGL)),
		DispOpt.OpenGL = 0;
	end

	if (~isfield(DispOpt,'ShowEdge') | isempty(DispOpt.ShowEdge)),
		DispOpt.ShowEdge = 1;
	end

	if (~isfield(DispOpt,'GraphType') | isempty(DispOpt.GraphType)),
		DispOpt.GraphType = 'Surf';
	end

	if (~isfield(DispOpt,'Shade') | isempty(DispOpt.Shade)),
		DispOpt.Shade = 'interp';
	end

	if (~isfield(DispOpt,'ColMap')),
		DispOpt.ColMap = [];
	end

	if(~isfield(DispOpt,'Zoom')),
		DispOpt.Zoom = '';
	end

	if(~isfield(DispOpt,'Dim') | isempty(DispOpt.Dim)),
		DispOpt.Dim = '2D';
	end

	if (isempty(DispOpt.ColMap)),
		DispOpt.ColMap = colormap;
		DispOpt.ColMapRange = [1 length(DispOpt.ColMap)];
	end

	if (isempty(DispOpt.ColMapRange)),
		DispOpt.ColMapRange = [1 length(DispOpt.ColMap)];
	end
	
	if (~isfield(DispOpt,'MaskValue')),
		DispOpt.MaskValue = [];
	end

	if (~isfield(DispOpt,'DispMode' ) | isempty(DispOpt.DispMode)),
		DispOpt.DispMode = 'New';
	end


	if (~isfield(DispOpt,'MaskColor' ) | isempty(DispOpt.MaskColor)),
		DispOpt.MaskColor = [0 0 0];
	end

	if (~isfield(DispOpt,'DataRange' ) | isempty(DispOpt.DataRange)),
		DispOpt.DataRange = [];
	end

	if (~isfield(DispOpt,'Grid') | isempty(DispOpt.Grid)),
		DispOpt.Grid = [1 1 1];
	end
		
	if (~isfield(DispOpt,'View')),
		DispOpt.View = [];
	end

%check for inconsistent input
	if (is_row(DispOpt.Grid) < 0 | length(DispOpt.Grid) ~= 3),
		err = ErrEval(FuncName,'Err_Bad .Grid vector size');	return
	end
	
	if (~eq_str(DispOpt.Shade,'flat') & ~eq_str(DispOpt.Shade,'none') & ~eq_str(DispOpt.Shade,'interp')),
		err = ErrEval(FuncName,'Err_Shade option specified is not supported.');	return
	end	
	
	if (eq_str(DispOpt.Shade,'flat') & DispOpt.OpenGL),
		err = ErrEval(FuncName,'Err_Cannot use flat shading with OpenGL option.');	return;
	end
	
	if (~eq_str(DispOpt.GraphType,'Mesh') & ~eq_str(DispOpt.GraphType,'Surf')),
		err = ErrEval(FuncName,'Err_Bad GraphType, only Mesh and Surf are allowed.');	return
	end
	
	if (~eq_str(DispOpt.DispMode,'New') & ~eq_str(DispOpt.DispMode,'Add')),
		err = ErrEval(FuncName,'Err_DispMode field can only be New or Add.');	return;	
	end

	if (~isempty(DispOpt.Zoom) & eq_str(DispOpt,'3D')),
		err = ErrEval('DispIVSurf','Err_You can not use zoom option with 3D displays'); return
	end

	if (~(size(ITvect,2) == 1 & DispOpt.TrueColor == 0) & ~(size(ITvect,2) == 3 & DispOpt.TrueColor == 1)),
		err = ErrEval(FuncName,'Err_ITvect size and .TrueColor specified do not match.');return
	end

	if (DispOpt.TrueColor & DispOpt.OpenGL),
		err = ErrEval(FuncName,'Err_Cannot use .TrueColor and OpenGL.');	return;
	end
	
   if (size(NodeList,2) == 1),
      NodeList = reshape (NodeList, 3, length(NodeList)/3); NodeList = NodeList';
   end
   if (size(FaceSetList,2) == 1),
      FaceSetList = reshape (FaceSetList, 3, length(FaceSetList)/3); FaceSetList = FaceSetList';
   end
	if (size(NodeList,2) ~= 4 & size(NodeList,2) ~= 3),
		err = ErrEval(FuncName,'Err_Bad Size for NodeList.');	return
	end	
   if (size(NodeList,2) == 4),
      fprintf(2,'Warning %s: First column is being trimmed, should not affect results ZSS Oct 21 04.\n');
      NodeList = NodeList(:,2:4);
   end
	
	if (size(FaceSetList,2) ~= 3),
		err = ErrEval(FuncName,'Err_Bad Size for FaceSetList.');	return
	end	
	
	%If ITvect is a 1x1 or a 1x3 value than you need to create a vector that's all constant
	if (size(ITvect,1) == 1),
		ITvect = ITvect .* ones (size(NodeList,1),size(ITvect,2));
	end
	
	%make sure that ITvect size matches that of the number of rows (Nodes) in NodeList
	Nit = size(ITvect,1);
	Nnd = size(NodeList,1);
	if (abs(Nit - Nnd) > 1),
		err= ErrEval(FuncName,'Err_Bad size match (more than one) between ITvect and NodeList\n\a');	return;
	elseif ((Nit - Nnd) == 1),
		%remove one point from ITvect
		ITvect = ITvect(1:Nit-1);
		ErrEval(FuncName,'Wrn_ITvect and Nodelist size missmatch, removed last point from ITvect.\n');
	elseif ((Nit - Nnd) == -1),
		%add a point to ITvect
		ITvect = [ITvect ; DispOpt.MaskValue];
		ErrEval(FuncName,'Wrn_ITvect and Nodelist size missmatch, add an empty point to ITvect.\n');
	end
	
	
	
%want special lighting ?
	DoLight = 0;

%Setup the figure handle
	if (FigHandle == 0),
		FigHndls.Fig(1) = figure;
	else
		FigHndls.Fig = FigHandle;		
	end

%decide on renderer
	if (DispOpt.OpenGL),
		if (DispOpt.verbose),	fprintf (1,'%s verbose: OpenGL mode, clearing axis.\n', FuncName);	end
		set(FigHndls.Fig(1),'Renderer','OpenGL');
		%You must clear the axis, or switching from one data set to the other
		%won't work right ???
		cla
	else
		if (DispOpt.verbose),	fprintf (1,'%s verbose: z-buffer mode.\n', FuncName);	end
		set(FigHndls.Fig(1),'Renderer','zbuffer');
	end

%If using indexed colors, scale ITvect to match color map
	if (~DispOpt.TrueColor), %We need to scale ITvect to fit the color map
				if (DispOpt.verbose),	fprintf (1,'%s verbose: Scaling ITvect to fit colormap\n', FuncName);	end
				%scale it to fit color map
				Opt.Format = 1;
				Opt.MaskValue = DispOpt.MaskValue;
				Opt.MaskColor = DispOpt.MaskColor;
				Opt.DataRange = DispOpt.DataRange;

				%extract the colormap you'll scale data to
				ColMap = DispOpt.ColMap(DispOpt.ColMapRange(1):DispOpt.ColMapRange(2),:);

				%scale data to fit the map
				[err,ITvect,ColMap] = ScaleToMap (ITvect,ColMap,Opt);

				%Now you must add the new color map (has a mask color added) to what's left over from the full map
				if (DispOpt.ColMapRange(1) > 1),
					NewFullMap = DispOpt.ColMap(1:DispOpt.ColMapRange(1)-1,:);
					NewFullMap = [NewFullMap;ColMap];
				else	NewFullMap = ColMap;		end

				if (DispOpt.ColMapRange(2) > size(DispOpt.ColMap,2)),
					NewFullMap2 = DispOpt.ColMap(DispOpt.ColMapRange(2)+1:size(DispOpt.ColMap,2),:);
					NewFullMap = [NewFullMap;NewFullMap2];
				end

				DispOpt.ColMap = NewFullMap;
	end			

%create the patches
	figure (FigHndls.Fig(1)); if (length (FigHndls.Fig) > 1), subplot ( FigHndls.Fig(2)); end

%decide if you should put this on top of old figure or not
	if (eq_str(DispOpt.DispMode,'New')),	
			hold off;
			cla;
			if (DispOpt.verbose),	fprintf (1,'%s verbose: Set hold to off\n', FuncName);	end
		else
			hold on;	
			if (DispOpt.verbose),	fprintf (1,'%s verbose: Set hold to on\n', FuncName);	end
	end

%display the data
	if (DispOpt.verbose),	fprintf (1,'%s verbose: Drawing Patches\n', FuncName);	end
	FigHndls.Patch = patch(	'Faces',FaceSetList,...
									'Vertices',NodeList,...
									'FaceVertexCData',ITvect);

%Make sure that The idexing is direct
	set (FigHndls.Patch,'CDataMapping','direct');

%decide on FaceColor based on whether it's a mesh or a surface
	if (eq_str (DispOpt.GraphType,'Surf')),
		set (FigHndls.Patch,'FaceColor',DispOpt.Shade);
		if (DispOpt.verbose),	fprintf (1,'%s verbose: Applied %s for Shade.\n', FuncName,DispOpt.Shade);	end
	elseif (eq_str (DispOpt.GraphType,'Mesh')),
		set (FigHndls.Patch,'FaceColor','none');
	end

%decide on EdgeColors
	if (DispOpt.ShowEdge),
		if (DispOpt.verbose),	fprintf (1,'%s verbose: Set Edgecolor to black\n', FuncName);	end
		set (FigHndls.Patch,'EdgeColor',[0 0 0]);
	else
		if (DispOpt.verbose),	fprintf (1,'%s verbose: Set Edgecolor to none\n', FuncName);	end
		set (FigHndls.Patch,'EdgeColor','none');
	end

%decide on the view based on the dimention of the display	

if (isempty(DispOpt.View)),	
	switch DispOpt.Dim,
		case '2D',
			if (DispOpt.verbose),	fprintf (1,'%s verbose: Set view to 2D\n', FuncName);	end
			view (2);
			axis xy;
		case '3D',
			if (DispOpt.verbose),	fprintf (1,'%s verbose: Set view to 3D\n', FuncName);	end
			view(-37.5,30);
			axis vis3d;
		otherwise,
			ErrEval('DispIVSurf','Wrn_DispOpt.Dim option not understood');
	end	
else
	view(DispOpt.View);
end
%Put up the colorbar and size it appropriately
if (~DispOpt.TrueColor & eq_str(DispOpt.DispMode,'New')),
	if (DispOpt.verbose),	fprintf (1,'%s verbose: Displaying colormap\n', FuncName);	end
	colormap (DispOpt.ColMap)
	colorbar;
	FigHndls.ColBar = colorbar;
	if (~isfield(DispOpt,'ColBarSize') | isempty(DispOpt.ColBarSize)),
		DispOpt.ColBarSize = get(FigHndls.ColBar,'Position');	end
	set (FigHndls.ColBar,'Position',DispOpt.ColBarSize);
else
	FigHndls.ColBar = [];
end

%Check for axes size and units
if (DispOpt.verbose),	fprintf (1,'%s verbose: Seting Axes properties\n', FuncName);	end
	
	FigHndls.Axes = gca;
	if (~isfield(DispOpt,'AxesSize') | isempty(DispOpt.AxesSize)),	
		DispOpt.AxesSize = get(FigHndls.Axes,'Position');	end
	if (~isfield(DispOpt,'AxesUnits') | isempty(DispOpt.AxesUnits)),	
		DispOpt.AxesUnits = get(FigHndls.Axes,'Units');	end
	set (FigHndls.Axes,'Units',DispOpt.AxesUnits);
	set (FigHndls.Axes,'Position',DispOpt.AxesSize);

	defx = 1;
	if (isfield(DispOpt,'Xlim')),
		if (~isempty(DispOpt.AxesSize)),
			DispOpt.Xlim = DispOpt.Xlim(:)';
			if (length(DispOpt.Xlim) ~= 2),
				err = ErrEval(FuncName,'Err_Bad size for DispOpt.Xlim'); return
			end
			set(FigHndls.Axes,'XLim',DispOpt.Xlim);
			defx = 0;
		end
	end
	defy = 1;
	if (isfield(DispOpt,'Ylim')),
		if (~isempty(DispOpt.AxesSize)),
			DispOpt.Ylim = DispOpt.Ylim(:)';
			if (length(DispOpt.Ylim) ~= 2),
				err = ErrEval(FuncName,'Err_Bad size for DispOpt.Ylim'); return
			end
			set(FigHndls.Axes,'YLim',DispOpt.Ylim);
			defy = 0;
		end
	end
	if (eq_str(DispOpt.Dim,'3D')),
		if (isfield(DispOpt,'Zlim')),
			if (~isempty(DispOpt.AxesSize)),
				DispOpt.Zlim = DispOpt.Zlim(:)';
				if (length(DispOpt.Zlim) ~= 2),
					err = ErrEval(FuncName,'Err_Bad size for DispOpt.Zlim'); return
				end
				set(FigHndls.Axes,'ZLim',DispOpt.Zlim);
			end
		end
	end
	
	if (defx & defy),	%make axis equal
		axis equal;
	end
	
	if (DispOpt.Grid(1))
		set (FigHndls.Axes,'XGrid','on');
	else	 set (FigHndls.Axes,'XGrid','off');	 end
	if (DispOpt.Grid(2))
		set (FigHndls.Axes,'YGrid','on');
	else	 set (FigHndls.Axes,'YGrid','off');	 end
	if (DispOpt.Grid(3))
		set (FigHndls.Axes,'ZGrid','on');
	else	 set (FigHndls.Axes,'ZGrid','off');	end

%put the Label on the axes
xlabel ('X');	ylabel('Y');	zlabel('Zee');

%Set the final Hold status to whatever is needed
	if (isfield(DispOpt,'Hold') & ~isempty(DispOpt.Hold)),	
		eval(['hold ' DispOpt.Hold ';']);	end

%Turn zooming on if needed
if (isfield(DispOpt,'Zoom') & ~isempty(DispOpt.Zoom)),	
	eval(['zoom ' DispOpt.Zoom ';']);	end
			
%set the lighting
			set(FigHndls.Patch,'FaceLighting', 'gouraud', 'AmbientStrength', 1);
			if (0),
			xl = get(FigHndls.Axes, 'Xlim');
			yl = get(FigHndls.Axes, 'Ylim');
			zl = get(FigHndls.Axes, 'Zlim');
			lightHandl1 = light('position', [10.*xl(1) mean(yl) mean(zl)], 'Style','infinite', 'Color', [0.5 0.5 0.5]);
			lightHandl2 = light('position', [10.*xl(2) mean(yl) mean(zl)], 'Style','infinite', 'Color', [0.5 0.5 0.5]);
			lightHandl3 = light('position', [mean(xl) 10.*yl(1) mean(zl)], 'Style','infinite', 'Color', [0.5 0.5 0.5]);
			lightHandl4 = light('position', [mean(xl) 10.*yl(2) mean(zl)], 'Style','infinite', 'Color', [0.5 0.5 0.5]);
			lightHandl5 = light('position', [mean(xl) mean(yl) 10.*zl(1)], 'Style','infinite', 'Color', [0.5 0.5 0.5]);
			lightHandl6 = light('position', [mean(xl) mean(yl) 10.*zl(2)], 'Style','infinite', 'Color', [0.5 0.5 0.5]);
			else %second method
				%camlight headlight
				camlight right
				camlight left
			end
if (DispOpt.verbose),	fprintf (1,'%s verbose: Done\n', FuncName);	end

err = 0;
return;
