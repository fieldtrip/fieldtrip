function varargout = plot2svg(param1,id,pixelfiletype)
%  Matlab to SVG converter
%  Prelinary version supporting 3D plots as well
%
%  Usage: plot2svg(filename,graphic handle,pixelfiletype)
%                  optional     optional     optional
%         or
%
%         plot2svg(figuresize,graphic handle,pixelfiletype)
%                   optional     optional      optional
%
%         pixelfiletype = 'png' (default), 'jpg'
%
%  Juerg Schwizer 23-Oct-2005
%  See http://www.zhinst.com/blogs/schwizer/ to get more informations
%
%  07.06.2005 - Bugfix axxindex (Index exceeds matrix dimensions)
%  19.09.2005 - Added possibility to select output format of pixel graphics
%  23.10.2005 - Bugfix cell array strings (added by Bill)
%               Handling of 'hggroups' and improved grouping of objects
%               Improved handling of pixel images (indexed and true color pictures)
%  23.10.2005 - Switched default pixelfromat to 'png'
%  07.11.2005 - Added handling of hidden axes for annotations (added by Bill)
%  03.12.2005 - Bugfix of viewBox to make Firefox 1.5 working
%  04.12.2005 - Improved handling of exponent values for log-plots
%               Improved markers
%  09.12.2005 - Bugfix '<' '>' '?' '"'
%  22.12.2005 - Implementation of preliminary 3D version
%               Clipping
%               Minor tick marks
%  22.01.2005 - Removed unused 'end'
%  29.10.2006 - Bugfix '°','±','µ','²','³','¼''½','¾','©''®'
%  17-04-2007 - Bugfix 'projection' in hggroup and hgtransform
%  27-01-2008 - Added Octave functionality (thanks to Jakob Malm)
%               Bugfixe cdatamapping (thanks to Tom)
%               Bugfix image data writing (thanks to Tom)
%               Patches includes now markers as well (needed for 'scatter'
%               plots (thanks to Phil)
%  04-02-2008 - Bugfix markers for Octave (thanks to Jakob Malm)
%  30-12-2008 - Bugfix image scaling and orientation
%               Bugfix correct backslash (thanks to Jason Merril)
%  20-06-2009 - Improvment of image handling (still some remaining issues)
%               Fix for -. line style (thanks to Ritesh Sood)
%  28-06-2009 - Improved depth sorting for patches and surface
%             - Bugfix patches
%             - Bugfix 3D axis handling
%  11-07-2009 - Support of FontWeight and FontAngle properties
%             - Improved markers (polygon instead of polyline for closed markers)
%             - Added character encoding entry to be fully SVG 1.1 conform
%  13-07-2009 - Support of rectangle for 2D
%             - Added preliminary support for SVG filters
%             - Added preliminary support for clipping with pathes
%             - Added preliminary support for turning axis tickmarks
%  18-07-2009 - Line style scaling with line width (will not match with png
%               output)
%             - Small optimizations for the text base line
%             - Bugfix text rotation versus shift
%             - Added more SVG filters
%             - Added checks for filter strings
%  21-07-2009 - Improved bounding box calculation for filters
%             - Bugfixes for text size / line distance
%             - Support of background box for text
%             - Correct bounding box for text objects
%  31-07-2009 - Improved support of filters
%             - Experimental support of animations
%  16-08-2009 - Argument checks for filters
%             - Rework of latex string handling
%             - 'sub' and 'super' workaround for Firefox and Inkscape
%  31-10-2009 - Bugfix for log axes (missing minor grid for some special
%               cases)
%  24-01-2010 - Bugfix nomy line #1102 (thanks to Pooya Jannaty)
%  17-02-2010 - Bugfix minor tickmarks for log axis scaling (thanks to
%               Harke Pera)
%             - Added more lex symbols
%  06-03-2010 - Automatic correction of illegal axis scalings by the user
%               (thanks to Juergen)
%             - Renamed plot2svg_beta to plot2svg
%  12-04-2010 - Improved Octave compatibility
%  05-05-2010 - Bugfix for ticklabels outside of the axis limits (thanks to
%               Ben Scandella)
%  30-10-2010 - Improved handling of empty cells for labels (thanks to 
%               Constantine)
%             - Improved HTML character coding (thanks to David Mack)
%             - Bugfix for last ')' (thanks to Jonathon Harding and Benjamin)
%             - Enabled scatter plots using hggroups
%             - Closing patches if they do not contain NaNs
%  10-11-2010 - Support of the 'Layer' keyword to but the grid on top of 
%               of the other axis content using 'top' (Many thanks to Justin
%               Ashmall)
%             - Tiny optimization of the grid display at axis borders
%  25-08-2011 - Fix for degree character (thanks to Manes Recheis)
%             - Fix for problems with dash-arrays in Inkscape (thanks to
%               Rüdiger Stirnberg)
%             - Modified shape of driangles (thanks to Rüdiger Stirnberg)
%  22-10-2011 - Removed versn as return value of function fileparts (thanks
%               to Andrew Scott)
%             - Fix for images (thanks to Roeland)
%  20-05-2012 - Added some security checks for empty data
%             - Fixed rotation for multiline text
%  25-08-2012 - Special handling of 1xn char arrays for tick labels
%               (thanks to David Plavcan)
%             - Fix for 'Index exceeds matrix dimensions' of axis labels
%               (thanks to Aslak Grinsted)
%             - Fix for another axis label problem (thanks to Ben Mitch)
%  15-09-2012 - Fix for linestyle none of rectangles (thanks to Andrew)
%             - Enabled scatter plot functionality

global PLOT2SVG_globals
global colorname
progversion='15-Sep-2012';
PLOT2SVG_globals.runningIdNumber = 0;
PLOT2SVG_globals.octave = false;
PLOT2SVG_globals.checkUserData = true;
PLOT2SVG_globals.ScreenPixelsPerInch = 90; % Default 90ppi
try
    PLOT2SVG_globals.ScreenPixelsPerInch = get(0, 'ScreenPixelsPerInch');
catch
    % Keep the default 90ppi
end
if nargout==1
    varargout={0};
end
disp(['   Matlab/Octave to SVG converter version ' progversion ', Juerg Schwizer (converter@bluewin.ch).'])
matversion=version;
if exist('OCTAVE_VERSION','builtin')
    PLOT2SVG_globals.octave = true;
    disp('   Info: PLOT2SVG runs in Octave mode.')
else
    if str2num(matversion(1))<6 % Check for matlab version and print warning if matlab version lower than version 6.0 (R.12)
        disp('   Warning: Future versions may no more support older versions than MATLAB R12.')
    end
end
if nargout > 1
    error('Function returns only one return value.')
end
if nargin<2 % Check if handle was included into function call, otherwise take current figure
    id=gcf;
end
if nargin==0
    if PLOT2SVG_globals.octave
        error('PLOT2SVG in Octave mode does not yet support a file menu. File name is needed during function call.')
    else
        [filename, pathname] = uiputfile( {'*.svg', 'SVG File (*.svg)'},'Save Figure as SVG File');
        if ~( isequal( filename, 0) || isequal( pathname, 0))    
            % yes. add backslash to path (if not already there)
            pathname = addBackSlash( pathname); 
            % check, if extension is allrigth
            if ( ~strcmpi( getFileExtension( filename), '.svg'))
                filename = [ filename, '.svg'];
            end
            finalname=[pathname filename];
        else
            disp('   Cancel button was pressed.')
            return
        end
    end
else
    if isnumeric(param1)
        if PLOT2SVG_globals.octave
            error('PLOT2SVG in Octave mode does not yet support a file menu. File name is needed during function call.')
        else
            [filename, pathname] = uiputfile( {'*.svg', 'SVG File (*.svg)'},'Save Figure as SVG File');  
            if ~( isequal( filename, 0) || isequal( pathname, 0))    
                % yes. add backslash to path (if not already there)
                pathname = addBackSlash( pathname); 
                % check, if ectension is allrigth
                if ( ~strcmpi( getFileExtension( filename), '.svg'))
                    filename = [ filename, '.svg'];
                end
                finalname=[pathname filename];
            else
                disp('   Cancel button was pressed.')
                return
            end
        end
    else
        finalname=param1;   
    end
end
% needed to see annotation axes
originalShowHiddenHandles = get(0, 'ShowHiddenHandles');
set(0, 'ShowHiddenHandles', 'on');
originalFigureUnits=get(id,'Units');
set(id,'Units','pixels');   % All data in the svg-file is saved in pixels
paperpos=get(id,'Position');
if ( nargin > 0)
    if isnumeric(param1)
        paperpos(3)=param1(1);
        paperpos(4)=param1(2);
    end
end
paperpos = paperpos * 90 / PLOT2SVG_globals.ScreenPixelsPerInch;
if (nargin < 3)
    PLOT2SVG_globals.pixelfiletype = 'png';
else
    PLOT2SVG_globals.pixelfiletype = pixelfiletype;
end
cmap=get(id,'Colormap');
colorname='';
for i=1:size(cmap,1)
    colorname(i,:)=sprintf('%02x%02x%02x',fix(cmap(i,1)*255),fix(cmap(i,2)*255),fix(cmap(i,3)*255));
end

% Open SVG-file
[pathstr,name] = fileparts(finalname);
%PLOT2SVG_globals.basefilename = fullfile(pathstr,name);
PLOT2SVG_globals.basefilepath = pathstr;
PLOT2SVG_globals.basefilename = name;
PLOT2SVG_globals.figurenumber = 1;
fid=fopen(finalname,'w');   % Create a new text file
fprintf(fid,'<?xml version="1.0" encoding="utf-8" standalone="no"?><!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n');    % Insert file header
fprintf(fid,'<svg preserveAspectRatio="xMinYMin meet" width="100%%" height="100%%" viewBox="0 0 %0.3f %0.3f" ',paperpos(3),paperpos(4));
fprintf(fid,' version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"');
%fprintf(fid,' onload="Init(evt)"');
fprintf(fid,'>\n');
fprintf(fid,'  <desc>Matlab Figure Converted by PLOT2SVG written by Juerg Schwizer</desc>\n');
%fprintf(fid,'  <script type="text/ecmascript" xlink:href="puzzle_script.js" />\n');
fprintf(fid,'  <g id="topgroup">\n');
group=1;
groups=[];
% Frame of figure
figcolor = searchcolor(id,get(id, 'Color'));
if (~ strcmp(figcolor, 'none'))
    % Draw rectangle in the background of the graphic frame to cover all
    % other graphic elements
    try % Octave does not have support for InvertHardcopy yet -- Jakob Malm
        if strcmp(get(id,'InvertHardcopy'),'on')
            fprintf(fid,'  <rect x="0" y="0" width="%0.3f" height="%0.3f" fill="#ffffff" stroke="none" />\n',paperpos(3),paperpos(4));
        else
            fprintf(fid,'  <rect x="0" y="0" width="%0.3f" height="%0.3f" fill="%s" stroke="none" />\n',paperpos(3),paperpos(4),figcolor);
        end
    catch
        fprintf(fid,'  <rect x="0" y="0" width="%0.3f" height="%0.3f" fill="%s" stroke="none" />\n',paperpos(3),paperpos(4),figcolor);
    end
end
% Search all axes
ax=get(id,'Children');
for j=length(ax):-1:1
    currenttype = get(ax(j),'Type');
    if strcmp(currenttype,'axes')
        group=group+1;
        groups=[groups group];
        group=axes2svg(fid,id,ax(j),group,paperpos);
    elseif strcmp(currenttype,'uicontrol')
        if strcmp(get(ax(j),'Visible'),'on')
            control2svg(fid,id,ax(j),group,paperpos);
        end
    elseif strcmp(currenttype, 'uicontextmenu') || ...
            strcmp(currenttype, 'uimenu') || ...
            strcmp(currenttype, 'hgjavacomponent') || ...
            strcmp(currenttype, 'uitoolbar')
        % ignore these types
    else
        disp(['   Warning: Unhandled main figure child type: ' currenttype]);
    end
end
fprintf(fid,'  </g>\n');
fprintf(fid,'</svg>\n');
fclose(fid);    % close text file
if nargout==1
    varargout={0};
end
set(id,'Units',originalFigureUnits);
set(0, 'ShowHiddenHandles', originalShowHiddenHandles);

function clippingIdString = clipping2svg(fid, id, ax, paperpos, axpos, projection, clippingIdString)
global PLOT2SVG_globals
if PLOT2SVG_globals.checkUserData && isstruct(get(id,'UserData'))
    struct_data = get(id,'UserData');
    if isfield(struct_data,'svg')
        if isfield(struct_data.svg,'ClippingPath')
            clip = struct_data.svg.ClippingPath;
            if ~isempty(clip)
                if size(clip, 2) ~=3
                    if size(clip, 2) ==2
                        clipx = clip(:, 1);
                        clipy = clip(:, 2);
                        clipz = zeros(size(clip,1),1);
                    else
                        error('The clipping vector has to be a nx3 or nx2 matrix.');
                    end
                else
                    clipx = clip(:, 1);
                    clipy = clip(:, 2);
                    clipz = clip(:, 3);
                end
                if strcmp(get(ax,'XScale'),'log')
                    clipx(find(clipx<=0)) = NaN;
                    clipx=log10(clipx);
                end
                if strcmp(get(ax,'YScale'),'log')
                    clipy(find(clipy<=0)) = NaN;
                    clipy=log10(clipy);
                end
                if strcmp(get(ax,'ZScale'),'log')
                    clipz(find(clipz<=0)) = NaN;
                    clipz=log10(clipz);
                end
                [x,y,z] = project(clipx,clipy,clipz,projection);
                x = (x*axpos(3)+axpos(1))*paperpos(3);
                y = (1-(y*axpos(4)+axpos(2)))*paperpos(4);
                clippingIdString = createId;
                fprintf(fid,'<clipPath id="%s">\n  <polygon fill="none" stroke="none" points="', clippingIdString);
                fprintf(fid,'%0.3f,%0.3f ',[x';y']);
                fprintf(fid,'"/>\n</clipPath>\n'); 
            end
        end
    end
end


function [angle, align] = improvedXLabel(id, angle, align)
global PLOT2SVG_globals
if PLOT2SVG_globals.checkUserData && isstruct(get(id,'UserData'))
    struct_data = get(id,'UserData');
    if isfield(struct_data,'svg')
        if isfield(struct_data.svg,'XTickLabelAngle')
            angle = struct_data.svg.XTickLabelAngle;
            align = 'Left';
        end
    end
end

function [angle, align] = improvedYLabel(id, angle, align)
global PLOT2SVG_globals
if PLOT2SVG_globals.checkUserData && isstruct(get(id,'UserData'))
    struct_data = get(id,'UserData');
    if isfield(struct_data,'svg')
        if isfield(struct_data.svg,'YTickLabelAngle')
            angle = struct_data.svg.YTickLabelAngle;
            align = 'Left';
        end
    end
end

function animation2svg(fid, id)
global PLOT2SVG_globals
if PLOT2SVG_globals.checkUserData && isstruct(get(id,'UserData'))
    struct_data = get(id,'UserData');
    if isfield(struct_data,'svg')
        if isfield(struct_data.svg,'Animation')
            animation = struct_data.svg.Animation;
            for i = 1:length(animation)
                if ~isfield(animation(i).SubAnimation, 'Type')
                    error(['Missing field ''Type'' for animation.']);
                end
                switch animation(i).SubAnimation.Type
                    case 'Opacity', type = 'opacity'; animationType = 0;
                    case 'Translate', type = 'translate'; animationType = 2;
                    case 'Scale', type = 'scale'; animationType = 2;
                    case 'Rotate', type = 'rotate'; animationType = 1;   
                    case 'skewX', type = 'skewX'; animationType = 1;
                    case 'skewY', type = 'skewY'; animationType = 1;
                    otherwise, error(['Unknown animation type ''' animation(i).SubAnimation.Type '''.']);
                end
                %fprintf(fid,'  <animate attributeType="XML" attributeName="%s" from="%0.3f" to="%0.3f" dur="%0.3fs" repeatCount="%s" />',...
                %   'opacity' , 0, 1, 5, 'indefinite');
                if animationType == 0
                    fprintf(fid,'  <animate attributeType="XML" attributeName="%s" dur="%0.3fs"', type, animation(i).SubAnimation.Duration);
                    fprintf(fid,' values="');
                    fprintf(fid,'%0.2f;', animation(i).SubAnimation.Value);
                    fprintf(fid,'" keyTimes="');
                    fprintf(fid,'%0.2f;', max(min(animation(i).SubAnimation.Key, 1), 0));
                    fprintf(fid,'" repeatCount="%s" calcMode="linear" />', 'indefinite');
                elseif animationType == 1
                    fprintf(fid,'  <animateTransform attributeName="transform" attributeType="XML" type="%s" dur="%0.3fs"', type, animation(i).SubAnimation.Duration);
                    fprintf(fid,' values="');
                    fprintf(fid,'%0.2f;', animation(i).SubAnimation.Value);
                    fprintf(fid,'" keyTimes="');
                    fprintf(fid,'%0.2f;', max(min(animation(i).SubAnimation.Key, 1), 0));
                    fprintf(fid,'" repeatCount="%s" calcMode="linear" additive="sum" />', 'indefinite');
                elseif animationType == 2
                    fprintf(fid,'  <animateTransform attributeName="transform" attributeType="XML" type="%s" dur="%0.3fs"', type, animation(i).SubAnimation.Duration);
                    fprintf(fid,' values="');
                    fprintf(fid,'%0.2f,%0.2f;', animation(i).SubAnimation.Value);
                    fprintf(fid,'" keyTimes="');
                    fprintf(fid,'%0.2f;', max(min(animation(i).SubAnimation.Key, 1), 0));
                    fprintf(fid,'" repeatCount="%s" calcMode="linear" additive="sum" />', 'indefinite');
                end
            end    
        end        
    end
end

function [filterString, boundingBox] = filter2svg(fid, id, boundingBoxAxes, boundingBoxElement)
global PLOT2SVG_globals
filterString = '';
boundingBox = boundingBoxAxes;
if PLOT2SVG_globals.checkUserData && isstruct(get(id,'UserData'))
    struct_data = get(id,'UserData');
    if isfield(struct_data,'svg')
        boundingBox = boundingBoxElement;
        absolute = true;
        offset = 0;
        if isfield(struct_data.svg,'BoundingBox')
            if isfield(struct_data.svg.BoundingBox, 'Type')
                switch struct_data.svg.BoundingBox.Type
                    case 'axes', boundingBox = boundingBoxAxes; absolute = true;
                    case 'element', boundingBox = boundingBoxElement; absolute = true;
                    case 'relative', boundingBox = boundingBoxElement; absolute = false;
                    otherwise
                        error(['Unknown bounding box type ''' struct_data.svg.BoundingBox.Type '''.']);
                end
            end
            if isfield(struct_data.svg.BoundingBox, 'Overlap')
                overlap = struct_data.svg.BoundingBox.Overlap;
                if absolute
                    boundingBox(1) = boundingBox(1) - overlap;
                    boundingBox(2) = boundingBox(2) - overlap;
                    boundingBox(3) = boundingBox(3) + 2 * overlap;
                    boundingBox(4) = boundingBox(4) + 2 * overlap;
                else
                    boundingBox(1) = boundingBox(1) - boundingBox(3) * overlap;
                    boundingBox(2) = boundingBox(2) - boundingBox(4) * overlap;
                    boundingBox(3) = boundingBox(3) + 2 * boundingBox(3) * overlap;
                    boundingBox(4) = boundingBox(4) + 2 * boundingBox(4) * overlap;    
                end
            end
            if isfield(struct_data.svg.BoundingBox, 'Visible') && strcmp(struct_data.svg.BoundingBox.Visible, 'on')
                % This functionality is very interesting for debugging of
                % bounding boxes of filters
                fprintf(fid,'<rect x="%0.3f" y="%0.3f" width="%0.3f" height="%0.3f" fill="none" stroke="#000000" stroke-dasharray="1,1" stroke-width="0.2pt" />\n', boundingBox(1), boundingBox(2), boundingBox(3), boundingBox(4));
            end
        end
        if isfield(struct_data.svg,'Filter')
            % Predefined filter sources. Additional filter sources will be
            % added later.
            predefinedSources = {'SourceGraphic','SourceAlpha','BackgroundImage','BackgroundAlpha','FillPaint','StrokePaint'};
            resultStrings = predefinedSources;
            filterId = createId;
            filterString  = ['filter="url(#' filterId ')"'];
            fprintf(fid,'<defs>\n');
            fprintf(fid,'  <filter x="%0.3f%%" y="%0.3f%%" width="%0.3f%%" height="%0.3f%%" id="%s">\n', 0, 0, 100, 100, filterId);
            %if absolute
            %    fprintf(fid,'  <filter x="%0.3f" y="%0.3f" width="%0.3f" height="%0.3f" filterUnits="userSpaceOnUse" id="%s">\n', boundingBox(1), boundingBox(2), boundingBox(3), boundingBox(4), filterId);
            %else
            %    fprintf(fid,'  <filter x="%0.3f%%" y="%0.3f%%" width="%0.3f%%" height="%0.3f%%" id="%s">\n', -(offset * 100), -(offset * 100), 100 + (offset * 200), 100 + (offset * 200), filterId);
            %    % Note: use default -10% for attribute x
            %    %       use default -10% for attribute y
            %    %       use default 120% for attribute width
            %    %       use default 120% for attribute height
            %end
            filter = struct_data.svg.Filter;
            for i = 1:length(filter)
                if isfield(filter(i).Subfilter, 'Type')
                    fprintf(fid,'  <%s', filter(i).Subfilter.Type);
                else
                    error(['Missing field ''Type'' for filter.'])
                end
                try
                    if isfield(filter(i).Subfilter, 'Position')
                        printAttributeArray(fid, {'x','y','width','height'}, filter(i).Subfilter, 'Position');
                    end
                    printAttributeString(fid, 'result', filter(i).Subfilter, 'Result');
                    % Add result string to the list in order to check the in
                    % strings.
                    resultStrings{length(resultStrings) + 1} = filter(i).Subfilter.Result;
                    % The strmatch below is a very inefficient search (Matlab limitation)
                    if ~isempty(strmatch(filter(i).Subfilter.Result, predefinedSources))
                        error('Usage of a predefined filter source as filter result string is not allowed.');
                    end
                    switch (filter(i).Subfilter.Type)
                        case 'feGaussianBlur'
                            printAttributeIn(fid, 'in', filter(i).Subfilter, 'Source', 'SourceGraphic', resultStrings);
                            printAttributeDouble(fid, 'stdDeviation', filter(i).Subfilter, 'Deviation');
                            fprintf(fid,' />\n');
                        case 'feImage'
                            printAttributeString(fid, 'xlink:href', filter(i).Subfilter, 'File');
                            printAttributeString(fid, 'preserveAspectRatio', filter(i).Subfilter, 'AspectRatio', 'xMidYMid meet');
                            fprintf(fid,' />\n');
                        case 'feComposite'    
                            printAttributeIn(fid, 'in', filter(i).Subfilter, 'Source1', 'SourceGraphic', resultStrings);
                            printAttributeIn(fid, 'in2', filter(i).Subfilter, 'Source2', 'SourceGraphic', resultStrings);
                            printAttributeList(fid, 'operator', filter(i).Subfilter, 'Operator', {'over','in','out','atop','xor','arithmetic'}, 'over'); % 'over' | 'in' | 'out' | 'atop' | 'xor' | 'arithmetic'
                            if isfield(filter(i).Subfilter, 'Operator') && strcmp(filter(i).Subfilter.Operator, 'arithmetic')
                                printAttributeArray(fid, {'k1','k2','k3','k4'}, filter(i).Subfilter, 'k');
                            end
                            fprintf(fid,' />\n');
                        case 'feSpecularLighting'
                            printAttributeDouble(fid, 'specularConstant', filter(i).Subfilter, 'SpecularConstant');
                            printAttributeDouble(fid, 'specularExponent', filter(i).Subfilter, 'SpecularExponent');
                            printAttributeDouble(fid, 'surfaceScale', filter(i).Subfilter, 'SurfaceScale');
                            fprintf(fid,' style="lighting-color:white"');
                            printAttributeIn(fid, 'in', filter(i).Subfilter, 'Source', 'SourceGraphic', resultStrings);
                            fprintf(fid,' >\n');
                            if isfield(filter(i).Subfilter, 'LightType')
                                fprintf(fid,' <%s', filter(i).Subfilter.LightType);
                                switch filter(i).Subfilter.LightType
                                    case 'feDistantLight'
                                        printAttributeDouble(fid, 'azimuth', filter(i).Subfilter, 'Azimuth');
                                        printAttributeDouble(fid, 'elevation', filter(i).Subfilter, 'Elevation');
                                    case 'fePointLight'
                                        printAttributeArray(fid, {'x','y','z'}, filter(i).Subfilter, 'Position');
                                    case 'feSpotLight'
                                        printAttributeArray(fid, {'x','y','z'}, filter(i).Subfilter, 'Position');
                                        printAttributeArray(fid, {'pointsAtX','pointsAtY','pointsAtZ'}, filter(i).Subfilter, 'PositionsAt');
                                        printAttributeDouble(fid, 'specularExponent', filter(i).Subfilter, 'LightSpecularExponent');
                                        printAttributeDouble(fid, 'limitingConeAngle', filter(i).Subfilter, 'LimitingConeAngle');
                                    otherwise, error(['Unknown light type ''' filter(i).Subfilter.LightType ''.']);
                                end
                                fprintf(fid,' />\n');
                            else
                                error('Missing field ''LightType''.');
                            end
                            fprintf(fid,'</%s>\n',filter(i).Subfilter.Type);
                        case 'feOffset'
                            printAttributeIn(fid, 'in', filter(i).Subfilter, 'Source', 'SourceGraphic', resultStrings);
                            printAttributeArray(fid, {'dx','dy'}, filter(i).Subfilter, 'Offset');
                            fprintf(fid,' />\n');
                        case 'feBlend'
                            printAttributeIn(fid, 'in', filter(i).Subfilter, 'Source1', 'SourceGraphic', resultStrings);
                            printAttributeIn(fid, 'in2', filter(i).Subfilter, 'Source2', 'SourceGraphic', resultStrings);
                            printAttributeList(fid, 'mode', filter(i).Subfilter, 'Mode', {'normal','multiply','screen','darken','lighten'}, 'normal');   % 'normal' | 'multiply' | 'screen' | 'darken' | 'lighten'
                            fprintf(fid,' />\n');
                        case 'feTurbulence'
                            printAttributeDouble(fid, 'baseFrequency', filter(i).Subfilter, 'BaseFrequency');
                            printAttributeDouble(fid, 'numOctaves', filter(i).Subfilter, 'NumOctaves');
                            printAttributeDouble(fid, 'seed', filter(i).Subfilter, 'Seed');
                            printAttributeList(fid, 'stitchTiles', filter(i).Subfilter, 'StitchTiles', {'stitch','noStitch'}, 'noStitch');   % stitch | noStitch
                            printAttributeList(fid, 'type', filter(i).Subfilter, 'TurbulenceType', {'fractalNoise','turbulence'}, 'turbulence');   % 'fractalNoise' | 'turbulence'
                            fprintf(fid,' />\n');
                        case 'feColorMatrix'
                            printAttributeIn(fid, 'in', filter(i).Subfilter, 'Source', 'SourceGraphic', resultStrings);
                            printAttributeList(fid, 'type', filter(i).Subfilter, 'ColorType', {'matrix','saturate','hueRotate','luminanceToAlpha'});   % 'matrix' | 'saturate' | 'hueRotate' | 'luminanceToAlpha'
                            if isfield(filter(i).Subfilter, 'ColorType') && strcmp(filter(i).Subfilter.ColorType, 'matrix')
                                if isfield(filter(i).Subfilter, 'Matrix') && (lenght(filter(i).Subfilter.Matrix) == 20)
                                    fprintf(fid,' values="');
                                    fprintf(fid,' %0.3f', filter(i).Subfilter.Matrix);
                                    fprintf(fid,'"');
                                else
                                    error('Field ''Matrix'' is missing or not a 5x4 matrix.');    
                                end
                            end
                            if isfield(filter(i).Subfilter, 'ColorType') && (strcmp(filter(i).Subfilter.ColorType, 'saturate') || strcmp(filter(i).Subfilter.ColorType, 'hueRotate'))
                                printAttributeDouble(fid, 'values', filter(i).Subfilter, 'Matrix');
                            end
                            fprintf(fid,' />\n');
                        case 'feFlood'    
                            printAttributeColor(fid, 'flood-color', filter(i).Subfilter, 'Color');
                            printAttributeDouble(fid, 'flood-opacity', filter(i).Subfilter, 'Opacity');
                            fprintf(fid,' />\n');
                        case 'feDisplacementMap'    
                            printAttributeIn(fid, 'in', filter(i).Subfilter, 'Source1', 'SourceGraphic', resultStrings);
                            printAttributeIn(fid, 'in2', filter(i).Subfilter, 'Source2', 'SourceGraphic', resultStrings);
                            printAttributeDouble(fid, 'scale', filter(i).Subfilter, 'Scale');
                            printAttributeList(fid, 'xChannelSelector', filter(i).Subfilter, 'xChannel', {'R','G','B','A'}, 'A');    % 'R' | 'G' | 'B' | 'A'
                            printAttributeList(fid, 'yChannelSelector', filter(i).Subfilter, 'yChannel', {'R','G','B','A'}, 'A');    % 'R' | 'G' | 'B' | 'A'
                            fprintf(fid,' />\n');
                        case 'feMerge'    
                            printAttributeIn(fid, 'in', filter(i).Subfilter, 'Source1', 'SourceGraphic', resultStrings);
                            printAttributeIn(fid, 'in2', filter(i).Subfilter, 'Source2', 'SourceGraphic', resultStrings);
                            printAttributeList(fid, 'mode', filter(i).Subfilter, 'Mode', {'normal','multiply','screen','darken','lighten'}); % 'normal' | 'multiply' | 'screen' | 'darken' | 'lighten'
                            fprintf(fid,' />\n');
                        case 'feMorphology'
                            printAttributeIn(fid, 'in', filter(i).Subfilter, 'Source', 'SourceGraphic', resultStrings);
                            printAttributeList(fid, 'operator', filter(i).Subfilter, 'Operator', {'erode','dilate'}); % 'erode' | 'dilate'
                            printAttributeDouble(fid, 'radius', filter(i).Subfilter, 'Radius');
                            fprintf(fid,' />\n');
                        case 'feTile'
                            printAttributeIn(fid, 'in', filter(i).Subfilter, 'Source', 'SourceGraphic', resultStrings);
                            fprintf(fid,' />\n');
                        case 'feDiffuseLighting'
                            printAttributeIn(fid, 'in', filter(i).Subfilter, 'Source', 'SourceGraphic', resultStrings);
                            fprintf(fid,' style="lighting-color:white"');
                            printAttributeDouble(fid, 'surfaceScale', filter(i).Subfilter, 'SurfaceScale');
                            printAttributeDouble(fid, 'diffuseConstant', filter(i).Subfilter, 'DiffuseConstant');
                            printAttributeDouble(fid, 'kernelUnitLength', filter(i).Subfilter, 'KernelUnitLength');
                            fprintf(fid,' />\n');
                            
                            %                     case 'feConvolveMatrix'
                            %                         
                            %                         printAttributeDouble(fid, 'order', filter(i).Subfilter);
                            % kernelMatrix
                            %                         printAttributeDouble(fid, 'divisor', filter(i).Subfilter, 1.0);
                            %                         printAttributeDouble(fid, 'bias', filter(i).Subfilter, 0);
                            % targetX
                            % targetY
                            %                         printAttributeString(fid, 'edgeMode', filter(i).Subfilter, 'EdgeMode', 'duplicate');
                            %                         fprintf(fid,' kernelUnitLength="1 1"');
                            %                         printAttributeString(fid, 'preserveAlpha', filter(i).Subfilter, 'PreserveAlpha', 'true');
                            %                         fprintf(fid,' />\n');
                        case {'feComponentTransfer','feConvolveMatrix'}
                            error('Filter not yet implemented.');
                        otherwise
                            error(['Unknown filter ''' filter(i).Subfilter.Type '''.']);
                    end
                catch
                    error([lasterr ' Error is caused by filter type ''' filter(i).Subfilter.Type '''.']);
                end
            end
            fprintf(fid,'  </filter>\n');
            fprintf(fid,'</defs>\n');
        end
    end
end

function printAttributeArray(fid, names, svgstruct, svgfield)
if isfield(svgstruct, svgfield)
    if isnumeric(svgstruct.(svgfield))
        if length(svgstruct.(svgfield)) ~= length(names)
            error(['Length mismatch for field ''' svgfield '''.'])    
        end
        for i = 1:length(names)
            fprintf(fid,' %s="%0.3f"', names{i}, svgstruct.(svgfield)(i));
        end
    else
        error(['Field ''' svgfield ''' must be numeric.']);
    end
else
    if nargin < 5
        error(['Missing field ''' svgfield '''.'])
    else
        for i = 1:length(names)
            fprintf(fid,' %s="%0.3f"', names{i}, default(i));
        end
    end
end

function printAttributeDouble(fid, name, svgstruct, svgfield, default)
if isfield(svgstruct, svgfield)
    if isnumeric(svgstruct.(svgfield))
        fprintf(fid,' %s="%0.3f"', name, svgstruct.(svgfield));
    else
        error(['Field ''' svgfield ''' must be numeric.']);
    end
else
    if nargin < 5
        error(['Missing field ''' svgfield '''.'])
    else
        fprintf(fid,' %s="%0.3f"', name, default);
    end
end

function printAttributeIn(fid, name, svgstruct, svgfield, default, resultStrings)
if isfield(svgstruct, svgfield)
    if ischar(svgstruct.(svgfield))
        % The strmatch below is a very inefficient search (Matlab limitation)
        if isempty(strmatch(svgstruct.(svgfield), resultStrings))
            error(['The source string ''' svgstruct.(svgfield) ''' was never a result string of a previous filter. Check for correct spelling.']);    
        else
            fprintf(fid,' %s="%s"', name, svgstruct.(svgfield));
        end    
    else
        error(['Field ''' svgfield ''' must be a string.']);
    end
else
    if nargin < 5
        error(['Missing field ''' svgfield '''.'])
    else
        fprintf(fid,' %s="%s"', name, default);
    end
end

function printAttributeString(fid, name, svgstruct, svgfield, default)
if isfield(svgstruct, svgfield)
    if ischar(svgstruct.(svgfield))
        fprintf(fid,' %s="%s"', name, svgstruct.(svgfield));
    else
        error(['Field ''' svgfield ''' must be a string.']);
    end
else
    if nargin < 5
        error(['Missing field ''' svgfield '''.'])
    else
        fprintf(fid,' %s="%s"', name, default);
    end
end

function printAttributeList(fid, name, svgstruct, svgfield, list, default)
if isfield(svgstruct, svgfield)
    if ischar(svgstruct.(svgfield))
        if isempty(strmatch(svgstruct.(svgfield), list))
            listString = strcat(list, ''' | ''');
            listString = [listString{:}];
            error(['Illegal string identifier ''' svgstruct.(svgfield) '''. Must be one out of the list: ''' listString(1:end-4) '.']);
        else
            fprintf(fid,' %s="%s"', name, svgstruct.(svgfield));
        end
    else
        error(['Field ''' svgfield ''' must be a string.']);
    end
else
    if nargin < 6
        error(['Missing field ''' svgfield '''.'])
    else
        fprintf(fid,' %s="%s"', name, default);
    end
end

function printAttributeColor(fid, name, svgstruct, svgfield, default)
if isfield(svgstruct, svgfield)
    if isnumeric(svgstruct.(svgfield))
        if length(svgstruct.(svgfield)) ~= 3
            error(['Color must be a 1x3 vector for field ''' svgfield '''.'])    
        else
            fprintf(fid,' %s="%s"', name, searchcolor(gca, svgstruct.(svgfield)));
        end
    else
        error(['Field ''' svgfield ''' must be a 1x3 vector.']);
    end
else
    if nargin < 5
        error(['Missing field ''' svgfield '''.'])
    else
        fprintf(fid,' %s="%s"', name, default);
    end
end


function frontTicks(fid, grouplabel, axpos, x, y, scolorname, linewidth, tick, index, edge_neighbours, c, valid_ticks, ticklength, tick_ratio, lim, drawBorder)
for k = 1:length(index)
    x_tick_end1 = interp1([0 1],[x(index(k)) x(edge_neighbours(index(k),c(1)))],ticklength*tick_ratio(c(3)),'linear','extrap');
    y_tick_end1 = interp1([0 1],[y(index(k)) y(edge_neighbours(index(k),c(1)))],ticklength*tick_ratio(c(3)),'linear','extrap');
    x_tick_end2 = interp1([0 1],[x(edge_neighbours(index(k),c(2))) x(edge_neighbours(edge_neighbours(index(k),c(2)),c(1)))],ticklength*tick_ratio(c(3)),'linear','extrap');
    y_tick_end2 = interp1([0 1],[y(edge_neighbours(index(k),c(2))) y(edge_neighbours(edge_neighbours(index(k),c(2)),c(1)))],ticklength*tick_ratio(c(3)),'linear','extrap');
    xg_line_start = interp1(lim,[x(index(k)) x(edge_neighbours(index(k),c(2)))],tick);
    yg_line_start = interp1(lim,[y(index(k)) y(edge_neighbours(index(k),c(2)))],tick);
    xg_line_end = interp1(lim,[x_tick_end1 x_tick_end2],tick);
    yg_line_end = interp1(lim,[y_tick_end1 y_tick_end2],tick);
    for i = valid_ticks
        line2svg(fid,grouplabel,axpos,[xg_line_start(i) xg_line_end(i)],[yg_line_start(i) yg_line_end(i)],scolorname,'-',linewidth)
    end
    if drawBorder
        line2svg(fid,grouplabel,axpos,[x(index(k)) x(edge_neighbours(index(k),c(2)))],[y(index(k)) y(edge_neighbours(index(k),c(2)))],scolorname,'-',linewidth)
    end
end

function gridLines(fid, grouplabel, axpos, x, y, scolorname, gridlinestyle, linewidth, axlim, axtick, axindex_inner, corners, c)
xg_line_start = interp1([axlim(1) axlim(2)],[x(corners(c(1))) x(corners(c(2)))], axtick);
yg_line_start = interp1([axlim(1) axlim(2)],[y(corners(c(1))) y(corners(c(2)))], axtick);
xg_line_end = interp1([axlim(1) axlim(2)],[x(corners(c(3))) x(corners(c(4)))], axtick);
yg_line_end = interp1([axlim(1) axlim(2)],[y(corners(c(3))) y(corners(c(4)))], axtick);
for i = axindex_inner
    line2svg(fid, grouplabel, axpos, [xg_line_start(i) xg_line_end(i)],[yg_line_start(i) yg_line_end(i)], scolorname, gridlinestyle, linewidth)
end

function minorGridLines(fid, grouplabel, axpos, x, y, scolorname, minor_gridlinestyle, linewidth, axlim, minor_axtick, corners, c)
xg_line_start = interp1([axlim(1) axlim(2)],[x(corners(c(1))) x(corners(c(2)))], minor_axtick);
yg_line_start = interp1([axlim(1) axlim(2)],[y(corners(c(1))) y(corners(c(2)))], minor_axtick);
xg_line_end = interp1([axlim(1) axlim(2)],[x(corners(c(3))) x(corners(c(4)))], minor_axtick);
yg_line_end = interp1([axlim(1) axlim(2)],[y(corners(c(3))) y(corners(c(4)))], minor_axtick);
for i = 1:length(xg_line_start)
    line2svg(fid, grouplabel, axpos, [xg_line_start(i) xg_line_end(i)],[yg_line_start(i) yg_line_end(i)], scolorname, minor_gridlinestyle, linewidth)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBFUNCTIONS %%%%%
% Create axis frame and insert all children of this axis frame
function group=axes2svg(fid,id,ax,group,paperpos)
global colorname
global PLOT2SVG_globals
originalAxesUnits=get(ax,'Units');
set(ax,'Units','normalized');
axpos=get(ax,'Position');
faces =    [1 2 4 3; 2 4 8 6; 3 4 8 7; 1 2 6 5; 1 5 7 3; 5 6 8 7];
%           x-y    ; y-z    ; x-z    ; y-z    ; x-z    ; x-y
corners(:,:,1) = [1 1 2 3 4; 2 1 3 2 4];
corners(:,:,2) = [2 2 4 6 8; 3 2 6 4 8];
corners(:,:,3) = [1 3 4 7 8; 3 3 7 4 8];
corners(:,:,4) = [1 1 2 5 6; 3 1 5 2 6];
corners(:,:,5) = [2 1 3 5 7; 3 1 5 3 7];
corners(:,:,6) = [1 5 6 7 8; 2 5 7 6 8];
edge_neighbours = [2 3 5; 1 4 6; 4 1 7; 3 2 8; 6 7 1; 5 8 2; 8 5 3; 7 6 4];
edge_opposite = [8 7 6 5 4 3 2 1];
nomx = [0 1 0 1 0 1 0 1];
nomy = [0 0 1 1 0 0 1 1];
nomz = [0 0 0 0 1 1 1 1];
[projection,edges] = get_projection(ax,id);
x = (edges(1,:)*axpos(3)+axpos(1))*paperpos(3);
y = (1-(edges(2,:)*axpos(4)+axpos(2)))*paperpos(4);    
% Depth Sort of view box edges 
[edge_z,edge_index]=sort(edges(3,:));
most_back_edge_index = edge_index(1);
% Back faces are plot box faces that are behind the plot (as seen by the
% view point)
back_faces = find(any(faces == most_back_edge_index,2));
front_faces = find(all(faces ~= most_back_edge_index,2));
groupax=group;
axlimx=get(ax,'XLim');
axlimy=get(ax,'YLim');
axlimz=get(ax,'ZLim');
axlimxori=axlimx;
axlimyori=axlimy;
axlimzori=axlimz;
if strcmp(get(ax,'XScale'),'log')
    axlimx=log10(axlimx);
    axlimx(find(isinf(axlimx)))=0;
end
if strcmp(get(ax,'YScale'),'log')
    axlimy=log10(axlimy);
    axlimy(find(isinf(axlimy)))=0;
end
if strcmp(get(ax,'ZScale'),'log')
    axlimz=log10(axlimz);
    axlimz(find(isinf(axlimz)))=0;
end
if strcmp(get(ax,'XDir'),'reverse')
    axlimx = fliplr(axlimx);
end
if strcmp(get(ax,'YDir'),'reverse')
    axlimy = fliplr(axlimy);
end
if strcmp(get(ax,'ZDir'),'reverse')
    axlimz = fliplr(axlimz);
end
axlimori = [axlimxori(1) axlimyori(1) axlimzori(1) axlimxori(2)-axlimxori(1) axlimyori(2)-axlimyori(1) axlimzori(2)-axlimzori(1)];
fprintf(fid,'  <g id ="%s">\n', createId);
axIdString = createId;
boundingBoxAxes = [min(x) min(y) max(x)-min(x) max(y)-min(y)];
fprintf(fid,'  <clipPath id="%s">\n',axIdString);
fprintf(fid,'    <rect x="%0.3f" y="%0.3f" width="%0.3f" height="%0.3f"/>\n',...
    boundingBoxAxes(1), boundingBoxAxes(2), boundingBoxAxes(3), boundingBoxAxes(4));
fprintf(fid,'  </clipPath>\n');
if strcmp(get(ax,'Visible'),'on')
    group=group+1;
    grouplabel=group;
    axxtick=get(ax,'XTick');
    axytick=get(ax,'YTick');
    axztick=get(ax,'ZTick');
    axlabelx=get(ax,'XTickLabel');
    axlabely=get(ax,'YTickLabel');
    axlabelz=get(ax,'ZTickLabel');
    % Workaround for Octave
    if PLOT2SVG_globals.octave
        if isempty(axlabelx)
            if strcmp(get(ax,'XScale'),'log')
                axlabelx = num2str(log10(axxtick)');
            else
                axlabelx = num2str(axxtick');
            end
        end
        if isempty(axlabely)
            if strcmp(get(ax,'YScale'),'log')
                axlabely = num2str(log10(axytick)');
            else
                axlabely = num2str(axytick');
            end
        end
        if isempty(axlabelz)
            if strcmp(get(ax,'ZScale'),'log')
                axlabelz = num2str(log10(axztick)');
            else
                axlabelz = num2str(axztick');
            end
        end
        if projection.xyplane
            axlabelz = [];
        end
    end
    gridlinestyle=get(ax,'GridLineStyle');
    minor_gridlinestyle=get(ax,'MinorGridLineStyle');
    try % Octave does not have 'TickLength' yet. --Jakob Malm
        both_ticklength = get(ax,'TickLength');
    catch
        both_ticklength = [ 0.01 0.025 ];
    end
    gridBehind = true; % Default setting
    try
        if strcmp(get(ax, 'Layer'), 'top') && projection.xyplane 
            gridBehind = false;
        end
    catch
        gridBehind = true;
    end
    if projection.xyplane
        ticklength = both_ticklength(1);
        xy_ratio = axpos(3)*paperpos(3)/ (axpos(4)*paperpos(4));
        if xy_ratio < 1
            tick_ratio = [1 1/xy_ratio 1];
        else
            tick_ratio = [xy_ratio 1 1];
        end
        if strcmp(get(ax,'TickDir'),'out')
            label_distance = -(0.02 + ticklength);
        else
            label_distance = -0.02;
        end
    else
        ticklength = both_ticklength(2);
        label_distance = -2*abs(ticklength);
        tick_ratio = [1 1 1];
    end
    linewidth = get(ax,'LineWidth');
    axxindex=find((axxtick >= axlimori(1)) & (axxtick <= (axlimori(1)+axlimori(4))));
    axyindex=find((axytick >= axlimori(2)) & (axytick <= (axlimori(2)+axlimori(5))));
    axzindex=find((axztick >= axlimori(3)) & (axztick <= (axlimori(3)+axlimori(6))));
    % remove sticks outside of the axes (-1 of legends)
    axxtick=axxtick(axxindex); 
    axytick=axytick(axyindex);
    axztick=axztick(axzindex);
    if length(axxtick) > 1
        minor_lin_sticks = (0.2:0.2:0.8)*(axxtick(2)-axxtick(1));
        minor_axxtick = [];
        for stick = [2*axxtick(1)-axxtick(2) axxtick]
            minor_axxtick = [minor_axxtick minor_lin_sticks + stick]; 
        end
        minor_axxtick = minor_axxtick(find(minor_axxtick > min(axlimx) & minor_axxtick < max(axlimx)));
    else
        minor_axxtick = [];
    end
    if length(axytick) > 1
        minor_lin_sticks = (0.2:0.2:0.8)*(axytick(2)-axytick(1));
        minor_axytick = [];
        for stick = [2*axytick(1)-axytick(2) axytick]
            minor_axytick = [minor_axytick minor_lin_sticks + stick]; 
        end
        minor_axytick = minor_axytick(find(minor_axytick > min(axlimy) & minor_axytick < max(axlimy)));
    else
        minor_axytick = [];
    end
    if length(axztick) > 1
        minor_lin_sticks = (0.2:0.2:0.8)*(axztick(2)-axztick(1));
        minor_axztick = [];
        for stick = [2*axztick(1)-axztick(2) axztick]
            minor_axztick = [minor_axztick minor_lin_sticks + stick]; 
        end
        minor_axztick = minor_axztick(find(minor_axztick > min(axlimz) & minor_axztick < max(axlimz)));
    else
        minor_axztick = [];
    end
    if strcmp(get(ax,'Box'),'on')
        axxindex_inner = find((axxtick > axlimori(1)) & (axxtick < (axlimori(1)+axlimori(4))));
        axyindex_inner = find((axytick > axlimori(2)) & (axytick < (axlimori(2)+axlimori(5))));
        axzindex_inner = find((axztick > axlimori(3)) & (axztick < (axlimori(3)+axlimori(6))));
    else
        axxindex_inner = find((axxtick >= axlimori(1)) & (axxtick <= (axlimori(1)+axlimori(4))));
        axyindex_inner = find((axytick >= axlimori(2)) & (axytick <= (axlimori(2)+axlimori(5))));
        axzindex_inner = find((axztick >= axlimori(3)) & (axztick <= (axlimori(3)+axlimori(6))));
    end
    minor_log_sticks = log10(0.2:0.1:0.9);
    if strcmp(get(ax,'TickDir'),'out')
        ticklength=-ticklength;
        valid_xsticks = 1:length(axxindex);
        valid_ysticks = 1:length(axyindex);
        valid_zsticks = 1:length(axzindex);
    else
        valid_xsticks = axxindex_inner;
        valid_ysticks = axyindex_inner;
        valid_zsticks = axzindex_inner;
    end
    if strcmp(get(ax,'XScale'),'log')
        axxtick = log10(get(ax,'XTick'));
        minor_axxtick = [];
        if ~isempty(axxtick)
            all_axxtick = axxtick(1):1:axxtick(end); 
            for stick = all_axxtick
                minor_axxtick = [minor_axxtick minor_log_sticks + stick]; 
            end
        end
        minor_axxtick = minor_axxtick(find(minor_axxtick > min(axlimx) & minor_axxtick < max(axlimx)));
    end
    if strcmp(get(ax,'YScale'),'log')
        axytick=log10(get(ax,'YTick'));
        minor_axytick = [];
        if ~isempty(axytick)
            all_axytick = axytick(1):1:axytick(end); 
            for stick = all_axytick
                minor_axytick = [minor_axytick minor_log_sticks + stick]; 
            end
        end
        minor_axytick = minor_axytick(find(minor_axytick > min(axlimy) & minor_axytick < max(axlimy)));
    end
    if strcmp(get(ax,'ZScale'),'log')
        axztick=log10(get(ax,'ZTick'));
        minor_axztick = [];
        if ~isempty(axztick)
            all_axztick = axztick(1):1:axztick(end); 
            for stick = all_axztick
                minor_axztick = [minor_axztick minor_log_sticks + stick]; 
            end
        end
        minor_axztick = minor_axztick(find(minor_axztick > min(axlimz) & minor_axztick < max(axlimz)));
    end
    % Draw back faces 
    linewidth=get(ax,'LineWidth');
    if ~strcmp(get(ax,'Color'),'none')
        background_color = searchcolor(id,get(ax,'Color'));
        background_opacity = 1;
    else
        background_color = '#000000';
        background_opacity = 0;
    end
    for p=1:size(back_faces)
        patch2svg(fid, group, axpos, x(faces(back_faces(p),:)), y(faces(back_faces(p),:)), background_color, '-', linewidth, 'none', background_opacity, 1.0, true)
    end
    for pindex = 1:size(back_faces)
        p = back_faces(pindex);
        for k = 1:size(corners,1)
            selectedCorners = squeeze(corners(k,:,p));
            switch corners(k,1,p)
                case 1 % x
                    % Draw x-grid
                    scolorname = searchcolor(id,get(ax,'XColor'));
                    if strcmp(get(ax,'XGrid'),'on') && gridBehind
                        if axlimx(1)~=axlimx(2)
                            gridLines(fid, grouplabel, axpos, x, y, scolorname, gridlinestyle, linewidth, axlimx, axxtick, axxindex_inner, selectedCorners, [2 3 4 5])
                            if strcmp(get(ax,'XTickMode'),'auto') && strcmp(get(ax,'XMinorGrid'),'on') && ~isempty(minor_axxtick)
                                minorGridLines(fid, grouplabel, axpos, x, y, scolorname, minor_gridlinestyle, linewidth, axlimx, minor_axxtick, selectedCorners, [2 3 4 5])                                
                            end
                        end
                    end
                    if projection.xyplane == false
                        if strcmp(get(ax,'Box'),'on')
                            line2svg(fid,grouplabel,axpos,[x(corners(k,2,p)) x(corners(k,3,p))],[y(corners(k,2,p)) y(corners(k,3,p))],scolorname,'-',linewidth);
                            line2svg(fid,grouplabel,axpos,[x(corners(k,4,p)) x(corners(k,5,p))],[y(corners(k,4,p)) y(corners(k,5,p))],scolorname,'-',linewidth);
                        else
                            if strcmp(get(ax,'XGrid'),'on')
                                line2svg(fid,grouplabel,axpos,[x(corners(k,2,p)) x(corners(k,3,p))],[y(corners(k,2,p)) y(corners(k,3,p))],scolorname,gridlinestyle,linewidth);
                                line2svg(fid,grouplabel,axpos,[x(corners(k,4,p)) x(corners(k,5,p))],[y(corners(k,4,p)) y(corners(k,5,p))],scolorname,gridlinestyle,linewidth);
                            end
                        end
                    end
                case 2 % y
                    % Draw y-grid
                    scolorname = searchcolor(id,get(ax,'YColor'));
                    if strcmp(get(ax,'YGrid'),'on') && gridBehind
                        if axlimy(1)~=axlimy(2)
                            gridLines(fid, grouplabel, axpos, x, y, scolorname, gridlinestyle, linewidth, axlimy, axytick, axyindex_inner, selectedCorners, [2 3 4 5])
                            if strcmp(get(ax,'YTickMode'),'auto') && strcmp(get(ax,'YMinorGrid'),'on') && ~isempty(minor_axytick)
                                minorGridLines(fid, grouplabel, axpos, x, y, scolorname, minor_gridlinestyle, linewidth, axlimy, minor_axytick, selectedCorners, [2 3 4 5])                                                                
                            end
                        end
                    end
                    if projection.xyplane == false
                        if strcmp(get(ax,'Box'),'on')
                            line2svg(fid,grouplabel,axpos,[x(corners(k,2,p)) x(corners(k,3,p))],[y(corners(k,2,p)) y(corners(k,3,p))],scolorname,'-',linewidth);
                            line2svg(fid,grouplabel,axpos,[x(corners(k,4,p)) x(corners(k,5,p))],[y(corners(k,4,p)) y(corners(k,5,p))],scolorname,'-',linewidth);
                        else
                            if strcmp(get(ax,'YGrid'),'on')
                                line2svg(fid,grouplabel,axpos,[x(corners(k,2,p)) x(corners(k,3,p))],[y(corners(k,2,p)) y(corners(k,3,p))],scolorname,gridlinestyle,linewidth);
                                line2svg(fid,grouplabel,axpos,[x(corners(k,4,p)) x(corners(k,5,p))],[y(corners(k,4,p)) y(corners(k,5,p))],scolorname,gridlinestyle,linewidth);
                            end
                        end
                    end
                case 3 % z
                    % Draw z-grid
                    scolorname = searchcolor(id,get(ax,'ZColor'));
                    if strcmp(get(ax,'ZGrid'),'on') && gridBehind
                        if axlimz(1)~=axlimz(2)
                            gridLines(fid, grouplabel, axpos, x, y, scolorname, gridlinestyle, linewidth, axlimz, axztick, axzindex_inner, selectedCorners, [2 3 4 5])
                            if strcmp(get(ax,'ZTickMode'),'auto') && strcmp(get(ax,'ZMinorGrid'),'on') && ~isempty(minor_axztick)
                                minorGridLines(fid, grouplabel, axpos, x, y, scolorname, minor_gridlinestyle, linewidth, axlimz, minor_axztick, selectedCorners, [2 3 4 5])                                                                
                            end
                        end
                    end
                    if projection.xyplane == false
                        if strcmp(get(ax,'Box'),'on')
                            line2svg(fid,grouplabel,axpos,[x(corners(k,2,p)) x(corners(k,3,p))],[y(corners(k,2,p)) y(corners(k,3,p))],scolorname,'-',linewidth);
                            line2svg(fid,grouplabel,axpos,[x(corners(k,4,p)) x(corners(k,5,p))],[y(corners(k,4,p)) y(corners(k,5,p))],scolorname,'-',linewidth);
                        else
                            if strcmp(get(ax,'ZGrid'),'on')
                                line2svg(fid,grouplabel,axpos,[x(corners(k,2,p)) x(corners(k,3,p))],[y(corners(k,2,p)) y(corners(k,3,p))],scolorname,gridlinestyle,linewidth);
                                line2svg(fid,grouplabel,axpos,[x(corners(k,4,p)) x(corners(k,5,p))],[y(corners(k,4,p)) y(corners(k,5,p))],scolorname,gridlinestyle,linewidth);
                            end
                        end
                    end
            end
        end
    end
end
fprintf(fid,'    <g>\n');
axchild=get(ax,'Children');
group = axchild2svg(fid,id,axIdString,ax,group,paperpos,axchild,axpos,groupax,projection,boundingBoxAxes);
fprintf(fid,'    </g>\n');
if strcmp(get(ax,'Visible'),'on')
    fprintf(fid,'    <g>\n');
    % Search axis for labeling
    if projection.xyplane
        [x_axis_point_value, x_axis_point_index_top] = min(y);
        [x_axis_point_value, x_axis_point_index_bottom] = max(y);
        if strcmp(get(ax,'Box'),'on')
            if strcmp(get(ax,'XAxisLocation'),'top')
                x_axis_point_index = [x_axis_point_index_top x_axis_point_index_bottom];
            else
                x_axis_point_index = [x_axis_point_index_bottom x_axis_point_index_top];
            end
        else
            if strcmp(get(ax,'XAxisLocation'),'top')
                x_axis_point_index = x_axis_point_index_top;
            else
                x_axis_point_index = x_axis_point_index_bottom;
            end
        end
        [y_axis_point_value, y_axis_point_index_left] = min(x);
        [y_axis_point_value, y_axis_point_index_right] = max(x);
        if strcmp(get(ax,'Box'),'on')
            if strcmp(get(ax,'YAxisLocation'),'right')
                y_axis_point_index = [y_axis_point_index_right y_axis_point_index_left];
            else
                y_axis_point_index = [y_axis_point_index_left y_axis_point_index_right];
            end
        else
            if strcmp(get(ax,'YAxisLocation'),'right')
                y_axis_point_index = y_axis_point_index_right;
            else
                y_axis_point_index = y_axis_point_index_left;
            end
        end
        [z_axis_point_value, z_axis_point_index] = min(x);   
    else
        [x_axis_point_value, x_axis_point_index] = max(y);
        [y_axis_value, y_axis_point_index] = max(y);
        [z_axis_point_value, z_axis_point_index] = min(x);   
    end
    % Draw grid
    for pindex = 1:size(front_faces)
        p = front_faces(pindex);
        for k = 1:size(corners,1)
            selectedCorners = squeeze(corners(k,:,p));
            switch corners(k,1,p)
                case 1 % x
                    % Draw x-grid
                    scolorname = searchcolor(id,get(ax,'XColor'));
                    if strcmp(get(ax,'XGrid'),'on') && gridBehind == false
                        if axlimx(1)~=axlimx(2)
                            gridLines(fid, grouplabel, axpos, x, y, scolorname, gridlinestyle, linewidth, axlimx, axxtick, axxindex_inner, selectedCorners, [2 3 4 5])
                            if strcmp(get(ax,'XTickMode'),'auto') && strcmp(get(ax,'XMinorGrid'),'on') && ~isempty(minor_axxtick)
                                minorGridLines(fid, grouplabel, axpos, x, y, scolorname, minor_gridlinestyle, linewidth, axlimx, minor_axxtick, selectedCorners, [2 3 4 5])
                            end
                        end
                    end
                case 2 % y
                    % Draw y-grid
                    scolorname = searchcolor(id,get(ax,'YColor'));
                    if strcmp(get(ax,'YGrid'),'on') && gridBehind == false
                        if axlimy(1)~=axlimy(2)
                            gridLines(fid, grouplabel, axpos, x, y, scolorname, gridlinestyle, linewidth, axlimy, axytick, axyindex_inner, selectedCorners, [2 3 4 5])
                            if strcmp(get(ax,'YTickMode'),'auto') && strcmp(get(ax,'YMinorGrid'),'on') && ~isempty(minor_axytick)
                                minorGridLines(fid, grouplabel, axpos, x, y, scolorname, minor_gridlinestyle, linewidth, axlimy, minor_axytick, selectedCorners, [2 3 4 5])                                                                
                            end
                        end
                    end
                case 3 % z
                    % Draw z-grid
                    scolorname = searchcolor(id,get(ax,'ZColor'));
                    if strcmp(get(ax,'ZGrid'),'on') && gridBehind == false
                        if axlimz(1)~=axlimz(2)
                            gridLines(fid, grouplabel, axpos, x, y, scolorname, gridlinestyle, linewidth, axlimz, axztick, axzindex_inner, selectedCorners, [2 3 4 5])
                            if strcmp(get(ax,'ZTickMode'),'auto') && strcmp(get(ax,'ZMinorGrid'),'on') && ~isempty(minor_axztick)
                                minorGridLines(fid, grouplabel, axpos, x, y, scolorname, minor_gridlinestyle, linewidth, axlimz, minor_axztick, selectedCorners, [2 3 4 5])                                                                
                            end
                        end
                    end
            end
        end
    end
    scolorname=searchcolor(id,get(ax,'XColor'));
    % Draw 'box' of x-axis
    if projection.xyplane == false
        if strcmp(get(ax,'Box'),'on')
            edge_line_index = [edge_opposite(most_back_edge_index) edge_neighbours(edge_opposite(most_back_edge_index),1)];
            line2svg(fid,grouplabel,axpos,x(edge_line_index),y(edge_line_index),scolorname,'-',linewidth)
        end
    end
    % Draw x-tick labels
    if (strcmp(get(ax,'XTickLabelMode'),'auto') && strcmp(get(ax,'XScale'),'log'))
        exponent = 1;
    else
        exponent = 0;
    end
    % Draw x-tick marks
    if (ticklength(1) ~= 0)
        if axlimx(1)~=axlimx(2)
            if (nomx(x_axis_point_index(1)))
                lim = [axlimx(2) axlimx(1)];    
            else
                lim = [axlimx(1) axlimx(2)];
            end
            x_label_end1 = interp1([0 1],[x(x_axis_point_index(1)) x(edge_neighbours(x_axis_point_index(1),2))],label_distance,'linear','extrap');
            y_label_end1 = interp1([0 1],[y(x_axis_point_index(1)) y(edge_neighbours(x_axis_point_index(1),2))],label_distance,'linear','extrap');
            x_label_end2 = interp1([0 1],[x(edge_neighbours(x_axis_point_index(1),1)) x(edge_neighbours(edge_neighbours(x_axis_point_index(1),1),2))],label_distance,'linear','extrap');
            y_label_end2 = interp1([0 1],[y(edge_neighbours(x_axis_point_index(1),1)) y(edge_neighbours(edge_neighbours(x_axis_point_index(1),1),2))],label_distance,'linear','extrap');
            xg_label_end = interp1(lim,[x_label_end1 x_label_end2],axxtick);
            yg_label_end = interp1(lim,[y_label_end1 y_label_end2],axxtick);            
            frontTicks(fid, grouplabel, axpos, x, y, scolorname, linewidth, ...
                axxtick, x_axis_point_index, edge_neighbours, [2 1 1], ...
                valid_xsticks,  ticklength, tick_ratio, lim, true);
            if strcmp(get(ax,'XTickMode'),'auto') && (strcmp(get(ax,'XMinorGrid'),'on') || strcmp(get(ax,'XScale'),'log')) && ~isempty(minor_axxtick)
                frontTicks(fid, grouplabel, axpos, x, y, scolorname, linewidth, ...
                    minor_axxtick, x_axis_point_index, edge_neighbours, [2 1 1], ...
                    1:length(minor_axxtick),  0.5 * ticklength, tick_ratio, lim, false);
            end
            if ~isempty(axlabelx) && ~(iscell(axlabelx) && all(cellfun(@isempty,axlabelx))) 
                if ischar(axlabelx) && size(axlabelx, 1) == 1
                    % Special handling of 1xn char arrays -> duplicate data
                    % for all ticks. Strange behavior but follows the
                    % behavior of Matlab
                    axlabelx = repmat(axlabelx, length(axxindex), 1);
                end
                % Note: 3D plot do not support the property XAxisLocation
                % setting 'top'.
                [angle, align] = improvedXLabel(ax, 0, 'Center');
                if strcmp(get(ax,'XAxisLocation'),'top') && (projection.xyplane == true)
                    for i = 1:length(axxindex)
                        label2svg(fid,grouplabel,axpos,ax,xg_label_end(i),yg_label_end(i),convertString(axlabelx(i,:)),align,angle,'bottom',1,paperpos,scolorname,exponent);
                    end
                else
                    for i = 1:length(axxindex)
                        label2svg(fid,grouplabel,axpos,ax,xg_label_end(i),yg_label_end(i),convertString(axlabelx(i,:)),align,angle,'top',1,paperpos,scolorname,exponent);
                    end
                end
            end
        end
    end
    scolorname=searchcolor(id,get(ax,'YColor'));
    % Draw 'box' of y-axis
    if projection.xyplane == false
        if strcmp(get(ax,'Box'),'on')
            edge_line_index = [edge_opposite(most_back_edge_index) edge_neighbours(edge_opposite(most_back_edge_index),2)];
            line2svg(fid,grouplabel,axpos,x(edge_line_index),y(edge_line_index),scolorname,'-',linewidth)
        end
    end
    % Draw y-tick labels
    if (strcmp(get(ax,'YTickLabelMode'),'auto') && strcmp(get(ax,'YScale'),'log'))
        exponent = 1;
    else
        exponent = 0;
    end
    % Draw y-tick marks
    if (ticklength(1) ~= 0)
        if axlimy(1)~=axlimy(2)
            if (nomy(y_axis_point_index(1)))
                lim = [axlimy(2) axlimy(1)];    
            else
                lim = [axlimy(1) axlimy(2)];
            end
            x_label_end1 = interp1([0 1],[x(y_axis_point_index(1)) x(edge_neighbours(y_axis_point_index(1),1))],label_distance,'linear','extrap');
            y_label_end1 = interp1([0 1],[y(y_axis_point_index(1)) y(edge_neighbours(y_axis_point_index(1),1))],label_distance,'linear','extrap');
            x_label_end2 = interp1([0 1],[x(edge_neighbours(y_axis_point_index(1),2)) x(edge_neighbours(edge_neighbours(y_axis_point_index(1),2),1))],label_distance,'linear','extrap');
            y_label_end2 = interp1([0 1],[y(edge_neighbours(y_axis_point_index(1),2)) y(edge_neighbours(edge_neighbours(y_axis_point_index(1),2),1))],label_distance,'linear','extrap');
            xg_label_end = interp1(lim,[x_label_end1 x_label_end2],axytick);
            yg_label_end = interp1(lim,[y_label_end1 y_label_end2],axytick);            
            frontTicks(fid, grouplabel, axpos, x, y, scolorname, linewidth, ...
                axytick, y_axis_point_index, edge_neighbours, [1 2 2], ...
                valid_ysticks,  ticklength, tick_ratio, lim, true);
            if strcmp(get(ax,'YTickMode'),'auto') && (strcmp(get(ax,'YMinorGrid'),'on') || strcmp(get(ax,'YScale'),'log')) && ~isempty(minor_axytick)
                frontTicks(fid, grouplabel, axpos, x, y, scolorname, linewidth, ...
                    minor_axytick, y_axis_point_index, edge_neighbours, [1 2 2], ...
                    1:length(minor_axytick), 0.5 * ticklength, tick_ratio, lim, false);
            end
            if ~isempty(axlabely) && ~(iscell(axlabely) && all(cellfun(@isempty,axlabely)))
                if ischar(axlabely) && size(axlabely, 1) == 1
                    % Special handling of 1xn char arrays -> duplicate data
                    % for all ticks. Strange behavior but follows the
                    % behavior of Matlab
                    axlabely = repmat(axlabely, length(axyindex), 1);
                end
                % Note: 3D plot do not support the property YAxisLocation
                % setting 'right'.
                if (projection.xyplane == true)
                    if strcmp(get(ax,'YAxisLocation'),'right')
                        [angle, align] = improvedYLabel(ax, 0, 'Left');
                        for i = 1:length(axyindex)
                            label2svg(fid,grouplabel,axpos,ax,xg_label_end(i),yg_label_end(i),convertString(axlabely(i,:)),align,angle,'middle',1,paperpos,scolorname,exponent);
                        end
                   else
                        [angle, align] = improvedYLabel(ax, 0, 'Right');
                        for i = 1:length(axyindex)
                            label2svg(fid,grouplabel,axpos,ax,xg_label_end(i),yg_label_end(i),convertString(axlabely(i,:)),align,angle,'middle',1,paperpos,scolorname,exponent);
                        end
                    end
                else
                    for i = 1:length(axyindex)
                        label2svg(fid,grouplabel,axpos,ax,xg_label_end(i),yg_label_end(i),convertString(axlabely(i,:)),'Center',0,'top',1,paperpos,scolorname,exponent);
                    end
                end
            end
        end
    end
    scolorname=searchcolor(id,get(ax,'ZColor'));
    % Draw 'box' of z-axis
    if projection.xyplane == false
        if strcmp(get(ax,'Box'),'on')
            edge_line_index = [edge_opposite(most_back_edge_index) edge_neighbours(edge_opposite(most_back_edge_index),3)];
            line2svg(fid,grouplabel,axpos,x(edge_line_index),y(edge_line_index),scolorname,'-',linewidth)
        end
    end
    if (strcmp(get(ax,'ZTickLabelMode'),'auto') && strcmp(get(ax,'ZScale'),'log'))
        exponent = 1;
    else
        exponent = 0;
    end
    % Draw z-tick marks
    if (ticklength(1) ~= 0)
        if axlimz(1)~=axlimz(2)
            if (nomz(z_axis_point_index(1)))
                lim = [axlimz(2) axlimz(1)];    
            else
                lim = [axlimz(1) axlimz(2)];
            end
            x_tick_end1 = interp1([0 1],[x(z_axis_point_index) x(edge_neighbours(z_axis_point_index,2))],ticklength*tick_ratio(3),'linear','extrap');
            y_tick_end1 = interp1([0 1],[y(z_axis_point_index) y(edge_neighbours(z_axis_point_index,2))],ticklength*tick_ratio(3),'linear','extrap');
            x_tick_end2 = interp1([0 1],[x(edge_neighbours(z_axis_point_index,3)) x(edge_neighbours(edge_neighbours(z_axis_point_index,3),2))],ticklength*tick_ratio(3),'linear','extrap');
            y_tick_end2 = interp1([0 1],[y(edge_neighbours(z_axis_point_index,3)) y(edge_neighbours(edge_neighbours(z_axis_point_index,3),2))],ticklength*tick_ratio(3),'linear','extrap');
            x_label_end1 = interp1([0 1],[x(z_axis_point_index) x(edge_neighbours(z_axis_point_index,2))],label_distance,'linear','extrap');
            y_label_end1 = interp1([0 1],[y(z_axis_point_index) y(edge_neighbours(z_axis_point_index,2))],label_distance,'linear','extrap');
            x_label_end2 = interp1([0 1],[x(edge_neighbours(z_axis_point_index,3)) x(edge_neighbours(edge_neighbours(z_axis_point_index,3),2))],label_distance,'linear','extrap');
            y_label_end2 = interp1([0 1],[y(edge_neighbours(z_axis_point_index,3)) y(edge_neighbours(edge_neighbours(z_axis_point_index,3),2))],label_distance,'linear','extrap');
            xg_line_start = interp1(lim,[x(z_axis_point_index) x(edge_neighbours(z_axis_point_index,3))],axztick);
            yg_line_start = interp1(lim,[y(z_axis_point_index) y(edge_neighbours(z_axis_point_index,3))],axztick);
            xg_line_end = interp1(lim,[x_tick_end1 x_tick_end2],axztick);
            yg_line_end = interp1(lim,[y_tick_end1 y_tick_end2],axztick);
            xg_label_end = interp1(lim,[x_label_end1 x_label_end2],axztick);
            yg_label_end = interp1(lim,[y_label_end1 y_label_end2],axztick);            
            for i = valid_zsticks
                line2svg(fid,grouplabel,axpos,[xg_line_start(i) xg_line_end(i)],[yg_line_start(i) yg_line_end(i)],scolorname,'-',linewidth)
            end
            line2svg(fid,grouplabel,axpos,[x(z_axis_point_index) x(edge_neighbours(z_axis_point_index,3))],[y(z_axis_point_index) y(edge_neighbours(z_axis_point_index,3))],scolorname,'-',linewidth)
            if ~isempty(axlabelz) && ~(iscell(axlabelz) && all(cellfun(@isempty,axlabelz)))
                if ischar(axlabelz) && size(axlabelz, 1) == 1
                    % Special handling of 1xn char arrays -> duplicate data
                    % for all ticks. Strange behavior but follows the
                    % behavior of Matlab
                    axlabelz = repmat(axlabelz, length(axzindex), 1);
                end
                for i = 1:length(axzindex)
                    label2svg(fid,grouplabel,axpos,ax,xg_label_end(i),yg_label_end(i),convertString(axlabelz(i,:)),'Right',0,'middle',1,paperpos,scolorname,exponent);
                end
            end
        end
    end
    exponent2svg(fid,groupax,axpos,paperpos,ax,axxtick,axytick,axztick)
    fprintf(fid,'    </g>\n');
end
fprintf(fid,'  </g>\n');
set(ax,'Units',originalAxesUnits);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% take any axis children and create objects for them
function group=axchild2svg(fid,id,axIdString,ax,group,paperpos,axchild,axpos,groupax,projection,boundingBoxAxes)
global colorname
global PLOT2SVG_globals
for i=length(axchild):-1:1
    if strcmp(get(axchild(i), 'Visible'), 'off')
        % do nothing
    elseif strcmp(get(axchild(i),'Type'),'line')
        scolorname=searchcolor(id,get(axchild(i),'Color'));
        linestyle=get(axchild(i),'LineStyle');
        linewidth=get(axchild(i),'LineWidth');
        marker=get(axchild(i),'Marker');
        markeredgecolor=get(axchild(i),'MarkerEdgeColor');
        if ischar(markeredgecolor)
            switch markeredgecolor
                case 'none',markeredgecolorname='none';
                otherwise,markeredgecolorname=scolorname;  % if markeredgecolorname is 'auto' or something else set the markeredgecolorname to the line color
            end    
        else    
            markeredgecolorname=searchcolor(id,markeredgecolor);
        end
        markerfacecolor=get(axchild(i),'MarkerFaceColor');
        if ischar(markerfacecolor)
            switch markerfacecolor
                case 'none',markerfacecolorname='none';
                otherwise,markerfacecolorname=scolorname;  % if markerfacecolorname is 'auto' or something else set the markerfacecolorname to the line color
            end
        else
            markerfacecolorname=searchcolor(id,markerfacecolor);
        end
        markersize=get(axchild(i),'MarkerSize')/1.5;
        linex = get(axchild(i),'XData');
        linex = linex(:)'; % Octave stores the data in a column vector
        if strcmp(get(ax,'XScale'),'log')
            linex(find(linex<=0)) = NaN;
            linex=log10(linex);
        end
        liney=get(axchild(i),'YData');
        liney = liney(:)'; % Octave stores the data in a column vector        
        if strcmp(get(ax,'YScale'),'log')
            liney(find(liney<=0)) = NaN;
            liney=log10(liney);
        end
        linez=get(axchild(i),'ZData');
        linez = linez(:)'; % Octave stores the data in a column vector        
        if isempty(linez)
            linez = zeros(size(linex));    
        end
        if strcmp(get(ax,'ZScale'),'log')
            linez(find(linez<=0)) = NaN;
            linez=log10(linez);
        end
        [x,y,z] = project(linex,liney,linez,projection);
        x = (x*axpos(3)+axpos(1))*paperpos(3);
        y = (1-(y*axpos(4)+axpos(2)))*paperpos(4);
        markerOverlap = 0;
        if ~strcmp(linestyle, 'none')
            markerOverlap = max(markerOverlap, convertunit(linewidth*0.5, 'points', 'pixels', axpos(4)));    
        end
        if ~strcmp(marker, 'none')
            markerOverlap = max(markerOverlap, convertunit(markersize, 'points', 'pixels', axpos(4)));    
        end
        boundingBoxElement = [min(x)-markerOverlap min(y)-markerOverlap max(x)-min(x)+2*markerOverlap max(y)-min(y)+2*markerOverlap];
        [filterString, boundingBox] = filter2svg(fid, axchild(i), boundingBoxAxes, boundingBoxElement);
        % put a line into a group with its markers
        if strcmp(get(axchild(i),'Clipping'),'on')
            clippingIdString = clipping2svg(fid, axchild(i), ax, paperpos, axpos, projection, axIdString);
            fprintf(fid,'<g id="%s" clip-path="url(#%s)" %s>\n', createId, clippingIdString, filterString);
        else
            fprintf(fid,'<g id="%s" %s>\n', createId, filterString);
        end
        if ~isempty(filterString)
            % Workaround for Inkscape filter bug
            fprintf(fid,'<rect x="%0.3f" y="%0.3f" width="%0.3f" height="%0.3f" fill="none" stroke="none" />\n', boundingBox(1), boundingBox(2), boundingBox(3), boundingBox(4));
        end
        line2svg(fid,groupax,axpos,x,y,scolorname,linestyle,linewidth)
        % put the markers into a subgroup of the lines
        fprintf(fid,'<g>\n');
        switch marker
            case 'none';
            case '.',group=group+1;circle2svg(fid,group,axpos,x,y,markersize*0.25,'none',markeredgecolorname,linewidth);
            case 'o',group=group+1;circle2svg(fid,group,axpos,x,y,markersize*0.75,markeredgecolorname,markerfacecolorname,linewidth);
            case '+',group=group+1;patch2svg(fid,group,axpos,x'*ones(1,5)+ones(length(linex),1)*[-1 1 NaN 0 0]*markersize,y'*ones(1,5)+ones(length(liney),1)*[0 0 NaN -1 1]*markersize,markeredgecolorname,'-',linewidth,markeredgecolorname, 1, 1, false);   
            case '*',group=group+1;patch2svg(fid,group,axpos,x'*ones(1,11)+ones(length(linex),1)*[-1 1 NaN 0 0 NaN -0.7 0.7 NaN -0.7 0.7]*markersize,y'*ones(1,11)+ones(length(liney),1)*[0 0 NaN -1 1 NaN 0.7 -0.7 NaN -0.7 0.7]*markersize,markeredgecolorname,'-',linewidth,markeredgecolorname, 1, 1, false);
            case 'x',group=group+1;patch2svg(fid,group,axpos,x'*ones(1,5)+ones(length(linex),1)*[-0.7 0.7 NaN -0.7 0.7]*markersize,y'*ones(1,5)+ones(length(liney),1)*[0.7 -0.7 NaN -0.7 0.7]*markersize,markeredgecolorname,'-',linewidth,markeredgecolorname, 1, 1, false);
            %% Octave keeps s, d, p and h in the HandleGraphics object, for the square, diamond, pentagram, and hexagram markers, respectively -- Jakob Malm
            case {'square', 's'},group=group+1;patch2svg(fid,group,axpos,x'*ones(1,5)+ones(length(linex),1)*[-1 -1 1 1 -1]*markersize,y'*ones(1,5)+ones(length(liney),1)*[-1 1 1 -1 -1]*markersize,markerfacecolorname,'-',linewidth,markeredgecolorname, 1, 1, true);
            case {'diamond', 'd'},group=group+1;patch2svg(fid,group,axpos,x'*ones(1,5)+ones(length(linex),1)*[-0.7071 0 0.7071 0 -0.7071]*markersize,y'*ones(1,5)+ones(length(liney),1)*[0 1 0 -1 0]*markersize,markerfacecolorname,'-',linewidth,markeredgecolorname, 1, 1, true);
            case {'pentagram', 'p'},group=group+1;patch2svg(fid,group,axpos,...
                    x'*ones(1,11)+ones(length(linex),1)*[0 0.1180 0.5 0.1910 0.3090 0 -0.3090 -0.1910 -0.5 -0.1180 0]*1.3*markersize,...
                    y'*ones(1,11)+ones(length(liney),1)*[-0.5257 -0.1625 -0.1625 0.0621 0.4253 0.2008 0.4253 0.0621 -0.1625 -0.1625 -0.5257]*1.3*markersize,markerfacecolorname,'-',linewidth,markeredgecolorname, 1, 1, true);
            case {'hexagram', 'h'},group=group+1;patch2svg(fid,group,axpos,...
                    x'*ones(1,13)+ones(length(linex),1)*[0 0.2309 0.6928 0.4619 0.6928 0.2309 0 -0.2309 -0.6928 -0.4619 -0.6928 -0.2309 0]*1*markersize,...
                    y'*ones(1,13)+ones(length(liney),1)*[0.8 0.4 0.4 0 -0.4 -0.4 -0.8 -0.4 -0.4 0 0.4 0.4 0.8]*1*markersize,markerfacecolorname,'-',linewidth,markeredgecolorname, 1, 1, true);    
            case '^',group=group+1;patch2svg(fid,group,axpos,x'*ones(1,4)+ones(length(linex),1)*[-1 1 0 -1]*markersize,y'*ones(1,4)+ones(length(liney),1)*[0.577 0.577 -1.155 0.577]*markersize,markerfacecolorname,'-',linewidth,markeredgecolorname, 1, 1, true);
            case 'v',group=group+1;patch2svg(fid,group,axpos,x'*ones(1,4)+ones(length(linex),1)*[-1 1 0 -1]*markersize,y'*ones(1,4)+ones(length(liney),1)*[-0.577 -0.577 1.155 -0.577]*markersize,markerfacecolorname,'-',linewidth,markeredgecolorname, 1, 1, true);
            case '<',group=group+1;patch2svg(fid,group,axpos,x'*ones(1,4)+ones(length(linex),1)*[0.577 0.577 -1.155 0.577]*markersize,y'*ones(1,4)+ones(length(liney),1)*[-1 1 0 -1]*markersize,markerfacecolorname,'-',linewidth,markeredgecolorname, 1, 1, true);
            case '>',group=group+1;patch2svg(fid,group,axpos,x'*ones(1,4)+ones(length(linex),1)*[-0.577 -0.577 1.155 -0.577]*markersize,y'*ones(1,4)+ones(length(liney),1)*[-1 1 0 -1]*markersize,markerfacecolorname,'-',linewidth,markeredgecolorname, 1, 1, true);
        end
        % close the marker group
        fprintf(fid,'</g>\n');
        animation2svg(fid, axchild(i));
        % close the line group
        fprintf(fid,'</g>\n');
    elseif strcmp(get(axchild(i),'Type'),'patch')
        flat_shading = 1;
        cmap=get(id,'Colormap');
        pointc = get(axchild(i),'FaceVertexCData');
        if isempty(pointc)
			% Workaround for octave
			pointc=get(axchild(i),'CData');
        end
        % Scale color if scaled color mapping is turned on
        if strcmp(get(axchild(i),'CDataMapping'),'scaled')
            clim=get(ax,'CLim');
            pointc=(pointc-clim(1))/(clim(2)-clim(1))*(size(cmap,1)-1)+1;
        end
        % Limit index to smallest or biggest color index
        pointc=max(pointc,1);
        pointc=min(pointc,size(cmap,1));
        if ~ischar(get(axchild(i),'FaceAlpha'))
            face_opacity = get(axchild(i),'FaceAlpha');
        else
            face_opacity = 1.0;
        end
        if ~ischar(get(axchild(i),'EdgeAlpha'))
            edge_opacity = get(axchild(i),'EdgeAlpha');
        else
            edge_opacity = 1.0;
        end
        linestyle = get(axchild(i),'LineStyle');
        linewidth = get(axchild(i),'LineWidth');
        marker = get(axchild(i),'Marker');
        markeredgecolor=get(axchild(i),'MarkerEdgeColor');
        markersize=get(axchild(i),'MarkerSize')/1.5;
        points=get(axchild(i),'Vertices')';
        if strcmp(get(ax,'XScale'),'log')
            points(1,:)=log10(points(1,:));
        end
        if strcmp(get(ax,'YScale'),'log')
            points(2,:)=log10(points(2,:));
        end
        % TODO LogZ
        if size(points,1)==2
            [x,y,z] = project(points(1,:),points(2,:),zeros(size(points(1,:))),projection);    
        else
            [x,y,z] = project(points(1,:),points(2,:),points(3,:),projection);
        end
        x = (x*axpos(3)+axpos(1))*paperpos(3);
        y = (1-(y*axpos(4)+axpos(2)))*paperpos(4);
        faces = get(axchild(i),'Faces');
        face_index = 1:size(faces,1);
        if size(points,1)==3;
            [z,face_index]=sort(sum(z(faces(:,:)),2));
            faces=faces(face_index,:);
        end
        markerOverlap = 0;
        if ~strcmp(linestyle, 'none')
            markerOverlap = max(markerOverlap, convertunit(linewidth*0.5, 'points', 'pixels', axpos(4)));    
        end
        if ~strcmp(marker, 'none')
            markerOverlap = max(markerOverlap, convertunit(markersize, 'points', 'pixels', axpos(4)));    
        end
        boundingBoxElement = [min(x)-markerOverlap min(y)-markerOverlap max(x)-min(x)+2*markerOverlap max(y)-min(y)+2*markerOverlap];
        [filterString, boundingBox] = filter2svg(fid, axchild(i), boundingBoxAxes, boundingBoxElement);
        if strcmp(get(axchild(i),'Clipping'),'on')
            clippingIdString = clipping2svg(fid, axchild(i), ax, paperpos, axpos, projection, axIdString);
            fprintf(fid,'<g id="%s" clip-path="url(#%s)" %s>\n', createId, clippingIdString, filterString);
        else
            fprintf(fid,'<g id="%s" %s>\n', createId, filterString);
        end
        if ~isempty(filterString)
            % Workaround for Inkscape filter bug
            fprintf(fid,'<rect x="%0.3f" y="%0.3f" width="%0.3f" height="%0.3f" fill="none" stroke="none" />\n', boundingBox(1), boundingBox(2), boundingBox(3), boundingBox(4));
        end
        for p = 1:size(faces,1)
            if ischar(get(axchild(i),'FaceColor'))
                if strcmp(get(axchild(i),'FaceColor'),'texturemap')
                    facecolorname='none';   % TO DO: texture map
                elseif strcmp(get(axchild(i),'FaceColor'),'none')
                    facecolorname='none';
                else
                    if size(pointc,1)==1
                        facecolor = pointc;    
                    elseif size(pointc,1)==size(faces,1)
                        if strcmp(get(axchild(i),'FaceColor'),'flat')
                            facecolor = pointc(face_index(p),:);
                        else
                            facecolor = pointc(face_index(p),:);
                            cdata = pointc(face_index(p),:);     % TO DO: color interpolation
                            flat_shading = 0;
                        end
                    elseif size(pointc,1)==size(points,2)
                        if strcmp(get(axchild(i),'FaceColor'),'flat')
                            facecolor = pointc(faces(p,1),:);
                        else
                            facecolor = pointc(faces(p,1));
                            cdata = pointc(faces(p,:),:);
                            flat_shading = 0;
                        end
                    else
                        error('Unsupported color handling for patches.');    
                    end
                    if ~isnan(facecolor)
                        if size(facecolor,2)==1
                            facecolorname = ['#' colorname(ceil(facecolor),:)];
                        else
                            if strcmp(get(axchild(i),'FaceColor'),'flat')  % Bugfix 27.01.2008
                                facecolorname = searchcolor(id,facecolor/64);
                            else
                                facecolorname = searchcolor(id,facecolor);    
                            end
                        end
                    else
                        facecolorname='none';
                    end
                end
            else
                facecolorname = searchcolor(id,get(axchild(i),'FaceColor'));       
            end
            if ischar(get(axchild(i),'EdgeColor'))
                if strcmp(get(axchild(i),'EdgeColor'),'none')
                    edgecolorname = 'none';
                else
                    if size(pointc,1)==1
                        edgecolor = pointc;    
                    elseif size(pointc,1)==size(faces,1)
                        edgecolor = pointc(p,:);
                    elseif size(pointc,1)==size(points,2)
                        if strcmp(get(axchild(i),'EdgeColor'),'flat')
                            edgecolor = pointc(faces(p,1));
                        else
                            edgecolor = pointc(faces(p,1));     % TO DO: color interpolation
                        end
                    else
                        error('Unsupported color handling for patches.');    
                    end
                    if ~isnan(edgecolor)
                        if size(edgecolor,2)==1
                            edgecolorname = ['#' colorname(ceil(edgecolor),:)];
                        else
                            if strcmp(get(axchild(i),'EdgeColor'),'flat')   % Bugfix 27.01.2008
                                edgecolorname = searchcolor(id,edgecolor/64);
                            else
                                edgecolorname = searchcolor(id,edgecolor);    
                            end
                        end
                    else
                        edgecolorname = 'none';
                    end
                end
            else
                edgecolorname = searchcolor(id,get(axchild(i),'EdgeColor'));       
            end
            % Close a patch if the coordinates do not contain NaNs 
            if any(isnan(x)) || any(isnan(y))
                closed = false;
            else
                closed = true;
            end
            if flat_shading
                patch2svg(fid, group, axpos, x(faces(p,:)), y(faces(p,:)), facecolorname, linestyle, linewidth, edgecolorname, face_opacity, edge_opacity, closed)
            else
                gouraud_patch2svg(fid, group, axpos, x(faces(p,:)), y(faces(p,:)), cdata, linestyle, linewidth, edgecolorname, face_opacity, edge_opacity, id)
            end
            if ~strcmp(marker, 'none')
                for q = 1:size(faces,2)
                    xmarker = x(faces(p,q));
                    ymarker = y(faces(p,q));
                    % put the markers into a subgroup of the lines
                    fprintf(fid,'<g>\n');
                    if ischar(markeredgecolor)
                        switch markeredgecolor
                            case 'none',markeredgecolorname='none';
                            otherwise,
                                % if markeredgecolorname is 'auto' or something
                                % else set the markeredgecolorname to the line color
                                markeredgecolorname = selectColor(axchild(i), id, p, q, points, pointc, ...
                                    colorname, faces, 'MarkerEdgeColor');
                        end
                    else
                        markeredgecolorname=searchcolor(id,markeredgecolor);
                    end
                    markerfacecolor=get(axchild(i),'MarkerFaceColor');
                    if ischar(markerfacecolor)
                        switch markerfacecolor
                            case 'none',markerfacecolorname='none';
                            otherwise,
                                markerfacecolorname = selectColor(axchild(i), id, p, q, points, pointc, ...
                                    colorname, faces,'MarkerFaceColor');
                        end
                    else
                        markerfacecolorname=searchcolor(id,markerfacecolor);
                    end
                    switch marker
                        case 'none';
                        case '.',group=group+1;circle2svg(fid,group,axpos,xmarker,ymarker,markersize*0.25,'none',markeredgecolorname,linewidth);
                        case 'o',group=group+1;circle2svg(fid,group,axpos,xmarker,ymarker,markersize*0.75,markeredgecolorname,markerfacecolorname,linewidth);
                        case '+',group=group+1;patch2svg(fid,group,axpos,xmarker'*ones(1,5)+ones(length(linex),1)*[-1 1 NaN 0 0]*markersize,ymarker'*ones(1,5)+ones(length(liney),1)*[0 0 NaN -1 1]*markersize,markeredgecolorname,'-',linewidth,markeredgecolorname, 1, 1, false);
                        case '*',group=group+1;patch2svg(fid,group,axpos,xmarker'*ones(1,11)+ones(length(linex),1)*[-1 1 NaN 0 0 NaN -0.7 0.7 NaN -0.7 0.7]*markersize,ymarker'*ones(1,11)+ones(length(liney),1)*[0 0 NaN -1 1 NaN 0.7 -0.7 NaN -0.7 0.7]*markersize,markeredgecolorname,'-',linewidth,markeredgecolorname, 1, 1, false);
                        case 'x',group=group+1;patch2svg(fid,group,axpos,xmarker'*ones(1,5)+ones(length(linex),1)*[-0.7 0.7 NaN -0.7 0.7]*markersize,ymarker'*ones(1,5)+ones(length(liney),1)*[0.7 -0.7 NaN -0.7 0.7]*markersize,markeredgecolorname,'-',linewidth,markeredgecolorname, 1, 1, false);
                            %% Octave keeps s, d, p and h in the HandleGraphics object, for the square, diamond, pentagram, and hexagram markers, respectively -- Jakob Malm
                        case {'square', 's'},group=group+1;patch2svg(fid,group,axpos,xmarker'*ones(1,5)+ones(length(linex),1)*[-1 -1 1 1 -1]*markersize,ymarker'*ones(1,5)+ones(length(liney),1)*[-1 1 1 -1 -1]*markersize,markerfacecolorname,'-',linewidth,markeredgecolorname, 1, 1, true);
                        case {'diamond', 'd'},group=group+1;patch2svg(fid,group,axpos,xmarker'*ones(1,5)+ones(length(linex),1)*[-0.7071 0 0.7071 0 -0.7071]*markersize,ymarker'*ones(1,5)+ones(length(liney),1)*[0 1 0 -1 0]*markersize,markerfacecolorname,'-',linewidth,markeredgecolorname, 1, 1, true);
                        case {'pentagram', 'p'},group=group+1;patch2svg(fid,group,axpos,...
                                xmarker'*ones(1,11)+ones(length(linex),1)*[0 0.1180 0.5 0.1910 0.3090 0 -0.3090 -0.1910 -0.5 -0.1180 0]*1.3*markersize,...
                                ymarker'*ones(1,11)+ones(length(liney),1)*[-0.5257 -0.1625 -0.1625 0.0621 0.4253 0.2008 0.4253 0.0621 -0.1625 -0.1625 -0.5257]*1.3*markersize,markerfacecolorname,'-',linewidth,markeredgecolorname, 1, 1, true);
                        case {'hexagram', 'h'},group=group+1;patch2svg(fid,group,axpos,...
                                xmarker'*ones(1,13)+ones(length(linex),1)*[0 0.2309 0.6928 0.4619 0.6928 0.2309 0 -0.2309 -0.6928 -0.4619 -0.6928 -0.2309 0]*1*markersize,...
                                ymarker'*ones(1,13)+ones(length(liney),1)*[0.8 0.4 0.4 0 -0.4 -0.4 -0.8 -0.4 -0.4 0 0.4 0.4 0.8]*1*markersize,markerfacecolorname,'-',linewidth,markeredgecolorname, 1, 1, true);
                        case '^',group=group+1;patch2svg(fid,group,axpos,xmarker'*ones(1,4)+ones(length(linex),1)*[-1 1 0 -1]*markersize,ymarker'*ones(1,4)+ones(length(liney),1)*[0.577 0.577 -1.155 0.577]*markersize,markerfacecolorname,'-',linewidth,markeredgecolorname, 1, 1, true);
                        case 'v',group=group+1;patch2svg(fid,group,axpos,xmarker'*ones(1,4)+ones(length(linex),1)*[-1 1 0 -1]*markersize,ymarker'*ones(1,4)+ones(length(liney),1)*[-0.577 -0.577 1.155 -0.577]*markersize,markerfacecolorname,'-',linewidth,markeredgecolorname, 1, 1, true);
                        case '<',group=group+1;patch2svg(fid,group,axpos,xmarker'*ones(1,4)+ones(length(linex),1)*[0.577 0.577 -1.155 0.577]*markersize,ymarker'*ones(1,4)+ones(length(liney),1)*[-1 1 0 -1]*markersize,markerfacecolorname,'-',linewidth,markeredgecolorname, 1, 1, true);
                        case '>',group=group+1;patch2svg(fid,group,axpos,xmarker'*ones(1,4)+ones(length(linex),1)*[-0.577 -0.577 1.155 -0.577]*markersize,ymarker'*ones(1,4)+ones(length(liney),1)*[-1 1 0 -1]*markersize,markerfacecolorname,'-',linewidth,markeredgecolorname, 1, 1, true);
                    end
                    % close the marker group
                    fprintf(fid,'</g>\n');
                end
            end
        end
        fprintf(fid,'</g>\n');
    elseif strcmp(get(axchild(i),'Type'),'surface')
        flat_shading = 1;
        cmap=get(id,'Colormap');
        [faces,points,pointc,alpha]=surface2patch(axchild(i));
        points=points';
        % Scale color if scaled color mapping is turned on
        if strcmp(get(axchild(i),'CDataMapping'),'scaled')
            clim=get(ax,'CLim');
            pointc=(pointc-clim(1))/(clim(2)-clim(1))*(size(cmap,1)-1)+1;
        end
        % Limit index to smallest or biggest color index
        pointc=max(pointc,1);
        if ~ischar(get(axchild(i),'FaceAlpha'))
            face_opacity = get(axchild(i),'FaceAlpha');
        elseif strcmp(get(axchild(i),'FaceAlpha'),'flat')
            face_opacity = alpha;
            switch get(axchild(i),'AlphaDataMapping')
                case {'direct'}
                    face_opacity = 1.0; % TODO
                case {'scaled'}
                    alim=get(ax,'ALim');
                    face_opacity=(face_opacity-alim(1))/(alim(2)-alim(1));
                case {'none'}
                    % Clip alpha data
                    face_opacity = min(1, face_opacity);
                    face_opacity = max(0, face_opacity);
                otherwise
                    error(['Unsupported AlphaDataMapping identifier ''' get(axchild(i),'AlphaDataMapping') '''.']);
            end
        else
            face_opacity = 1.0;
        end
        if ~ischar(get(axchild(i),'EdgeAlpha'))
            edge_opacity = get(axchild(i),'EdgeAlpha');
        else
            edge_opacity = 1.0;
        end
        pointc=min(pointc,size(cmap,1));
        linestyle=get(axchild(i),'LineStyle');
        linewidth=get(axchild(i),'LineWidth');
        if strcmp(get(ax,'XScale'),'log')
            points(1,:)=log10(points(1,:));
        end
        if strcmp(get(ax,'YScale'),'log')
            points(2,:)=log10(points(2,:));
        end
        if size(points,1)==3
            if strcmp(get(ax,'ZScale'),'log')
                points(3,:)=log10(points(3,:));
            end   
        end
        if size(points,1)==3
            [x,y,z] = project(points(1,:),points(2,:),points(3,:),projection); 
        else
            [x,y,z] = project(points(1,:),points(2,:),zeros(size(points(1,:))),projection); 
        end
        x = (x*axpos(3)+axpos(1))*paperpos(3);
        y = (1-(y*axpos(4)+axpos(2)))*paperpos(4);
        face_index = 1:size(faces,1);
        if size(points,1)==3;
            [z,face_index]=sort(sum(z(faces(:,:)),2));
            faces=faces(face_index,:);
        end
        markerOverlap = 0;
        if ~strcmp(linestyle, 'none')
            markerOverlap = max(markerOverlap, convertunit(linewidth*0.5, 'points', 'pixels', axpos(4)));    
        end
        boundingBoxElement = [min(x)-markerOverlap min(y)-markerOverlap max(x)-min(x)+2*markerOverlap max(y)-min(y)+2*markerOverlap];
        [filterString, boundingBox] = filter2svg(fid, axchild(i), boundingBoxAxes, boundingBoxElement);
        if strcmp(get(axchild(i),'Clipping'),'on')
            clippingIdString = clipping2svg(fid, axchild(i), ax, paperpos, axpos, projection, axIdString);
            fprintf(fid,'<g clip-path="url(#%s)" id="%s" %s>\n',clippingIdString, createId, filterString);
        else
            fprintf(fid,'<g id="%s" %s>\n', createId, filterString);
        end
        if ~isempty(filterString)
            % Workaround for Inkscape filter bug
            fprintf(fid,'<rect x="%0.3f" y="%0.3f" width="%0.3f" height="%0.3f" fill="none" stroke="none" />\n', boundingBox(1), boundingBox(2), boundingBox(3), boundingBox(4));
        end
        for p=1:size(faces,1)
            if ischar(get(axchild(i),'FaceColor'))
                if strcmp(get(axchild(i),'FaceColor'),'texturemap')
                    facecolorname='none';   % TO DO: texture map
                elseif strcmp(get(axchild(i),'FaceColor'),'none')
                    facecolorname='none';
                else
                    if size(pointc,1)==1
                        facecolor = pointc;    
                    elseif size(pointc,1)==size(faces,1)
                        facecolor = pointc(face_index(p),:);
                    elseif size(pointc,1)==size(points,2)
                        if strcmp(get(axchild(i),'FaceColor'),'flat')
                            facecolor = pointc(faces(p,1));
                        else
                            facecolor = pointc(faces(p,1));
                            cdata = pointc(faces(p,:));
                            flat_shading = 0;
                        end
                    else
                        error('Unsupported color handling for patches.');    
                    end
                    if ~isnan(facecolor)
                        if size(facecolor,2)==1
                            facecolorname = ['#' colorname(ceil(facecolor),:)];
                        else
                            facecolorname = searchcolor(id,facecolor);    
                        end
                    else
                        facecolorname='none';
                    end
                end
            else
                facecolorname = searchcolor(id,get(axchild(i),'FaceColor'));       
            end
            if size(face_opacity,1)==1
                face_opacity_value = face_opacity;
            elseif size(face_opacity,1)==size(faces,1)
                face_opacity_value = face_opacity(p,:);
            elseif size(face_opacity,1)==size(points,2)
                face_opacity_value = face_opacity(faces(p,1));
            else
                error('Unsupported face alpha value handling for patches.');
            end
            if ischar(get(axchild(i),'EdgeColor'))
                if strcmp(get(axchild(i),'EdgeColor'),'none')
                    edgecolorname = 'none';
                else
                    if size(pointc,1)==1
                        edgecolor = pointc;    
                    elseif size(pointc,1)==size(faces,1)
                        edgecolor = pointc(p,:);
                    elseif size(pointc,1)==size(points,2)
                        if strcmp(get(axchild(i),'EdgeColor'),'flat')
                            edgecolor = pointc(faces(p,1));
                        else
                            edgecolor = pointc(faces(p,1));     % TO DO: color interpolation
                        end
                    else
                        error('Unsupported color handling for patches.');    
                    end
                    if ~isnan(edgecolor)
                        if size(edgecolor,2)==1
                            edgecolorname = ['#' colorname(ceil(edgecolor),:)];
                        else
                            edgecolorname = searchcolor(id,edgecolor);    
                        end
                    else
                        edgecolorname = 'none';
                    end
                end
            else
                edgecolorname = searchcolor(id,get(axchild(i),'EdgeColor'));       
            end
            if flat_shading
                patch2svg(fid, group, axpos, x(faces(p,:)), y(faces(p,:)), facecolorname, linestyle, linewidth, edgecolorname,  face_opacity_value, edge_opacity, false)
            else
                gouraud_patch2svg(fid, group, axpos, x(faces(p,:)), y(faces(p,:)), cdata, linestyle, linewidth, edgecolorname, face_opacity_value, edge_opacity,id)
            end
        end
        fprintf(fid,'</g>\n');
    elseif strcmp(get(axchild(i),'Type'),'rectangle')
        scolorname = searchcolor(id,get(axchild(i),'EdgeColor'));
        fcolorname = searchcolor(id,get(axchild(i),'FaceColor'));
        linewidth = get(axchild(i),'LineWidth');
        position = get(axchild(i),'Position');
        posx = [position(1) position(1)+position(3)];
        if strcmp(get(ax,'XScale'),'log')
            posx(find(posx <= 0)) = NaN;
            posx=log10(posx);
        end
        posy = [position(2) position(2)+position(4)];
        if strcmp(get(ax,'YScale'),'log')
            posy(find(posy <= 0)) = NaN;
            posy=log10(posy);
        end
        posz=[0 0];
        linestyle = get(axchild(i),'LineStyle');
        if strcmp(linestyle,'none');
            scolorname = 'none';
        end
        pattern = lineStyle2svg(linestyle, linewidth);
        [x,y,z] = project(posx,posy,posz,projection);
        x = (x*axpos(3)+axpos(1))*paperpos(3);
        y = (1-(y*axpos(4)+axpos(2)))*paperpos(4);
        rect = [min(x) min(y) max(x)-min(x) max(y)-min(y)];
        curvature = get(axchild(i),'Curvature');
        curvature(1) = curvature(1)*rect(3)*0.5;
        curvature(2) = curvature(2)*rect(4)*0.5;
        markerOverlap = 0;
        if ~strcmp(linestyle, 'none')
            markerOverlap = max(markerOverlap, convertunit(linewidth*0.5, 'points', 'pixels', axpos(4)));    
        end
        boundingBoxElement = rect + [-markerOverlap -markerOverlap 2*markerOverlap 2*markerOverlap];
        [filterString, boundingBox] = filter2svg(fid, axchild(i), boundingBoxAxes, boundingBoxElement);
        % put a rectangle into a group with its markers
        if strcmp(get(axchild(i),'Clipping'),'on')
            clippingIdString = clipping2svg(fid, axchild(i), ax, paperpos, axpos, projection, axIdString);
            fprintf(fid,'<g id="%s" clip-path="url(#%s)" %s>\n', createId, clippingIdString, filterString);
        else
            fprintf(fid,'<g id="%s" %s>\n', createId, filterString);
        end
        if ~isempty(filterString)
            % Workaround for Inkscape filter bug
            fprintf(fid,'<rect x="%0.3f" y="%0.3f" width="%0.3f" height="%0.3f" fill="none" stroke="none" />\n', boundingBox(1), boundingBox(2), boundingBox(3), boundingBox(4));
        end
        fprintf(fid,'<rect x="%0.3f" y="%0.3f" width="%0.3f" height="%0.3f" rx="%0.3f" ry="%0.3f" fill="%s" stroke="%s" stroke-width="%0.1fpt" %s />\n',...
            rect(1), rect(2), rect(3), rect(4), curvature(1), curvature(2), fcolorname, scolorname, linewidth, pattern);
        % close the rectangle group
        fprintf(fid,'</g>\n');
    elseif strcmp(get(axchild(i),'Type'),'text')
		if PLOT2SVG_globals.octave
			extent = [0 0 0 0];
		else
			extent = get(axchild(i),'Extent');
		end
        margin = get(axchild(i),'Margin');
        facecolor = get(axchild(i),'BackgroundColor');
        edgecolor = get(axchild(i),'EdgeColor');
        linewidth = get(axchild(i),'LineWidth');
        linestyle = get(axchild(i),'LineStyle');
        if ischar(facecolor)
            if ~strcmp(facecolor,'none')
                error('Illegal face color for text.');    
            else
                facecolorname = 'none';
            end
        else
            facecolorname = searchcolor(id, facecolor);       
        end
        if ischar(edgecolor)
            if ~strcmp(edgecolor,'none')
                error('Illegal edge color for text.');
            else
                edgecolorname = 'none';
            end
        else
            edgecolorname = searchcolor(id, edgecolor);       
        end
        extentx = [extent(1) extent(1)+extent(3)];
        extenty = [extent(2) extent(2)+extent(4)];
        extentz = [0 0];
        [x,y,z] = project(extentx,extenty,extentz,projection);
        x = (x*axpos(3)+axpos(1))*paperpos(3);
        y = (1-(y*axpos(4)+axpos(2)))*paperpos(4);
        box = [min(x)-margin min(y)-margin max(x)-min(x)+2*margin max(y)-min(y)+2*margin];
        markerOverlap = 0;
        if ~strcmp(linestyle, 'none')
            markerOverlap = max(markerOverlap, convertunit(linewidth*0.5, 'points', 'pixels', axpos(4)));    
        end
        boundingBoxElement = [min(x)-markerOverlap min(y)-markerOverlap max(x)-min(x)+2*markerOverlap max(y)-min(y)+2*markerOverlap];
        [filterString, boundingBox] = filter2svg(fid, axchild(i), boundingBoxAxes, boundingBoxElement);
        if strcmp(get(axchild(i),'Clipping'),'on') && ~PLOT2SVG_globals.octave
            clippingIdString = clipping2svg(fid, axchild(i), ax, paperpos, axpos, projection, axIdString);
            fprintf(fid,'<g id="%s" clip-path="url(#%s)" %s>\n', createId, clippingIdString, filterString);
        else
            fprintf(fid,'<g id="%s" %s>\n', createId, filterString);
        end
        if ~isempty(filterString)
            % Workaround for Inkscape filter bug
            fprintf(fid,'<rect x="%0.3f" y="%0.3f" width="%0.3f" height="%0.3f" fill="none" stroke="none" />\n', boundingBox(1), boundingBox(2), boundingBox(3), boundingBox(4));
        end
        if ~strcmp(edgecolorname, 'none') || ~strcmp(facecolorname, 'none')
            pattern = lineStyle2svg(linestyle, linewidth);
            fprintf(fid,'<rect x="%0.3f" y="%0.3f" width="%0.3f" height="%0.3f" fill="%s" stroke="%s" stroke-width="%0.1fpt" %s />\n', ...
                box(1), box(2), box(3), box(4), facecolorname, edgecolorname, linewidth, pattern);
        end
        text2svg(fid,1,axpos,paperpos,axchild(i),ax,projection)
        fprintf(fid,'</g>\n');
    elseif strcmp(get(axchild(i),'Type'),'image')
        cmap=get(id,'Colormap');
        pointx=get(axchild(i),'XData');
        pointy=get(axchild(i),'YData');
        % If the XData is a vector we only use start and stop for the image
        if (size(pointx, 1) > 2) || (size(pointx, 2) > 1)
            pointx = [pointx(1) pointx(end)];   
        end
        if (size(pointy, 1) > 2) || (size(pointy, 2) > 1) 
            pointy = [pointy(1) pointy(end)];   
        end
        if (size(pointx, 1) > 1) && (size(pointy, 1) > 1)
            [x,y,z] = project(pointx,zeros(size(pointx)),zeros(size(pointx)),projection);
        else
            [x,y_dummy,z] = project(pointx,zeros(size(pointx)),zeros(size(pointx)),projection);
            [x_dummy,y,z] = project(zeros(size(pointy)),pointy,zeros(size(pointy)),projection);
        end
        pointc=get(axchild(i),'CData');
        %pointcclass = class(pointc);  % Bugfix proposed by Tom
        if strcmp(get(axchild(i),'CDataMapping'),'scaled')
            clim=get(ax,'CLim');
            pointc=(pointc-clim(1))/(clim(2)-clim(1))*(size(cmap,1) - 1) + 1; % Bugfix proposed by Tom
            %pointcclass = 'double'; % since range is now [0->size(cmap,1)-1]  % Bugfix proposed by Tom
        end
        data_aspect_ratio = get(ax,'DataAspectRatio');
        if length(x) == 2
            if size(pointc, 2) == 1
                halfwidthx = abs(x(2) - x(1)) * data_aspect_ratio(1);
            else
                halfwidthx = abs(x(2) - x(1))/(size(pointc,2) - 1);   
            end
        else
            halfwidthx = data_aspect_ratio(1);
        end
        if length(y) == 2
            if size(pointc, 1) == 1
                halfwidthy = abs(y(2) - y(1)) * data_aspect_ratio(2);
            else
                halfwidthy = abs(y(2) - y(1))/(size(pointc,1) - 1);   
            end
        else
            halfwidthy = data_aspect_ratio(2);
        end
        if length(pointx) > 1
            if xor(strcmp(get(ax,'XDir'),'reverse'), pointx(1) > pointx(2))
                if ndims(pointc) < 3
                    pointc=fliplr(pointc);
                elseif ndims(pointc) == 3
                    for j = 1:size(pointc,3)
                        pointc(:,:,j)=fliplr(pointc(:,:,j));
                    end
                else
                    error('Invalid number of dimensions of data.');
                end
            end
        end
        if length(pointy) > 1
            if xor(strcmp(get(ax,'YDir'),'reverse'), pointy(1) > pointy(2))
                if ndims(pointc) < 3
                    pointc=flipud(pointc);
                elseif ndims(pointc) == 3
                    for j = [1:size(pointc,3)]
                        pointc(:,:,j)=flipud(pointc(:,:,j));
                    end
                else
                    error('Invalid number of dimensions of data.');
                end
            end
        end
        % pointc = cast(pointc,pointcclass);  % Bugfix proposed by Tom
        % Function 'cast' is not supported by old Matlab versions
        if (~isa(pointc, 'double') && ~isa(pointc, 'single'))
            if strcmp(get(axchild(i),'CDataMapping'),'scaled')
                pointc = double(pointc);
            else
                pointc = double(pointc) + 1;
            end
        end
        if ndims(pointc) ~= 3
            pointc = max(min(round(double(pointc)),size(cmap,1)),1);
        end
        CameraUpVector = get(ax,'CameraUpVector');
        filename = [PLOT2SVG_globals.basefilename sprintf('%03d',PLOT2SVG_globals.figurenumber) '.' PLOT2SVG_globals.pixelfiletype];
        PLOT2SVG_globals.figurenumber = PLOT2SVG_globals.figurenumber + 1;
        if isempty(PLOT2SVG_globals.basefilepath)
            current_path = pwd;
        else
            current_path = PLOT2SVG_globals.basefilepath;
        end
        if exist(fullfile(current_path,filename),'file')
            lastwarn('');
            delete(filename);
            if strcmp(lastwarn,'File not found or permission denied.')
                error('Cannot write image file. Make sure that no image is opened in an other program.')    
            end
        end
        if ndims(pointc) < 3
            pointc = flipud(pointc);
        elseif ndims(pointc) == 3
            for j = 1:size(pointc,3)
                pointc(:,:,j)=flipud(pointc(:,:,j));
            end
        else
            error('Invalid number of dimensions of data.');
        end
        if ndims(pointc) == 3
            % pointc is not indexed
            imwrite(pointc,fullfile(PLOT2SVG_globals.basefilepath,filename),PLOT2SVG_globals.pixelfiletype);
        else
            % pointc is probably indexed
            if PLOT2SVG_globals.octave
				pointc = max(2, pointc);
            end
            imwrite(pointc,cmap,fullfile(PLOT2SVG_globals.basefilepath,filename),PLOT2SVG_globals.pixelfiletype);
        end
            lx=(size(pointc,2)*halfwidthx)*axpos(3)*paperpos(3);
        	ly=(size(pointc,1)*halfwidthy)*axpos(4)*paperpos(4);
        if strcmp(get(ax,'DataAspectRatioMode'),'manual')
            pointsx=((min(x) - halfwidthx/2)*axpos(3)+axpos(1))*paperpos(3);
            pointsy=(1-((max(y) + halfwidthy/2)*axpos(4)+axpos(2)))*paperpos(4);
        else
            pointsx=axpos(1)*paperpos(3);
            pointsy=(1-(axpos(4)+axpos(2)))*paperpos(4);
        end
        [filterString, boundingBox] = filter2svg(fid, axchild(i), boundingBoxAxes, boundingBoxAxes);
        if strcmp(get(axchild(i),'Clipping'),'on')
            clippingIdString = clipping2svg(fid, axchild(i), ax, paperpos, axpos, projection, axIdString);
            fprintf(fid,'<g id="%s" clip-path="url(#%s)" %s>\n', createId, clippingIdString, filterString);
            if ~isempty(filterString)
                % Workaround for Inkscape filter bug
                fprintf(fid,'<rect x="%0.3f" y="%0.3f" width="%0.3f" height="%0.3f" fill="none" stroke="none" />\n', boundingBox(1), boundingBox(2), boundingBox(3), boundingBox(4));
            end
            fprintf(fid,'<image x="%0.3f" y="%0.3f" width="%0.3f" height="%0.3f" image-rendering="optimizeSpeed" preserveAspectRatio="none" xlink:href="%s" />\n', pointsx, pointsy, lx, ly, filename);
            fprintf(fid,'</g>\n');
        else
            fprintf(fid,'<g id="%s" %s>\n', createId, filterString);
            if ~isempty(filterString)
                % Workaround for Inkscape filter bug
                fprintf(fid,'<rect x="%0.3f" y="%0.3f" width="%0.3f" height="%0.3f" fill="none" stroke="none" />\n', boundingBox(1), boundingBox(2), boundingBox(3), boundingBox(4));
            end
            fprintf(fid,'<image x="%0.3f" y="%0.3f" width="%0.3f" height="%0.3f" image-rendering="optimizeSpeed" preserveAspectRatio="none" xlink:href="%s" />\n', pointsx, pointsy, lx, ly, filename);
            fprintf(fid,'</g>\n');
        end
    elseif strcmp(get(axchild(i),'Type'), 'hggroup')
        % handle group types (like error bars)
        % FIXME: they are not yet perfectly handled, there are more options
        % that are not used
        [filterString, boundingBox] = filter2svg(fid, axchild(i), boundingBoxAxes, boundingBoxAxes);
        if strcmp(get(axchild(i),'Clipping'),'on')
            clippingIdString = clipping2svg(fid, axchild(i), ax, paperpos, axpos, projection, axIdString);
            fprintf(fid,'<g id="%s" clip-path="url(#%s)" %s>\n', createId, clippingIdString, filterString);
        else
            fprintf(fid, '<g id="%s" %s>', createId, filterString);
        end
        if ~isempty(filterString)
            % Workaround for Inkscape filter bug
            fprintf(fid,'<rect x="%0.3f" y="%0.3f" width="%0.3f" height="%0.3f" fill="none" stroke="none" />\n', boundingBox(1), boundingBox(2), boundingBox(3), boundingBox(4));
        end
        group=axchild2svg(fid,id,axIdString,ax,group,paperpos,get(axchild(i), 'Children'),axpos,groupax,projection,boundingBoxAxes);
        fprintf(fid, '</g>');
    elseif strcmp(get(axchild(i),'Type'), 'hgtransform')
        if strcmpi(get(axchild(i), 'Visible'), 'on')
            [filterString, boundingBox] = filter2svg(fid, axchild(i), boundingBoxAxes, boundingBoxAxes);
            if strcmp(get(axchild(i),'Clipping'),'on')
                clippingIdString = clipping2svg(fid, axchild(i), ax, paperpos, axpos, projection, axIdString);
                fprintf(fid,'<g id="%s" clip-path="url(#%s)" %s>\n', createId, clippingIdString, filterString);
            else
                fprintf(fid, '<g id="%s" %s>', createId, filterString);
            end
            if ~isempty(filterString)
                % Workaround for Inkscape filter bug
                fprintf(fid,'<rect x="%0.3f" y="%0.3f" width="%0.3f" height="%0.3f" fill="none" stroke="none" />\n', boundingBox(1), boundingBox(2), boundingBox(3), boundingBox(4));
            end
            group=axchild2svg(fid,id,axIdString,ax,group,paperpos,get(axchild(i), 'Children'),axpos,groupax,projection,boundingBoxAxes);
            fprintf(fid, '</g>');
        end
    else
        disp(['   Warning: Unhandled child type: ' get(axchild(i),'Type')]);
    end
end

function result = selectColor(axchild, id, p, q, points, pointc, colorname, faces, type)
if size(pointc, 1) == 1
    color = pointc;
elseif size(pointc, 1) == size(faces, 1)
    color = pointc(p,:);
elseif size(pointc, 1) == size(points, 2)
    if strcmp(get(axchild, type), 'flat')
        color = pointc(faces(p, q));
    else
        color = pointc(faces(p, q));     % TO DO: color interpolation
    end
else
    error('Unsupported color handling for patches.');
end
if ~isnan(color)
    if size(color, 2) == 1
        result = ['#' colorname(ceil(color), :)];
    else
        if strcmp(get(axchild, type), 'flat')   % Bugfix 27.01.2008
            result = searchcolor(id, color / 64);
        else
            result = searchcolor(id, color);
        end
    end
else
    result = 'none';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a patch (filled area)
function patch2svg(fid,group,axpos,xtot,ytot,scolorname,style,width, edgecolorname, face_opacity, edge_opacity, closed)
if closed
    type = 'polygon';
else
    type = 'polyline';
end
pattern = lineStyle2svg(style, width);
if strcmp(style, 'none')
    edge_opacity = 0.0;
end
for i = 1:size(xtot, 1)
    x = xtot(i, :);
    y = ytot(i, :);
    if all(~isnan(x)) && all(~isnan(y))
        for j = 1:20000:length(x)
            xx = x(j:min(length(x), j + 19999));
            yy = y(j:min(length(y), j + 19999));
            if ~strcmp(edgecolorname,'none') || ~strcmp(scolorname,'none')
                fprintf(fid,'      <%s fill="%s" fill-opacity="%0.2f" stroke="%s" stroke-width="%0.1fpt" stroke-opacity="%0.2f" %s points="',...
                    type, scolorname, face_opacity, edgecolorname, width, edge_opacity, pattern);
                fprintf(fid,'%0.3f,%0.3f ',[xx;yy]);
                fprintf(fid,'"/>\n');        
            end
        end
    else
        parts = find(isnan(x) + isnan(y));
        if ~isempty(parts) && (parts(1) ~= 1)
            parts = [0 parts];
        end
        if parts(length(parts)) ~= length(x)
            parts = [parts length(x) + 1];
        end
        for j = 1:(length(parts) - 1)
            xx = x((parts(j)+1):(parts(j+1)-1));
            yy = y((parts(j)+1):(parts(j+1)-1));
            if ~strcmp(edgecolorname,'none') || ~strcmp(scolorname,'none')
                if ~isempty(xx)
                    fprintf(fid,'      <%s fill="%s" fill-opacity="%0.2f" stroke="%s" stroke-width="%0.1fpt" stroke-opacity="%0.2f" %s points="',...
                        type, scolorname, face_opacity, edgecolorname, width, edge_opacity, pattern);
                    fprintf(fid,'%0.3f,%0.3f ',[xx;yy]);
                    fprintf(fid,'"/>\n');        
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a patch (filled area)
function gouraud_patch2svg(fid,group,axpos,xtot,ytot,cdata,style,width, edgecolorname, face_opacity, edge_opacity,id)
global colorname
pattern = lineStyle2svg(style, width);
if strcmp(style, 'none')
    edge_opacity = 0.0;
end
for i=1:size(xtot,1)
    x = xtot(i,:);
    y = ytot(i,:);
    if (any(isnan(x)) || any(isnan(y)))
        fprintf('Warning: Found NaN in Gouraud patch.\n')
    else
        % If there are more than 2 edges always 3 edges are taken togehter
        % to form a triangle
        if length(x) > 2
            for j = 3:length(x)
                coord = [x([1 j-1 j]);y([1 j-1 j])];
                face_color = cdata(1,:);
                face_color2 = cdata(j-1,:);
                face_color3 = cdata(j,:);
                delta = coord(:,3)-coord(:,2);
                if det([delta (coord(:,1)-coord(:,2))]) ~= 0
                    if ~isnan(face_color)
                        IDstring1 = createId;
                        IDstring2 = createId;
                        if size(face_color2,2)==1
                            face_color_name2 = ['#' colorname(ceil(face_color2),:)];
                        else
                            face_color_name2 = searchcolor(id,face_color2/64);    
                        end
                        if size(face_color3,2)==1
                            face_color_name3 = ['#' colorname(ceil(face_color3),:)];
                        else
                            face_color_name3 = searchcolor(id,face_color3/64);    
                        end
                        grad_end=(delta)*(delta'*(coord(:,1)-coord(:,2)))/(delta'*delta) + coord(:,2);
                        if size(face_color,2)==1
                            face_color_name = ['#' colorname(ceil(face_color),:)];
                        else
                            face_color_name = searchcolor(id,face_color/64);    
                        end
                        fprintf(fid,'<defs>\n');
                        fprintf(fid,'<linearGradient id="%s" gradientUnits="userSpaceOnUse" x1="%0.3f" y1="%0.3f" x2="%0.3f" y2="%0.3f">\n',...
                            IDstring1, coord(1,2), coord(2,2), coord(1,3), coord(2,3));
                        fprintf(fid,'<stop offset="0" stop-color="%s" stop-opacity="1"/>\n',face_color_name2);
                        fprintf(fid,'<stop offset="1" stop-color="%s" stop-opacity="1"/>\n',face_color_name3);
                        fprintf(fid,'</linearGradient>\n');
                        fprintf(fid,'<linearGradient id="%s" gradientUnits="userSpaceOnUse" x1="%0.3f" y1="%0.3f" x2="%0.3f" y2="%0.3f">\n',...
                            IDstring2, coord(1,1), coord(2,1), grad_end(1), grad_end(2));
                        fprintf(fid,'<stop offset="0" stop-color="%s" stop-opacity="1"/>\n',face_color_name);
                        fprintf(fid,'<stop offset="1" stop-color="%s" stop-opacity="0"/>\n',face_color_name);
                        fprintf(fid,'</linearGradient>\n');
                        fprintf(fid,'</defs>\n');    
                        % Open group
                        temp_string = sprintf('%0.3f,%0.3f ',coord);
                        fprintf(fid,'<g opacity="%0.2f">\n',face_opacity);
                        fprintf(fid,'<polyline fill="url(#%s)" stroke="none" points="%s"/>\n',IDstring1,temp_string);
                        fprintf(fid,'<polyline fill="url(#%s)" stroke="none" points="%s"/>\n',IDstring2,temp_string);
                        % Close group
                        fprintf(fid,'</g>\n');
                    end
                end
            end
        end
        % Last we draw the line around the patch
        if ~strcmp(edgecolorname,'none')
            fprintf(fid,'<polygon fill="none" stroke="%s" stroke-width="%0.1fpt" stroke-opacity="%0.2f" %s points="',...
                edgecolorname, width, edge_opacity, pattern);
            fprintf(fid,'%0.3f,%0.3f ',[x;y]);
            fprintf(fid,'"/>\n');        
        end
    end
end

function patternString = lineStyle2svg(lineStyle, lineWidth)
% Create the string for the line style
% Note: On Matlab the line style is not identical on screen and in a png.
%       We will try to match with the screen format.
scaling = max(1, lineWidth * 0.4);
switch lineStyle
    case '--', patternString = sprintf('stroke-dasharray="%0.1f,%0.1f"', 8*scaling, 2*scaling);
    case ':', patternString = sprintf('stroke-dasharray="%0.1f,%0.1f"', 2*scaling, 2*scaling);
    case '-.', patternString = sprintf('stroke-dasharray="%0.1f,%0.1f,%0.1f,%0.1f"', 8*scaling, 2*scaling, 2*scaling, 2*scaling);
    otherwise, patternString = 'stroke-dasharray="none"';   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a line segment
% this algorthm was optimized for large segement counts
function line2svg(fid, group, axpos, x, y, scolorname, style, width)
if ~strcmp(style,'none')
    pattern = lineStyle2svg(style, width);
    if (isnan(x) == zeros(size(x)) & isnan(y) == zeros(size(y)))
        for j = 1:20000:length(x)
            xx = x(j:min(length(x), j + 19999));
            yy = y(j:min(length(y), j + 19999));
            fprintf(fid,'      <polyline fill="none" stroke="%s" stroke-width="%0.1fpt" %s points="', scolorname, width, pattern);
            fprintf(fid,'%0.3f,%0.3f ',[xx;yy]);
            fprintf(fid,'"/>\n');
        end
    else
        parts = find(isnan(x) + isnan(y));
        if ~isempty(parts) && (parts(1) ~= 1)
            parts=[0 parts];
        end
        if parts(length(parts)) ~= length(x)
            parts = [parts length(x) + 1];
        end
        for j = 1:(length(parts) - 1)
            xx = x((parts(j) + 1):(parts(j + 1) - 1));
            yy = y((parts(j) + 1):(parts(j + 1) - 1));
            if ~isempty(xx)
                fprintf(fid,'      <polyline fill="none" stroke="%s" stroke-width="%0.1fpt" %s points="', scolorname, width, pattern);
                fprintf(fid,'%0.3f,%0.3f ', [xx;yy]);
                fprintf(fid,'"/>\n');
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a circle
function circle2svg(fid,group,axpos,x,y,radius,markeredgecolorname,markerfacecolorname,width)
for j = 1:length(x)
    if ~(isnan(x(j)) || isnan(y(j)))
        if ~strcmp(markeredgecolorname,'none') || ~strcmp(markerfacecolorname,'none')
            fprintf(fid,'<circle cx="%0.3f" cy="%0.3f" r="%0.3f" fill="%s" stroke="%s" stroke-width="%0.1fpt" />\n',x(j),y(j),radius,markerfacecolorname,markeredgecolorname,width);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function control2svg(fid,id,ax,group,paperpos)
global PLOT2SVG_globals
set(ax,'Units','pixels');
pos=get(ax,'Position');
pict=getframe(id,pos);
if isempty(pict.colormap)
    pict.colormap=colormap;
end
filename = [PLOT2SVG_globals.basefilename sprintf('%03d',PLOT2SVG_globals.figurenumber) '.' PLOT2SVG_globals.pixelfiletype];
PLOT2SVG_globals.figurenumber = PLOT2SVG_globals.figurenumber + 1;
if isempty(PLOT2SVG_globals.basefilepath)
    current_path = pwd;
else
    current_path = PLOT2SVG_globals.basefilepath;
end
if exist(fullfile(current_path,filename),'file')
    lastwarn('');
    delete(filename);
    if strcmp(lastwarn,'File not found or permission denied.')
        error('Cannot write image file. Make sure that no image is opened in an other program.')    
    end
end
imwrite(pict.cdata,fullfile(PLOT2SVG_globals.basefilepath,filename),PLOT2SVG_globals.pixelfiletype);
set(ax,'Units','normalized');
posNorm=get(ax,'Position');
posInches(1)=posNorm(1)*paperpos(3);
posInches(2)=posNorm(2)*paperpos(4);
posInches(3)=posNorm(3)*paperpos(3);
posInches(4)=posNorm(4)*paperpos(4);
lx = posInches(3);
ly = posInches(4);
pointsx = posInches(1);
pointsy = paperpos(4)-posInches(2)-posInches(4);
fprintf(fid,'<image x="%0.3f" y="%0.3f" width="%0.3f" height="%0.3f" image-rendering="optimizeSpeed" xlink:href="%s" />\n', pointsx, pointsy, lx, ly, filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a text in the axis frame
% the position of the text has to be adapted to the axis scaling
function text2svg(fid,group,axpos,paperpos,id,ax,projection)
global PLOT2SVG_globals;
originalTextUnits=get(id,'Units');
if PLOT2SVG_globals.octave
	set(id,'Units','data');
else
	set(id,'Units','Data');
end
textpos=get(id,'Position');
if PLOT2SVG_globals.octave
    xlim = get(ax, 'XLim');
    ylim = get(ax, 'YLim');
    zlim = get(ax, 'ZLim');
    if get(ax, 'XLabel') == id
        textpos = textpos + [mean(xlim) ylim(1)-diff(ylim)./axpos(4)*0.06 0];
    elseif get(ax, 'YLabel') == id
        textpos = textpos + [xlim(1)-diff(xlim)./axpos(3)*0.06 mean(ylim) 0];
    elseif get(ax, 'ZLabel') == id
        if projection.xyplane
            return;
        end
        textpos = textpos + [xlim(1)-diff(xlim)./axpos(3)*0.06 0 mean(zlim)];
    elseif get(ax, 'Title') == id
        textpos = textpos + [mean(xlim) ylim(2)+diff(ylim)./axpos(4)*0.01 0];
    end
end
textfontsize = get(id,'FontSize');
fontsize = convertunit(get(id,'FontSize'),get(id,'FontUnits'),'points', axpos(4));   % convert fontsize to inches
paperposOriginal = get(gcf,'Position');
font_color = searchcolor(id,get(id,'Color'));
if strcmp(get(ax,'XScale'),'log')
    textpos(1) = log10(textpos(1));
end
if strcmp(get(ax,'YScale'),'log')
    textpos(2) = log10(textpos(2));
end
if strcmp(get(ax,'ZScale'),'log')
    textpos(3) = log10(textpos(3));
end
[x,y,z] = project(textpos(1), textpos(2), textpos(3), projection);
x = (x * axpos(3) + axpos(1)) * paperpos(3);
y = (1 - (y * axpos(4) + axpos(2))) * paperpos(4);
textvalign = get(id,'VerticalAlignment');
textalign = get(id,'HorizontalAlignment');
texttext = get(id,'String');
textrot = get(id,'Rotation');
dx = sin(textrot * pi / 180) * (fontsize * 1.25 * 1.2);
dy = cos(textrot * pi / 180) * (fontsize * 1.25 * 1.2);
lines = max(size(get(id,'String'),1),1);
if size(texttext,2)~=0
    j = 1;
    for i = 0:1:(lines - 1)
        if iscell(texttext)
            label2svg(fid, group, axpos, id, x + i * dx, y + i * dy, convertString(texttext{j}), textalign, textrot, textvalign, lines, paperpos, font_color, 0)
        else
            label2svg(fid, group, axpos, id, x + i * dx, y + i * dy, convertString(texttext(j,:)), textalign, textrot, textvalign, lines, paperpos, font_color, 0)
        end
        j = j + 1;   
    end
else
    label2svg(fid,group,axpos,id,x,y,'',textalign,textrot,textvalign,lines,paperpos,font_color,0)
end
set(id,'Units',originalTextUnits);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adds the exponents to the axis thickmarks if needed
% MATLAB itself offers no information about this exponent scaling
% the exponent have therefore to be extracted from the thickmarks
function exponent2svg(fid,group,axpos,paperpos,ax,axxtick,axytick,axztick)
global PLOT2SVG_globals
if strcmp(get(ax,'XTickLabelMode'),'auto') && strcmp(get(ax,'XScale'),'linear')
    fontsize=convertunit(get(ax,'FontSize'),get(ax,'FontUnits'),'points', axpos(4));   % convert fontsize to inches
    font_color=searchcolor(ax,get(ax,'XColor'));
    if PLOT2SVG_globals.octave
        % Octave stores XTickLabel in a cell array, which does not work nicely with str2num. --Jakob Malm
        axlabelx = get(ax, 'XTickLabel');
        numlabels = zeros(length(axlabelx), 1);
        for ix = 1:length (axlabelx)
            numlabels(ix) = str2num(axlabelx{ix});
        end
    else
        numlabels = str2num(get(ax,'XTickLabel'));
    end
    labelpos = axxtick;%get(ax,'XTick');
    numlabels = numlabels(:);
    labelpos = labelpos(:);
    indexnz = find(labelpos ~= 0);
    if (~isempty(indexnz) && ~isempty(numlabels))
        ratio=numlabels(indexnz)./labelpos(indexnz);
        if round(log10(ratio(1))) ~= 0 && ratio(1) ~= 0
            exptext = sprintf('&#215; 10<tspan font-size="%0.1fpt" dy="%0.1fpt">%g</tspan>', 0.7*fontsize, -0.7*fontsize,-log10(ratio(1)));
            label2svg(fid,group,axpos,ax,(axpos(1)+axpos(3))*paperpos(3),(1-axpos(2))*paperpos(4)+3*fontsize,exptext,'right',0,'top',1,paperpos,font_color,0)           
        end
    end
end
if strcmp(get(ax,'YTickLabelMode'),'auto') && strcmp(get(ax,'YScale'),'linear')
    fontsize=convertunit(get(ax,'FontSize'),get(ax,'FontUnits'),'points', axpos(4));
    font_color=searchcolor(ax,get(ax,'YColor'));
    if PLOT2SVG_globals.octave
        % Octave stores YTickLabel in a cell array, which does not work nicely with str2num. --Jakob Malm
        axlabely = get(ax, 'YTickLabel');
        numlabels = zeros(length(axlabely), 1);
        for ix = 1:length(axlabely)
            numlabels(ix) = str2num(axlabely{ix});
        end        
    else
        numlabels = str2num(get(ax,'YTickLabel'));
    end
    labelpos = axytick;%get(ax,'YTick');
    numlabels = numlabels(:);
    labelpos = labelpos(:);
    indexnz = find(labelpos ~= 0);
    if (~isempty(indexnz) && ~isempty(numlabels))
        ratio = numlabels(indexnz)./labelpos(indexnz);
        if round(log10(ratio(1))) ~= 0 && ratio(1) ~= 0
            exptext = sprintf('&#215; 10<tspan font-size="%0.1fpt" dy="%0.1fpt">%g</tspan>',0.7*fontsize,-0.7*fontsize,-log10(ratio(1)));
            label2svg(fid,group,axpos,ax,axpos(1)*paperpos(3),(1-(axpos(2)+axpos(4)))*paperpos(4)-0.5*fontsize,exptext,'left',0,'bottom',1,paperpos,font_color,0)           
        end
    end
end
if strcmp(get(ax,'ZTickLabelMode'),'auto') && strcmp(get(ax,'ZScale'),'linear')
    fontsize=convertunit(get(ax,'FontSize'),get(ax,'FontUnits'),'points', axpos(4));
    font_color=searchcolor(ax,get(ax,'ZColor'));
    if PLOT2SVG_globals.octave
        % Octave stores ZTickLabel in a cell array, which does not work nicely with str2num. --Jakob Malm
        axlabelz = get (ax, 'ZTickLabel');
        numlabels = zeros(length(axlabelz), 1);
        for ix = 1:length(axlabelz)
            numlabels(ix) = str2num(axlabelz{ix});
        end
    else
        numlabels = str2num(get(ax,'ZTickLabel'));
    end
    labelpos = axztick;%get(ax,'ZTick');
    numlabels = numlabels(:);
    labelpos = labelpos(:);
    indexnz = find(labelpos ~= 0);
    if (~isempty(indexnz) && ~isempty(numlabels))
        ratio = numlabels(indexnz)./labelpos(indexnz);
        if round(log10(ratio(1))) ~= 0 && ratio(1) ~= 0
            exptext = sprintf('&#215; 10<tspan font-size="%0.1fpt" dy="%0.1fpt">%g</tspan>',0.7*fontsize,-0.7*fontsize,-log10(ratio(1)));
            label2svg(fid,group,axpos,ax,axpos(1)*paperpos(3),(1-(axpos(2)+axpos(4)))*paperpos(4)-0.5*fontsize,exptext,'left',0,'top',1,paperpos,font_color,0)           
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a label in the figure
% former versions of FrameMaker supported the commands FDY and FDX to shift the text
% this commands were replaced by a shift parameter that is normed by the font size
function label2svg(fid,group,axpos,id,x,y,tex,align,angle,valign,lines,paperpos,font_color,exponent)
if isempty(tex)
    return;
end
textfontname=get(id,'FontName');
if strcmp(textfontname, '*')
    textfontname = 'Arial';
end
set(id,'FontUnits','points');
textfontsize=get(id,'FontSize');
if isfield(get(id),'Interpreter')
    if strcmp(get(id,'Interpreter'),'tex')
        latex=1;
    elseif strcmp(get(id,'Interpreter'),'latex')
        latex=1;
    else
        latex=0;
    end
else
    latex=1;
end
fontsize=convertunit(get(id,'FontSize'),get(id,'FontUnits'),'points', axpos(4));   % convert fontsize to inches
paperposOriginal=get(gcf,'Position');
fontsize=fontsize*paperpos(4)/paperposOriginal(4);
textfontsize=textfontsize*paperpos(4)/paperposOriginal(4);
fontweight = get(id,'FontWeight');
switch lower(fontweight)
    case 'bold', fweight = ' font-weight="bold"';
    case 'light', fweight = ' font-weight="lighter"';
    case 'demi', fweight = ' font-weight="lighter"';
    otherwise, fweight = '';   % default 'normal'
end
fontangle = get(id,'FontAngle');
switch lower(fontangle)
    case 'italic', fangle = ' font-style="italic"';
    case 'oblique', fangle = ' font-style="oblique"';
otherwise, fangle = '';  % default 'normal'
end
% Note: The attribute 'alignment-baseline' cannot be used as it is often
% badly supported. Therfore, we use shifts.
switch lower(valign)
     case 'top',shift=fontsize*1.18;
     case 'cap',shift=fontsize*0.95;
     case 'middle',shift = -((lines-1)/2*fontsize*1.25*1.2) + fontsize * 0.45;
     case 'bottom',shift = -((lines-1)*fontsize*1.25*1.2) + fontsize * -0.25;
     otherwise,shift=0;
end
switch lower(align)
    case 'right', anchor = 'end'; 
    case 'center',anchor = 'middle';
    otherwise,anchor = 'start';
end
if iscellstr(tex)
    tex = strvcat(tex);
elseif ~ ischar(tex)
    error('Invalid character type');
end    
if latex==1
    tex=strrep(tex,'$','');
    %if ~exponent
    %  try
    %        tex = texlabel(tex);
    %    catch
    %        fprintf('  Warning: Error during conversion to a latex string.');
    %    end
    %end
    tex=strrep(tex,'\alpha','{&#945;}');
    tex=strrep(tex,'\beta','{&#946;}');
    tex=strrep(tex,'\gamma','{&#947;}');
    tex=strrep(tex,'\delta','{&#948;}');
    tex=strrep(tex,'\epsilon','{&#949;}');
    tex=strrep(tex,'\zeta','{&#950;}');
    tex=strrep(tex,'\eta','{&#951;}');
    tex=strrep(tex,'\theta','{&#952;}');
    tex=strrep(tex,'\vartheta','{&#977;}');
    tex=strrep(tex,'\iota','{&#953;}');
    tex=strrep(tex,'\kappa','{&#954;}');
    tex=strrep(tex,'\lambda','{&#955;}');
    tex=strrep(tex,'\mu','{&#181;}');
    tex=strrep(tex,'\nu','{&#957;}');
    tex=strrep(tex,'\xi','{&#958;}');
    tex=strrep(tex,'\pi','{&#960;}');
    tex=strrep(tex,'\rho','{&#961;}');
    tex=strrep(tex,'\sigma','{&#963;}');
    tex=strrep(tex,'\varsigma','{&#962;}');
    tex=strrep(tex,'\tau','{&#964;}');
    tex=strrep(tex,'\upsilon','{&#965;}');
    tex=strrep(tex,'\phi','{&#966;}');
    tex=strrep(tex,'\chi','{&#967;}');
    tex=strrep(tex,'\psi','{&#968;}');
    tex=strrep(tex,'\omega','{&#969;}');
    tex=strrep(tex,'\Gamma','{&#915;}');
    tex=strrep(tex,'\Delta','{&#916;}');
    tex=strrep(tex,'\Theta','{&#920;}');
    tex=strrep(tex,'\Lambda','{&#923;}');
    tex=strrep(tex,'\Xi','{&#926;}');
    tex=strrep(tex,'\Pi','{&#928;}');
    tex=strrep(tex,'\Sigma','{&#931;}');
    tex=strrep(tex,'\Tau','{&#932;}');
    tex=strrep(tex,'\Upsilon','{&#933;}');
    tex=strrep(tex,'\Phi','{&#934;}');
    tex=strrep(tex,'\Psi','{&#936;}');
    tex=strrep(tex,'\Omega','{&#937;}');
    tex=strrep(tex,'\infty','{&#8734;}');
    tex=strrep(tex,'\pm','{&#177;}');
    tex=strrep(tex,'\Im','{&#8465;}');
    tex=strrep(tex,'\Re','{&#8476;}');
    tex=strrep(tex,'\approx','{&#8773;}');
    tex=strrep(tex,'\leq','{&#8804;}');
    tex=strrep(tex,'\geq','{&#8805;}');
    tex=strrep(tex,'\times','{&#215;}');
    tex=strrep(tex,'\leftrightarrow','{&#8596;}');
    tex=strrep(tex,'\leftarrow','{&#8592;}');
    tex=strrep(tex,'\uparrow','{&#8593;}');
    tex=strrep(tex,'\rightarrow','{&#8594;}');
    tex=strrep(tex,'\downarrow','{&#8595;}');
    tex=strrep(tex,'\circ','{&#186;}');
    tex=strrep(tex,'\propto','{&#8733;}');
    tex=strrep(tex,'\partial','{&#8706;}');
    tex=strrep(tex,'\bullet','{&#8226;}');
    tex=strrep(tex,'\div','{&#247;}');

	tex=strrep(tex,'\sum','{&#8721;}');	
	tex=strrep(tex,'\ast','{&#8727;}');	
	tex=strrep(tex,'\sqrt','{&#8730;}');	
    tex=strrep(tex,'\angle','{&#8736;}');	
	tex=strrep(tex,'\wedge','{&#8743;}');	
	tex=strrep(tex,'\land','{&#8743;}');	
	tex=strrep(tex,'\vee','{&#8744;}');	
	tex=strrep(tex,'\lor','{&#8744;}');	
	tex=strrep(tex,'\cap','{&#8745;}');	
	tex=strrep(tex,'\cup','{&#8746;}');	
	tex=strrep(tex,'\int','{&#8747;}');	
%&there4;&#8756;	
	tex=strrep(tex,'\sim','{&#8764;}');	

	tex=strrep(tex,'\forall','{&#8704;}');	
	tex=strrep(tex,'\partial','{&#8706;}');	
	tex=strrep(tex,'\exists','{&#8707;}');	
	tex=strrep(tex,'\emptyset','{&#8709;}');	
	tex=strrep(tex,'\nabla','{&#8711;}');	
	tex=strrep(tex,'\in','{&#8712;}');	
	tex=strrep(tex,'\notin','{&#8713;}');	
	tex=strrep(tex,'\ni','{&#8715;}');	
	tex=strrep(tex,'\prod','{&#8719;}');	

	tex=strrep(tex,'\cong','{&#8773;}');	
	tex=strrep(tex,'\approx','{&#8776;}');	
	tex=strrep(tex,'\neq','{&#8800;}');	
	tex=strrep(tex,'\equiv','{&#8801;}');	
	tex=strrep(tex,'\leq','{&#8804;}');	
	tex=strrep(tex,'\geq','{&#8805;}');	
	tex=strrep(tex,'\subset','{&#8834;}');	
	tex=strrep(tex,'\supset','{&#8835;}');	
%&nsub;&#8836;	
	tex=strrep(tex,'\subseteq','{&#8838;}');	
	tex=strrep(tex,'\supseteq','{&#8839;}');	
    tex=strrep(tex,'\oplus','{&#8853;}');	
    tex=strrep(tex,'\otimes','{&#8855;}');	
    tex=strrep(tex,'\bot','{&#8869;}');	
    tex=strrep(tex,'\cdot','{&#8901;}');	
	tex=strrep(tex,'\bullet','{&#8226;}');	
	tex=strrep(tex,'\ldots','{&#8230;}');	
	tex=strrep(tex,'\prime','{&#8242;}');	
%	&#8243; double prime	
%	&#8254; oline	
    
    
    tex=strrep(tex,'\\','{&#92;}');
    tex=strrep(tex,'\{','{&#123;}');
    tex=strrep(tex,'\}','{&#125;}');
    tex=strrep(tex,'\_','{&#95;}');
    tex=strrep(tex,'\^','{&#94;}');
    
    %fprintf('%s\n', tex);
    tex=latex2svg(tex, textfontname, textfontsize, 0);
end
if isempty(tex)
    return;
end
if exponent
    tex=sprintf('10<tspan font-size="%0.1fpt" dy="%0.1fpt">%s</tspan>', 0.7*textfontsize, -0.7*textfontsize, tex);
    shift = shift + 0.4*fontsize;   % Small correction to make it look nicer
end
% Note: Obviously, Matlab is using font sizes that are rounded to decimal
% pt sizes. This may cause problems for very small figures. But we have to
% follow due to the background size.
%fprintf('%s\n', tex);
fprintf(fid,'  <g transform="translate(%0.3f,%0.3f)">\n', x - shift * sin(-angle/180*pi), y + shift * cos(-angle/180*pi));
fprintf(fid,'    <g transform="rotate(%0.1f)">\n',-angle);
fprintf(fid,'      <text x="%0.3f" y="%0.3f" font-family="%s" text-anchor="%s" font-size="%0.0fpt"%s%s fill="%s" >', 0, 0, textfontname, anchor, textfontsize, fweight, fangle, font_color);
fprintf(fid,'%s',tex);
fprintf(fid,'</text>\n'); 
fprintf(fid,'    </g>\n');
fprintf(fid,'  </g>\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% converts LATEX strings into SVG strings
function returnvalue = latex2svg(StringText, font, size, style)
returnvalue = StringText;
try
if ~isempty(StringText)
    bracket = find(StringText == '{' | StringText == '}');
    bracketCounter = zeros(1,length(StringText));
    bracketCounter(StringText == '{') = 1;
    bracketCounter(StringText == '}') = -1;
    bracketCounter = cumsum(bracketCounter);
    if bracketCounter(end) ~= 0
        fprintf(['Warning: Number of open and closed braces is not equal. Latex string ''' StringText ''' will not be converted.\n']);
    elseif any(bracketCounter < 0)
        fprintf(['Warning: Found a closed brace without a previously opened brace. Latex string ''' StringText ''' will not be converted.\n']);
    else
        if isempty(bracket)
            if ~isempty(find(StringText == '^' | StringText == '_' | StringText == '\' ))
                returnvalue = ['<tspan>' singleLatex2svg(StringText, size) '</tspan>'];    
                % Clean up empty tspan elements
                % More could be done here, but with huge effort to make it
                % match all special cases.
                returnvalue = strrep(returnvalue, '></tspan><tspan>', '>');
            else
                returnvalue = StringText;    
            end
        else
            returnvalue = '<tspan>';
            lastValidCharacter = 1;
            localSize = size;
            restoreOffset = [];
            restoreSize = [];
            for i = 1:length(bracket)
                lastValidCharacterOffset = 1;
                if StringText(bracket(i)) == '{'
                    % Found '{'
                    removeCharacters = 1;
                    localOffset = 0;
                    restoreOffset(bracketCounter(bracket(i))) = 0;
                    restoreSize(bracketCounter(bracket(i))) = localSize;
                    if (bracket(i) > 1)
                        if StringText(bracket(i) - 1) == '_'
                            localOffset = -0.4 * localSize;
                            localSize = 0.7 * localSize;
                            restoreOffset(bracketCounter(bracket(i))) = -localOffset;
                            removeCharacters = 2;
                        elseif StringText(bracket(i) - 1) == '^'
                            localOffset = 0.5 * localSize;
                            localSize = 0.7 * localSize;
                            restoreOffset(bracketCounter(bracket(i))) = -localOffset;
                            removeCharacters = 2;
                        end
                    end
                    returnvalue = [returnvalue singleLatex2svg(StringText(lastValidCharacter:bracket(i) - removeCharacters), restoreSize(bracketCounter(bracket(i)))) '</tspan><tspan'];
                    if localOffset ~= 0
                        returnvalue = [returnvalue ' dy="' num2str(-localOffset, '%0.1f') 'pt"'];
                    end
                    if localSize ~= restoreSize(bracketCounter(bracket(i)))
                        returnvalue = [returnvalue ' font-size="' num2str(localSize, '%0.0f') 'pt"'];                        
                    end
                    returnvalue = [returnvalue '>'];
                else
                    % Found '}'
                    returnvalue = [returnvalue singleLatex2svg(StringText(lastValidCharacter:bracket(i) - 1), localSize) '</tspan><tspan'];
                    if restoreOffset(bracketCounter(bracket(i) - 1)) ~= 0
                        returnvalue = [returnvalue ' dy="' num2str(-restoreOffset(bracketCounter(bracket(i) - 1)), '%0.1f') 'pt"'];
                    end
                    if restoreSize(bracketCounter(bracket(i) - 1)) ~= localSize
                        localSize = restoreSize(bracketCounter(bracket(i) - 1));
                        returnvalue = [returnvalue ' font-size="' num2str(localSize, '%0.0f') 'pt"'];
                    end
                    returnvalue = [returnvalue '>'];
                end
                lastValidCharacter = bracket(i) + lastValidCharacterOffset;
            end
            if lastValidCharacter <= length(StringText)
                returnvalue = [returnvalue singleLatex2svg(StringText(lastValidCharacter:end), localSize)];
            end
            returnvalue = [returnvalue '</tspan>'];
            % Clean up empty tspan elements
            % More could be done here, but with huge effort to make it
            % match all special cases.
            returnvalue = strrep(returnvalue, '></tspan><tspan>', '>');
        end
    end
end
catch
    fprintf(['Warning: Error ''' lasterr ''' occurred during conversion. Latex string ''' StringText ''' will not be converted.\n']);
end

function StringText = singleLatex2svg(StringText, size)
index = find(StringText == '_' | StringText == '^');
if ~isempty(index)
    if index(end) == length(StringText)
        % Remove orphan '_' or '^'
        index = index(1:end-1);
    end
    for i = length(index):-1:1
        if StringText(index(i)) == '_'
            localOffset = 0.4 * size;
            localSize = 0.7 * size;
            StringText = [StringText(1:index(i)-1) ...
                    '</tspan><tspan dy="' num2str(localOffset, '%0.1f') 'pt" font-size="' num2str(localSize, '%0.0f') 'pt">' ...
                    StringText(index(i)+1) ...
                    '</tspan><tspan dy="' num2str(-localOffset, '%0.1f') 'pt" font-size="' num2str(size, '%0.0f') 'pt">' ...
                    StringText(index(i)+2:end)];
        else
            localOffset = 0.5 * size;
            localSize = 0.7 * size;
            StringText = [StringText(1:index(i)-1) ...
                    '</tspan><tspan dy="' num2str(-localOffset, '%0.1f') 'pt" font-size="' num2str(localSize, '%0.0f') 'pt">' ...
                    StringText(index(i)+1) ...
                    '</tspan><tspan dy="' num2str(localOffset, '%0.1f') 'pt" font-size="' num2str(size, '%0.0f') 'pt">' ...
                    StringText(index(i)+2:end)];
        end
    end
end
if ~isempty(strfind(StringText, '\bf'))
    StringText = strrep(StringText, '\bf', '</tspan><tspan font-weight="bold">');
end
if ~isempty(strfind(StringText, '\it'))
    StringText = strrep(StringText, '\it', '</tspan><tspan font-style="italic">');
end
if ~isempty(strfind(StringText, '\sl'))
    StringText = strrep(StringText, '\sl', '</tspan><tspan font-style="oblique">');
end
if ~isempty(strfind(StringText, '\rm'))
    StringText = strrep(StringText, '\rm', '</tspan><tspan font-style="normal" font-weight="normal">');
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function name=searchcolor(id,value)
if ischar(value)
    name = value;
else
    name=sprintf('#%02x%02x%02x',fix(value(1)*255),fix(value(2)*255),fix(value(3)*255));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rvalue = convertunit(value, from, to, parentheight)
% From SVG 1.1. Specification:
% "1pt" equals "1.25px" (and therefore 1.25 user units)
% "1pc" equals "15px" (and therefore 15 user units)
% "1mm" would be "3.543307px" (3.543307 user units)
% "1cm" equals "35.43307px" (and therefore 35.43307 user units)
% "1in" equals "90px" (and therefore 90 user units)
if nargin < 4
    parentheight = 1.25;    % Default
end
switch lower(from)  % convert from input unit to points
    case 'pixels', rvalue = value * 0.8;
    case 'points', rvalue = value;
    case 'centimeters', rvalue = value / 2.54*72;
    case 'inches', rvalue = value * 72; % 72 points = 1 inch
    case 'normalized', rvalue = value * (parentheight * 0.8);
    otherwise, error(['Unknown unit ' from '.']);
end
switch lower(to)    % convert from points to specified unit
    case 'pixels', rvalue = rvalue * 1.25;
    case 'points'; % do nothing
    case 'centimeters', rvalue = rvalue * 2.54 / 72;
    case 'inches', rvalue = rvalue / 72;    % 72 points = 1 inch
    case 'normalized', rvalue = value / (parentheight * 0.8);
    otherwise, error(['Unknown unit ' to '.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function strString=addBackSlash( strSlash)
% adds a backslash at the last position of the string (if not already there)
if ( strSlash(end) ~= filesep)
    strString = [ strSlash filesep];
else
    strString = strSlash;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function strExt=getFileExtension( strFileName)
% returns the file extension of a filename
[path, name, strExt] = fileparts( strFileName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function StringText=convertString(StringText)
if iscell(StringText) % Octave stores some strings in cell arrays. --Jakob Malm
    StringText = StringText{1};
end
if ~isempty(StringText)
    StringText=strrep(StringText,'&','&amp;');  % Do not change sequence !!
    StringText=strrep(StringText,'\\','\');
    StringText=strrep(StringText,'<','&lt;');
    StringText=strrep(StringText,'>','&gt;');
    StringText=strrep(StringText,'"','&quot;');
    % Workaround for Firefox and Inkscape
    StringText=strrep(StringText,'°','&#176;');
    %StringText=strrep(StringText,'°','&deg;');
    StringText=strrep(StringText,'±','&plusmn;');
    StringText=strrep(StringText,'µ','&micro;');
    StringText=strrep(StringText,'²','&sup2;');
    StringText=strrep(StringText,'³','&sup3;');
    StringText=strrep(StringText,'¼','&frac14;');
    StringText=strrep(StringText,'½','&frac12;');
    StringText=strrep(StringText,'¾','&frac34;');
    StringText=strrep(StringText,'©','&copy;');
    StringText=strrep(StringText,'®','&reg;');
    if any(StringText > 190)
        StringText=strrep(StringText,'¿','&#191;');
        StringText=strrep(StringText,'À','&#192;');
        StringText=strrep(StringText,'Á','&#193;');
        StringText=strrep(StringText,'Â','&#194;');
        StringText=strrep(StringText,'Ã','&#195;');
        StringText=strrep(StringText,'Ä','&#196;');
        StringText=strrep(StringText,'Å','&#197;');
        StringText=strrep(StringText,'Æ','&#198;');
        StringText=strrep(StringText,'Ç','&#199;');
        StringText=strrep(StringText,'È','&#200;');
        StringText=strrep(StringText,'É','&#201;');
        StringText=strrep(StringText,'Ê','&#202;');
        StringText=strrep(StringText,'Ë','&#203;');
        StringText=strrep(StringText,'Ì','&#204;');
        StringText=strrep(StringText,'Í','&#205;');
        StringText=strrep(StringText,'Î','&#206;');
        StringText=strrep(StringText,'Ï','&#207;');
        StringText=strrep(StringText,'Ð','&#208;');
        StringText=strrep(StringText,'Ñ','&#209;');
        StringText=strrep(StringText,'Ò','&#210;');
        StringText=strrep(StringText,'Ó','&#211;');
        StringText=strrep(StringText,'Ô','&#212;');
        StringText=strrep(StringText,'Õ','&#213;');
        StringText=strrep(StringText,'Ö','&#214;');
        StringText=strrep(StringText,'×','&#215;');
        StringText=strrep(StringText,'Ø','&#216;');
        StringText=strrep(StringText,'Ù','&#217;');
        StringText=strrep(StringText,'Ú','&#218;');
        StringText=strrep(StringText,'Û','&#219;');
        StringText=strrep(StringText,'Ü','&#220;');
        StringText=strrep(StringText,'Ý','&#221;');
        StringText=strrep(StringText,'Þ','&#222;');
        StringText=strrep(StringText,'ß','&#223;');
        StringText=strrep(StringText,'à','&#224;');
        StringText=strrep(StringText,'á','&#225;');
        StringText=strrep(StringText,'â','&#226;');
        StringText=strrep(StringText,'ã','&#227;');
        StringText=strrep(StringText,'ä','&#228;');
        StringText=strrep(StringText,'å','&#229;');
        StringText=strrep(StringText,'æ','&#230;');
        StringText=strrep(StringText,'ç','&#231;');
        StringText=strrep(StringText,'è','&#232;');
        StringText=strrep(StringText,'é','&#233;');
        StringText=strrep(StringText,'ê','&#234;');
        StringText=strrep(StringText,'ë','&#235;');
        StringText=strrep(StringText,'ì','&#236;');
        StringText=strrep(StringText,'í','&#237;');
        StringText=strrep(StringText,'î','&#238;');
        StringText=strrep(StringText,'ï','&#239;');
        StringText=strrep(StringText,'ð','&#240;');
        StringText=strrep(StringText,'ñ','&#241;');
        StringText=strrep(StringText,'ò','&#242;');
        StringText=strrep(StringText,'ó','&#243;');
        StringText=strrep(StringText,'ô','&#244;');
        StringText=strrep(StringText,'õ','&#245;');
        StringText=strrep(StringText,'ö','&#246;');
        StringText=strrep(StringText,'÷','&#247;');
        StringText=strrep(StringText,'ø','&#248;');
        StringText=strrep(StringText,'ù','&#249;');
        StringText=strrep(StringText,'ú','&#250;');
        StringText=strrep(StringText,'û','&#251;');
        StringText=strrep(StringText,'ü','&#252;');
        StringText=strrep(StringText,'ý','&#253;');
        StringText=strrep(StringText,'þ','&#254;');
        StringText=strrep(StringText,'ÿ','&#255;');
    end
    StringText=deblank(StringText);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IdString = createId
global PLOT2SVG_globals
IdString = ['ID' sprintf('%06d',PLOT2SVG_globals.runningIdNumber)];
PLOT2SVG_globals.runningIdNumber = PLOT2SVG_globals.runningIdNumber + 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [projection, edges] = get_projection(ax,id)
global PLOT2SVG_globals
xc = get(ax,'CameraTarget');
phi = get(ax,'CameraViewAngle');
vi = get(ax,'View');
xi = get(ax,'XLim');
yi = get(ax,'YLim');
zi = get(ax,'ZLim');
projection.aspect_scaling = get(ax,'DataAspectRatio');
if strcmp(get(ax,'XScale'),'log')
    if strcmp(get(ax,'XLimMode'),'manual') && any(get(ax,'XLim') == 0)
        % Fix illegal scalings set by the user
        % -> replace all 0 with automatic calculated values (child limits)
        xlimM = get(ax,'XLim');
        set(ax,'XLimMode','auto');
        xlimA = get(ax,'XLim');
        xlimM(xlimM == 0) = xlimA(xlimM == 0);
        set(ax,'XLimMode','manual');
        set(ax,'XLim', xlimM);
    end
    xi = log10(get(ax,'XLim'));
end
if strcmp(get(ax,'YScale'),'log')
    if strcmp(get(ax,'YLimMode'),'manual') && any(get(ax,'YLim') == 0)
        % Fix illegal scalings set by the user
        % -> replace all 0 with automatic calculated values (child limits)
        ylimM = get(ax,'YLim');
        set(ax,'YLimMode','auto');
        ylimA = get(ax,'YLim');
        ylimM(ylimM == 0) = ylimA(ylimM == 0);
        set(ax,'YLimMode','manual');
        set(ax,'YLim', ylimM);
    end
    yi = log10(get(ax,'YLim'));
end
if strcmp(get(ax,'ZScale'),'log')
    if strcmp(get(ax,'ZLimMode'),'manual') && any(get(ax,'ZLim') == 0)
        % Fix illegal scalings set by the user
        % -> replace all 0 with automatic calculated values (child limits)
        zlimM = get(ax,'ZLim');
        set(ax,'ZLimMode','auto');
        zlimA = get(ax,'ZLim');
        zlimM(zlimM == 0) = zlimA(zlimM == 0);
        set(ax,'ZLimMode','manual');
        set(ax,'ZLim', zlimM);
    end
    zi = log10(get(ax,'ZLim'));
end
projection.xi = xi;
projection.yi = yi;
projection.zi = zi;
xc(1) = (xc(1) - xi(1))/(xi(2)-xi(1));
xc(2) = (xc(2) - yi(1))/(yi(2)-yi(1));
xc(3) = (xc(3) - zi(1))/(zi(2)-zi(1));
if strcmp(get(ax,'XScale'),'log')
    x = [xi(1) xi(2) xi(1) xi(2) xi(1) xi(2) xi(1) xi(2)] - log10(projection.aspect_scaling(1));
else
    x = [xi(1) xi(2) xi(1) xi(2) xi(1) xi(2) xi(1) xi(2)]/projection.aspect_scaling(1);
end
if strcmp(get(ax,'YScale'),'log')
    y = [yi(1) yi(1) yi(2) yi(2) yi(1) yi(1) yi(2) yi(2)] - log10(projection.aspect_scaling(2));
else
    y = [yi(1) yi(1) yi(2) yi(2) yi(1) yi(1) yi(2) yi(2)]/projection.aspect_scaling(2);
end
if strcmp(get(ax,'ZScale'),'log')
    z = [zi(1) zi(1) zi(1) zi(1) zi(2) zi(2) zi(2) zi(2)] - log10(projection.aspect_scaling(3));
else
    z = [zi(1) zi(1) zi(1) zi(1) zi(2) zi(2) zi(2) zi(2)]/projection.aspect_scaling(3);    
end
if PLOT2SVG_globals.octave
        projection.A = get(ax,'x_ViewTransform');
        projection.A(3,:) = -projection.A(3,:);
        projection.A(1:3,4) = 0;
else
    if strcmp(get(ax,'Projection'),'orthographic')
        projection.A = viewmtx(vi(1),vi(2));
    else
        projection.A = viewmtx(vi(1),vi(2),phi,xc);
    end
end
if (vi(1) == 0) && (mod(vi(2),90) == 0)
    projection.xyplane = true;
else
    projection.xyplane = false;
end
axpos = get(ax,'Position');
figpos = get(id,'Position');
[m,n] = size(x);
x4d = [x(:),y(:),z(:),ones(m*n,1)]';
x2d = projection.A*x4d;
x2 = zeros(m,n); y2 = zeros(m,n); z2 = zeros(m,n);
x2(:) = x2d(1,:)./x2d(4,:);
y2(:) = x2d(2,:)./x2d(4,:);
projection.ax = ax;
projection.xrange = max(x2) - min(x2);
projection.yrange = max(y2) - min(y2);
projection.xoffset = (max(x2) + min(x2))/2;
projection.yoffset = (max(y2) + min(y2))/2;
if (strcmp(get(ax,'PlotBoxAspectRatioMode'),'manual') || strcmp(get(ax,'DataAspectRatioMode'),'manual'))
      if (projection.xrange*axpos(4)*figpos(4) < projection.yrange*axpos(3)*figpos(3))
          projection.xrange = projection.yrange*axpos(3)*figpos(3)/axpos(4)/figpos(4);
      else
          projection.yrange = projection.xrange*axpos(4)*figpos(4)/axpos(3)/figpos(3);
      end
end
x2(:) = (x2d(1,:)./x2d(4,:) - projection.xoffset)/projection.xrange + 0.5;
y2(:) = (x2d(2,:)./x2d(4,:) - projection.yoffset)/projection.yrange + 0.5;
z2(:) =  x2d(3,:);
edges = [x2; y2; z2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x2,y2,z2] = project(x,y,z,projection)
[m,n] = size(x);
if strcmp(get(projection.ax,'XDir'),'reverse')
    xi = projection.xi;
    x = (1 - (x - xi(1)) / (xi(2) - xi(1))) * (xi(2) - xi(1)) + xi(1);
end
if strcmp(get(projection.ax,'YDir'),'reverse')
    yi = projection.yi;
    y = (1 - (y - yi(1)) / (yi(2) - yi(1))) * (yi(2) - yi(1)) + yi(1);
end
if strcmp(get(projection.ax,'XScale'),'log')
    x = x - log10(projection.aspect_scaling(1));
else
    x = x/projection.aspect_scaling(1);
end
if strcmp(get(projection.ax,'YScale'),'log')
    y = y - log10(projection.aspect_scaling(2));
else
    y = y/projection.aspect_scaling(2);
end
if strcmp(get(projection.ax,'ZScale'),'log')
    z = z - log10(projection.aspect_scaling(3));
else
    z = z/projection.aspect_scaling(3);
end
x4d = [x(:), y(:), z(:), ones(m*n,1)]';
x2d = projection.A*x4d;
x2 = zeros(m,n); y2 = zeros(m,n); z2 = zeros(m,n);
x2(:) = (x2d(1,:)./x2d(4,:) - projection.xoffset)/projection.xrange + 0.5;
y2(:) = (x2d(2,:)./x2d(4,:) - projection.yoffset)/projection.yrange + 0.5;
z2(:) =  x2d(3,:);
%x = [0 1 0 1 0 1 0 1];
%y = [0 0 1 1 0 0 1 1];
%z = [0 0 0 0 1 1 1 1];


function [f, v, fvc, fva] = surface2patch(s)
x = get(s, 'xdata');
y = get(s, 'ydata');
z = get(s, 'zdata');
c = get(s, 'cdata');
a = get(s, 'AlphaData');
if ~isempty(x) && ~isequal(size(x),size(z))
    x = repmat(x(:)',size(z,1),1);
end
if ~isempty(y) && ~isequal(size(y),size(z))
  y = repmat(y(:),1,size(z,2));
end
[m n]= size(z);
if isempty(x)
    [x y] = meshgrid(1:n, 1:m);
end
[cm cn cp] = size(c);
[am an ap] = size(a);
%if cm==(m-1) & cn==(n-1)
%    cmode = 'f';
%elseif cm==m & cn==n
%    cmode = 'v';
%else  
%    cmode = '';
%end
v = [x(:) y(:) z(:)];
q = [1:m*n-m-1]';
q(m:m:end) = [];
fvc = reshape(c, [cm*cn cp]);
fva = reshape(a, [am*an ap]);
f = [q q+m q+m+1 q+1];
