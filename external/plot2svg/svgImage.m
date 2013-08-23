function svgImage(s, file, aspectRatio, result)
% Adds a feImage SVG filter
% PRELIMINARY IMPLEMENTATION (Parameters may change)
%
% svgImage(s, file, aspectRatio, result)
% Parameters:
%   s : Array of plot object handles
%   file : Pixel graphics file name (png or jpeg) with extension.
%   aspectRatio: 'none' -> scale to bounding box limits
%                'xMinYMin meet', 'xMinYMin slice', 'xMidYMid meet', ...
%                -> see SVG 1.1 specification
%   result : String that identifies the filter result for following filter
%            stages.   
for i = 1:length(s)
    userdata = get(s(i),'UserData');
    if isfield(userdata, 'svg') && isfield(userdata.svg, 'Filter')
        next = length(userdata.svg.Filter) + 1;
    else
        next = 1;
    end
    userdata.svg.Filter(next).Subfilter.Type = 'feImage';
    userdata.svg.Filter(next).Subfilter.File = file;
    userdata.svg.Filter(next).Subfilter.AspectRatio = aspectRatio;
    userdata.svg.Filter(next).Subfilter.Result = result;
    set(s(i),'UserData', userdata);
end
