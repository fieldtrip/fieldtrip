function svgBoundingBox(s, type, overlap, visible)
% Configures the bounding box of a SVG filter
% PRELIMINARY IMPLEMENTATION (Parameters may change)
%
% svgBoundingBox(s, type, overlap, visible)
% Parameters:
%   s : Array of plot object handles
%   type : [axes, element, relative]
%          Sets the filter bounding box to cover the axis reagion (axes), the
%          element extension (element or relative). Axes gives usually the
%          best results but may be slower.
%   overlap : Many filters need an overlap to work correctly.
%             Typical values for type 'axes' and 'element' -> 10
%             Typical values for type 'relative' -> 0.1
%   visible : Debugging functionality to see the bounding box used for an
%             object
for i = 1:length(s)
    userdata = get(s(i),'UserData');
    userdata.svg.BoundingBox.Visible = visible;    % Useful for debugging of bounding box for filters
    userdata.svg.BoundingBox.Type = type;          % [axes, element, relative]
    userdata.svg.BoundingBox.Overlap = overlap;
    set(s(i),'UserData', userdata);
end
