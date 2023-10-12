function svgOffset(s, source, offset, result)
% Adds a feOffset SVG filter
% PRELIMINARY IMPLEMENTATION (Parameters may change)
%
% svgOffset(s, source, offset, result)
% Parameters:
%   s : Array of plot object handles
%   source : Any previous defined filter result string, 'SourceGraphic',
%            or 'SourceAlpha'.
%   offset : Offset value [x y]
%   result : String that identifies the filter result for following filter
%            stages.   
for i = 1:length(s)
    userdata = get(s(i),'UserData');
    if isfield(userdata, 'svg') && isfield(userdata.svg, 'Filter')
        next = length(userdata.svg.Filter) + 1;
    else
        next = 1;
    end
    userdata.svg.Filter(next).Subfilter.Type = 'feOffset';
    userdata.svg.Filter(next).Subfilter.Source = source;
    userdata.svg.Filter(next).Subfilter.Offset = offset;
    userdata.svg.Filter(next).Subfilter.Result = result;
    set(s(i),'UserData', userdata);
end
