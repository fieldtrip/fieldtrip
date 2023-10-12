function svgLuminanceToAlpha(s, source, result)
% Adds a feColorMatrix SVG filter that maps luminance to alpha values
% PRELIMINARY IMPLEMENTATION (Parameters may change)
%
% svgLuminanceToAlpha(s, source, result)
% Parameters:
%   s : Array of plot object handles
%   source : Any previous defined filter result string, 'SourceGraphic',
%            or 'SourceAlpha'.
%   result : String that identifies the filter result for following filter
%            stages.
for i = 1:length(s)
    userdata = get(s(i),'UserData');
    if isfield(userdata, 'svg') && isfield(userdata.svg, 'Filter')
        next = length(userdata.svg.Filter) + 1;
    else
        next = 1;
    end
    userdata.svg.Filter(next).Subfilter.Type = 'feColorMatrix';
    userdata.svg.Filter(next).Subfilter.Source = source;
    userdata.svg.Filter(next).Subfilter.Result = result;
    userdata.svg.Filter(next).Subfilter.ColorType = 'luminanceToAlpha';
    set(s(i),'UserData', userdata);
end
