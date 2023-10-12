function svgGaussianBlur(s, source, deviation, result)
% Adds a feGaussianBlur SVG filter
% PRELIMINARY IMPLEMENTATION (Parameters may change)
%
% svgGaussianBlur(s, source, deviation, result)
% Parameters:
%   s : Array of plot object handles
%   source : Any previous defined filter result string, 'SourceGraphic',
%            or 'SourceAlpha'.
%   deviation : Blur strength
%   result : String that identifies the filter result for following filter
%            stages.   
for i = 1:length(s)
    userdata = get(s(i),'UserData');
    if isfield(userdata, 'svg') && isfield(userdata.svg, 'Filter')
        next = length(userdata.svg.Filter) + 1;
    else
        next = 1;
    end
    userdata.svg.Filter(next).Subfilter.Type = 'feGaussianBlur';
    userdata.svg.Filter(next).Subfilter.Deviation = deviation;
    userdata.svg.Filter(next).Subfilter.Source = source;
    userdata.svg.Filter(next).Subfilter.Result = result;
    set(s(i),'UserData', userdata);
end
