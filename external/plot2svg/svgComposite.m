function svgComposite(s, source1, source2, operator, result, k)
% Adds a feComposite SVG filter
% PRELIMINARY IMPLEMENTATION (Parameters may change)
%
% svgComposite(s, source1, source2, operator, result)
% Parameters:
%   s : Array of plot object handles
%   source1 : Any previous defined filter result string, 'SourceGraphic',
%             or 'SourceAlpha'.
%   source2 : Any previous defined filter result string, 'SourceGraphic',
%             or 'SourceAlpha'.
%   operator : Operator 'over','in','out','atop','xor','arithmetic'
%              -> see SVG 1.1 specification. 'arithmetic' is not yet
%              supported.
%   result : String that identifies the filter result for following filter
%            stages.   
for i = 1:length(s)
    userdata = get(s(i),'UserData');
    if isfield(userdata, 'svg') && isfield(userdata.svg, 'Filter')
        next = length(userdata.svg.Filter) + 1;
    else
        next = 1;
    end
    userdata.svg.Filter(next).Subfilter.Type = 'feComposite';
    userdata.svg.Filter(next).Subfilter.Source1 = source1;
    userdata.svg.Filter(next).Subfilter.Source2 = source2;
    userdata.svg.Filter(next).Subfilter.Operator = operator;
    if nargin > 5
        userdata.svg.Filter(next).Subfilter.k = k;
    end
    userdata.svg.Filter(next).Subfilter.Result = result;
    set(s(i),'UserData', userdata);
end