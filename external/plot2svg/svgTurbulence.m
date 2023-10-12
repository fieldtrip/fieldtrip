function svgTurbulence(s, frequency, octaves, seed, stitch, type, result)
% Adds a feTurbulence SVG filter
% PRELIMINARY IMPLEMENTATION (Parameters may change)
%
% svgTurbulence(s, frequency, octaves, seed, stitch, type, result)
% Parameters:
%   s : Array of plot object handles
%   frequency : Base frequency, typical 0.05
%   octaves : Octaves, typical 2
%   seed : Seed value for random generator, typical 0..255
%   stitch : Stitch tiles [stitch, noStitch]
%   type : Turbulence type [fractalNoise, turbulence]
%   result : String that identifies the filter result for following filter
%            stages.   
for i = 1:length(s)
    userdata = get(s(i),'UserData');
    if isfield(userdata, 'svg') && isfield(userdata.svg, 'Filter')
        next = length(userdata.svg.Filter) + 1;
    else
        next = 1;
    end
    userdata.svg.Filter(next).Subfilter.Type = 'feTurbulence';
    userdata.svg.Filter(next).Subfilter.BaseFrequency = frequency;
    userdata.svg.Filter(next).Subfilter.NumOctaves = octaves;
    userdata.svg.Filter(next).Subfilter.Seed = seed;
    userdata.svg.Filter(next).Subfilter.StitchTiles = tiles;
    userdata.svg.Filter(next).Subfilter.TurbulenceType = type;
    userdata.svg.Filter(next).Subfilter.Result = result;
    set(s(i),'UserData', userdata);
end
