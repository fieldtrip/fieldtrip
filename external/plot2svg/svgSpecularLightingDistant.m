function svgSpecularLightingDistant(s, source, specularConstant, specularExponent, surfaceScale, azimuth, elevation, result)
% Adds a feSpecularLighting SVG filter with distant light source
% PRELIMINARY IMPLEMENTATION (Parameters may change)
%
% svgSpecularLightingDistant(s, source, specularConstant, specularExponent, surfaceScale, azimuth, elevation, result)
% Parameters:
%   s : Array of plot object handles
%   source : Any previous defined filter result string, 'SourceGraphic',
%            or 'SourceAlpha'.
%   specularConstant : Specular constant
%   specularExponent : Specular exponent
%   surfaceScale : Surface scaling factor
%   azimuth : Light azimuth angle [deg], typical 225.
%   elevation : Light elevation angle [deg], typical 45.
%   result : String that identifies the filter result for following filter
%            stages.
for i = 1:length(s)
    userdata = get(s(i),'UserData');
    if isfield(userdata, 'svg') && isfield(userdata.svg, 'Filter')
        next = length(userdata.svg.Filter) + 1;
    else
        next = 1;
    end
    userdata.svg.Filter(next).Subfilter.Type = 'feSpecularLighting';
    userdata.svg.Filter(next).Subfilter.Source = source;
    userdata.svg.Filter(next).Subfilter.Result = result;
    userdata.svg.Filter(next).Subfilter.SpecularConstant = specularConstant; 
    userdata.svg.Filter(next).Subfilter.SpecularExponent = specularExponent;
    userdata.svg.Filter(next).Subfilter.SurfaceScale = surfaceScale;
    userdata.svg.Filter(next).Subfilter.LightType = 'feDistantLight';
    userdata.svg.Filter(next).Subfilter.Azimuth = azimuth;
    userdata.svg.Filter(next).Subfilter.Elevation = elevation;
    set(s(i),'UserData', userdata);
end
