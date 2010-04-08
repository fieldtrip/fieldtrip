function [X,Y,Width,Height,Lbl] = read_lay(layoutname);

% READ_LAY reads an electrode or gradiometer layout file
% Layout files are used for topoplotting and multiplotting.
%
% Use as
%   [X, Y, Width, Height, Lbl] = read_lay(layoutname)

if ~exist(layoutname, 'file')
  error(sprintf('could not open layout file: %s', layoutname));
end

% discard the channel number
[chNum,X,Y,Width,Height,Lbl] = textread(layoutname,'%f %f %f %f %f %q');
