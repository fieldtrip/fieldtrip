function [atlas, varargout] = handle_atlas_input(atlas, varargin)
% HANDLE_ATLAS_INPUT handles user input to specify volumetric atlases in some coordinate. It
% does two things: (1) call FT_READ_ATLAS to read the atlas from file, if it is specified as a
% string input, and (2) if the optional second data input argument is provided, and it has a
% coordsys and/or unit field, checks the coordinate systems and units of the atlas and the
% input against each other.
%
% This code was taken from ft_sourceplot to avoid duplication upon adding similar functionality
% to ft_sourceplot_interactive.

if ischar(atlas)
  % initialize the atlas
  [~, f, ~] = fileparts(atlas);
  fprintf(['reading ', f, ' atlas coordinates and labels\n']);
  atlas = ft_read_atlas(atlas);
end

% ensure that the atlas is formatted properly
if nargin > 1
  data_hasunit = isfield(varargin{1}, 'unit');
  data_hascoordsys = isfield(varargin{1}, 'coordsys');
  atlas = ft_checkdata(atlas, 'hasunit', data_hasunit, 'hascoordsys', data_hascoordsys);

  if data_hasunit
    % ensure that the units are consistent, convert the units if required
    atlas = ft_convert_units(atlas, varargin{1}.unit);
  end
  if data_hascoordsys
    atlas = fixcoordsys(atlas);
    for k = 1:numel(varargin)
      if ~isfield(varargin{k}, 'coordsys')
        continue;
      end
      % ensure that the coordinate systems match
      varargin{k} = fixcoordsys(varargin{k});
      assert(isequal(varargin{k}.coordsys, atlas.coordsys), 'coordinate systems do not match');
    end
  end
end

varargout = varargin;

end