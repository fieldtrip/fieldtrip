function [source] = fixinside(source, opt);

% FIXINSIDE ensures that the region of interest (which is indicated by the
% field "inside") is consistently defined for source structures and volume
% structures. Furthermore, it solves backward compatibility problems.
%
% Use as
%   [source] = fixinside(source, 'logical');
% or
%   [source] = fixinside(source, 'index');

% Copyright (C) 2006, Robert Oostenveld
%
% $Log: fixinside.m,v $
% Revision 1.2  2009/10/01 12:31:52  jansch
% enforce inside vector to be column
%
% Revision 1.1  2009/07/02 08:04:36  roboos
% moved fixdimord and fixinside from private to public
%
% Revision 1.3  2007/07/17 10:34:34  roboos
% prevent display of source structure on screen by adding ;
%
% Revision 1.2  2007/05/29 16:08:13  ingnie
% also automatic construction of inside from source.dim (with roboos)
%
% Revision 1.1  2006/03/30 12:24:34  roboos
% Implemented private/fixinside, which facilitates consistent
% handling of source/volume data. Improved documentation. Fixed some
% bugs related to inconsistent handling of ROIs (i.e. inside/outside)
%


if nargin<2
  opt = 'logical';
end

if ~isfield(source, 'inside')
  if isfield(source, 'pos')
    % assume that all positions are inside the region of interest
    source.inside  = [1:size(source.pos,1)]';
    source.outside = [];
  elseif isfield(source, 'dim')
    source.inside  = [1:prod(source.dim)]';
    source.outside = [];
  end
end

if ~isfield(source, 'inside')
  % nothing to do
  return;
end

% determine the format
if isa(source.inside, 'logical')
  logicalfmt = 1;
elseif all(source.inside(:)==0 | source.inside(:)==1)
  source.inside = logical(source.inside);
  logicalfmt = 1;
else
  logicalfmt = 0;
end

if ~logicalfmt && strcmp(opt, 'logical')
  % convert to a logical array
  if ~isfield(source, 'outside')
    source.outside = [];
  end
  inside(source.inside)  = (1==1);  % true
  inside(source.outside) = (1==0);  % false
  source.inside = inside(:);
  if isfield(source, 'outside')
    source = rmfield(source, 'outside');
  end
elseif logicalfmt && strcmp(opt, 'index')
  % convert to a vectors with indices
  tmp = source.inside;
  source.inside  = find( tmp(:));
  source.outside = find(~tmp(:));
else
  % nothing to do
end
