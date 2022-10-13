function [labelx, labely, labelz] = coordsys2label(coordsys, format, both)

% COORDSYS2LABEL returns the labels for the three axes, given the symbolic
% string representation of the coordinate system.
%
% Use as
%   [labelx, labely, labelz] = coordsys2label(coordsys, format, both)
%
% The scalar argument 'format' results in return values like these
%   1) 'right'
%   2) 'the right'
%   3) '+X (right)'
%
% The boolean argument 'both' results in return values like these
%   0) 'right'              i.e. only the direction that it is pointing to
%   1) {'left' 'right'}     i.e. both the directions that it is pointing from and to
%
% See also FT_DETERMINE_COORDSYS, FT_PLOT_AXES, FT_HEADCOORDINATES, SETVIEWPOINT

% Copyright (C) 2017-2021, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$


% FIXME this function could also return a label for the origin
% see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=3304

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(coordsys) && ~strcmp(coordsys, 'unknown')
  
  % the coordsys consists of three letters for the direction of the positive axes
  % or of a string that relates to external software, an atlas or a template
  if length(coordsys)==3 && length(intersect(coordsys, 'rlapis'))==3
    axis = 'XYZ';
    label = cell(1,3);
    for i=1:3
      switch coordsys(i)
        case 'l'
          label{i} = {['-' axis(i) ' (right)'],      ['+' axis(i) ' (left)']};
        case 'r'
          label{i} = {['-' axis(i) ' (left)'],       ['+' axis(i) ' (right)']};
        case 'i'
          label{i} = {['-' axis(i) ' (superior)'],   ['+' axis(i) ' (inferior)']};
        case 's'
          label{i} = {['-' axis(i) ' (inferior)'],   ['+' axis(i) ' (superior)']};
        case 'a'
          label{i} = {['-' axis(i) ' (posterior)'],  ['+' axis(i) ' (anterior)']};
        case 'p'
          label{i} = {['-' axis(i) ' (anterior)'],   ['+' axis(i) ' (posterior)']};
        otherwise
          ft_error('incorrect letter in the coordsys');
      end % switch
    end % for each of the three axes
    
    labelx = label{1};
    labely = label{2};
    labelz = label{3};
    
  else
    switch lower(coordsys)
      case {'ras' 'scanras' 'nifti' 'neuromag' 'itab' 'acpc' 'spm' 'mni' 'tal'}
        % the nifti coordinate system is defined according to https://doi.org/10.1016/j.jneumeth.2016.03.001 and https://github.com/fieldtrip/website/pull/444
        % but see also https://brainder.org/2012/09/23/the-nifti-file-format/ which shows that it can be more complex than that
        labelx = {'-X (left)'      '+X (right)'   };
        labely = {'-Y (posterior)' '+Y (anterior)'};
        labelz = {'-Z (inferior)'  '+Z (superior)'};
      case {'als' 'ctf' '4d' 'bti' 'eeglab'}
        labelx = {'-X (posterior)' '+X (anterior)'};
        labely = {'-Y (right)'     '+Y (left)'};
        labelz = {'-Z (inferior)'  '+Z (superior)'};
      case {'lps' 'scanlps' 'dicom'}
        labelx = {'-X (right)'     '+X (left)'};
        labely = {'-Y (anterior)'  '+Y (posterior)'};
        labelz = {'-Z (inferior)'  '+Z (superior)'};
      case {'rsp' 'paxinos'}
        labelx = {'-X (left)'      '+X (right)'};
        labely = {'-Y (inferior)'  '+Y (superior)'};
        labelz = {'-Z (anterior)'  '+Z (posterior)'};
      otherwise
        % the coordinate system is unknown
        ft_warning('cannot determine the labels for the axes of the "%s" coordinate system', coordsys);
        labelx = {'-X (unknown)' '+X (unknown)'};
        labely = {'-Y (unknown)' '+Y (unknown)'};
        labelz = {'-Z (unknown)' '+Z (unknown)'};
    end
    
  end % if it matches with 'rlasif'
  
else
  % the coordinate system is not specified
  labelx = {'-X (unknown)' '+X (unknown)'};
  labely = {'-Y (unknown)' '+Y (unknown)'};
  labelz = {'-Z (unknown)' '+Z (unknown)'};
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch format
  case 1
    % remove the -, +, X, Y, Z , (, ) and space
    labelx{1} = trim(labelx{1});
    labelx{2} = trim(labelx{2});
    labely{1} = trim(labely{1});
    labely{2} = trim(labely{2});
    labelz{1} = trim(labelz{1});
    labelz{2} = trim(labelz{2});
    
  case 2
    % remove the -, +, X, Y, Z , (, ) and space
    labelx{1} = trim(labelx{1});
    labelx{2} = trim(labelx{2});
    labely{1} = trim(labely{1});
    labely{2} = trim(labely{2});
    labelz{1} = trim(labelz{1});
    labelz{2} = trim(labelz{2});
    % insert 'the' for left and right
    if any(strcmp(labelx{1}, {'left',  'right'})), labelx{1} = ['the ' labelx{1}]; end
    if any(strcmp(labelx{2}, {'left',  'right'})), labelx{2} = ['the ' labelx{2}]; end
    if any(strcmp(labely{1}, {'left',  'right'})), labely{1} = ['the ' labely{1}]; end
    if any(strcmp(labely{2}, {'left',  'right'})), labely{2} = ['the ' labely{2}]; end
    if any(strcmp(labelz{1}, {'left',  'right'})), labelz{1} = ['the ' labelz{1}]; end
    if any(strcmp(labelz{2}, {'left',  'right'})), labelz{2} = ['the ' labelz{2}]; end
    
  case 3
    % keep it as it is
    
  otherwise
    ft_error('unsupported option')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch both
  case false
    % only keep the direction that is it pointing towards
    labelx = labelx{2};
    labely = labely{2};
    labelz = labelz{2};
  case true
    % keep it as it is
  otherwise
    ft_error('unsupported option')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = trim(str)
str(str=='+') = [];
str(str=='-') = [];
str(str=='(') = [];
str(str==')') = [];
str(str=='X') = [];
str(str=='Y') = [];
str(str=='Z') = [];
str(str==' ') = [];
