function grad = itab2grad(header_info)

% ITAB2GRAD converts the original Chieti ITAB header structure into a gradiometer
% definition that is compatible with FieldTrip forward and inverse computations
%
% See also READ_HEADER

% Copyright (C) 2009, Robert Oostenveld, Donders Institute for Brain, Cognition and Behaviour
% Copyright (C) 2009, Stefania Della Penna, ITAB, University Chiety, Italy
%
% $Log: itab2grad.m,v $
% Revision 1.1  2009/10/16 07:31:18  roboos
% renamed chieti into itab for consistency with other formats
%
% Revision 1.1  2009/10/13 10:11:50  roboos
% first implementation, based on code from Stefania
%

grad = struct;
for i=1:header_info.nmagch
  grad.label{i}   = header_info.ch(i).label;
  grad.pnt(i,1:3) = [header_info.ch(i).pos(1).r_s.comp];
  grad.ori(i,1:3) = [header_info.ch(i).pos(1).u_s.comp];
end
grad.unit  = 'mm';
grad.tra   = eye(header_info.nmagch);
grad.label = grad.label(:); % should be column vector
