function [res] = fiff_pick_info(info,sel)
%
% [res] = fiff_pick_info(info,sel)
%
% Pick desired channels from measurement info
%
% res       - Info modified according to sel
% info      - The original data
% sel       - List of channels to select
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%
%
%   Revision 1.2  2006/04/23 15:29:40  msh
%   Added MGH to the copyright
%
%   Revision 1.1  2006/04/22 10:59:30  msh
%   Added fiff_pick_info
%
%

me='MNE:fiff_pick_info';

if nargin == 1
    sel = [];
elseif nargin ~= 2
    error(me,'Incorrect number of arguments');
end

res = info;
if isempty(sel)
    return;
end


res.chs      = res.chs(sel);
res.ch_names = res.ch_names(sel);
res.nchan    = length(sel);

return;

end
