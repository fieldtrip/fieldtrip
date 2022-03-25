function h=setfpos(x,y)
%setfpos Set figure position.
%	
%
% h = setfpos(X,Y)
%   X - Horizontal length.
%   Y - Vertical height.
% h is the figure handle to the figure just handled
% Sets current figure lower left position to X by Y pixels (measured from
% lower left corner of screen. On image, max pixels are about
% 1264, 980
%
%  adapted from function setfsize.m
%
%				Ziad Saad Sept 27 97


z = get(gcf,'position');
set(gcf,'position',[x y z(3) z(4)]);
h = gcf;
return;
