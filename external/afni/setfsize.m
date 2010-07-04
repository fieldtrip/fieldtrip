function h=setfsize(wx,wy)
% h=setfsize (wx,wy)
%
% Set figure position.	
%
% h = setfsize(X,Y)
%   X - Horizontal length.
%   Y - Vertical height.
% h is the figure handle to the figure just handled
% Sets current figure lower left position to X by Y pixels (measured from
% lower left corner of screen. On image, max pixels are about
% 1280, 1024, you can use  get (0,'ScreenSize') to get the screen size
%
%  adapted from function setfsize.m
%
%				Ziad Saad Sept 27 97


z = get(gcf,'position');
set(gcf,'position',[z(1) z(2) wx wy]);
h = gcf;
return;
