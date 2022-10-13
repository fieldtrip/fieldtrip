function ShowCmap(M, h)
% ShowCmap(M, h)
% a function to display and examine a colormap
% M (nx3): color map
% h: figure handle
%
% see also readXcol, ROIcmap, rgbdectohex, and ScaleToMap
%
   figure(h); clf
   subplot (211);
   colormap(M);
   image ([1:1:size(M,1)]);
   str = sprintf('%d colors color map. Pick colors with mouse\nHit "enter" to quit', size(M,1));
   title(str, 'fontsize',14);
   drawnow;

   i = 0;
   subplot (269);cla

   x1 = 1;
   while (~isempty(x1)),
      [x1,y] = ginput (1);
      if (~isempty(x1)),
         x1 = floor(x1(1)-0.5)+1;
         subplot (269);
         addsquare([0 i], [2.5 1+i], M(x1,:)); hold on
         plot (-0.2, i+0.5, 'k*');
         axis ([-1 3 -1 11]);
         str = sprintf ('Col %d: %.3g %.3g %.3g', x1, M(x1,1), M(x1,2), M(x1,3));
         ht = text (3, 0.3+i, 0, str, 'fontsize',14, 'color', M(x1,:));
         title (str);
         i = rem(i +1, 10);
      end
   end
return;
