function h=similaritygraph(coord,S,lim,lcolor)
%function h=similaritygraph(coord,[S],[lim],[lcolor])
%
%PURPOSE
%
%To visualize a weighted graph
% 
%INPUT 
%
%[An argument in brackets is optional. If it isn't  given or it's
% an empty matrix/string, the function will use a default value.] 
%
% coord    (Mx2 matrix) coordinates
% [S]      (MxM matrix) similarity matrix: if this not given, no edges
%            are drawn, only vertices (black dots).
% [lim]     (vector) 1xK of monotonically increasing values between
%            0 and 1 sets limits for line colors: default [0.5 0.7 0.9]
% [lcolor] (Kx3 matrix) of RGB colors that define line colors
%            for each interval ]lim(1),lim(2)], ]lim(2),lim(3)],
%            ..., ]lim(end-1),lim(end)]  default: shades of red
%
%OUTPUT
%
% h         (struct) various graphic handles and other info 
%             (see details)
%
%DETAILS
%
%Draws a weighted graph where the black points are vertices and
%the red shaded lines between them are edges. The level is
%determined by dividing values of S (the weights) into bins between the
% <= lim(1)                 not drawn 
% ]lim(1),lim(2)]           colormap(1,:)    default: light red
%  ...                        ...               ...
% ]lim(end-1),lim(end)]     colormap(end,:)  default: bright red 
% >lim(2)                   not drawn
%
%NOTE The function always adds to a plot (turns 'hold on' temporarily).
%
%Output is a structure of handles to graphic objects etc
%
% h.graph   (vector) graphic handles to all edges (lines)
% h.vertex  (vector) graphic handles to all vertices (markers)
% h.example (vector) contains one example (graphic handle) of each
%             different line type (gray level); needed for legend 
% h.text    (cell array of strings) contains labels for lines in h.example
%             legend(h.example,h.text) sets a proper legend for the
%             graph line colors
%
%If there are no values in ith bin, the shade is left out, and
%h.graph(i)=[], h.example(i)=NaN, and h.text(i)='' will be set.
%
%USED IN
% icassoGraph

%COPYRIGHT NOTICE
%This function is a part of Icasso software library
%Copyright (C) 2003-2005 Johan Himberg
%
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

isHold=ishold;
hold on;

if nargin<3|isempty(lim),
  lim=[.5 .75 .9 1];
end

if nargin<2,
  S=zeros(size(coord,1));
end

lim=lim(:)';

if any(diff(lim)<=0),
  error('Limits must be monotonically increasing');
end

h.edge=[];
h.vertex=[];
h.example=[];
h.text=[];

% Change this to alter vertex color and line attrib. and line color
% scale
MARKERSIZE=3;
MARKERCOLOR=[0 0 0];
LINEWIDTH=1;

Nclass=length(lim)-1;
N=size(S,1);

if all(S(:)<lim(1)),
  isLines=0;
else
  isLines=1;
end

if isLines,
  if nargin<4|isempty(lcolor)
    lcolor=redscale(Nclass+1);
    lcolor=redscale(Nclass+1,0.9);
    lcolor(1:end-1,:)=lcolor(2:end,:);
  elseif size(lcolor,1)~=Nclass,
    error('Wrong number of colors.');
  end
  
  for i=1:Nclass, 
    S_=S;
    S_(S_<=lim(i) | S_>lim(i+1))=0; 
    links(:,:,i)=S_;
    hold on;
    [dummy,dummy,hd]=som_grid(S_,[N 1],'coord',coord,'marker','none', ...
			      'linecolor',lcolor(i,:),'linewidth',i.*LINEWIDTH);
    if ~isempty(hd),
      h.example(i)=hd(1); h.edge{i}=hd;
      h.text{i}=[sprintf('%0.2f',lim(i)) '<s_{ij}\leq ' sprintf('%0.2f',lim(i+1))]; 
    else
      % No lines in this slot
      h.example(i)=NaN; h.edge{i}=[];
      h.text{i}='';
    end
  end
  
  [dummy,h.vertex]=som_grid(S_,[N 1],'coord',coord,'marker','o','line','none',...
			    'markercolor',MARKERCOLOR,'markersize',MARKERSIZE);
else
  
  [dummy,h.vertex]=som_grid('rect',[N 1],'coord',coord,'marker','o','line','none',...
			    'markercolor',MARKERCOLOR,'markersize', ...
			    MARKERSIZE); 
  h.edge=[]; 
  h.example=h.vertex(1); 
  lim=[]; 
  h.text{1}='all s_{ij} below threshold: no lines shown';
end

if length(MARKERSIZE)==N,
  set(h.node,'markerfacecolor','none');
end

if ~isHold,
  hold off;
end

if nargout==0,
  clear h;
end
