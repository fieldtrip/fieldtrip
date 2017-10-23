function [grootte,ordening] = cluster_ward(realdist,optie,absoluut,aantal)

% CLUSTER_WARD Hierarchical cluster algorithm
%
% [grootte,ordening] = cluster(realdist,optie,absoluut,aantal)
%
% input: distance matrix, clustering optie, absolute or relative scale,
%        number of required clusters
% output: weights of cluster, ordening of dendogram
% three options: single-link, complete-link, and ward's algorithm
% produces dendogram
%
% See also WARD_DISTANCE

% Copyright (C) 1998, Tom Heskens

if nargin < 2,
  optie = 'ward';
end
if nargin < 3,
  absoluut = 1;
end
if nargin < 4,
  aantal = 2;
end

N = length(realdist);

if (optie(1:4) == 'sing'),
   a1 = 0.5*ones(1,N);
   a2 = 0.5*ones(1,N);
   a3 = zeros(1,N);
   a4 = -0.5*ones(1,N);
   ward = 0;
elseif (optie(1:4) == 'comp'),
   a1 = 0.5*ones(1,N);
   a2 = 0.5*ones(1,N);
   a3 = zeros(1,N);
   a4 = 0.5*ones(1,N);
   ward = 0;
elseif (optie(1:4) == 'ward'),
   ward = 1;
   a4 = zeros(1,N);
end



gew = ones(N,N);
maxie = 1000*max(max(realdist));
for i=1:N,
for j=1:i,
    realdist(i,j) = maxie;
end
end
level = zeros(N,1);
samen = zeros(N,2);
for i=1:N,
    [d1,flap] = min(realdist);
    [level(i),c2] = min(d1);
    c1 = flap(c2);
    samen(i,:) = [c1,c2];
    if (i>1),
       gew(:,i) = gew(:,i-1);
    end
    if (ward),
       n3 = gew(:,i);
       n1 = gew(c1,i)*ones(size(n3));
       n2 = gew(c2,i)*ones(size(n3));
       noemer = n1+n2+n3;
       a1 = (n1+n3)./noemer;
       a2 = (n2+n3)./noemer;
       a3 = -n3./noemer;
    end
    for j=1:c1-1,
        realdist(j,c1) = a1(j)*realdist(j,c1)+a2(j)*realdist(j,c2)+...
           a3(j)*realdist(c1,c2)+a4(j)*abs(realdist(j,c1)-realdist(j,c2));
        realdist(j,c2) = maxie;
    end
    for j=c1+1:c2-1,
        realdist(c1,j) = a1(j)*realdist(c1,j)+a2(j)*realdist(j,c2)+...
           a3(j)*realdist(c1,c2)+a4(j)*abs(realdist(c1,j)-realdist(j,c2));
        realdist(j,c2) = maxie;
    end
    for j=c2+1:N,
        realdist(c1,j) = a1(j)*realdist(c1,j)+a2(j)*realdist(c2,j)+...
           a3(j)*realdist(c1,c2)+a4(j)*abs(realdist(c1,j)-realdist(c2,j));
        realdist(c2,j) = maxie;
    end
    realdist(c1,c2) = maxie;
    gew(c1,i) = gew(c1,i)+gew(c2,i);
    gew(c2,i) = 0;
end

if (absoluut==0) level = 1:N; end

dummy = N:-1:1;
samen = samen(dummy,:);
level = level(dummy);
gew = [gew(:,dummy),ones(N,1)];
from = samen(:,1);
off = samen(:,2);

lijn = zeros(1,N);
laatst = zeros(1,N);
lijn(1) = N/2;
laatst(1) = 1.02*level(2);

groepnaam = 1:N;
grootte = zeros(size(groepnaam));
hold off
plot(laatst(1),lijn(1),'.');
%axis('off');
hold on
for i=2:N,
    oldy = lijn(from(i));
    oldx = laatst(from(i));
    newx = level(i);
    if (i-1 == aantal),
       groepnaam = find(lijn > 0);
       grootte = gew(groepnaam,i);
    end
    gew1 = gew(from(i),i+1);
    gew2 = gew(off(i),i+1);
    [d1,d2] = max([gew1,gew2]);
    teken = 3-2*d2;
    newy1 = oldy-teken*gew2/2;
    newy2 = oldy+teken*gew1/2;
    plot([oldx,newx],[oldy,oldy]);
    plot([newx,newx],[newy1,newy2]);

    if (gew1==1),
       plot([0,newx],[newy1,newy1]);
    end
    if (gew2==1),
       plot([0,newx],[newy2,newy2]);
    end

    laatst(from(i)) = newx;
    laatst(off(i)) = newx;
    lijn(from(i)) = newy1;
    lijn(off(i)) = newy2;

end
hold off

[dummy,ordening] = sort(lijn);
for i=1:length(groepnaam),
    qqq(i) = find(ordening == groepnaam(i));
end
[dummy,ind] = sort(qqq);
grootte = grootte(ind)';

set(gca,'YTick',[1:N]-0.5);
clear string
flap = sprintf('%%%dd',floor(log10(N))+1);
for i=1:N,
  string(i,:) = sprintf(flap,ordening(i));
end
set(gca,'YTickLabel',string)

