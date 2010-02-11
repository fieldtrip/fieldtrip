function slr_view_data(t, x, dim, w)
% View the distribution of data in feature space with projection of two
% dimensional plane specified by 'dim'. 
% 
%  -- Example
% > slr_view_data(t, x)
% View features value of first two dimension of "x".
% > slr_view_data(t, x, dim)
% View features value of two dimension of "x" projected onto the dimension
% specified 'dim'.
% > slr_view_data(t, x, dim, w)
% % View features value of two dimension of "x" projected onto the
% dimension specified 'dim' and the boundary specified by "w". 
%
% 2006/09/20 OY
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.

if nargin < 3
    dim1 = 1;
    dim2 = 2;
else
    dim1 = dim(1);
    dim2 = dim(2);
end
    
t = label2num(t);

ix1 = find(t == 1);
ix2 = find(t == 2);

plot(x(ix1,dim1),x(ix1,dim2), 'bx');
hold on;
plot(x(ix2,dim1),x(ix2,dim2), 'rd');

if nargin == 4
    w1 = w(dim1);
    w2 = w(dim2);
    c = w(end);
    
    minx1 = min(x(:,dim1))*1.0;
    maxx1 = max(x(:,dim1))*1.0;
    minx2 = min(x(:,dim2))*1.0;
    maxx2 = max(x(:,dim2))*1.0;
    
    [X1, X2] = meshgrid([minx1:(maxx1-minx1)/10:maxx1],[minx2:(maxx2-minx2)/10:maxx2]);
    Z = w1.*X1+w2.*X2+c;
    
    hold on;
    [col,h]=contour(X1,X2,Z,[0,0]);
    set(h, 'linewidth', 2);
end