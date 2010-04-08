function [X] = slr_make_kernel(x, mode, xcenter, R)
% Make a explanatry matrix using kernels
% 
% [X] = slr_make_kernel(x, mode, xcenter, R)
% 
% -- Input
% x : data      (N*D)
% xcenter  : center of kernel (Nref*D)
% R : width of Gaussian kernel
%
% -- Output
% X : Explanatory matrix 
%
% 2007/10/03 Okito Yamashita 
% * Linear Kernel
% 2005/12/04 Okito Yamashita
%
%     [ K(x(1), xc(1)), K(x(1), xc(2)), ... K(x(1), xc(Nref)) ]
%     [ K(x(2), xc(1)), K(x(2), xc(2)), ... K(x(2), xc(Nref)) ]
% K = [                        :                              ]
%     [                        :                              ]
%     [ K(x(N), xc(1)), K(x(N), xc(2)), ... K(x(N), xc(Nref)) ]
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.

if nargin < 2
    error('mode must be specified.')
end

if nargin < 3 | isempty(xcenter), % for training
    xcenter = x;
end

N = size(x,1);
D = size(x,2);
Ncenter = size(xcenter,1);
K = zeros(N, Ncenter);


switch mode
    case 'linear',
       for i = 1 : Ncenter
           xref = xcenter(i,:);
           K(:,i) = xref * x';
       end
       X = K;    

    case 'Gaussian',

        if nargin < 4
            error('-- You must specify the width parameter of Gaussian kernel');
        end

        for i = 1 : Ncenter
            xref = xcenter(i,:);
            d = x - repmat(xref,[N,1]);
            dd = sum(d.^2,2);
            K(:,i) = dd'/D;  % normalization
        end

        X = exp(-R*K);  %Gaussian kernel

    otherwise,


end
