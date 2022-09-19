function y = cauchyPdf(x,scale,location)
% Cauchy probability density function 
% 
% INPUT 
% x = variable
% scale = Scale of the distribution [1].
% location = x-value where the peak occurs [0]
% OUTPUT
% y = density

if nargin <3
    location =0;
    if nargin <2
        scale = 1;
    end
end

y =1/(pi.*scale).*(scale.^2)./((x-location).^2+scale.^2);
