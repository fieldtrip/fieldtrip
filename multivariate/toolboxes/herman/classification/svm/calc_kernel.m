
function K = calc_kernel(data1,data2,kernel_type,kernel_param)

% CALC_KERNEL implements kernel calculations
% 
% INPUT
%       data1, data2                  - two data sets for which kernel product is to be calculated
%       kernel_type, kernel_param (p) - kernel settings
%             'linear'           - x1*x2 (p=1 but it does not matter) 
%             'poly'             - (x1*x2+1)^p
%             'gaussian'         - exp(norm(x1-x2)/(2*p^2))  
%             'rbf'              - allows for inhomogeneous (p is a vector then) weighing of feature components by sigmas in Gaussian kernel (the same as 'gaussian' for scalar p)
%             'htrbf'            - heavy tailed rbf
%
%
% USE
%      it is an auxiliary function used within SVMMETHOD object (see svmmethod.m)

% Pawel Herman, 2009

if nargin<4
    kernel_param = 1;
end

if nargin<3
    kernel_type='linear';
end


switch lower(kernel_type)
    case 'poly'

        if ~isscalar(kernel_param)
            error('kernel_param should be a scalar for kernel poly');
        end

        ps= data1 * data2';

        if kernel_param > 1
            K =(ps+1).^kernel_param;
        else
            error('kernel_param should be greater than 1 for kernel poly');
        end;

    case 'rbf'   %handles both homogenous (see gaussian) and inhomogenous cases
        
        d = length(kernel_param);
        if d==1          
            kernel_param=ones(1,d)*kernel_param;
        elseif size(kernel_param,1)==d && size(kernel_param,2)==1
            kernel_param = kernel_param';
        else
            error('Kernel_param should be a scalar or a vector');
        end

        sigmas = diag(1./kernel_param.^2);
        ps = data1 * sigmas * data2';
        [n,m]=size(ps);
        norm1 = sum(data1.^2*sigmas,2);
        norm2 = sum(data2.^2*sigmas,2);
        ps = -2*ps + repmat(norm1,1,m) + repmat(norm2',n,1) ;

        K = exp(-ps/2);

    case 'htrbf'    % heavy tailed RBF 
        
        if length(kernel_param)~=2
            error('There are supposed to be two parameters for heavy tailed RBF kernel (htrbf)');
        end
        b=kernel_param(2);
        a=kernel_param(1);
        n2 = size(data2,1);
        n1 = size(data1,1);
        for i=1:n2
            ps(:,i) = sum( abs((data1.^a - ones(n1,1)*data2(i,:).^a)).^b   ,2);
        end;

        K = exp(-ps);

    case 'gaussian'    
        if length(kernel_param)~=1
            error('There is supposed to be one parameter for gaussian kernel (gaussian->homogeneous rbf)');
        end
        a = kernel_param;
        n2 = size(data2,1);
        n1 = size(data1,1);
        for i=1:n2
            ps(:,i) = sum( abs((data1 - ones(n1,1)*data2(i,:))).^2 ,2)./kernel_param.^2;
        end;

        K = exp(-ps/2);

    case 'linear'    
        
        K = data1 * data2';
end




