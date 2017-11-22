function out = spm_dartel_resids(job)
% Generate residuals in a form suitable for generating a Fisher kernel
% FORMAT spm_dartel_residuals(job)
% job.flowfields
% job.images
% job.template
% job.K
%
% The aim is to obtain better pattern recognition through using
% Fisher kernels.  See  Bishop's PRML or the work of Jaakkola and
% Haussler for more information.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_dartel_resids.m 5668 2013-10-03 18:34:18Z guillaume $


PG = job.template{1};
PI = job.images;
PU = job.flowfields;
K  = job.K;
fwhm = job.fwhm;

N  = nifti(PG);
dm = size(N.dat);
n  = dm(4);


spm_progress_bar('Init',dm(3),'Generating Covariances','Planes complete');
% Load mean values from template
mu = zeros([dm(1:3),n+1],'single');
mu(:,:,:,end)=1;
for i=1:dm(4),
    mu(:,:,:,i)   = N.dat(:,:,:,i);
    mu(:,:,:,end) = mu(:,:,:,end)-mu(:,:,:,i); % background
end
mu = min(max(mu,0),1);
mu = (mu+0.05)/(1+(n+1)*0.05); % Note the fudge factor


% Generate weighting images derived from a Gaussian approximation
% to the multinomial distribution.
wt = zeros([dm(1:3),n,n+1],'single');
cv = zeros([dm(1:2),n+1,n+1]);
I  = eye(n+1);
for k=1:dm(3),

    % Compute voxel-wise covariance matrix for multinomial
    % distribution.  This is the inverse of the Hessian matrix
    % in Equation 4.110 of Bishop's PRML book.
    for i=1:n+1,
        for j=1:n+1,
            cv(:,:,i,j) = (I(i,j)-mu(:,:,k,i)).*mu(:,:,k,j);
        end
    end

    % For rotating on to plane
    R = null(ones(n+1))';

    % Convert to a form suited to building a Fisher kernel.
    % See the workings at the bottom of the file in order to
    % understand better.
    % Note that the dimensionality is reduced because all
    % points lie on e.g., x1+x2+x3=1, where x1>=0, x2>=0, x3>=0
    for i=1:dm(2),
        for j=1:dm(1),
            [u,s]=svd(reshape(cv(j,i,:),n+1,n+1));
            s   = diag(s);
            s   = diag(s(1:n).^(-0.5));
            u   = u(:,1:n);
            tmp = R*u*s*u';
            wt(j,i,k,:,:) = reshape(tmp,1,1,1,n,n+1);
        end
    end
    spm_progress_bar('Set',k);
end
spm_progress_bar('Clear');


% Re-load the means, because the means used for generating the
% covariance matrices were modified slightly (with a fudge factor).
mu(:,:,:,end)=1;
for i=1:n,
    mu(:,:,:,i)   = N.dat(:,:,:,i);
    mu(:,:,:,end) = mu(:,:,:,end)-mu(:,:,:,i);
end
mu = min(max(mu,0),1);



% Now generate the residual images
spm_progress_bar('Init',numel(PI{1}),'Creating Resids');

for i=1:numel(PI{1}),

    % Derive deformations and Jacobian determinants from the
    % flow fields.
    NU = nifti(PU{i});
    [pth,nam,ext,num] = spm_fileparts(NU.dat.fname);
    fprintf('%s: def',nam); drawnow;
    u      = single(squeeze(NU.dat(:,:,:,1,:)));
    [y,dt] = dartel3('Exp',u,[K -1 1]);
    y1 = double(y(:,:,:,1));
    y2 = double(y(:,:,:,2));
    y3 = double(y(:,:,:,3));
    clear y

    % Determinant for Jacobian transformation of variables
    % Note that the square root of this value is used so that
    % the integral of the variance is conserved after warping,
    % rather than the usual mean.
    % On second thoughts, maybe it should not be square rooted,
    % as this makes more sense if the residuals are subsequently
    % smoothed.
    % dt = sqrt(dt);
 
    f  = zeros([NU.dat.dim(1:3),n+1],'single');
 
    [pth,nam,ext,num] = spm_fileparts(PI{1}{i});
    NI = nifti(fullfile(pth,[nam ext]));
    NO = NI;
    NO.dat.fname     = fullfile('.',['resid_' nam '.nii']);
    NO.dat.scl_slope = 1.0;
    NO.dat.scl_inter = 0.0;
    NO.dat.dtype     = 'float32-le';
    NO.dat.dim       = [NU.dat.dim(1:3) 1 n];
    NO.mat           = NU.mat;
    NO.mat0          = NU.mat;
    NO.mat_intent    = 'Aligned';
    NO.mat0_intent   = 'Aligned';
    NO.descrip = 'Normalised residuals';
    create(NO);


    % Compute the warped tissue probabilities
    f(:,:,:,end) = 1; % Background
    for j=1:n,
        fprintf(' %d', j); drawnow;

        [pth,nam,ext,num] = spm_fileparts(PI{j}{i});
        NI = nifti(fullfile(pth,[nam ext]));
        if sum(sum((NI.mat  - NU.mat ).^2)) < 0.0001 && ...
           sum(sum((NI.mat0 - NU.mat0).^2)) < 0.0001,
            ty1 = y1;
            ty2 = y2;
            ty3 = y3;
        else
            mat = NI.mat;
            if ~isempty(NI.extras) && isstruct(NI.extras) && isfield(NI.extras,'mat'),
                mat1 = NI.extras.mat;
                if size(mat1,3) >= j && sum(sum(mat1(:,:,j).^2)) ~=0,
                    mat = mat1;
                end;
            end;
            M   = mat\NU.mat0;
            ty1 = M(1,1)*y1 + M(1,2)*y2 + M(1,3)*y3 + M(1,4);
            ty2 = M(2,1)*y1 + M(2,2)*y2 + M(2,3)*y3 + M(2,4);
            ty3 = M(3,1)*y1 + M(3,2)*y2 + M(3,3)*y3 + M(3,4);
        end;
        spl_param  = [3 3 3  1 1 1];
        cf         = spm_bsplinc(NI.dat(:,:,:,1,1),spl_param);
        f(:,:,:,j) = spm_bsplins(cf,ty1,ty2,ty3,spl_param);
        f(:,:,:,end) = f(:,:,:,end) - f(:,:,:,j);
        clear cf ty1 ty2 ty3
    end;
    clear y1 y2 y3


    % Residuals
    f = f - mu;

    % Scale in order to obtain approximate Mahalinobis distances
    % for use in generating Fisher kernel
    for j=1:n,
        res = zeros(dm(1:3));
        for k=1:n+1,
            res = res + wt(:,:,:,j,k).*f(:,:,:,k);
            fprintf('.'); drawnow;
        end
        res = res.*dt; % Jacobian transform
        if fwhm>0,
            vx = sqrt(sum(NO.mat(1:3,1:3).^2));
            spm_smooth(res,res,fwhm/vx); % Note the abuse of MATLAB
        end
        NO.dat(:,:,:,1,j) = res; 
    end

    fprintf('\n');
    spm_progress_bar('Set',i);
end
spm_progress_bar('Clear');

PI  = job.images{1};
out.files = cell(numel(PI),1);
for i=1:numel(PI),
    [pth,nam,ext] = fileparts(PI{i});
    fname         = fullfile(pwd,['resid_' nam '.nii']);
    out.files{i} = fname;
end





if false,
    % Workings in order to check the code


    %------------------------------------------------------------
    % Check covariance matrices for building Fisher kernel.
    %------------------------------------------------------------
    % Random "tissue classes"
    n = 10000;
   %x = double(rand(n,1)>0.9);
    x = rand(n,4).*repmat([0.5 0.5 0.5 0.5],n,1);
    [mx,ind] = max(x,[],2);
    x = zeros(size(x));
    for i=1:n, x(i,ind(i)) = 1; end;


    % Means
    mn = sum(x)/size(x,1)

    % Empirically determined covariance matrix
    res=x-repmat(mn,n,1);
    vr = res'*res/n;


    % Covariance matrix from means
    cv = zeros(4);
    I  = eye(4);
    for i=1:4,
        for j=1:4,
            cv(i,j)=(I(i,j)-mn(i))*mn(j);
        end
    end
    cv



    %------------------------------------------------------------
    % More workings related to the use of Fisher kernels
    % (Jaakkola and Haussler, 1999).  In this example, it is
    % for a binomial distribution
    %------------------------------------------------------------
    % Eq. 6.32 of PRML (numerically derived)
    g  = inline(['((log(mu+1e-5)*x + log(1-mu-1e-5)*(1-x)) '...
                '- (log(mu-1e-5)*x + log(1-mu+1e-5)*(1-x)))/2e-5'],'mu','x');

    for it=1:4,
        x1 = double(rem(it-1,2));
        x2 = double(rem(floor((it-1)/2),2));

        % A range f mean
        MU=0.05:0.01:0.95;
        for i=1:numel(MU),
            mu     = MU(i);

            % Fisher information matrix (Eq. 6.34)
            F      = mu*g(mu,1)^2 + (1-mu)*g(mu,0)^2;

            % Eq. 6.33
            sc1(i) = g(mu,x1)*inv(F)*g(mu,x2);

            % The aproach used for binomial distributions
            vr     = (1-mu)*mu;
            sc2(i) = (x1-mu)*inv(vr)*(x2-mu);
        end
        subplot(2,2,it);
        plot(MU,[sc1',sc2'])
        xlabel('\mu')
        ylabel('k(x,x'')')
        title(['x=' num2str(x1) ', x''=' num2str(x2)]);
    end


end

