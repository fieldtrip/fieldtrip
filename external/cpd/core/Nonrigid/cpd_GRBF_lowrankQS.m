function [Q,S]=cpd_GRBF_lowrankQS(Y, beta, numeig, eigfgt);

[M,D]=size(Y);
hsigma=sqrt(2)*beta;

OPTS.issym=1;
OPTS.isreal=1;
OPTS.disp=0;

% if we do not use FGT we can construct affinity matrix G and find the
% first eigenvectors/values directly
if ~eigfgt
   G=cpd_G(Y,Y,beta);
   [Q,S]=eigs(G,numeig,'lm',OPTS);
   return; % stop here
end 

%%% FGT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if we use FGT than we can find eigenvectors without constructing G,
% all we need to give is the matrix vector product Gx, which can be
% implemented through FGT.

e          = 8;      % Ratio of far field (default e = 10)
K          = round(min([sqrt(M) 100])); % Number of centers (default K = sqrt(Nx))
p          = 6;      % Order of truncation (default p = 8)


[Q,S]=eigs(@grbf,M,numeig,'lm',OPTS);


    function y=grbf(x,beta) % return y=Gx, without explicitelly constructing G
        [xc , A_k] = fgt_model(Y' , x', hsigma, e,K,p);
        y = fgt_predict(Y' , xc , A_k , hsigma,e);
    end

end