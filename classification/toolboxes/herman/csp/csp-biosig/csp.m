function [V,d] = csp(ECM,k_pat,Mode)
% CSP computes common spatial patterns
% 	this version supports multiple classes using a One-vs-Rest scheme
%
% [V] = csp(ECM)
%
% ECM(k,:,:) is the extended covariance matrices (see COVM) for class k  
% V	each column represents one CSP component. 
%
% REFERENCES: 
% [1] Koles ZJ, Soong AC.
% 	EEG source localization: implementing the spatio-temporal decomposition approach.
% 	Electroencephalogr Clin Neurophysiol. 1998 Nov;107(5):343-52
% [2] Ramoser, H.; Muller-Gerking, J.; Pfurtscheller, G.;
% 	Optimal spatial filtering of single trial EEG during imagined hand movement
%	Rehabilitation Engineering, IEEE Transactions on [see also IEEE Trans. on Neural Systems and Rehabilitation]
%	Volume 8,  Issue 4,  Dec. 2000 Page(s):441 - 446 
% [3] Dornhege, G.; Blankertz, B.; Curio, G.; Muller, K.-R.;
%    	Boosting bit rates in noninvasive EEG single-trial classifications by feature combination and multiclass paradigms
% 	Biomedical Engineering, IEEE Transactions on
% 	Volume 51,  Issue 6,  June 2004 Page(s):993 - 1002
%	Digital Object Identifier 10.1109/TBME.2004.827088 
% [4] Lemm, S.; Blankertz, B.; Curio, G.; Muller, K.-R.;
%    	Spatio-Spectral Filters for Improving the Classification of Single Trial EEG
% 	Biomedical Engineering, IEEE Transactions on
% 	Volume 52,  Issue 9,  Sept. 2005 Page(s):1541 - 1548
%	Digital Object Identifier 10.1109/TBME.2005.851521 

%	$Id: csp.m,v 1.2 2008/10/02 13:28:44 schloegl Exp $
%	Copyright (C) 2007 by Alois Schloegl <a.schloegl@ieee.org>
%	This is part of the BIOSIG-toolbox http://biosig.sf.net/

sz = size(ECM); 

for k=1:sz(1),
  	[mu,sd,COV(k,:,:),xc,N,R2]=decovm(squeeze(ECM(k,:,:)));
end; 


V = repmat(NaN,sz(2)-1,2*k_pat*sz(1));
d = V(1,:);

if 0,

elseif strcmpi(Mode,'CSP0');  
	% common diagonalization 
	[P,D] = eig(squeeze(sum(COV,1))); 
	P = diag(sqrt(1./diag(D)))*P';

    % class specific components
	%V = repmat(NaN,sz(2)-1,2*sz(1));
	for k = 1:sz(1), 
  		C = P * squeeze(COV(k,:,:)) * P';
  		[R,d1]  = eig((C+C')/2);
		[d1,ix] = sort(diag(d1)); 
	  	V(:,2*k_pat*k+[1-2*k_pat:0]) = P'*R(:,ix([1:k_pat,end-k_pat+1:end]));
	  	d(1,2*k_pat*k+[1-2*k_pat:0]) = d1(ix([1:k_pat,end-k_pat+1:end]))';
	end; 

elseif strcmpi(Mode,'CSP3');  
	% do actual CSP calculation as generalized eigenvalues
	R = permute(COV,[2,3,1]);
	for k = 1:sz(1), 
		[W,D] = eig(R(:,:,k),sum(R,3));
		V(:,2*k_pat*k+[1-2*k_pat:0]) = W(:,[1:k_pat,end-k_pat+1:end]);
        d(1,2*k_pat*k+[1-2*k_pat:0]) = diag(D([1:k_pat,end-k_pat+1:end],[1:k_pat,end-k_pat+1:end]));
	end;
end; 


