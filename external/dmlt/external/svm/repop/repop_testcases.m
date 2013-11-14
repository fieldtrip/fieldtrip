function []=repop_testcases(testType)
%
% This file contains lots of test-cases to test the performance of the repop
% files vs. the matlab built-ins.
%
% N.B. there appears to be a bug in MATLAB when comparing mixed
% complex/real + double/single values in a max/min
% 
% Copyright 2006-     by Jason D.R. Farquhar (jdrf@zepler.org)
% Permission is granted for anyone to copy, use, or modify this
% software and accompanying documents for any uncommercial
% purposes, provided this copyright notice is retained, and note is
% made of any changes that have been made. This software and
% documents are distributed without any warranty, express or
% implied
if ( nargin<1 || isempty(testType) ) testType={'acc','timing'}; end;

if ( ~isempty(strmatch('acc',testType)) ) 
fprintf('-------------------  Accuracy tests -------------------\n');

X=complex(randn(10,100),randn(10,100)); 
Y=complex(randn(size(X)),randn(size(X)));
fprintf('\n****************\n Double Real X, Double Real Y\n******************\n')
accuracyTests(real(X),real(Y),'dRdR')
fprintf('\n****************\n Double Complex X, Double Real Y\n******************\n')
accuracyTests(X,real(Y),'dCdR')
fprintf('\n****************\n Double Real X, Double Complex Y\n******************\n')
accuracyTests(real(X),Y,'dRdC')
fprintf('\n****************\n Double Complex X, Double Complex Y\n******************\n')
accuracyTests(X,Y,'dCdC')

fprintf('\n****************\n Double Real X, Single Real Y\n******************\n')
accuracyTests(real((X)),real(single(Y)),'dRsR')
fprintf('\n****************\n Double Complex X, Single Real Y\n******************\n')
accuracyTests((X),real(single(Y)),'dCsR')
fprintf('\n****************\n Double Real X, Single Complex Y\n******************\n')
accuracyTests((real(X)),single(Y),'dRsC')
fprintf('\n****************\n Double Complex X, Single Complex Y\n******************\n')
accuracyTests((X),single(Y),'dCsC')

fprintf('\n****************\n Single Real X, Double Real Y\n******************\n')
accuracyTests(real(single(X)),real((Y)),'sRdR')
fprintf('\n****************\n Single Complex X, Double Real Y\n******************\n')
accuracyTests(single(X),real((Y)),'sCdR')
fprintf('\n****************\n Single Real X, Double  Complex Y\n******************\n')
accuracyTests(single(real(X)),(Y),'sRdC')
fprintf('\n****************\n Single Complex X,Double Complex Y\n******************\n')
accuracyTests(single(X),(Y),'sCdC')

fprintf('\n****************\n Single Real X, Single Real Y\n******************\n')
accuracyTests(real(single(X)),real(single(Y)),'sRsR')
fprintf('\n****************\n Single Complex X, Single Real Y\n******************\n')
accuracyTests(single(X),real(single(Y)),'sCsR')
fprintf('\n****************\n Single Real X, Single Complex Y\n******************\n')
accuracyTests(single(real(X)),single(Y),'sRsC')
fprintf('\n****************\n Single Complex X, Single Complex Y\n******************\n')
accuracyTests(single(X),single(Y),'sCsC')

fprintf('All tests passed\n');
end

if ( ~isempty(strmatch('timing',testType)) )
fprintf('-------------------  Timing tests -------------------\n');

X=complex(randn(100,1000),randn(100,1000)); 
Y=complex(randn(size(X)),randn(size(X)));
timingTests(real(X),real(Y),'[100x1000] RR');
timingTests(X,Y,'[100x1000] CC');

end

return;

function []=accuracyTests(X,Y,str)
% PLUS
unitTest([str ' Matx + col Vec'],X,'+',Y(:,1),repop(X,'+',Y(:,1)),X+repmat(Y(:,1),[1,size(X,2)]));  
unitTest([str ' Matx + row Vec'],X,'+',Y(1,:),repop(X,'+',Y(1,:)),X+repmat(Y(1,:),[size(X,1),1]));
unitTest([str ' Matx + Matx(:,1:2)'],X,'+',Y(:,1:2),repop(X,'+',Y(:,1:2),'m'),X+repmat(Y(:,1:2),[1,size(X,2)/2]));
unitTest([str ' Matx + Matx'],X,'+',Y(:,:),repop(X,'+',Y(:,:)),X+repmat(Y(:,:),[1,1])); 
% TIMES
unitTest([str ' Matx * col Vec'],X,'*',Y(:,1),repop(X,'*',Y(:,1)),X.*repmat(Y(:,1),[1,size(X,2)]));
unitTest([str ' Matx * row Vec'],X,'*',Y(1,:),repop(X,'*',Y(1,:)),X.*repmat(Y(1,:),[size(X,1),1]));
unitTest([str ' Matx * Matx(:,1:2)'],X,'*',Y(:,1:2),repop(X,'*',Y(:,1:2),'m'),X.*repmat(Y(:,1:2),[1,size(X,2)/2]));
unitTest([str ' Matx * Matx'],X,'*',Y(:,:),repop(X,'*',Y(:,:)),X.*repmat(Y(:,:),[1,1]));  

% other operations
unitTest([str ' Matx - row vec'],X,'-',Y(:,1),repop(X,'-',Y(:,1)),X-repmat(Y(:,1),[1,size(X,2)])); 
unitTest([str ' Matx / row vec'],X,'/',Y(:,1),repop(X,'/',Y(:,1)),X./repmat(Y(:,1),[1,size(X,2)]));
unitTest([str ' Matx \ row vec'],X,'\',Y(:,1),repop(X,'\',Y(:,1)),X.\repmat(Y(:,1),[1,size(X,2)]));
unitTest([str ' Matx ^ row vec'],X,'^',Y(:,1),repop(X,'^',Y(:,1)),X.^repmat(Y(:,1),[1,size(X,2)]),1e-5);
unitTest([str ' Matx == row vec'],X,'==',Y(:,1),repop(X,'==',Y(:,1)),X==repmat(Y(:,1),[1,size(X,2)]));
unitTest([str ' Matx ~= row vec'],X,'~=',Y(:,1),repop(X,'~=',Y(:,1)),X~=repmat(Y(:,1),[1,size(X,2)]));
% N.B. repop tests with complex inputs use the magnitude of the input!
tX=X; tY=Y; if( ~isreal(X) | ~isreal(Y) ) tX=abs(X); tY=abs(Y); end;
unitTest([str ' Matx < row vec'],X,'<',Y(:,1),repop(X,'<',Y(:,1)),tX<repmat(tY(:,1),[1,size(X,2)]));
unitTest([str ' Matx > row vec'],X,'>',Y(:,1),repop(X,'>',Y(:,1)),tX>repmat(tY(:,1),[1,size(X,2)]));
unitTest([str ' Matx <= row vec'],X,'<=',Y(:,1),repop(X,'<=',Y(:,1)),tX<=repmat(tY(:,1),[1,size(X,2)]));
unitTest([str ' Matx >= row vec'],X,'>=',Y(:,1),repop(X,'>=',Y(:,1)),tX>=repmat(tY(:,1),[1,size(X,2)]));
%unitTest([str ' min Matx, row vec'],repop('min',X,Y(:,1)),min(X,repmat(Y(:,1),[1, size(X,2)])));
%unitTest([str ' max Matx, row vec'],repop('max',X,Y(:,1)),max(X,repmat(Y(:,1),[1,size(X,2)])));
%return;
%function []=inplaceaccuracyTests(X,Y,str)

return;

% Inplace operations test
% PLUS
% N.B. need the Z(1)=Z(1); to force to make a "deep" copy, i.e. not just
% equal pointers
Z=X;unitTest([str ' Matx + col Vec (inplace)'],Z,'+',Y(:,1),repop(Z,'+',Y(:,1),'i'),X+repmat(Y(:,1),[1,size(X,2)]));  
Z=X;unitTest([str ' Matx + row Vec (inplace)'],Z,'+',Y(1,:),repop(Z,'+',Y(1,:),'i'),X+repmat(Y(1,:),[size(X,1),1]));
Z=X;unitTest([str ' Matx + Matx(:,1:2) (inplace)'],Z,'+',Y(:,1:2),repop(Z,'+',Y(:,1:2),'mi'),X+repmat(Y(:,1:2),[1,size(X,2)/2]));
Z=X;unitTest([str ' Matx + Matx (inplace)'],Z,'+',Y(:,:),repop(Z,'+',Y(:,:),'i'),X+repmat(Y(:,:),[1,1])); 
% TIMES
Z=X;unitTest([str ' Matx * col Vec (inplace)'],Z,'*',Y(:,1),repop(Z,'*',Y(:,1),'i'),X.*repmat(Y(:,1),[1,size(X,2)]));
Z=X;unitTest([str ' Matx * row Vec (inplace)'],Z,'*',Y(1,:),repop(Z,'*',Y(1,:),'i'),X.*repmat(Y(1,:),[size(X,1),1]));
Z=X;unitTest([str ' Matx * Matx(:,1:2) (inplace)'],Z,'*',Y(:,1:2),repop(Z,'*',Y(:,1:2),'mi'),X.*repmat(Y(:,1:2),[1,size(X,2)/2]));
Z=X;unitTest([str ' Matx * Matx (inplace)'],Z,'*',Y(:,:),repop(Z,'*',Y(:,:),'i'),X.*repmat(Y(:,:),[1,1]));  
% other operations
Z=X;unitTest([str ' Matx - row vec (inplace)'],Z,'-',Y(:,1),repop(Z,'-',Y(:,1),'i'),X-repmat(Y(:,1),[1,size(X,2)])); 
Z=X;unitTest([str ' Matx \ row vec (inplace)'],Z,'\',Y(:,1),repop(Z,'\',Y(:,1),'i'),X.\repmat(Y(:,1),[1,size(X,2)]));
Z=X;unitTest([str ' Matx / row vec (inplace)'],Z,'/',Y(:,1),repop(Z,'/',Y(:,1),'i'),X./repmat(Y(:,1),[1,size(X,2)]));
Z=X;unitTest([str ' Matx ^ row vec (inplace)'],Z,'^',Y(:,1),repop(Z,'^',Y(:,1),'i'),X.^repmat(Y(:,1),[1,size(X,2)]),1e-5);
%unitTest([str ' min Matx, row vec (inplace)'],repop('min',X,Y(:,1),'i'),min(X,repmat(Y(:,1),[1, size(X,2)])));
%unitTest([str ' max Matx, row vec (inplace)'],repop('max',X,Y(:,1),'i'),max(X,repmat(Y(:,1),[1,size(X,2)])));


function []=timingTests(X,Y,str)
%PLUS
fprintf('%s Matx + Scalar\n',str);
tic; for i=1:1000; Z=repop(X,'+',10); end;
fprintf('%30s %gs\n','repop',toc/1000);
tic; Z=X;Z(1)=Z(1);for i=1:1000; Z=repop(Z,'+',10,'i'); end;
fprintf('%30s %gs\n','repop (inplace)',toc/1000);
tic; for i=1:1000; T=X+10;end;
fprintf('%30s %gs\n','MATLAB',toc/1000);

fprintf('%s Matx + col vec\n',str);
tic, for i=1:1000; Z=repop(X,'+',Y(:,1)); end;
fprintf('%30s %gs\n','repop',toc/1000); % = .05  / .01
tic, Z=X;Z(1)=Z(1); for i=1:1000; Z=repop(Z,'+',Y(:,1),'i'); end;
fprintf('%30s %gs\n','repop (inplace)',toc/1000); % = .05  / .01
tic, for i=1:1000; Z=X+repmat(Y(:,1),1,size(X,2));end;
fprintf('%30s %gs\n','MATLAB',toc/1000); % = .05  / .01

fprintf('%s Matx + row vec\n',str);
tic, for i=1:1000; Z=repop(X,'+',Y(1,:)); end;
fprintf('%30s %gs\n','repop',toc/1000); % = .05  / .01
tic, Z=X;Z(1)=Z(1); for i=1:1000; Z=repop(Z,'+',Y(1,:),'i'); end;
fprintf('%30s %gs\n','repop (inplace)',toc/1000); % = .05  / .01
tic, for i=1:1000; Z=X+repmat(Y(1,:),size(X,1),1);end;
fprintf('%30s %gs\n','MATLAB',toc/1000); % = .05  / .01

fprintf('%s Matx + Matx(:,1:2)\n',str);
tic; for i=1:1000; Z=repop(X,'+',Y(:,1:2),'m'); end;
fprintf('%30s %gs\n','repop',toc/1000); % = .05  / .01   
tic, Z=X;Z(1)=Z(1); for i=1:1000; Z=repop(Z,'+',Y(:,1:2),'im'); end;
fprintf('%30s %gs\n','repop (inplace)',toc/1000); % = .05  / .01
tic; for i=1:1000; T=X+repmat(Y(:,1:2),[1,size(X,2)/2]);end;
fprintf('%30s %gs\n','MATLAB',toc/1000); % = .05  / .01


%TIMES
fprintf('%s Matx * col Vec\n',str);
tic; for i=1:1000; Z=repop(X,'*',Y(:,1)); end;
fprintf('%30s %gs\n','repop',toc/1000);
tic; Z=X;Z(1)=Z(1); for i=1:1000; Z=repop(Z,'*',Y(:,1),'i'); end;
fprintf('%30s %gs\n','repop (inplace)',toc/1000);
tic; for i=1:1000; Z=X.*repmat(Y(:,1),1,size(X,2));end;
fprintf('%30s %gs\n','MATLAB',toc/1000);
tic; for i=1:1000; Z=spdiags(Y(:,1),0,size(X,1),size(X,1))*X;end;
fprintf('%30s %gs\n','MATLAB (spdiags)',toc/1000);

fprintf('%s Matx * row Vec\n',str);
tic; for i=1:1000; Z=repop(X,'*',Y(1,:)); end;
fprintf('%30s %gs\n','repop',toc/1000);
tic; Z=X;Z(1)=Z(1); for i=1:1000; Z=repop(Z,'*',Y(1,:),'i'); end;
fprintf('%30s %gs\n','repop (inplace)',toc/1000);
tic; for i=1:1000; Z=X.*repmat(Y(:,1),1,size(X,2));end;
fprintf('%30s %gs\n','MATLAB',toc/1000);
tic; for i=1:1000; Z=X*spdiags(Y(1,:)',0,size(X,2),size(X,2));end;
fprintf('%30s %gs\n','MATLAB (spdiags)',toc/1000);

fprintf('%s Matx * Matx(:,1:2)\n',str);
tic; for i=1:1000; Z=repop(X,'*',Y); end;
fprintf('%30s %gs\n','repop',toc/1000);
tic; Z=X;Z(1)=Z(1); for i=1:1000; Z=repop(Z,'*',Y,'i'); end;
fprintf('%30s %gs\n','repop (inplace)',toc/1000);
tic; for i=1:1000; Z=X.*repmat(Y(:,1:2),[1,size(X,2)/2]);end;
fprintf('%30s %gs\n','MATLAB',toc/1000);

return;

% simple function to check the accuracy of a test and report the result
function [testRes,trueRes,diff]=unitTest(testStr,A,op,B,testRes,trueRes,tol)
global LOGFILE;
if ( ~isempty(LOGFILE) ) % write tests and result to disc
   writeMxInfo(LOGFILE,A);
   fprintf(LOGFILE,'%s\n',op);
   writeMxInfo(LOGFILE,B);
   fprintf(LOGFILE,'=\n');
   writeMxInfo(LOGFILE,trueRes);
   fprintf(LOGFILE,'\n');
end
if ( nargin < 7 ) 
   if ( isa(trueRes,'double') ) tol=1e-11; 
   elseif ( isa(trueRes,'single') ) tol=1e-5; 
   elseif ( isa(trueRes,'integer') ) 
      warning('Integer inputs!'); tol=1;       
   elseif ( isa(trueRes,'logical') ) tol=0;
   end
end;
diff=abs(testRes-trueRes)./max(1,abs(testRes+trueRes));
fprintf('%45s = %0.3g ',testStr,max(diff(:)));
if ( max(diff(:)) > tol ) 
   if ( exist('mimage') )
      mimage(squeeze(testRes),squeeze(trueRes),squeeze(diff))
   end
   warning([testStr ': failed!']);
   fprintf('Type return to continue\n'); keyboard;
else
   fprintf('Passed \n');
end
