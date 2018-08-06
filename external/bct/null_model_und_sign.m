function [W0,R] = null_model_und_sign(W,bin_swaps,wei_freq)
%NULL_MODEL_UND_SIGN     Random graphs with preserved weight, degree and
%                        strength distributions
%
%   W0 = null_model_und_sign(W);
%   W0 = null_model_und_sign(W,bin_swaps);
%   W0 = null_model_und_sign(W,bin_swaps,wei_freq);
%   [W0 R] = null_model_und_sign(W,bin_swaps,wei_freq);
%
%   This function randomizes an undirected network with positive and
%   negative weights, while preserving the degree and strength
%   distributions. This function calls randmio_und_signed.m
%
%   Inputs: W,          Undirected weighted connection matrix
%           bin_swaps,  Average number of swaps of each edge in binary randomization.
%                           bin_swap=5 is the default (each edge rewired 5 times)
%                           bin_swap=0 implies no binary randomization
%           wei_freq,   Frequency of weight sorting in weighted randomization
%                           wei_freq must be in the range of: 0 < wei_freq <= 1
%                           wei_freq=1 implies that weights are resorted at each step
%                                   (default in older [<2011] versions of MATLAB)
%                           wei_freq=0.1 implies that weights are sorted at each 10th step
%                                   (faster, default in newer versions of Matlab)
%
%   Output:     W0,     Randomized weighted connection matrix
%               R,      Correlation coefficient between strength sequences
%                           of input and output connection matrices
%
%   Notes:
%       The value of bin_swaps is ignored when binary topology is fully
%   connected (e.g. when the network has no negative weights).
%       Randomization may be better (and execution time will be slower) for
%   higher values of bin_swaps and wei_freq. Higher values of bin_swaps may
%   enable a more random binary organization, and higher values of wei_freq
%   may enable a more accurate conservation of strength sequences.
%       R are the correlation coefficients between positive and negative
%   strength sequences of input and output connection matrices and are
%   used to evaluate the accuracy with which strengths were preserved. Note
%   that correlation coefficients may be a rough measure of
%   strength-sequence accuracy and one could implement more formal tests
%   (such as the Kolmogorov-Smirnov test) if desired.
%
%   Example usage:
%
%   %Create random weights matrix
%
%   W=tril(randn(100),-1); W=W+W.';
%
%   %Compute one instance of null model (slow execution time):
%   %bin_swaps=5,   rewire each binary edge 5 times on average
%   %wei_freq=1,    sort all edges at every step
%
%   tic; [W0_slow R_slow]=null_model_und_sign(W,5,1); R_slow, toc
%
%   R_slow =
%       0.9720    0.9746
%   Elapsed time is 2.112715 seconds.
%
%   %Compute another instance of of null model (fast execution time):
%   %bin_swaps=5,   rewire each binary edge 5 times on average
%   %wei_freq=0.1,  sort all edges at every 10th step (10=1/0.1)
%
%   tic; [W0_fast R_fast]=null_model_und_sign(W,5,0.1); R_fast, toc
%
%   R_fast =
%       0.9623    0.9789
%   Elapsed time is 0.584797 seconds.
%
%
%   Reference: Rubinov and Sporns (2011) Neuroimage 56:2068-79
%
%
%   2011-2015, Mika Rubinov, U Cambridge

%   Modification History
%   Mar 2011: Original
%   Sep 2012: Edge-sorting acceleration
%   Dec 2015: Enforce preservation of negative degrees in sparse
%             networks with negative weights (thanks to Andrew Zalesky).

%#ok<*ASGLU>

if ~exist('bin_swaps','var')
    bin_swaps=5;
end
if ~exist('wei_freq','var')
    if nargin('randperm')==1
        wei_freq=1;
    else
        wei_freq=0.1;
    end
end

if wei_freq<=0 || wei_freq>1
    error('wei_freq must be in the range of: 0 < wei_freq <= 1.')
end
if wei_freq && wei_freq<1 && nargin('randperm')==1
    warning('wei_freq may only equal 1 in older (<2011) versions of MATLAB.')
    wei_freq=1;
end

n=size(W,1);                                                %number of nodes
W(1:n+1:end)=0;                                             %clear diagonal
Ap = W>0;                                                   %positive adjacency matrix
An = W<0;                                                   %negative adjacency matrix

if nnz(Ap)<(n*(n-1))                                        %if Ap is not full
    W_r  = randmio_und_signed(W,bin_swaps);
    Ap_r = W_r>0;
    An_r = W_r<0;
else
    Ap_r = Ap;
    An_r = An;
end

W0=zeros(n);                                                %null model network
for s=[1 -1]
    switch s                                                %switch sign (positive/negative)
        case 1
            S=sum(W.*Ap,2);                                 %positive strength
            Wv=sort(W(triu(Ap)));                           %sorted weights vector
            [I,J]=find(triu(Ap_r));                         %weights indices
            Lij=n*(J-1)+I;                                  %linear weights indices
        case -1
            S=sum(-W.*An,2);                                %negative strength
            Wv=sort(-W(triu(An)));                          %sorted weights vector
            [I,J]=find(triu(An_r));                         %weights indices
            Lij=n*(J-1)+I;                                  %linear weights indices
    end
    
    P=(S*S.');                                              %expected weights matrix
    
    if wei_freq==1
        for m=numel(Wv):-1:1                                %iteratively explore all weights
            [dum,Oind]=sort(P(Lij));                        %get indices of Lij that sort P
            r=ceil(rand*m);
            o=Oind(r);                                      %choose random index of sorted expected weight
            W0(Lij(o)) = s*Wv(r);                           %assign corresponding sorted weight at this index
            
            f = 1 - Wv(r)/S(I(o));                          %readjust expected weight probabilities for node I(o)
            P(I(o),:) = P(I(o),:)*f;                        %[1 - Wv(r)/S(I(o)) = (S(I(o)) - Wv(r))/S(I(o))]
            P(:,I(o)) = P(:,I(o))*f;
            f = 1 - Wv(r)/S(J(o));                          %readjust expected weight probabilities for node J(o)
            P(J(o),:) = P(J(o),:)*f;                        %[1 - Wv(r)/S(J(o)) = (S(J(o)) - Wv(r))/S(J(o))]
            P(:,J(o)) = P(:,J(o))*f;
            
            S([I(o) J(o)]) = S([I(o) J(o)])-Wv(r);          %readjust strengths of nodes I(o) and J(o)
            Lij(o)=[];                                      %remove current index from further consideration
            I(o)=[];
            J(o)=[];
            Wv(r)=[];                                       %remove current weight from further consideration
        end
    else
        wei_period = round(1/wei_freq);                     %convert frequency to period
        for m=numel(Wv):-wei_period:1                       %iteratively explore at the given period
            [dum,Oind]=sort(P(Lij));                        %get indices of Lij that sort P
            R=randperm(m,min(m,wei_period)).';
            O = Oind(R);
            W0(Lij(O)) = s*Wv(R);                           %assign corresponding sorted weight at this index
            
            WA = accumarray([I(O);J(O)],Wv([R;R]),[n,1]);   %cumulative weight
            IJu = any(WA,2);
            F = 1-WA(IJu)./S(IJu);
            F = F(:,ones(1,n));                             %readjust expected weight probabilities for node I(o)
            P(IJu,:) = P(IJu,:).*F;                         %[1 - Wv(r)/S(I(o)) = (S(I(o)) - Wv(r))/S(I(o))]
            P(:,IJu) = P(:,IJu).*F.';
            S(IJu) = S(IJu)-WA(IJu);                      	%re-adjust strengths of nodes I(o) and J(o)
            
            O=Oind(R);
            Lij(O)=[];                                      %remove current index from further consideration
            I(O)=[];
            J(O)=[];
            Wv(R)=[];                                       %remove current weight from further consideration
        end
    end
end
W0=W0+W0.';

rpos=corrcoef(sum( W.*(W>0)),sum( W0.*(W0>0)));
rneg=corrcoef(sum(-W.*(W<0)),sum(-W0.*(W0<0)));
R=[rpos(2) rneg(2)];
