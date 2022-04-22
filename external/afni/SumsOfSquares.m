function [err, fstat, intensity_new, dfterm_new, dfdenom, tnames_new, LC] = SumsOfSquares(y, NF, FL, ntot, nterms, Qd, s, sindices, dfbothSS, modw, modwo, tnames, dfterm, dfe, dsgn, N_Brik, Contr)
%
%   [err,] = ss.m ()
%
%Purpose:
%
%
%
%Input Parameters:
%
%	 y: values from all factor combinations at one voxel in vector format
%   ntot: total number of files = length of vector y
%   tername: term names such as
%   termlist:
%   dmat: design matrix in the converted regression model
%   cmat: constraints matrix
%
%Output Parameters:
%   err : 0 No Problem
%       : 1  Problems
%   fstat: a vector of F values for those terms at a voxel
%   intensity_new: intensity for those terms at a voxel, which are defined as the sqare root of MS terms
%   dfterm: vector of degree of freedom for those terms (numerators) at each voxel (but the same for all voxels).
%           This is why I don't differentiate them among the voxels in anova.m
%   dfdenom: vector of denominator
%
%Key Terms:
%
%More Info :
%
%
%
%
%     Author : Gang Chen
%     Date : Tue Mar 23 14:05:05 EST 2004
%     SSCC/NIMH/ National Institutes of Health, Bethesda MD 20892


%Define the function name for easy referencing
FuncName = 'ss.m';

%Debug Flag
DBG = 1;

%initailize return variables
err = 1;

dfboth = dfbothSS;
y0 = y;  % y gets mutated down below. Keep the original y so that the contrasts can be done later at the end.

% Fit the full model and compute the residual sum of squares
mu =  mean(y);
y = y - mu;         % two passes to improve accuracy
mu2 = mean(y);
%mu = mu + mu2;       % THIS CAN BE DELETED FOR MY CASE since it is not used from this point!!!!!
y = y - mu2;
sst = sum(y.^2);    % SSTO

[ssx, y] = dofit(y, Qd);      % y: predicted value; ssx:
sse = sst - ssx;   %SS(Error) = SS(Total) - SS(Rest)

% Fit each model, get its residual SS and d.f.
ssw(1:nterms) = ssx;    % for full model we already know the results

for j=length(sindices):-1:1
   % Find the next model index to fit
   k = sindices(j);

   % Look in unsorted arrays to see if we have already fit this model
   if j>nterms
      k0 = k+nterms;
   else
      k0 = k;
   end
   if dfboth(k0)~=-1
      continue
   end

   % Find the model with this index
   if (j > nterms)
      thismod = modwo(k, :);
   else
      thismod = modw(k, :);
   end

   % Fit this model
		
%   [ssx0,dfx0] = QRfit(X, y, C);
   ssx0 = dofit(y, s(j).Qdt);

   % Use these results for each term that requires them
   mod0 = repmat(thismod, nterms, 1);
   k = find(all(modw == mod0, 2));
   ssw(k) = ssx0;
   k = find(all(modwo == mod0, 2));
   sswo(k) = ssx0;
end
clear mod0


% Compute the sum of squares attributed to each term
ssterm = ssw - sswo;
ssterm(dfterm==0) = 0;    % make this exact


% ========================================================
% !!!The ones for dfterm are the same for all voxles: is there a way to move them out of this routine?
% Also the switches should stay in PreProc.m instead of being here!
%=========================================================

switch NF
   case 1,
      fstat = repmat(0, size(ssterm));  % only ss term: SSTR -- treatment sum of squares

   case 2,
      fstat = repmat(0, size(ssterm));

   case 3,
      switch dsgn
	      case 1,
	         fstat = repmat(0, [1 N_Brik]);
			
	      case 2,  % 7 terms: 1 (A); 2 (B); 3 (C); 4 (AB), 5 (AC), 6 BC, 7 ABC				
	         fstat = repmat(0, [1 N_Brik]);   % same as fstat = repmat(0, [1 N_Brik]);
				
	      case {3, 4,} % BXC(A): only 5 terms - 1 (A); 2 (B); 3 C(A); 4 (AB), 5 BC(A)
	         							
	         ssterm(3) = ssterm(3) + ssterm(5);  % SSC(A) = SSC + SSAC
	         dfterm(3) = dfterm(3) + dfterm(5);        %  DF C(A) = DF C + DF AC
		
	         ssterm(6) = ssterm(6) + ssterm(7);   % SSBC(A) = SSBC + SSABC
	         dfterm(6) = dfterm(6) + dfterm(7);
				
				fstat = repmat(0, [1 N_Brik]); 	
								
      end  % switch dsgn	

   case 4,

      switch dsgn
         % Order of the 15 terms: 1 (A); 2 (B); 3 (C); 4 (D); 5 (AB); 6 (AC); 7 (AD); 8 (BC); 9 (BD); 10 (CD)
	      %        11 (ABC); 12 (ABD); 13 (ACD); 14 (BCD); 15 (ABCD)
	      case {1, 2,}
	         fstat = repmat(0, [1 N_Brik]);
	      case {3, 4,}	% only 11 terms in nesting case without AD, ABD, ACD, and ABCD: 1 (A); 2 (B); 3 (C); 4 (D); 5 (AB); 6 (AC);
			% 7 (BC); 8 (BD); 9 (CD); 10 (ABC); 11 (BCD); 	      			
            ssterm(4) = ssterm(4) + ssterm(7);       %SSD(A) = SSD + SSAD
	         dfterm(4) = dfterm(4) + dfterm(7);   %DF(D(A)) = DF(D) + DF(AD)

            ssterm(9) = ssterm(9) + ssterm(12);      %SSBD(A) = SSB + SSABD
            dfterm(9) = dfterm(9) + dfterm(12);   %DF(BD(A)) = DF(BD) + DF(ABD)
	
            ssterm(10) = ssterm(10) + ssterm(13);    %SSCD(A) = SSC + SSACD
            dfterm(10) = dfterm(10) + dfterm(13);   %DF(CD(A)) = DF(CD) + DF(ACD)
	
            ssterm(14) = ssterm(14) + ssterm(15);    %SSBCD(A) = SSBCD + SSABCD
            dfterm(14) = dfterm(14) + dfterm(15);   %DF(BCD(A)) = DF(BCD)+ DF(ABCD)
	
            fstat = repmat(0, [1 N_Brik]);
				
         case 5,
			% There are only 9 terms: 1 (A); 2 (B); 3 (C); 4 (D); 5 (AB); 6 (AC); 7 (BC); 8 (CD); 9 (ABC)
						
            ssterm(4) = ssterm(4) + ssterm(7) + ssterm(9) + ssterm(12);       %SSD(AB) = SSD + SSAD + SSBD + SSABD
            dfterm(4) = dfterm(4) + dfterm(7) + dfterm(9) + dfterm(12);   %DF(D(AB)) = (DF(A)+ DF(B) + DF(AB) + 1)DF(D)	
			
            ssterm(10) = ssterm(10) + ssterm(13) + ssterm(14) + ssterm(15);       %SSCD(AB) = SSCD + SSACD + SSBCD + SSABCD
            dfterm(10) = dfterm(10) + dfterm(13) + dfterm(14) + dfterm(15);   %DF(CD(AB)) = (DF(A)+ DF(B) + DF(AB) + 1)DF(D)	
			
            fstat = repmat(0, [1 N_Brik]);						
      end
		
	case 5,
	
	   switch dsgn   %Order of the 2^5 - 1 = 31 terms: 1 (A); 2 (B); 3 (C); 4 (D); 5 (E); 6 (AB); 7 (AC); 8 (AD); 9 (AE);
		                % 10 (BC); 11 (BD); 12 (BE); 13 (CD); 14 (CE); 15 (DE) 16 (ABC); 17 (ABD); 18 (ABE); 19 (ACD); 20 (ACE); 21 (ADE);
							 % 22 (BCD); 23 (BCE); 24 (BDE); 25 (CDE); 26 (ABCD); 27 (ABCE); 28 (ABDE) 29 (ACDE); 30 (BCDE); 31 (ABCDE)
		
		case {1, 2,} 		
		   fstat = repmat(0, [1 N_Brik]);   % pure cross design
		
		case {3,4},   % only 23 terms
		
			ssterm(5) = ssterm(5) + ssterm(9);   %SSE(A) = SSE + SSAE
	      dfterm(5) = dfterm(5) + dfterm(9);   %DF(E(A)) = DF(E) + DF(AE)
			
			ssterm(12) = ssterm(12) + ssterm(18);   %SSBE(A) = SSBE + SSABE
	      dfterm(12) = dfterm(12) + dfterm(18);   %DF(BE(A)) = DF(BE) + DF(ABE)
			
			ssterm(14) = ssterm(14) + ssterm(20);   %SSCE(A) = SSCE + SSACE
	      dfterm(14) = dfterm(14) + dfterm(20);   %DF(CE(A)) = DF(CE) + DF(ACE)
			
			ssterm(15) = ssterm(15) + ssterm(21);   %SSDE(A) = SSDE + SSADE
	      dfterm(15) = dfterm(15) + dfterm(21);   %DF(DE(A)) = DF(DE) + DF(ADE)
			
			ssterm(23) = ssterm(23) + ssterm(27);   %SSBCE(A) = SSBCE + SSABCE
	      dfterm(23) = dfterm(23) + dfterm(27);   %DF(BCE(A)) = DF(BCE) + DF(ABCE)
			
			ssterm(24) = ssterm(24) + ssterm(28);   %SSBDE(A) = SSBDE + SSABDE
	      dfterm(24) = dfterm(24) + dfterm(28);   %DF(BDE(A)) = DF(BDE) + DF(ABDE)
			
			ssterm(25) = ssterm(25) + ssterm(29);   %SSCDE(A) = SSCDE + SSACDE
	      dfterm(25) = dfterm(25) + dfterm(29);   %DF(CDE(A)) = DF(CDE) + DF(ACDE)
			
			ssterm(30) = ssterm(30) + ssterm(31);   %SSBCDE(A) = SSBCDE + SSABCDE
	      dfterm(30) = dfterm(30) + dfterm(31);   %DF(BCDE(A)) = DF(BCDE) + DF(ABCDE)
			
			fstat = repmat(0, [1 N_Brik]);
		end % switch dsgn	
			
end  % switch NF

% Compute the mean square for each term
msterm = ssterm ./ max(1, dfterm');
mse = sse * (dfe>0) / max(1, dfe);
intensity = sqrt(msterm);   %Intensity for each term

% Equated computed and expected mean squares, then solve for
%           variance component estimates


switch NF
   case 1,
      dfterm_new = dfterm; tnames_new = tnames;	msterm_new = msterm; intensity_new = intensity;
      msdenom = repmat(mse, size(msterm)); dfdenom = repmat(dfe, size(msterm));

   case 2,  %totally 3 terms: A, B, abd AXB.
      dfterm_new = dfterm; tnames_new = tnames;	msterm_new = msterm; intensity_new = intensity;
		
      switch dsgn  % Allocate denominator and its d.f. for each F ratio
      case 1,
         msdenom = repmat(mse, size(msterm));
	      dfdenom = repmat(dfe, size(msterm));
      case 2,	
   	   msdenom = [msterm(3), mse, mse];
	      dfdenom = [dfterm(3), dfe, dfe];	
      case 3,
	      msdenom = [msterm(3), msterm(3), mse];
	      dfdenom = [dfterm(3), dfterm(3), dfe];				
   end   % Close swtich dsgn
		
   case 3,
   switch dsgn
      case 1,
	      dfterm_new = dfterm; tnames_new = tnames;	msterm_new = msterm; intensity_new = intensity;
	      msdenom = repmat(mse, size(msterm)); dfdenom = repmat(dfe, size(msterm));
			
      case 2,   % 7 terms: 1 (A); 2 (B); 3 (C); 4 (AB), 5 (AC), 6 BC, 7 ABC
	      dfterm_new = dfterm; tnames_new = tnames;	msterm_new = msterm; intensity_new = intensity;
	      msdenom = [msterm(5), msterm(6), mse, msterm(7), mse, mse, mse];
	      dfdenom = [dfterm(5), dfterm(6), dfe, dfterm(7), dfe, dfe, dfe];
			
      case 3,  % 5 terms: 1 A; 2 B; 3 C(A); 4 (AB), 5  BC(A)
         msterm_new = [msterm(1:4), msterm(6)];   % Throw out those four which do not exist for nesting: AC and ABC.
	      intensity_new = [intensity(1:4), intensity(6)];    % Throw out those four which do not exist for nesting.
	      dfterm_new = [dfterm(1:4)', dfterm(6)'];
	      tnames_new = [tnames(1:4); tnames(6)];     % Only preserve those valid for nesting. Semicolon for coloumn catenation	
         msdenom = [msterm(3), msterm(6), mse, msterm(6), mse,0,0];  % pad 2 extra 0's for potential complaints
	      dfdenom = [dfterm(3), dfterm(6), dfe, dfterm(6), dfe,0,0];
				
      case 4,  % 5 terms: 1 A; 2 B; 3 C(A); 4 (AB), 5  BC(A)
         msterm_new = [msterm(1:4), msterm(6)];   % Throw out those four which do not exist for nesting: AC and ABC.
	      intensity_new = [intensity(1:4), intensity(6)];    % Throw out those four which do not exist for nesting.
	      dfterm_new = [dfterm(1:4)', dfterm(6)'];
	      tnames_new = [tnames(1:4); tnames(6)];     % Only preserve those valid for nesting. Semicolon for coloumn catenation	
         msdenom = [msterm(4), mse, msterm(6), mse, mse,0,0];  % 2 extra 0's are for error-prone problem down below for contrasts
	      dfdenom = [dfterm(4), dfe, dfterm(6), dfe, dfe,0,0];	
   end			

   case  4,
      switch dsgn     %Order of the 2^4 - 1 = 15 terms: 1 (A); 2 (B); 3 (C); 4 (D); 5 (AB); 6 (AC); 7 (AD); 8 (BC); 9 (BD); 10 (CD)
	                   % 11 (ABC); 12 (ABD); 13 (ACD); 14 (BCD); 15 (ABCD)
	   case 1,
	      dfterm_new = dfterm; tnames_new = tnames;	msterm_new = msterm; intensity_new = intensity;
	      msdenom = repmat(mse, size(msterm)); dfdenom = repmat(dfe, size(msterm));
	   case 2,
	      dfterm_new = dfterm; tnames_new = tnames;	msterm_new = msterm; intensity_new = intensity;				
	      msdenom = [msterm(7), msterm(9), msterm(10), mse, msterm(12), msterm(13), mse, msterm(14), mse, mse, msterm(15), mse, mse, mse, mse];
	      dfdenom = [dfterm(7), dfterm(9), dfterm(10), dfe, dfterm(12), dfterm(13), dfe, dfterm(14), dfe, dfe, dfterm(15), dfe, dfe, dfe, dfe];
	   case 3,   % only 11 terms in nesting case without AD, ABD, ACD, and ABCD: 1 (A); 2 (B); 3 (C); 4 (D); 5 (AB); 6 (AC);
			% 7 (BC); 8 (BD); 9 (CD); 10 (ABC); 11 (BCD);
         msterm_new = [msterm(1:6), msterm(8:11), msterm(14)];   % Throw out those four which do not exist for nesting: AD, ABD, ACD, and ABCD.
	      intensity_new = [intensity(1:6), intensity(8:11), intensity(14)];    % Throw out those four which do not exist for nesting.
	      dfterm_new = [dfterm(1:6)', dfterm(8:11)', dfterm(14)'];
	      tnames_new = [tnames(1:6); tnames(8:11); tnames(14)];     % Only preserve those valid for nesting. Semicolon for coloumn catenation	
         msdenom = [msterm(4), msterm(9), msterm(10), mse, msterm(9), msterm(10), msterm(14), mse, mse, msterm(14), mse];
	      dfdenom = [dfterm(4), dfterm(9), dfterm(10), dfe, dfterm(9), dfterm(10), dfterm(14), dfe, dfe, dfterm(14), dfe];
	   case 4,
	      msterm_new = [msterm(1:6), msterm(8:11), msterm(14)];   % Throw out those four which do not exist for nesting: AD, ABD, ACD, and ABCD.
	      intensity_new = [intensity(1:6), intensity(8:11), intensity(14)];    % Throw out those four which do not exist for nesting.
	      dfterm_new = [dfterm(1:6)', dfterm(8:11)', dfterm(14)'];
	      tnames_new = [tnames(1:6); tnames(8:11); tnames(14)];     % Only preserve those valid for nesting. Semicolon for coloumn catenation	
	      msdenom = [msterm(6), msterm(8), mse, msterm(10), msterm(11), mse, mse, msterm(14), mse, mse, mse];  % denominator MS
   	   dfdenom = [dfterm(6), dfterm(8), dfe, dfterm(10), dfterm(11), dfe, dfe, dfterm(14), dfe, dfe, dfe];  % denominator DF
	   case 5, % only 9 terms in nesting case without AD, BD, ABD, ACD, BCD and ABCD: 1 (A); 2 (B); 3 (C); 4 (D); 5 (AB); 6 (AC); 7 (BC); 8 (CD); 9 (ABC)
	      msterm_new = [msterm(1:6), msterm(8), msterm(10:11)];   % Throw out those four which do not exist for nesting: AD, BD, ABD, ACD, BCD and ABCD.
	      intensity_new = [intensity(1:6), intensity(8), intensity(10:11)];    % Throw out those 6 which do not exist for nesting.
	      dfterm_new = [dfterm(1:6)', dfterm(8), dfterm(10:11)'];
	      tnames_new = [tnames(1:6); tnames(8); tnames(10:11)];     % Only preserve those valid for nesting. Semicolon for coloumn catenation	
	      msdenom = [msterm(4), msterm(4), msterm(10), mse, msterm(4), msterm(10), msterm(10), mse, msterm(10),0,0,0,0,0,0];
	    % denominator MS: 6 extra 0's are for error-prone problem down below for contrasts
   	   dfdenom = [dfterm(4), dfterm(4), dfterm(10), dfe, dfterm(4), dfterm(10), dfterm(10), dfe, dfterm(10),0,0,0,0,0,0];  % denominator DF				
   end	% Close switch dsgn
	
	case  5,
      switch dsgn     %Order of the 2^5 - 1 = 31 terms: 1 (A); 2 (B); 3 (C); 4 (D); 5 (E); 6 (AB); 7 (AC); 8 (AD); 9 (AE);
		                % 10 (BC); 11 (BD); 12 (BE); 13 (CD); 14 (CE); 15 (DE) 16 (ABC); 17 (ABD); 18 (ABE); 19 (ACD); 20 (ACE); 21 (ADE);
							 % 22 (BCD); 23 (BCE); 24 (BDE); 25 (CDE); 26 (ABCD); 27 (ABCE); 28 (ABDE) 29 (ACDE); 30 (BCDE); 31 (ABCDE)
	   case 1,  % 31 terms
	      dfterm_new = dfterm; tnames_new = tnames;	msterm_new = msterm; intensity_new = intensity;
	      msdenom = repmat(mse, size(msterm)); dfdenom = repmat(dfe, size(msterm));
	   case 2,   % 31 terms
	      dfterm_new = dfterm; tnames_new = tnames;	msterm_new = msterm; intensity_new = intensity;				
	      msdenom = [msterm(9), msterm(12), msterm(14), msterm(15), mse, msterm(18), msterm(20), msterm(21), mse, msterm(23), ...
			   msterm(24), mse, msterm(25), mse, mse, msterm(27), msterm(28), mse, msterm(29), mse, mse, msterm(30), mse, mse, mse, ...
				msterm(31), mse, mse, mse, mse, mse];
	      dfdenom = [dfterm(9), dfterm(12), dfterm(14), dfterm(15), dfe, dfterm(18), dfterm(20), dfterm(21), dfe, dfterm(23), ...
			   dfterm(24), dfe, dfterm(25), dfe, dfe, dfterm(27), dfterm(28), dfe, dfterm(29), dfe, dfe, dfterm(30), dfe, dfe, dfe, ...
				dfterm(31), dfe, dfe, dfe, dfe, dfe];
		case 3,  % 23 terms
		   msterm_new = [msterm(1:8), msterm(10:17), msterm(19), msterm(22:26), msterm(30)];
			intensity_new = [intensity(1:8), intensity(10:17), intensity(19), intensity(22:26), intensity(30)];
			dfterm_new = [dfterm(1:8)', dfterm(10:17)', dfterm(19), dfterm(22:26)', dfterm(30)];
			tnames_new = [tnames(1:8); tnames(10:17); tnames(19); tnames(22:26); tnames(30)];
			msdenom = [msterm(5), msterm(12), msterm(14), msterm(15), mse, msterm(12), msterm(14), msterm(15), msterm(23), msterm(24), ...
			   mse, msterm(25), mse, mse, msterm(23), msterm(24), msterm(25), msterm(30), mse, mse, mse, msterm(30), mse];
			dfdenom = [dfterm(5), dfterm(5), dfterm(5), dfterm(5), dfe, dfterm(12), dfterm(14), dfterm(15), dfterm(23), dfterm(24), ...
			   dfe, dfterm(25), dfe, dfe, dfterm(23), dfterm(24), dfterm(25), dfterm(30),dfe, dfe, dfe, dfterm(30), dfe];  % denominator DF
		
		case 4,  % 23 terms
		   msterm_new = [msterm(1:8), msterm(10:17), msterm(19), msterm(22:26), msterm(30)];
			intensity_new = [intensity(1:8), intensity(10:17), intensity(19), intensity(22:26), intensity(30)];
			dfterm_new = [dfterm(1:8)', dfterm(10:17)', dfterm(19), dfterm(22:26)', dfterm(30)];
			tnames_new = [tnames(1:8); tnames(10:17); tnames(19); tnames(22:26); tnames(30)];
			msdenom = [msterm(8), msterm(11), msterm(13), mse, msterm(15), msterm(17), msterm(19), mse, msterm(22), mse, ...
			   msterm(24), mse, msterm(25), mse, msterm(26), mse, mse, mse, msterm(30), mse, mse, mse, mse];
			dfdenom = [dfterm(8), dfterm(11), dfterm(13), dfe, dfterm(15), dfterm(17), dfterm(19), dfe, dfterm(23), dfe, ...
			   dfterm(24), dfe, dfterm(25), dfe, dfterm(26), dfe, dfe, dfe, dfterm(30), dfe, dfe, dfe, dfe];  % denominator DF		
		
		end
	
end  % Close swtich NF


% Compute an F statistic for each term

t = (msdenom>0);  %Only calculate it when denominator is greater than 0
fstat(t) = msterm_new(t) ./ msdenom(t);

%===================================================
% Contrast tests
if (Contr.do == 0),
   LC = 0;  % assign this so that output argument LC would not be empty
else


if (NF == 1),
% Get t values for contrast tests within each factor (1st order)

% dmat(:, num_col0+1:num_col1+num_col0) is the matrix for 1st order contrasts
   if (Contr.ord1.tot > 0),
      for (i = 1:1:Contr.ord1.tot),
         LC.t1(i).t = 0;  % initializtion in case it is assigned later on due to denominator of 0			
         LC.t1(i).value = Contr.ord1.cnt(i).vec * y0;   % intensity for this 1st order contrast
         tmp = msdenom(Contr.ord1.cnt(i).idx1)*Contr.ord1.cnt(i).scalar;
         if (tmp > 0), LC.t1(i).t = LC.t1(i).value/sqrt(tmp); end
      end
   end
end %if (NF == 1)

if (NF == 2),
% Get t values for contrast tests within each factor (1st order)

% dmat(:, num_col0+1:num_col1+num_col0) is the matrix for 1st order contrasts
   if (Contr.ord1.tot > 0),
      for (i = 1:1:Contr.ord1.tot),
         LC.t1(i).t = 0;  % initializtion in case it is assigned later on due to denominator of 0			
	      LC.t1(i).value = Contr.ord1.cnt(i).vec * y0;   % intensity for this 1st order contrast
	      tmp = msdenom(Contr.ord1.cnt(i).idx1)*Contr.ord1.cnt(i).scalar;
	      if (tmp > 0), LC.t1(i).t = LC.t1(i).value/sqrt(tmp); end
      end
   end
	
	if (Contr.ord2.tot > 0),    % 7 terms: 1 (A); 2 (B); 3 (AB)
      for (i = 1:1:Contr.ord2.tot),
         LC.t2(i).t = 0;  % initializtion in case it is assigned later on due to denominator of 0
	      LC.t2(i).value = Contr.ord2.cnt(i).vec * y0;   % intensity for this 2nd order contrast
	      tmp = msdenom(Contr.ord2.cnt(i).idx2)*Contr.ord2.cnt(i).scalar;
	      if (tmp > 0), LC.t2(i).t = LC.t2(i).value/sqrt(tmp); end
      end %for (i = 1:1:Contr.ord2.tot),
   end	%if (Contr.ord2.tot > 0),
	
end %if (NF == 2)

if (NF == 3),
% Get t values for contrast tests within each factor (1st order)

% dmat(:, num_col0+1:num_col1+num_col0) is the matrix for 1st order contrasts
   if (Contr.ord1.tot > 0),
      for (i = 1:1:Contr.ord1.tot),
         LC.t1(i).t = 0;  % initializtion in case it is assigned later on due to denominator of 0			
	      LC.t1(i).value = Contr.ord1.cnt(i).vec * y0;   % intensity for this 1st order contrast
	      tmp = msdenom(Contr.ord1.cnt(i).idx1)*Contr.ord1.cnt(i).scalar;
	      if (tmp > 0), LC.t1(i).t = LC.t1(i).value/sqrt(tmp); end
      end
   end	

% Get t values for contrast tests (2nd order)

% dmat(:, num_col0+num_col1+1:num_col1+num_col0+num_col2) is the matrix for 2nd order contrasts
   if (Contr.ord2.tot > 0),    % 7 terms: 1 (A); 2 (B); 3 (C); 4 (AB), 5 (AC), 6 BC, 7 ABC
      for (i = 1:1:Contr.ord2.tot),
         LC.t2(i).t = 0;  % initializtion in case it is assigned later on due to denominator of 0
	      switch Contr.ord2.cnt(i).idx1
	         case 1,
	            switch Contr.ord2.cnt(i).idx2
		            case 2, what = msdenom(4);   % for MSAB
						case 3, what = msdenom(5) * (dsgn == 1 | dsgn == 2) + msdenom(3) * (dsgn == 4);   % MSAC
						% For design type 4 -- BXC(A): the denominator for this contrast C(A) is MSBC(A), the one for main effect of C(A)!
	            end	
	         case 2,
	         if (Contr.ord2.cnt(i).idx2 == 3),
		         what = msdenom(6) * (dsgn == 1 | dsgn == 2) + msdenom(5)* (dsgn == 3 | dsgn == 4);  % MSBC
	         else fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause; end		
	      end %switch Contr.ord2.cnt(i).idx1
	      LC.t2(i).value = Contr.ord2.cnt(i).vec * y0;   % intensity for this 2nd order contrast
	      tmp = what*Contr.ord2.cnt(i).scalar;
	      if (tmp > 0), LC.t2(i).t = LC.t2(i).value/sqrt(tmp); end
      end %for (i = 1:1:Contr.ord2.tot),
   end	%if (Contr.ord2.tot > 0),

% dmat(:, num_col0+num_col1+1:num_col1+num_col0+num_col2) is the matrix for 2nd order contrasts
   if (Contr.ord3.tot > 0),
      for (i = 1:1:Contr.ord3.tot),
         LC.t3(i).t = 0;  % initializtion in case it is assigned later on due to denominator of 0
	 switch Contr.ord3.cnt(i).idx1
	    case 1,
	       switch Contr.ord3.cnt(i).idx2
		  case 2,
		     if (Contr.ord3.cnt(i).idx3 == 3),
			     what = msdenom(7) * (dsgn == 1 | dsgn == 2);   % MSABC
		     else fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause; end
		  case 3, fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;	
	       end	
	    case 2, fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;
	    case 3,   % Less likely occur
	       fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;
	 end
	 LC.t3(i).value = Contr.ord3.cnt(i).vec * y0;   % intensity for this 2nd order contrast
	 tmp = what*Contr.ord3.cnt(i).scalar;
	 if (tmp > 0), LC.t3(i).t = LC.t3(i).value/sqrt(tmp); end
      end %for (i = 1:1:Contr.ord3.tot),
   end
end % If (NF == 3)

if (NF == 4),
% Get t values for contrast tests within each factor (1st order)

% dmat(:, num_col0+1:num_col1+num_col0) is the matrix for 1st order contrasts
   if (Contr.ord1.tot > 0),
      for (i = 1:1:Contr.ord1.tot),
         LC.t1(i).t = 0;  % initializtion in case it is assigned later on due to denominator of 0			
	      LC.t1(i).value = Contr.ord1.cnt(i).vec * y0;   % intensity for this 1st order contrast
	      tmp = msdenom(Contr.ord1.cnt(i).idx1)*Contr.ord1.cnt(i).scalar;
	      if (tmp > 0), LC.t1(i).t = LC.t1(i).value/sqrt(tmp); end
      end
   end	

% Get t values for contrast tests within each factor (2nd order)
% dmat(:, num_col0+num_col1+1:num_col1+num_col0+num_col2) is the matrix for 2nd order contrasts
   if (Contr.ord2.tot > 0),
      for (i = 1:1:Contr.ord2.tot),
         LC.t2(i).t = 0;  % initializtion in case it is assigned later on due to denominator of 0
%		switch dsgn
         switch Contr.ord2.cnt(i).idx1
	         case 1,
	         switch Contr.ord2.cnt(i).idx2
		         case 2, what = msdenom(5);   % MSAB
		         case 3, what = msdenom(6);   % MSAC
		         case 4, what = msdenom(7) * (dsgn == 1 | dsgn == 2);   % MSAD
	         end	
	         case 2,
	         switch Contr.ord2.cnt(i).idx2
	            case 3, what = msdenom(8) * (dsgn == 1 | dsgn == 2) + msdenom(7)* (dsgn == 3 | dsgn == 4 | dsgn == 5);  % MSBC
	            case 4, what = msdenom(9) * (dsgn == 1 | dsgn == 2) + msdenom(8)* (dsgn == 3 | dsgn == 4);  % Less likely to occur: MSBD	
	         end		
	         case 3,   % Less likely occur
	         switch Contr.ord2.cnt(i).idx2
		         case 4, what = msdenom(10)* (dsgn == 1 | dsgn == 2) + msdenom(9)* (dsgn == 3 | dsgn == 4) + msdenom(7)*(dsgn == 5);  % Less likely occur: MSCD		
	         end
	      end		
%		end
		% Can i remove the following loop by doing it in a matrix operation fashion?
		%for (j = 1:1:Contr2.cnt(i).NT),  % each term of this contrast
		%   LC2(i).value = LC2(i).value + Contr2.cnt(i).coef(j) * dmat(:, Contr2.cnt(i).code(j).pos)' * y0;
			%scalar = scalar  + Contr2.cnt(i).coef(j) * Contr2.cnt(i).coef(j);
		%end
		
         LC.t2(i).value = Contr.ord2.cnt(i).vec * y0;   % intensity for this 2nd order contrast
	      tmp = what*Contr.ord2.cnt(i).scalar;
	      if (tmp > 0), LC.t2(i).t = LC.t2(i).value/sqrt(tmp); end
      end
   end	

   % dmat(:, num_col0+num_col1+1:num_col1+num_col0+num_col2) is the matrix for 2nd order contrasts
   if (Contr.ord3.tot > 0),
      for (i = 1:1:Contr.ord3.tot),
         LC.t3(i).t = 0;  % initializtion in case it is assigned later on due to denominator of 0
%		if (dsgn == 3),
	      switch Contr.ord3.cnt(i).idx1
	
				case 1,
	         switch Contr.ord3.cnt(i).idx2
		
					case 2,
		         switch Contr.ord3.cnt(i).idx3
			         case 3, what = msdenom(11) * (dsgn == 1 | dsgn == 2) + msdenom(10) * (dsgn == 3 | dsgn == 4) + msdenom(9) * (dsgn == 5);   % MSABC
			         case 4, what = msdenom(12) * (dsgn == 1 | dsgn == 2); %  MSABD not exist for (dsgn == 3 | dsgn == 4 | dsgn == 5)
		         end
		
					case 3,
		         if (Contr.ord3.cnt(i).idx3 == 4), what = msdenom(13) * (dsgn == 1 | dsgn == 2);  % MSACD
		         else fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause; end
		
					case 4, fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;	
	         end	
	
			   case 2,
	         switch Contr.ord3.cnt(i).idx2
		         case 3,
		         if (Contr.ord3.cnt(i).idx3 == 4), what = msdenom(14) * (dsgn == 1 | dsgn == 2) + msdenom(11) * (dsgn == 3 | dsgn == 4);  % MSBCD
		         else fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;end
		         case 4,
		         fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;	
	       end				
	    case 3,   % Less likely occur
	       fprintf('\nSomething is wrong in the contrast coding!\n');
	       fprintf(2,'Halted: Ctrl+c to exit'); pause;
	    case 4,
	       fprintf('\nSomething is wrong in the contrast coding!\n');
	       fprintf(2,'Halted: Ctrl+c to exit'); pause;	
	 end
%		end
	 LC.t3(i).value = Contr.ord3.cnt(i).vec * y0;   % intensity for this 3rd order contrast
	 tmp = what*Contr.ord3.cnt(i).scalar;
	 if (tmp > 0), LC.t3(i).t = LC.t3(i).value/sqrt(tmp); end
      end
   end	
end % if (NF == 4)


if (NF == 5),
% Get t values for contrast tests within each factor (1st order)

% dmat(:, num_col0+1:num_col1+num_col0) is the matrix for 1st order contrasts
   if (Contr.ord1.tot > 0),
      for (i = 1:1:Contr.ord1.tot),
         LC.t1(i).t = 0;  % initializtion in case it is assigned later on due to denominator of 0			
	      LC.t1(i).value = Contr.ord1.cnt(i).vec * y0;   % intensity for this 1st order contrast
	      tmp = msdenom(Contr.ord1.cnt(i).idx1)*Contr.ord1.cnt(i).scalar;
	      if (tmp > 0), LC.t1(i).t = LC.t1(i).value/sqrt(tmp); end
      end
   end	

% Get t values for contrast tests within each factor (2nd order)
% dmat(:, num_col0+num_col1+1:num_col1+num_col0+num_col2) is the matrix for 2nd order contrasts
%Order of the 2^5 - 1 = 31 terms: 1 (A); 2 (B); 3 (C); 4 (D); 5 (E); 6 (AB); 7 (AC); 8 (AD); 9 (AE);
		                % 10 (BC); 11 (BD); 12 (BE); 13 (CD); 14 (CE); 15 (DE) 16 (ABC); 17 (ABD); 18 (ABE); 19 (ACD); 20 (ACE); 21 (ADE);
							 % 22 (BCD); 23 (BCE); 24 (BDE); 25 (CDE); 26 (ABCD); 27 (ABCE); 28 (ABDE) 29 (ACDE); 30 (BCDE); 31 (ABCDE)

   if (Contr.ord2.tot > 0),
      for (i = 1:1:Contr.ord2.tot),
         LC.t2(i).t = 0;  % initializtion in case it is assigned later on due to denominator of 0
%		switch dsgn
         switch Contr.ord2.cnt(i).idx1
	         case 1,
	         switch Contr.ord2.cnt(i).idx2
		         case 2, what = msdenom(6);   % denominator for MSAB
		         case 3, what = msdenom(7);   % MSAC
		         case 4, what = msdenom(8);   % MSAD
					case 5, what = msdenom(9) * (dsgn == 1 | dsgn == 2);   % MSAE
	         end	
	         case 2,
	         switch Contr.ord2.cnt(i).idx2
	            case 3, what = msdenom(10) * (dsgn == 1 | dsgn == 2) + msdenom(9) * (dsgn == 3 | dsgn == 4);  % MSBC
	            case 4, what = msdenom(11) * (dsgn == 1 | dsgn == 2) + msdenom(10) * (dsgn == 3 | dsgn == 4);  % Less likely to occur: MSBD	
					case 5, what = msdenom(12) * (dsgn == 1 | dsgn == 2) + msdenom(11) * (dsgn == 3 | dsgn == 4);  % Less likely to occur: MSBE
	         end		
	         case 3,   % Less likely occur
	         switch Contr.ord2.cnt(i).idx2
		         case 4, what = msdenom(13)* (dsgn == 1 | dsgn == 2) + msdenom(12) * (dsgn == 3 | dsgn == 4);  % Less likely occur: MSCD
					case 5, what = msdenom(14)* (dsgn == 1 | dsgn == 2) + msdenom(13) * (dsgn == 3 | dsgn == 4);  % Less likely occur: MSCE	
	         end
				case 4,   % Less likely occur
	         switch Contr.ord2.cnt(i).idx2
					case 5, what = msdenom(15)* (dsgn == 1 | dsgn == 2) + msdenom(14) * (dsgn == 3 | dsgn == 4);  % Less likely occur: MSDE	
	         end
				
	      end		
%		end
		% Can i remove the following loop by doing it in a matrix operation fashion?
		%for (j = 1:1:Contr2.cnt(i).NT),  % each term of this contrast
		%   LC2(i).value = LC2(i).value + Contr2.cnt(i).coef(j) * dmat(:, Contr2.cnt(i).code(j).pos)' * y0;
			%scalar = scalar  + Contr2.cnt(i).coef(j) * Contr2.cnt(i).coef(j);
		%end
		
         LC.t2(i).value = Contr.ord2.cnt(i).vec * y0;   % intensity for this 2nd order contrast
	      tmp = what*Contr.ord2.cnt(i).scalar;
	      if (tmp > 0), LC.t2(i).t = LC.t2(i).value/sqrt(tmp); end
      end
   end	

   % dmat(:, num_col0+num_col1+1:num_col1+num_col0+num_col2) is the matrix for 2nd order contrasts
   if (Contr.ord3.tot > 0),
      for (i = 1:1:Contr.ord3.tot),
         LC.t3(i).t = 0;  % initializtion in case it is assigned later on due to denominator of 0
%		if (dsgn == 3),
	      switch Contr.ord3.cnt(i).idx1
	
				case 1,
	         switch Contr.ord3.cnt(i).idx2
		
					case 2,
		         switch Contr.ord3.cnt(i).idx3
			         case 3, what = msdenom(16) * (dsgn == 1 | dsgn == 2) + msdenom(15) * (dsgn == 3 | dsgn == 4);   % MSABC
			         case 4, what = msdenom(17) * (dsgn == 1 | dsgn == 2) + msdenom(16) * (dsgn == 3 | dsgn == 4); %  MSABD
						case 5, what = msdenom(18) * (dsgn == 1 | dsgn == 2); %  MSABE
		         end   % switch Contr.ord3.cnt(i).idx3
		
					case 3,
		         switch Contr.ord3.cnt(i).idx3
					   case 4, what = msdenom(19) * (dsgn == 1 | dsgn == 2) + msdenom(17) * (dsgn == 3 | dsgn == 4); %  MSACD
						case 5, what = msdenom(20) * (dsgn == 1 | dsgn == 2); %  MSACE
		         end  % switch Contr.ord3.cnt(i).idx3
		
					case 4,
					if (Contr.ord3.cnt(i).idx3 == 5), what = msdenom(21) * (dsgn == 1 | dsgn == 2); % ADE
					else fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;	
	            end	
	         end % switch Contr.ord3.cnt(i).idx2
			
				case 2,
	         switch Contr.ord3.cnt(i).idx2
		         case 3,
		         switch Contr.ord3.cnt(i).idx3
					   case 4, what = msdenom(22) * (dsgn == 1 | dsgn == 2) + msdenom(18) * (dsgn == 3 | dsgn == 4); %  MSBCD
						case 5, what = msdenom(23) * (dsgn == 1 | dsgn == 2) + msdenom(19) * (dsgn == 3 | dsgn == 4); %  MSBCE
					end % switch Contr.ord3.cnt(i).idx3	
		         case 4,
					if (Contr.ord3.cnt(i).idx3 == 5),
					if (dsgn == 1 | dsgn == 2),
					   what = msdenom(24);
					elseif (dsgn == 3)
					   what = msdenom(20);	
					end	
%					what = msdenom(24) * (dsgn == 1 | dsgn == 2) + msdenom(20) * (dsgn == 3); % BDE
		         else fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;	end
	         end		% switch Contr.ord3.cnt(i).idx2		
	         case 3,
	            if (Contr.ord3.cnt(i).idx2 == 4 & Contr.ord3.cnt(i).idx3 == 5),
			         if (dsgn == 1 | dsgn == 2),
						   what = msdenom(25);
						elseif (dsgn == 3),
						   what = msdenom(21);
						end		
%						what = msdenom(25) * (dsgn == 1 | dsgn == 2) + msdenom(21) * (dsgn == 3); %  MSCDE
			      else 	
			         fprintf('\nSomething is wrong in the contrast coding!\n');
	               fprintf(2,'Halted: Ctrl+c to exit'); pause;
			      end
	         case 4,
	            fprintf('\nSomething is wrong in the contrast coding!\n');
	            fprintf(2,'Halted: Ctrl+c to exit'); pause;	
	      end  % switch Contr.ord3.cnt(i).idx1
%		end
	      LC.t3(i).value = Contr.ord3.cnt(i).vec * y0;   % intensity for this 3rd order contrast
	      tmp = what*Contr.ord3.cnt(i).scalar;
	      if (tmp > 0), LC.t3(i).t = LC.t3(i).value/sqrt(tmp); end
      end % for (i = 1:1:Contr.ord3.tot)
   end	% if (Contr.ord3.tot > 0)

   if (Contr.ord4.tot > 0),
      for (i = 1:1:Contr.ord4.tot),
         LC.t4(i).t = 0;  % initializtion in case it is assigned later on due to denominator of 0
	      switch Contr.ord4.cnt(i).idx1
	
				case 1,
	         switch Contr.ord4.cnt(i).idx2		
					case 2,
		         switch Contr.ord4.cnt(i).idx3
			         case 3,
						switch Contr.ord4.cnt(i).idx4
						   case 4,
							if (dsgn == 1 | dsgn == 2),
							   what = msdenom(26);
							elseif (dsgn == 3),
							   what = msdenom(22);
							end							
%							what = msdenom(26) * (dsgn == 1 | dsgn == 2) + msdenom(22) * (dsgn == 3); %  MSABCD
							case 5,
							if (dsgn == 1 | dsgn == 2),
							   what = msdenom(27);
							end							
%							what = msdenom(27) * (dsgn == 1 | dsgn == 2); %  MSABCE
						end % switch Contr.ord4.cnt(i).idx4	
			         case 4,
						if (Contr.ord4.cnt(i).idx4 == 5),
						if (dsgn == 1 | dsgn == 2),
						   what = msdenom(28);
						end	
%						what = msdenom(28) * (dsgn == 1 | dsgn == 2); % ABDE
						else fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;	
	               end						
						case 5, fprintf('\nSomething is wrong in the contrast coding!\n');
	                  fprintf(2,'Halted: Ctrl+c to exit'); pause;	
	            end  % switch Contr.ord4.cnt(i).idx3
	
					case 3,
		         switch Contr.ord4.cnt(i).idx3
					   case 4,
						if (Contr.ord4.cnt(i).idx4 == 5),
						if (dsgn == 1 | dsgn == 2),
						   what = msdenom(29);
						end						
%						what = msdenom(29) * (dsgn == 1 | dsgn == 2); % ACDE
						else fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;	
	               end
						case 5, fprintf('\nSomething is wrong in the contrast coding!\n');
	                  fprintf(2,'Halted: Ctrl+c to exit'); pause;	
	            end  % switch Contr.ord4.cnt(i).idx3
		
					case 4,
					   fprintf('\nSomething is wrong in the contrast coding!\n');
	               fprintf(2,'Halted: Ctrl+c to exit'); pause;
	         end % switch Contr.ord4.cnt(i).idx2
			
				case 2,
	         switch Contr.ord4.cnt(i).idx2
		         case 3,
		         switch Contr.ord4.cnt(i).idx3
					   case 4,
						if (Contr.ord4.cnt(i).idx4 == 5),
						if (dsgn == 1 | dsgn == 2),
						   what = msdenom(30);
						end
%						what = msdenom(30) * (dsgn == 1 | dsgn == 2) + msdenom(23) * (dsgn == 3); % BCDE
						else fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;	
	               end
						case 5, fprintf('\nSomething is wrong in the contrast coding!\n');
	                  fprintf(2,'Halted: Ctrl+c to exit'); pause;
					end	% switch Contr.ord4.cnt(i).idx3	
		         case 4,
					fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;	
	         end		% switch Contr.ord4.cnt(i).idx2		
	         case {3,4,5,}
	            fprintf('\nSomething is wrong in the contrast coding!\n');
	            fprintf(2,'Halted: Ctrl+c to exit'); pause;			
	      end  % switch Contr.ord4.cnt(i).idx1
	      LC.t4(i).value = Contr.ord4.cnt(i).vec * y0;   % intensity for this 3rd order contrast
	      tmp = what*Contr.ord4.cnt(i).scalar;
	      if (tmp > 0), LC.t4(i).t = LC.t4(i).value/sqrt(tmp); end
      end % for (i = 1:1:Contr.ord4.tot)
   end	% if (Contr.ord4.tot > 0)	
	
end % if (NF == 5)


end %if (Contr.do == 0)

err = 0;
return;
