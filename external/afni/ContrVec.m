function [err, Contr] = ContrVec(order, n, NF, group, dmat, Contr, FL, num_col)
%
%   [err,] = ContrVec.m ()
%
%Purpose:
%
%   Obtain a vector for each contrast so that effect size is calculated in SumsOfSquares.m
%
%Input Parameters:
%
%
%
%Output Parameters:
%   err : 0 No Problem
%       : 1  Problems
%
%
%
%Key Terms:
%
%More Info :
%
%
%
%
%     Author : Gang Chen
%     Date : Wed Dec 1 10:15:29 EST 2004
%     SSCC/NIMH/ National Institutes of Health, Bethesda MD 20892


%Define the function name for easy referencing
FuncName = 'ContrVec.m';

%Debug Flag
DBG = 1;

%initialize return variables
err = 1;

switch order

case 1,   % dmat(:, num_col0+1:num_col0+num_col(1)) is the matrix for 1st order contrasts

   for (i = 1:1:Contr.ord1.tot),
      Contr.ord1.cnt(i).vec = zeros(1, n);  % n = ntot in anovan.m: total number of files, and no. of arrows in dmat
		Contr.ord1.cnt(i).scalar = 0;
      for (j = 1:1:Contr.ord1.cnt(i).NT),   % for each term in a contrast
         first = 0;
         Contr.ord1.cnt(i).code(j).pos = 0;			
			shift = 1;   % Grand mean
	      for (k = 1:1:NF),
            if (Contr.ord1.cnt(i).code(j).str(k) >= 'a' & Contr.ord1.cnt(i).code(j).str(k) <= 'z'),
				   tmpv = 10 + Contr.ord1.cnt(i).code(j).str(k) - 'a';
				else tmpv = Contr.ord1.cnt(i).code(j).str(k) - '0'; end  % tmpv = level #
						
				if (tmpv ~= 0),
	            first = first + 1;  % the first non-zero index backward! For 1st-roder, it is the only nonzero index.
			      if (first == 1),
			         Contr.ord1.cnt(i).idx1 = k;
%				      tmp = FL(k).N_level;
	            elseif (first > 1),
	               fprintf('\nError in contrast coding: more than one non-zero index in 1st order constrasts!\n');
	               fprintf(2,'Halted: Ctrl+c to exit'); pause;
			      end  % if (first == 1)
	         end   % if (tmpv ~= 0)
         end   % for (k = 1:1:NF)
			if (Contr.ord1.cnt(i).idx1 > 1),
			   for (k = 1:1:Contr.ord1.cnt(i).idx1-1),
			      shift = shift + FL(k).N_level;
				end					
			end
			if (Contr.ord1.cnt(i).code(j).str(Contr.ord1.cnt(i).idx1) >= 'a' & Contr.ord1.cnt(i).code(j).str(Contr.ord1.cnt(i).idx1) <= 'z'),
			   tmpv = 10 + Contr.ord1.cnt(i).code(j).str(Contr.ord1.cnt(i).idx1) - 'a';
			else tmpv = Contr.ord1.cnt(i).code(j).str(Contr.ord1.cnt(i).idx1) - '0'; end
					
         Contr.ord1.cnt(i).code(j).pos = shift + tmpv;
			
%        tmp = FL(Contr.ord1.cnt(i).idx1).N_level/n;
%         Contr.ord1.cnt(i).scalar = tmp * Contr.ord1.cnt(i).coef*Contr.ord1.cnt(i).coef';		
%         Contr.ord1.cnt(i).vec = Contr.ord1.cnt(i).vec + Contr.ord1.cnt(i).coef(j) * dmat(:, Contr.ord1.cnt(i).code(j).pos)';
			
			count=0;   % The part should be valid for both balanced and unbalanced designs
%			tmpstr = char (group{Contr.ord1.cnt(i).idx1});
         for (ll = 1:1:n) % Should have a more decent way to do this in a matrix fashion instead of a lenthy loop?
            tmpstr = char (group{Contr.ord1.cnt(i).idx1}(ll));
				if strcmp(tmpstr, FL(Contr.ord1.cnt(i).idx1).level(tmpv).expr)
	            count = count + 1;   % total number of occurences at tmpv-th level
	         end
         end
			Contr.ord1.cnt(i).vec = Contr.ord1.cnt(i).vec + Contr.ord1.cnt(i).coef(j) * dmat(:, Contr.ord1.cnt(i).code(j).pos)'/count;			
			Contr.ord1.cnt(i).scalar = Contr.ord1.cnt(i).scalar + Contr.ord1.cnt(i).coef(j)*Contr.ord1.cnt(i).coef(j)/count;    % For variance coef: j-th term in i-th 1st-order contrast.				
					
      end	% for (j = 1:1:Contr.ord1.cnt(i).NT)
%      Contr.ord1.cnt(i).vec = Contr.ord1.cnt(i).vec*tmp;   % assuming balanced design!!!
   end   %for (i = 1:1:Contr1.tot)

case 2,  % dmat(:, num_col0+num_col(1)+1:num_col0+num_col(1)+num_col(2)) is the matrix for 2nd order contrasts

   % when the 1st factor is collapsed
   %   shift1 = (FL(1).N_level*(FL(2).N_level + FL(3).N_level + FL(4).N_level))*(str2num(Contr2.cnt(i).code(1).str(1)) == 0);
   % when the 2nd factor is collapsed
   %   shift2 = (FL(2).N_level*(FL(3).N_level + FL(4).N_level))*(str2num(Contr2.cnt(i).code(1).str(2)) == 0);
   % when the 3rd factor is collapsed
   %   shift3 = (FL(3).N_level*FL(4).N_level)*(str2num(Contr2.cnt(i).code(1).str(3)) == 0);
   %
   % when the 4th factor is collapsed
   %   shift4 = FL(4).N_level*(str2num(Contr2.cnt(i).code(1).str(4)) == 0);
   %
   %   Contr2.cnt(i).shift = num_col0 + num_col(1) + 1 + shift1 + shift2 + shift3 + shift4;  % position shifting for the two collapsed indices

   shift = 1 + num_col(1);   % Grand mean plus factor means
   for (i = 1:1:Contr.ord2.tot),
	
   % for the two non-zero indices
	Contr.ord2.cnt(i).vec = zeros(1, n);  % n = ntot in anovan.m: total number of files, and no. of arrows in dmat
	Contr.ord2.cnt(i).scalar = 0;
	for (j = 1:1:Contr.ord2.cnt(i).NT),   % For each term in this 2nd-order contrast
      first = 0;
	   Contr.ord2.cnt(i).code(j).pos = 0;
		for (k = 1:1:NF),
         if (Contr.ord2.cnt(i).code(j).str(k) >= 'a' & Contr.ord2.cnt(i).code(j).str(k) <= 'z'),
			   tmpv = 10 + Contr.ord2.cnt(i).code(j).str(k) - 'a';
			else tmpv = Contr.ord2.cnt(i).code(j).str(k) - '0'; end
						
			if (tmpv ~= 0),
            first = first + 1;  % the first non-zero index backward!
	        	if (first == 1),
		         Contr.ord2.cnt(i).idx1 = k;
%       			tmp = FL(k).N_level;
		      elseif (first == 2), Contr.ord2.cnt(i).idx2 = k; end
		   end
      end  %for (k = 1:1:NF)
		sec_fctr = Contr.ord2.cnt(i).idx2;
      switch Contr.ord2.cnt(i).idx1
      case 1,   % 2nd-order contrasts: Ax
%         sec_fctr = Contr.ord2.cnt(i).idx2;
			if (Contr.ord2.cnt(i).code(j).str(1) >= 'a' & Contr.ord2.cnt(i).code(j).str(1) <= 'z'),
			   tmpv1 = 10 + Contr.ord2.cnt(i).code(j).str(1) - 'a';
			else tmpv1 = Contr.ord2.cnt(i).code(j).str(1) - '0'; end
			if (Contr.ord2.cnt(i).code(j).str(sec_fctr) >= 'a' & Contr.ord2.cnt(i).code(j).str(sec_fctr) <= 'z'),
			   tmpv2 = 10 + Contr.ord2.cnt(i).code(j).str(sec_fctr) - 'a';
         else tmpv2 = Contr.ord2.cnt(i).code(j).str(sec_fctr) - '0'; end % tmpv = level #	
			
			Contr.ord2.cnt(i).code(j).pos = shift + (tmpv1 - 1) * FL(sec_fctr).N_level + tmpv2;			
			if (sec_fctr > 2),
			for (ll = 2:1:(sec_fctr-1)),  % shift over AB, AC, ...
			   Contr.ord2.cnt(i).code(j).pos = Contr.ord2.cnt(i).code(j).pos + FL(1).N_level*FL(ll).N_level;
			end % for (ll = 2:1:(sec_fctr-1))
			end % if (sec_fctr > 2)
						
%			switch Contr.ord2.cnt(i).idx2
%         case 2,
%            if (Contr.ord2.cnt(i).code(j).str(1) >= 'a' & Contr.ord2.cnt(i).code(j).str(1) <= 'z'),
%				   tmpv1 = 10 + Contr.ord2.cnt(i).code(j).str(1) - 'a';
%            else tmpv1 = Contr.ord2.cnt(i).code(j).str(1) - '0'; end
%				if (Contr.ord2.cnt(i).code(j).str(2) >= 'a' & Contr.ord2.cnt(i).code(j).str(2) <= 'z'),
%				   tmpv2 = 10 + Contr.ord2.cnt(i).code(j).str(2) - 'a';
%            else tmpv2 = Contr.ord2.cnt(i).code(j).str(2) - '0'; end % tmpv = level #
									
%				Contr.ord2.cnt(i).code(j).pos = shift + (tmpv1 - 1) * FL(2).N_level + tmpv2; % AB
%         case 3,
%            if (Contr.ord2.cnt(i).code(j).str(1) >= 'a' & Contr.ord2.cnt(i).code(j).str(1) <= 'z'),
%				   tmpv1 = 10 + Contr.ord2.cnt(i).code(j).str(1) - 'a';
%            else tmpv1 = Contr.ord2.cnt(i).code(j).str(1) - '0'; end
%				if (Contr.ord2.cnt(i).code(j).str(3) >= 'a' & Contr.ord2.cnt(i).code(j).str(3) <= 'z'),
%				   tmpv2 = 10 + Contr.ord2.cnt(i).code(j).str(3) - 'a';
%            else tmpv2 = Contr.ord2.cnt(i).code(j).str(3) - '0'; end									
									
%				Contr.ord2.cnt(i).code(j).pos = shift + FL(1).N_level*FL(2).N_level+ ...
%               (tmpv1 - 1) * FL(3).N_level + tmpv2; % AC	

%			case 4,   % AD
%			   if (Contr.ord2.cnt(i).code(j).str(1) >= 'a' & Contr.ord2.cnt(i).code(j).str(1) <= 'z'),
%				   tmpv1 = 10 + Contr.ord2.cnt(i).code(j).str(1) - 'a';
%            else tmpv1 = Contr.ord2.cnt(i).code(j).str(1) - '0'; end
%				if (Contr.ord2.cnt(i).code(j).str(4) >= 'a' & Contr.ord2.cnt(i).code(j).str(4) <= 'z'),
%				   tmpv2 = 10 + Contr.ord2.cnt(i).code(j).str(4) - 'a';
%            else tmpv2 = Contr.ord2.cnt(i).code(j).str(4) - '0'; end
				
%				Contr.ord2.cnt(i).code(j).pos = shift + FL(1).N_level*FL(2).N_level+ FL(1).N_level*FL(3).N_level+...
%               (tmpv1 - 1) * FL(4).N_level + tmpv2; % AD
%	      end				
      case 2,  % 2nd-order contrasts Bx
%         sec_fctr = Contr.ord2.cnt(i).idx2;
			if (Contr.ord2.cnt(i).code(j).str(2) >= 'a' & Contr.ord2.cnt(i).code(j).str(2) <= 'z')
			   tmpv1 = 10 + Contr.ord2.cnt(i).code(j).str(2) - 'a';
		   else tmpv1 = Contr.ord2.cnt(i).code(j).str(2) - '0'; end
			if (Contr.ord2.cnt(i).code(j).str(sec_fctr) >= 'a' & Contr.ord2.cnt(i).code(j).str(sec_fctr) <= 'z'),
			   tmpv2 = 10 + Contr.ord2.cnt(i).code(j).str(sec_fctr) - 'a';
		   else tmpv2 = Contr.ord2.cnt(i).code(j).str(sec_fctr) - '0'; end	
			
			Contr.ord2.cnt(i).code(j).pos = shift + (tmpv1 - 1) * FL(sec_fctr).N_level + tmpv2;
			for (ll = 2:1:NF), % shift over all those Ax's
			   Contr.ord2.cnt(i).code(j).pos = Contr.ord2.cnt(i).code(j).pos + FL(1).N_level*FL(ll).N_level;
			end
						
			if (sec_fctr > 3),
			for (ll = 3:1:(sec_fctr-1)),  % shift over Bx's
			   Contr.ord2.cnt(i).code(j).pos = Contr.ord2.cnt(i).code(j).pos + FL(2).N_level*FL(ll).N_level;
			end % for (ll = 3:1:(sec_fctr-1))
			end % if (sec_fctr > 3)				
			
%			if (Contr.ord2.cnt(i).idx2 == 3),
%            if (Contr.ord2.cnt(i).code(j).str(2) >= 'a' & Contr.ord2.cnt(i).code(j).str(2) <= 'z'),
%   			   tmpv1 = 10 + Contr.ord2.cnt(i).code(j).str(2) - 'a';
%		      else tmpv1 = Contr.ord2.cnt(i).code(j).str(2) - '0'; end
%				if (Contr.ord2.cnt(i).code(j).str(3) >= 'a' & Contr.ord2.cnt(i).code(j).str(3) <= 'z'),
%				   tmpv2 = 10 + Contr.ord2.cnt(i).code(j).str(3) - 'a';
%		      else tmpv2 = Contr.ord2.cnt(i).code(j).str(3) - '0'; end	
									
%				if (NF == 3), Contr.ord2.cnt(i).code(j).pos = shift + FL(1).N_level*(FL(2).N_level + FL(3).N_level) + ...
%		          (tmpv1 - 1) * FL(3).N_level + tmpv2; % BC
%				else Contr.ord2.cnt(i).code(j).pos = shift + FL(1).N_level*(FL(2).N_level + FL(3).N_level + FL(4).N_level) + ...
%		          (tmpv1 - 1) * FL(3).N_level + tmpv2; % BC
%				end	 	
%		      tmp = FL(2).N_level*FL(3).N_level/n;
%			else fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause; end

		case 3,  % 2nd-order contrasts Cx
		   if (Contr.ord2.cnt(i).code(j).str(3) >= 'a' & Contr.ord2.cnt(i).code(j).str(3) <= 'z')
			   tmpv1 = 10 + Contr.ord2.cnt(i).code(j).str(3) - 'a';
		   else tmpv1 = Contr.ord2.cnt(i).code(j).str(3) - '0'; end
			if (Contr.ord2.cnt(i).code(j).str(sec_fctr) >= 'a' & Contr.ord2.cnt(i).code(j).str(sec_fctr) <= 'z'),
			   tmpv2 = 10 + Contr.ord2.cnt(i).code(j).str(sec_fctr) - 'a';
		   else tmpv2 = Contr.ord2.cnt(i).code(j).str(sec_fctr) - '0'; end	
			
			Contr.ord2.cnt(i).code(j).pos = shift + (tmpv1 - 1) * FL(sec_fctr).N_level + tmpv2;
			for (ll = 2:1:NF), % shift over all those Ax's
			   Contr.ord2.cnt(i).code(j).pos = Contr.ord2.cnt(i).code(j).pos + FL(1).N_level*FL(ll).N_level;
			end
			
			for (ll = 3:1:NF), % shift over all those Bx's
			   Contr.ord2.cnt(i).code(j).pos = Contr.ord2.cnt(i).code(j).pos + FL(2).N_level*FL(ll).N_level;
			end
						
			if (sec_fctr > 4),
			for (ll = 4:1:(sec_fctr-1)),  % shift over Bx's
			   Contr.ord2.cnt(i).code(j).pos = Contr.ord2.cnt(i).code(j).pos + FL(3).N_level*FL(ll).N_level;
			end % for (ll = 4:1:(sec_fctr-1))
			end % if (sec_fctr > 4)
			
		case 4,  % 2nd-order contrasts Dx
		   if (Contr.ord2.cnt(i).code(j).str(4) >= 'a' & Contr.ord2.cnt(i).code(j).str(4) <= 'z')
			   tmpv1 = 10 + Contr.ord2.cnt(i).code(j).str(4) - 'a';
		   else tmpv1 = Contr.ord2.cnt(i).code(j).str(4) - '0'; end
			if (Contr.ord2.cnt(i).code(j).str(sec_fctr) >= 'a' & Contr.ord2.cnt(i).code(j).str(sec_fctr) <= 'z'),
			   tmpv2 = 10 + Contr.ord2.cnt(i).code(j).str(sec_fctr) - 'a';
		   else tmpv2 = Contr.ord2.cnt(i).code(j).str(sec_fctr) - '0'; end	
			
			Contr.ord2.cnt(i).code(j).pos = shift + (tmpv1 - 1) * FL(sec_fctr).N_level + tmpv2;
			for (ll = 2:1:NF), % shift over all those Ax's
			   Contr.ord2.cnt(i).code(j).pos = Contr.ord2.cnt(i).code(j).pos + FL(1).N_level*FL(ll).N_level;
			end
			
			for (ll = 3:1:NF), % shift over all those Bx's
			   Contr.ord2.cnt(i).code(j).pos = Contr.ord2.cnt(i).code(j).pos + FL(2).N_level*FL(ll).N_level;
			end
			
			for (ll = 4:1:NF), % shift over all those Cx's
			   Contr.ord2.cnt(i).code(j).pos = Contr.ord2.cnt(i).code(j).pos + FL(3).N_level*FL(ll).N_level;
			end
						
			if (sec_fctr > 5),
			for (ll = 5:1:(sec_fctr-1)),  % shift over Bx's
			   Contr.ord2.cnt(i).code(j).pos = Contr.ord2.cnt(i).code(j).pos + FL(4).N_level*FL(ll).N_level;
			end % for (ll = 5:1:(sec_fctr-1))
			end % if (sec_fctr > 5)
				
			
%			fprintf('\nSomething is wrong in the contrast coding!\n');
%		   fprintf(2,'Halted: Ctrl+c to exit'); pause;			
		end  %switch Contr.ord2.cnt(i).idx1
		
		count = 0;
		for (ll = 1:1:n)
			tmpstr1 = char (group{Contr.ord2.cnt(i).idx1}(ll));
			tmpstr2 = char (group{Contr.ord2.cnt(i).idx2}(ll));			
			if strcmp(tmpstr1, FL(Contr.ord2.cnt(i).idx1).level(tmpv1).expr) & strcmp(tmpstr2, FL(Contr.ord2.cnt(i).idx2).level(tmpv2).expr)
			   count = count + 1;   % total number of occurences at tmpv-th level
			end	
		end
		Contr.ord2.cnt(i).vec = Contr.ord2.cnt(i).vec + Contr.ord2.cnt(i).coef(j) * dmat(:, Contr.ord2.cnt(i).code(j).pos)'/count;			
		Contr.ord2.cnt(i).scalar = Contr.ord2.cnt(i).scalar + Contr.ord2.cnt(i).coef(j)*Contr.ord2.cnt(i).coef(j)/count;    % For variance coef: j-th term in i-th 1st-order contrast.		

%		Contr.ord2.cnt(i).scalar = tmp * Contr.ord2.cnt(i).coef*Contr.ord2.cnt(i).coef';		
%		Contr.ord2.cnt(i).vec = Contr.ord2.cnt(i).vec + Contr.ord2.cnt(i).coef(j) * dmat(:, Contr.ord2.cnt(i).code(j).pos)';			
   end	 %for (j = 1:1:Contr.ord2.cnt(i).NT),
%	Contr.ord2.cnt(i).vec = Contr.ord2.cnt(i).vec*tmp;   % assuming balanced design!!!
   end   %for (i = 1:1:Contr2.tot)
			
case 3,

   shift = 1 + num_col(1) + num_col(2);   % skip grand mean, factor means and 2nd-roder terms
   for (i = 1:1:Contr.ord3.tot),
   % for the  non-zero indices
	   Contr.ord3.cnt(i).vec = zeros(1, n);  % n = ntot in anovan.m: total number of files, and no. of arrows in dmat
	   Contr.ord3.cnt(i).scalar = 0;
		for (j = 1:1:Contr.ord3.cnt(i).NT),   % for each term in a contrast
         first = 0;
		   Contr.ord3.cnt(i).code(j).pos = 0;
			
		   for (k = 1:1:NF),
            if (Contr.ord3.cnt(i).code(j).str(k) >= 'a' & Contr.ord3.cnt(i).code(j).str(k) <= 'z'),
				   tmpv = 10 + Contr.ord3.cnt(i).code(j).str(k) - 'a';
			   else tmpv = Contr.ord3.cnt(i).code(j).str(k) - '0'; end
			
			   if (tmpv ~= 0),
	            first = first + 1;  % the first non-zero index backward!
		         % if (first == 1), tmp = str2num(Contr3.cnt(i).code(j).str(k)) - 1;	  % (level# - 1) for the first collapsed index	   		
				   % elseif (first == 2), Contr3.cnt(i).code(j).shift = Contr3.cnt(i).code(j).shift + tmp * FL(k).N_level str2num(Contr3.cnt(i).code(j).str(k);
				   if (first == 1),
				      Contr.ord3.cnt(i).idx1 = k;
					   %Contr3.cnt(i).code(j).pos = Contr3.cnt(i).shift + str2num(Contr3.cnt(i).code(j).str(k)) - 1;
%					   tmp = FL(k).N_level;
		         elseif (first == 2),
				      Contr.ord3.cnt(i).idx2 = k;
					   %Contr3.cnt(i).code(j).pos = Contr3.cnt(i).code(j).pos + (str2num(Contr3.cnt(i).code(j).str(k)) - 1) * tmp;
		            %else fprintf('\nError in contrast coding: more than two non-zero indices!\n');
		            %   fprintf(2,'Halted: Ctrl+c to exit'); pause;
				   end
				   if (first == 3), Contr.ord3.cnt(i).idx3 = k; end
		      end
	      end  % for (k = 1:1:NF)
		   switch Contr.ord3.cnt(i).idx1
		   case 1,
			   switch Contr.ord3.cnt(i).idx2
				   case 2,
					   switch Contr.ord3.cnt(i).idx3
						   case 3, % ABC
							
							if (Contr.ord3.cnt(i).code(j).str(1) >= 'a' & Contr.ord3.cnt(i).code(j).str(1) <= 'z'),
								tmpv1 = 10 + Contr.ord3.cnt(i).code(j).str(1) - 'a';
					      else tmpv1 = Contr.ord3.cnt(i).code(j).str(1) - '0'; end
							if (Contr.ord3.cnt(i).code(j).str(2) >= 'a' & Contr.ord3.cnt(i).code(j).str(2) <= 'z'),
							   tmpv2 = 10 + Contr.ord3.cnt(i).code(j).str(2) - 'a';
					      else tmpv2 = Contr.ord3.cnt(i).code(j).str(2) - '0'; end
							if (Contr.ord3.cnt(i).code(j).str(3) >= 'a' & Contr.ord3.cnt(i).code(j).str(3) <= 'z'),
							   tmpv3 = 10 + Contr.ord3.cnt(i).code(j).str(3) - 'a';
					      else tmpv3 = Contr.ord3.cnt(i).code(j).str(3) - '0'; end
							
							Contr.ord3.cnt(i).code(j).pos = shift + (tmpv1 - 1) * FL(2).N_level * FL(3).N_level + ...
					         (tmpv2 - 1) * FL(3).N_level + tmpv3; % ABC
%							tmp = FL(1).N_level*FL(2).N_level*FL(3).N_level/n;
							
							case 4, % ABD
							
							if (Contr.ord3.cnt(i).code(j).str(1) >= 'a' & Contr.ord3.cnt(i).code(j).str(1) <= 'z'),
								tmpv1 = 10 + Contr.ord3.cnt(i).code(j).str(1) - 'a';
					      else tmpv1 = Contr.ord3.cnt(i).code(j).str(1) - '0'; end
							if (Contr.ord3.cnt(i).code(j).str(2) >= 'a' & Contr.ord3.cnt(i).code(j).str(2) <= 'z'),
							   tmpv2 = 10 + Contr.ord3.cnt(i).code(j).str(2) - 'a';
					      else tmpv2 = Contr.ord3.cnt(i).code(j).str(2) - '0'; end
							if (Contr.ord3.cnt(i).code(j).str(4) >= 'a' & Contr.ord3.cnt(i).code(j).str(4) <= 'z'),
							   tmpv3 = 10 + Contr.ord3.cnt(i).code(j).str(4) - 'a';
					      else tmpv3 = Contr.ord3.cnt(i).code(j).str(4) - '0'; end
							
							% shift the columns for ABC and part of ABD
							Contr.ord3.cnt(i).code(j).pos = shift + FL(1).N_level * FL(2).N_level * FL(3).N_level + ...
							   (tmpv1 - 1) * FL(2).N_level * FL(4).N_level + (tmpv2 - 1) * FL(4).N_level + tmpv3; % ABD
%							   tmp = FL(1).N_level*FL(2).N_level*FL(4).N_level/n;

							case 5, % ABE
							
							if (Contr.ord3.cnt(i).code(j).str(1) >= 'a' & Contr.ord3.cnt(i).code(j).str(1) <= 'z'),
								tmpv1 = 10 + Contr.ord3.cnt(i).code(j).str(1) - 'a';
					      else tmpv1 = Contr.ord3.cnt(i).code(j).str(1) - '0'; end
							if (Contr.ord3.cnt(i).code(j).str(2) >= 'a' & Contr.ord3.cnt(i).code(j).str(2) <= 'z'),
							   tmpv2 = 10 + Contr.ord3.cnt(i).code(j).str(2) - 'a';
					      else tmpv2 = Contr.ord3.cnt(i).code(j).str(2) - '0'; end
							if (Contr.ord3.cnt(i).code(j).str(5) >= 'a' & Contr.ord3.cnt(i).code(j).str(5) <= 'z'),
							   tmpv3 = 10 + Contr.ord3.cnt(i).code(j).str(5) - 'a';
					      else tmpv3 = Contr.ord3.cnt(i).code(j).str(5) - '0'; end
							
							% shift the columns for ABC, ABD and part of ABE
							Contr.ord3.cnt(i).code(j).pos = shift + FL(1).N_level * FL(2).N_level * FL(3).N_level + ...
							   + FL(1).N_level * FL(2).N_level * FL(4).N_level + ...
								(tmpv1 - 1) * FL(2).N_level * FL(5).N_level + (tmpv2 - 1) * FL(5).N_level + tmpv3; % ABE

						end   % switch Contr.ord3.cnt(i).idx3
					case 3,
					   if (Contr.ord3.cnt(i).idx3 == 4),  % ACD
						   if (Contr.ord3.cnt(i).code(j).str(1) >= 'a' & Contr.ord3.cnt(i).code(j).str(1) <= 'z'),
								tmpv1 = 10 + Contr.ord3.cnt(i).code(j).str(1) - 'a';
					      else tmpv1 = Contr.ord3.cnt(i).code(j).str(1) - '0'; end
							if (Contr.ord3.cnt(i).code(j).str(3) >= 'a' & Contr.ord3.cnt(i).code(j).str(3) <= 'z'),
							   tmpv2 = 10 + Contr.ord3.cnt(i).code(j).str(3) - 'a';
					      else tmpv2 = Contr.ord3.cnt(i).code(j).str(3) - '0'; end
							if (Contr.ord3.cnt(i).code(j).str(4) >= 'a' & Contr.ord3.cnt(i).code(j).str(4) <= 'z'),
							   tmpv3 = 10 + Contr.ord3.cnt(i).code(j).str(4) - 'a';
					      else tmpv3 = Contr.ord3.cnt(i).code(j).str(4) - '0'; end
							
							Contr.ord3.cnt(i).code(j).pos = shift + FL(1).N_level * FL(2).N_level * FL(3).N_level + FL(1).N_level * FL(2).N_level * FL(4).N_level + ...
							   FL(1).N_level * FL(2).N_level * FL(5).N_level + (tmpv1 - 1) * FL(3).N_level * FL(4).N_level + (tmpv2 - 1) * FL(4).N_level + tmpv3; % ACD
%							tmp = FL(1).N_level*FL(3).N_level*FL(4).N_level/n;
						else fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;		
						end
					case 4, fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;					
				end				
			case 2,
			   switch Contr.ord3.cnt(i).idx2
				   case 3,
					   switch Contr.ord3.cnt(i).idx3
						
							case 4,  % BCD
							if (Contr.ord3.cnt(i).code(j).str(2) >= 'a' & Contr.ord3.cnt(i).code(j).str(2) <= 'z'),
								tmpv1 = 10 + Contr.ord3.cnt(i).code(j).str(2) - 'a';
					      else tmpv1 = Contr.ord3.cnt(i).code(j).str(2) - '0'; end
							if (Contr.ord3.cnt(i).code(j).str(3) >= 'a' & Contr.ord3.cnt(i).code(j).str(3) <= 'z'),
							   tmpv2 = 10 + Contr.ord3.cnt(i).code(j).str(3) - 'a';
					      else tmpv2 = Contr.ord3.cnt(i).code(j).str(3) - '0'; end
							if (Contr.ord3.cnt(i).code(j).str(4) >= 'a' & Contr.ord3.cnt(i).code(j).str(4) <= 'z'),
							   tmpv3 = 10 + Contr.ord3.cnt(i).code(j).str(4) - 'a';
					      else tmpv3 = Contr.ord3.cnt(i).code(j).str(4) - '0'; end
							
							Contr.ord3.cnt(i).code(j).pos = shift + FL(1).N_level*FL(2).N_level*FL(3).N_level + FL(1).N_level*FL(2).N_level*FL(4).N_level + ...
						      FL(1).N_level * FL(2).N_level * FL(5).N_level + FL(1).N_level*FL(3).N_level*FL(4).N_level + FL(1).N_level * FL(3).N_level * FL(5).N_level + ...
								FL(1).N_level * FL(4).N_level * FL(5).N_level + (tmpv1 - 1) * FL(3).N_level * FL(4).N_level + (tmpv2 - 1)* FL(4).N_level + tmpv3;
								
							case 5, % BCE
							if (Contr.ord3.cnt(i).code(j).str(2) >= 'a' & Contr.ord3.cnt(i).code(j).str(2) <= 'z'),
								tmpv1 = 10 + Contr.ord3.cnt(i).code(j).str(2) - 'a';
					      else tmpv1 = Contr.ord3.cnt(i).code(j).str(2) - '0'; end
							if (Contr.ord3.cnt(i).code(j).str(3) >= 'a' & Contr.ord3.cnt(i).code(j).str(3) <= 'z'),
							   tmpv2 = 10 + Contr.ord3.cnt(i).code(j).str(3) - 'a';
					      else tmpv2 = Contr.ord3.cnt(i).code(j).str(3) - '0'; end
							if (Contr.ord3.cnt(i).code(j).str(5) >= 'a' & Contr.ord3.cnt(i).code(j).str(5) <= 'z'),
							   tmpv3 = 10 + Contr.ord3.cnt(i).code(j).str(5) - 'a';
					      else tmpv3 = Contr.ord3.cnt(i).code(j).str(5) - '0'; end
							
							Contr.ord3.cnt(i).code(j).pos = shift + FL(1).N_level*FL(2).N_level*FL(3).N_level + FL(1).N_level*FL(2).N_level*FL(4).N_level + ...
						      FL(1).N_level*FL(2).N_level*FL(5).N_level + FL(1).N_level*FL(3).N_level*FL(4).N_level + ...
								FL(1).N_level*FL(3).N_level*FL(5).N_level + FL(1).N_level*FL(4).N_level*FL(5).N_level + FL(2).N_level * FL(3).N_level * FL(4).N_level + ...
								(tmpv1 - 1) * FL(3).N_level * FL(5).N_level + (tmpv2 - 1)* FL(5).N_level + tmpv3;	
						
					end	% switch Contr.ord3.cnt(i).idx3
					
					case 4, % Do I want BDE here?
					   if (Contr.ord3.cnt(i).idx3 == 5),
						   if (Contr.ord3.cnt(i).code(j).str(2) >= 'a' & Contr.ord3.cnt(i).code(j).str(2) <= 'z'),
								tmpv1 = 10 + Contr.ord3.cnt(i).code(j).str(2) - 'a';
					      else tmpv1 = Contr.ord3.cnt(i).code(j).str(2) - '0'; end
							if (Contr.ord3.cnt(i).code(j).str(4) >= 'a' & Contr.ord3.cnt(i).code(j).str(4) <= 'z'),
							   tmpv2 = 10 + Contr.ord3.cnt(i).code(j).str(4) - 'a';
					      else tmpv2 = Contr.ord3.cnt(i).code(j).str(4) - '0'; end
							if (Contr.ord3.cnt(i).code(j).str(5) >= 'a' & Contr.ord3.cnt(i).code(j).str(4) <= 'z'),
							   tmpv3 = 10 + Contr.ord3.cnt(i).code(j).str(5) - 'a';
					      else tmpv3 = Contr.ord3.cnt(i).code(j).str(5) - '0'; end
							
							Contr.ord3.cnt(i).code(j).pos = shift + FL(1).N_level*FL(2).N_level*FL(3).N_level + FL(1).N_level*FL(2).N_level*FL(4).N_level + ...
						         FL(1).N_level*FL(2).N_level*FL(5).N_level + FL(1).N_level*FL(3).N_level*FL(4).N_level + ...
								   FL(1).N_level*FL(3).N_level*FL(5).N_level + FL(1).N_level*FL(4).N_level*FL(5).N_level + ...
									FL(2).N_level*FL(3).N_level*FL(4).N_level + FL(2).N_level*FL(3).N_level*FL(5).N_level + ...
									(tmpv1 - 1) * FL(4).N_level * FL(5).N_level + (tmpv2 - 1)* FL(5).N_level + tmpv3;
						else fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;	end
					
				end	% switch Contr.ord3.cnt(i).idx2
			case 3,  % Do I want CDE here?
			   if (Contr.ord3.cnt(i).idx1 == 3 & Contr.ord3.cnt(i).idx2 == 4 & Contr.ord3.cnt(i).idx3 == 5),
				   if (Contr.ord3.cnt(i).code(j).str(3) >= 'a' & Contr.ord3.cnt(i).code(j).str(3) <= 'z'),
						tmpv1 = 10 + Contr.ord3.cnt(i).code(j).str(3) - 'a';
					else tmpv1 = Contr.ord3.cnt(i).code(j).str(2) - '0'; end
					if (Contr.ord3.cnt(i).code(j).str(4) >= 'a' & Contr.ord3.cnt(i).code(j).str(4) <= 'z'),
					   tmpv2 = 10 + Contr.ord3.cnt(i).code(j).str(4) - 'a';
					else tmpv2 = Contr.ord3.cnt(i).code(j).str(4) - '0'; end
					if (Contr.ord3.cnt(i).code(j).str(5) >= 'a' & Contr.ord3.cnt(i).code(j).str(4) <= 'z'),
					   tmpv3 = 10 + Contr.ord3.cnt(i).code(j).str(5) - 'a';
					else tmpv3 = Contr.ord3.cnt(i).code(j).str(5) - '0'; end
							
					Contr.ord3.cnt(i).code(j).pos = shift + FL(1).N_level*FL(2).N_level*FL(3).N_level + FL(1).N_level*FL(2).N_level*FL(4).N_level + ...
					   FL(1).N_level*FL(2).N_level*FL(5).N_level + FL(1).N_level*FL(3).N_level*FL(4).N_level + ...
					   FL(1).N_level*FL(3).N_level*FL(5).N_level + FL(1).N_level*FL(4).N_level*FL(5).N_level + ...
						FL(2).N_level*FL(3).N_level*FL(4).N_level + FL(2).N_level*FL(3).N_level*FL(5).N_level + ...
							+ FL(2).N_level*FL(4).N_level*FL(5).N_level + (tmpv1 - 1) * FL(4).N_level * FL(5).N_level + (tmpv2 - 1)* FL(5).N_level + tmpv3;
				else fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;	end
								
			case 4,
			   fprintf('\nSomething is wrong in the contrast coding!\n');
		      fprintf(2,'Halted: Ctrl+c to exit'); pause;			
		   end   % switch Contr.ord3.cnt(i).idx1
			
			count = 0;
		   for (ll = 1:1:n)  % Should have a more decent way to do this in a matrix fashion instead of a lenthy loop?
		      tmpstr1 = char (group{Contr.ord3.cnt(i).idx1}(ll));
				tmpstr2 = char (group{Contr.ord3.cnt(i).idx2}(ll));
				tmpstr3 = char (group{Contr.ord3.cnt(i).idx3}(ll));
				if strcmp(tmpstr1, FL(Contr.ord3.cnt(i).idx1).level(tmpv1).expr) & strcmp(tmpstr2, FL(Contr.ord3.cnt(i).idx2).level(tmpv2).expr) ...
				   & strcmp(tmpstr3, FL(Contr.ord3.cnt(i).idx3).level(tmpv3).expr)
			      count = count + 1;   % total number of occurences at tmpv-th level
			   end	
			end
		   Contr.ord3.cnt(i).vec = Contr.ord3.cnt(i).vec + Contr.ord3.cnt(i).coef(j) * dmat(:, Contr.ord3.cnt(i).code(j).pos)'/count;			
		   Contr.ord3.cnt(i).scalar = Contr.ord3.cnt(i).scalar + Contr.ord3.cnt(i).coef(j)*Contr.ord3.cnt(i).coef(j)/count;    % For variance coef: j-th term in i-th 3rd-order contrast.	
			
%		   Contr.ord3.cnt(i).scalar = tmp * Contr.ord3.cnt(i).coef*Contr.ord3.cnt(i).coef';		
%		   Contr.ord3.cnt(i).vec = Contr.ord3.cnt(i).vec + Contr.ord3.cnt(i).coef(j) * dmat(:, Contr.ord3.cnt(i).code(j).pos)';			
      end	% for (j = 1:1:Contr.ord3.cnt(i).NT)
%	   Contr.ord3.cnt(i).vec = Contr.ord3.cnt(i).vec*tmp;   % assuming balanced design!!!
   end   %for (i = 1:1:Contr3.tot)

case 4,

   shift = 1 + num_col(1) + num_col(2) + num_col(3);   % skip grand mean, factor means and 2nd-roder terms
   for (i = 1:1:Contr.ord4.tot),
   % for the  non-zero indices
	   Contr.ord4.cnt(i).vec = zeros(1, n);  % n = ntot in anovan.m: total number of input files, and no. of arrows in dmat
	   Contr.ord4.cnt(i).scalar = 0;
		for (j = 1:1:Contr.ord4.cnt(i).NT),   % for each term in a contrast
         first = 0;
		   Contr.ord4.cnt(i).code(j).pos = 0;
			
		   for (k = 1:1:NF),
            if (Contr.ord4.cnt(i).code(j).str(k) >= 'a' & Contr.ord4.cnt(i).code(j).str(k) <= 'z'),
				   tmpv = 10 + Contr.ord4.cnt(i).code(j).str(k) - 'a';
			   else tmpv = Contr.ord4.cnt(i).code(j).str(k) - '0'; end
			
			   if (tmpv ~= 0),
	            first = first + 1;  % the first non-zero index
					% Maybe I should replace all idx with idx() in *.m and Contr.ord4.cnt(i).idx(first) = k when first ~= 0 later?
				   switch first
					case 1,
				      Contr.ord4.cnt(i).idx1 = k;
		         case 2,
				      Contr.ord4.cnt(i).idx2 = k;
				   case 3,
				      Contr.ord4.cnt(i).idx3 = k;
					case 4,
					   Contr.ord4.cnt(i).idx4 = k;
					end		
		      end  % if (tmpv ~= 0)
	      end  % for (k = 1:1:NF)
		   switch Contr.ord4.cnt(i).idx1
		   case 1,
			   switch Contr.ord4.cnt(i).idx2
				case 2,
				   switch Contr.ord4.cnt(i).idx3
					case 3,
						switch Contr.ord4.cnt(i).idx4
						case 4,	% ABCD				
							
							if (Contr.ord4.cnt(i).code(j).str(1) >= 'a' & Contr.ord4.cnt(i).code(j).str(1) <= 'z'),
								tmpv1 = 10 + Contr.ord4.cnt(i).code(j).str(1) - 'a';
					      else tmpv1 = Contr.ord4.cnt(i).code(j).str(1) - '0'; end
							if (Contr.ord4.cnt(i).code(j).str(2) >= 'a' & Contr.ord4.cnt(i).code(j).str(2) <= 'z'),
							   tmpv2 = 10 + Contr.ord4.cnt(i).code(j).str(2) - 'a';
					      else tmpv2 = Contr.ord4.cnt(i).code(j).str(2) - '0'; end
							if (Contr.ord4.cnt(i).code(j).str(3) >= 'a' & Contr.ord4.cnt(i).code(j).str(3) <= 'z'),
							   tmpv3 = 10 + Contr.ord4.cnt(i).code(j).str(3) - 'a';
					      else tmpv3 = Contr.ord4.cnt(i).code(j).str(3) - '0'; end
							if (Contr.ord4.cnt(i).code(j).str(4) >= 'a' & Contr.ord4.cnt(i).code(j).str(4) <= 'z'),
							   tmpv4 = 10 + Contr.ord4.cnt(i).code(j).str(4) - 'a';
					      else tmpv4 = Contr.ord4.cnt(i).code(j).str(4) - '0'; end
							
							Contr.ord4.cnt(i).code(j).pos = shift + (tmpv1 - 1) * FL(2).N_level * FL(3).N_level * FL(4).N_level + ...
					         (tmpv2 - 1) * FL(3).N_level * FL(4).N_level + (tmpv3 - 1) * FL(4).N_level + tmpv4;

						case 5, % ABCE
						   if (Contr.ord4.cnt(i).code(j).str(1) >= 'a' & Contr.ord4.cnt(i).code(j).str(1) <= 'z'),
								tmpv1 = 10 + Contr.ord4.cnt(i).code(j).str(1) - 'a';
					      else tmpv1 = Contr.ord4.cnt(i).code(j).str(1) - '0'; end
							if (Contr.ord4.cnt(i).code(j).str(2) >= 'a' & Contr.ord4.cnt(i).code(j).str(2) <= 'z'),
							   tmpv2 = 10 + Contr.ord4.cnt(i).code(j).str(2) - 'a';
					      else tmpv2 = Contr.ord4.cnt(i).code(j).str(2) - '0'; end
							if (Contr.ord4.cnt(i).code(j).str(3) >= 'a' & Contr.ord4.cnt(i).code(j).str(3) <= 'z'),
							   tmpv3 = 10 + Contr.ord4.cnt(i).code(j).str(3) - 'a';
					      else tmpv3 = Contr.ord4.cnt(i).code(j).str(3) - '0'; end
							if (Contr.ord4.cnt(i).code(j).str(5) >= 'a' & Contr.ord4.cnt(i).code(j).str(5) <= 'z'),
							   tmpv4 = 10 + Contr.ord4.cnt(i).code(j).str(5) - 'a';
					      else tmpv4 = Contr.ord4.cnt(i).code(j).str(5) - '0'; end							
							
							Contr.ord4.cnt(i).code(j).pos = shift + FL(1).N_level * FL(2).N_level * FL(3).N_level * FL(4).N_level + ...
					         (tmpv1 - 1) * FL(2).N_level * FL(3).N_level * FL(5).N_level + (tmpv2 - 1) * FL(3).N_level * FL(5).N_level + ...
								(tmpv3 - 1) * FL(5).N_level + tmpv4; % ABCD
						end  % switch Contr.ord4.cnt(i).idx4
					case 4,
					   if (Contr.ord4.cnt(i).idx4 == 5)
						   if (Contr.ord4.cnt(i).code(j).str(1) >= 'a' & Contr.ord4.cnt(i).code(j).str(1) <= 'z'),
								tmpv1 = 10 + Contr.ord4.cnt(i).code(j).str(1) - 'a';
					      else tmpv1 = Contr.ord4.cnt(i).code(j).str(1) - '0'; end
							if (Contr.ord4.cnt(i).code(j).str(2) >= 'a' & Contr.ord4.cnt(i).code(j).str(2) <= 'z'),
							   tmpv2 = 10 + Contr.ord4.cnt(i).code(j).str(2) - 'a';
					      else tmpv2 = Contr.ord4.cnt(i).code(j).str(2) - '0'; end
							if (Contr.ord4.cnt(i).code(j).str(4) >= 'a' & Contr.ord4.cnt(i).code(j).str(4) <= 'z'),
							   tmpv3 = 10 + Contr.ord4.cnt(i).code(j).str(4) - 'a';
					      else tmpv3 = Contr.ord4.cnt(i).code(j).str(4) - '0'; end
							if (Contr.ord4.cnt(i).code(j).str(5) >= 'a' & Contr.ord4.cnt(i).code(j).str(5) <= 'z'),
							   tmpv4 = 10 + Contr.ord4.cnt(i).code(j).str(5) - 'a';
					      else tmpv4 = Contr.ord4.cnt(i).code(j).str(5) - '0'; end							
							
							Contr.ord4.cnt(i).code(j).pos = shift + FL(1).N_level * FL(2).N_level * FL(3).N_level * FL(4).N_level + ...
					         FL(1).N_level * FL(2).N_level * FL(3).N_level * FL(5).N_level + ...
								(tmpv1 - 1) * FL(2).N_level * FL(4).N_level * FL(5).N_level + (tmpv2 - 1) * FL(4).N_level * FL(5).N_level + ...
								(tmpv3 - 1) * FL(5).N_level + tmpv4; % ABCD
						else fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;	end  % if (Contr.ord4.cnt(i).idx4 == 5)
					end	% switch Contr.ord4.cnt(i).idx3
				case 3,
				   if (Contr.ord4.cnt(i).idx3 == 4 &  Contr.ord4.cnt(i).idx4 == 5), % ACDE
						   if (Contr.ord4.cnt(i).code(j).str(1) >= 'a' & Contr.ord4.cnt(i).code(j).str(1) <= 'z'),
								tmpv1 = 10 + Contr.ord4.cnt(i).code(j).str(1) - 'a';
					      else tmpv1 = Contr.ord4.cnt(i).code(j).str(1) - '0'; end
							if (Contr.ord4.cnt(i).code(j).str(3) >= 'a' & Contr.ord4.cnt(i).code(j).str(3) <= 'z'),
							   tmpv2 = 10 + Contr.ord4.cnt(i).code(j).str(3) - 'a';
					      else tmpv2 = Contr.ord4.cnt(i).code(j).str(3) - '0'; end
							if (Contr.ord4.cnt(i).code(j).str(4) >= 'a' & Contr.ord4.cnt(i).code(j).str(4) <= 'z'),
							   tmpv3 = 10 + Contr.ord4.cnt(i).code(j).str(4) - 'a';
					      else tmpv3 = Contr.ord4.cnt(i).code(j).str(4) - '0'; end
							if (Contr.ord4.cnt(i).code(j).str(5) >= 'a' & Contr.ord4.cnt(i).code(j).str(5) <= 'z'),
							   tmpv4 = 10 + Contr.ord4.cnt(i).code(j).str(5) - 'a';
					      else tmpv4 = Contr.ord4.cnt(i).code(j).str(5) - '0'; end							
							
							Contr.ord4.cnt(i).code(j).pos = shift + FL(1).N_level * FL(2).N_level * FL(3).N_level * FL(4).N_level + ...
					         FL(1).N_level * FL(2).N_level * FL(3).N_level * FL(5).N_level + FL(1).N_level * FL(2).N_level * FL(4).N_level * FL(5).N_level + ...
								(tmpv1 - 1) * FL(3).N_level * FL(4).N_level * FL(5).N_level + (tmpv2 - 1) * FL(4).N_level * FL(5).N_level + ...
								(tmpv3 - 1) * FL(5).N_level + tmpv4; % ABCD
					else fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;	end			
				end  % switch Contr.ord4.cnt(i).idx2
			case 2,
			   if (Contr.ord4.cnt(i).idx2 == 3 & Contr.ord4.cnt(i).idx3 == 4 & Contr.ord4.cnt(i).idx4 == 5),  % BCDE
						   if (Contr.ord4.cnt(i).code(j).str(2) >= 'a' & Contr.ord4.cnt(i).code(j).str(2) <= 'z'),
								tmpv1 = 10 + Contr.ord4.cnt(i).code(j).str(2) - 'a';
					      else tmpv1 = Contr.ord4.cnt(i).code(j).str(2) - '0'; end
							if (Contr.ord4.cnt(i).code(j).str(3) >= 'a' & Contr.ord4.cnt(i).code(j).str(3) <= 'z'),
							   tmpv2 = 10 + Contr.ord4.cnt(i).code(j).str(3) - 'a';
					      else tmpv2 = Contr.ord4.cnt(i).code(j).str(3) - '0'; end
							if (Contr.ord4.cnt(i).code(j).str(4) >= 'a' & Contr.ord4.cnt(i).code(j).str(4) <= 'z'),
							   tmpv3 = 10 + Contr.ord4.cnt(i).code(j).str(4) - 'a';
					      else tmpv3 = Contr.ord4.cnt(i).code(j).str(4) - '0'; end
							if (Contr.ord4.cnt(i).code(j).str(5) >= 'a' & Contr.ord4.cnt(i).code(j).str(5) <= 'z'),
							   tmpv4 = 10 + Contr.ord4.cnt(i).code(j).str(5) - 'a';
					      else tmpv4 = Contr.ord4.cnt(i).code(j).str(5) - '0'; end							
							
							Contr.ord4.cnt(i).code(j).pos = shift + FL(1).N_level * FL(2).N_level * FL(3).N_level * FL(4).N_level + ...
					         FL(1).N_level * FL(2).N_level * FL(3).N_level * FL(5).N_level + FL(1).N_level * FL(2).N_level * FL(4).N_level * FL(5).N_level + ...
								FL(1).N_level * FL(3).N_level * FL(4).N_level * FL(5).N_level + (tmpv1 - 1) * FL(3).N_level * FL(4).N_level * FL(5).N_level + ...
								(tmpv2 - 1) * FL(4).N_level * FL(5).N_level + (tmpv3 - 1) * FL(5).N_level + tmpv4; % ABCD
				else fprintf('\nSomething is wrong in the contrast coding!\n'); fprintf(2,'Halted: Ctrl+c to exit'); pause;	end								
			case {3,4,}
			   fprintf('\nSomething is wrong in the contrast coding!\n');
		      fprintf(2,'Halted: Ctrl+c to exit'); pause;			
		   end   % switch Contr.ord4.cnt(i).idx1
			
			count = 0;
		   for (ll = 1:1:n)  % Should have a more decent way to do this in a matrix fashion instead of a lenthy loop?
		      tmpstr1 = char (group{Contr.ord4.cnt(i).idx1}(ll));
				tmpstr2 = char (group{Contr.ord4.cnt(i).idx2}(ll));
				tmpstr3 = char (group{Contr.ord4.cnt(i).idx3}(ll));
				tmpstr4 = char (group{Contr.ord4.cnt(i).idx4}(ll));
				if strcmp(tmpstr1, FL(Contr.ord4.cnt(i).idx1).level(tmpv1).expr) & strcmp(tmpstr2, FL(Contr.ord4.cnt(i).idx2).level(tmpv2).expr) ...
				   & strcmp(tmpstr3, FL(Contr.ord4.cnt(i).idx3).level(tmpv3).expr) & strcmp(tmpstr4, FL(Contr.ord4.cnt(i).idx4).level(tmpv4).expr)
			      count = count + 1;   % total number of occurences at tmpv-th level
			   end	
			end % for (ll = 1:1:n)
		   Contr.ord4.cnt(i).vec = Contr.ord4.cnt(i).vec + Contr.ord4.cnt(i).coef(j) * dmat(:, Contr.ord4.cnt(i).code(j).pos)'/count;			
		   Contr.ord4.cnt(i).scalar = Contr.ord4.cnt(i).scalar + Contr.ord4.cnt(i).coef(j)*Contr.ord4.cnt(i).coef(j)/count;    % For variance coef: j-th term in i-th 4th-order contrast.				
      end	% for (j = 1:1:Contr.ord4.cnt(i).NT)
   end   %for (i = 1:1:Contr4.tot)

end % switch order

err = 0;
return;
