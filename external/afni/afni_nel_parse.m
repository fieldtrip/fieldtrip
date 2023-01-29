function [nel] = afni_nel_parse(nellcell)
   %should really check that this is nel and not ngr ...

   %assuming we have nel

   %get the header part
   cnel = char(nellcell); clear nellcell;
   ncnel = length(cnel);

   %break up into head and data
      expr = '(?<head><(\w+).*?>)(?<data>.*)<';
         %using expr = '(?<head><(\w+).*?>)(?<data>.*<.*>)';
         %would get all of the data chunk including its closing tag
         %but then I'd have to strip it...
      [pp] = regexp(cnel,expr,'names');

   %Now get the header fields into structs
   hh = regexp(pp.head, '(?<lhs>\w*)="(?<rhs>\S*)"','names');
   nel = cell2struct({hh(:).rhs}, {hh(:).lhs},2); clear hh;


   tt = nel.ni_type;
   nel.vec_typ=zeros(1,1000); nn=0;
   while (~isempty(tt)),
      [ttt,tt] = strtok(tt, ',');
      %look for the N*type syntax
         N = 1;
         [tttt,ttt] = strtok(ttt,'*');
         Ntttt = str2double(tttt);
         if (~isnan(Ntttt)),  %have a number, get the type
            N = Ntttt;
            tttt = strtok(ttt,'*');
         end
         nel.vec_typ(nn+1:1:nn+N) = afni_ni_rowtype_name_to_code(tttt);
         nn = nn + N;
   end

   %interpret some things
   nel.name = char(strtok(regexp(pp.head,'<\w+', 'match'),'<'));
   nel.vec_typ = nel.vec_typ(1:nn);
   nel.vec_len = str2num(nel.ni_dimen);
   nel.vec_num = length(nel.vec_typ);

   %Now parse the data
   if (~afni_ni_is_numeric_type(nel.vec_typ)),
      fprintf(2,'Data not all numeric, will not parse it');
      nel.data = pp.data;
   else
      nel = afni_nel_parse_data(nel, pp.data);
   end

   clear pp

   return;

function [nel] = afni_nel_parse_data(nel, data)
      nel.data = reshape(sscanf(data,'%f'), nel.vec_num, nel.vec_len)';
   return;
