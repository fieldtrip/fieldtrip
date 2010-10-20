K=32;
order = 4;

B = zeros(K,1+order);
A = zeros(K,1+order);
B(1,1) = 1;
A(1,1) = 1;

for k=2:K; 
  [B(k,:),A(k,:)] = cheby2(order, 20, 0.5/k); 
end

fmtString = '  {%12.8f';
for i=1:order
   fmtString = [fmtString ',%12.8f'];
end
fmtString = [fmtString '}'];

fprintf(1,'double chebyCoefsB[%i][%i] = {\n', K, 1+order);
for k=1:K;
   fprintf(1,fmtString, B(k,:));
   if k==K
      fprintf(1,'\n};\n');
   else
      fprintf(1,',\n');
   end 
end
        

fprintf(1,'double chebyCoefsA[%i][%i] = {\n', K, 1+order);
for k=1:K;
   fprintf(1,fmtString, A(k,:));
   if k==K
      fprintf(1,'\n};\n');
   else
      fprintf(1,',\n');
   end 
end