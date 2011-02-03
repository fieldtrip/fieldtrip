function l = logdet(M,errorDet)

[R,p] = chol(M);
if p ~= 0
   l = errorDet;
else
   l = 2*sum(log(diag(R)));
end;

