function A = make_A(M,x1,x2,x3,dG1,dG2,dG3,wt,lkp)
% Matrix of rate of change of weighted difference w.r.t. parameter changes

% copy-pasted from spm_realign (vers 8)

p0 = [0 0 0  0 0 0  1 1 1  0 0 0];
A  = zeros(numel(x1),length(lkp));
for i=1:length(lkp)
    pt         = p0;
    pt(lkp(i)) = pt(i)+1e-6;
    [y1,y2,y3] = coords(pt,M,M,x1,x2,x3);
    tmp        = sum([y1-x1 y2-x2 y3-x3].*[dG1 dG2 dG3],2)/(-1e-6);
    if ~isempty(wt), A(:,i) = tmp.*wt;
    else A(:,i) = tmp; end
end
