function [P1, Pt1, PX, L, sigma2, st]=cpd_Pfast(X,T,sigma2, outliers, sigma2_init, fgt)

switch fgt,
    case 1
        if sigma2<0.05, sigma2=0.05; end;
        [P1, Pt1, PX, L]=cpd_P_FGT(X, T, sigma2, outliers, sigma2_init); st='(FGT)';

    case 2
        if (sigma2 > 0.015*sigma2_init) % FGT sqrt(2/(N+M)
            [P1, Pt1, PX, L]=cpd_P_FGT(X, T, sigma2, outliers, sigma2_init); st='(FGT)';
        else
            % quite FGT, switch to the truncated kernel approximation
            [P1,Pt1, PX, L]=cpd_Pappmex(X,T, sigma2 ,outliers,1e-3);  st='(Truncated)';
        end
end