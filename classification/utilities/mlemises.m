function params = mlemises(phases) %Fisher p.88

    mu = angle(sum(exp(i*phases)));
    R = abs(mean(exp(i*phases)));
    %since we have typically a small number of trials we use the following
    %estimation of k

    if R < 0.53 % besseli approximation
        k = 2*R + R^3 + 5*(R^5)/6;
    elseif (R>=0.53&&R<0.85)
        k = -0.4 + 1.39*R + 0.43/(1-R);
    elseif R>=0.85
        k = 1/(R^3-4*(R^2)+3*R);
    end

    n = length(phases);
    if k<2
        k = k-2/(n*k);
        if k<0,k = 0;
        end
    end
    if k>2
        k = (n-1)^3*k/(n^3+n);
    end

    params = [mu k];
end