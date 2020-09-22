function a = computePrior(pvi, experienced, party, parms)
%
%  Compute prior mean on intercept given fundementals and precomputed parms.coefs
%
    x = ones(size(pvi,1), 6);
    x(:, 2) = pvi;
    republican = (party==1);
    democratic = (party==-1);
    x(:, 3) = experienced;
    x(:, 4) = democratic;
    x(:, 5) = republican;
    x(:, 6) = pvi*republican;
    a = x*parms.coefs;
end