function a = computePrior(pvi, experienced, republican, parms)
    x = ones(size(pvi,1), 6);
    x(:, 2) = pvi;
    democratic = 1 - republican;
    x(:, 3) = experienced;
    x(:, 4) = democratic;
    x(:, 5) = republican;
    x(:, 6) = pvi*republican;
    a = x*parms.coefs;
end