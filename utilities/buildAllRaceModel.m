function [bestnlZ, besthyp] = buildAllRaceModel(data, years, states, parms)
    meanmask = [true, false, false, false];
    covmask = [false, true, true, false];
    firmmask = [false, false, false, true];
    mc = {@meanConst};
    ml = {@meanLinear};
    md = {@meanDiscrete, parms.nfirm};
    mb = {@logsqrtbinom};
    mml = {@meanMask, meanmask, ml};
    mmc = {@meanMask, meanmask, mc};
    mmd = {@meanMask, firmmask, md};
    meanfunc = {@meanSum, {mml, mmc, mmd}};
    cm = {@covMaterniso, 3};
    cdb = {@covDiag, mb};
    cdf = {@covDiag, md};
    cmf = {@covMask, {firmmask, cdf}};
    cmd = {@covMask, {covmask, cdb}};
    cmm = {@covMask, {meanmask, cm}};
    covfunc = {@covSum, {cmm, cmd, cmf}};
    likfunc = @likGauss;
    inffunc = @infExact;

    % setting hyperpriors
    % using linear model as hyperpriors for mean prior
    
    mu_ml = 0; sigma_ml = 0.05/250;
    mu_mc = 0.5; sigma_mc = 0.05;    

    pg_ml = {@priorGauss, mu_ml, sigma_ml^2};
    pg_mc = {@priorGauss, mu_mc, sigma_mc^2};
    prior.mean = {pg_ml, pg_mc};
    
    for i=1:parms.nfirm
        mu_md = 0; sigma_md = 0.05;
        pg_md = {@priorGauss, mu_md, sigma_md^2};
        prior.mean{2+i} = pg_md;
    end

    % length scale: 30 days in average,
    % no less than 10 days, no more than 90 days
    % output scale: binomial variance
    mu_ls = log(30); sigma_ls = log(90/30)/2;
    mu_os = log(1/40); sigma_os = log(2);
    
    pg_ls = {@priorGauss, mu_ls, sigma_ls^2};
    pg_os = {@priorGauss, mu_os, sigma_os^2};
    prior.cov = {pg_ls, pg_os};
    
    for i=1:parms.nfirm
        mu_vf = log(0.1); sigma_vf = log(2)/2;
        pg_vf = {@priorGauss, mu_vf, sigma_vf^2};
        prior.cov{2+i} = pg_vf;
    end

    % total variance: CI for unknown sigma (chi-squared)
    mu_lik = log(0.05); sigma_lik = log(2)/2;
        
    pg_lik = {@priorGauss, mu_lik, sigma_lik^2};
    prior.lik = {pg_lik};
    
    [xs, ys] = getCellArrayData(data, years, states);
   
    im = {@infPrior, inffunc, prior};
    par = {meanfunc,covfunc,likfunc, xs, ys};
    mfun = @minimize_v2;
    
    bestnlZ = 0;
    for i=1:5
        hyp = sample_prior(prior);
        hyp = feval(mfun, hyp,  @gp_likelihood_independent, -100, im, par{:});
        [nlZ, ~] = gp_likelihood_independent(hyp, im, par{:});
        if nlZ < bestnlZ, bestnlZ = nlZ; besthyp = hyp; end
    end

    % we do not care about ymu, ys2
    % [~, ~, fmu, fs2] = gp(besthyp, inffunc, meanfunc, covfunc, likfunc, xs, ys, xstar);
end