function [meanfunc, covfunc, likfunc, inffunc, prior] = model(parms)

    meanmask = [true, false, false];
    covmask = [false, true, true];
    mc = {@meanConst};
    ml = {@meanLinear};
    mml = {@meanMask, meanmask, ml};
    mmc = {@meanMask, meanmask, mc};

    % try no firm
    meanfunc = {@meanSum, {mml, mmc}};
    cm = {@covMaterniso, 3};
    mb = {@logsqrtbinom};
    cdb = {@covDiag, mb};
    cmd = {@covMask, {covmask, cdb}};
    cmm = {@covMask, {meanmask, cm}};

    cs = {@covLINiso};
    ci = {@covConst};
    cms = {@covMask, {meanmask, cs}};

    covfunc = {@covSum, {cmm, cmd, cms, ci}};
    likfunc = @likGauss;
    inffunc = @infExact;

    mu_ml = 0.0; sigma_ml = 0.002;
    mu_mc = 0.5; sigma_mc = 0.1;
    prior.slope = [mu_ml, sigma_ml];
    prior.intercept = [mu_mc, sigma_mc];
%     pc_ml = {@priorClamped};
%     pc_mc = {@priorClamped};
%     prior.mean = {pc_ml, pc_mc};
% 
% 
%     mu_ls = log(30); sigma_ls = log(90/30)/2;
%     mu_os = log(1/40); sigma_os = log(2);
%     pg_ls = {@priorGauss, mu_ls, sigma_ls^2};
%     pg_os = {@priorGauss, mu_os, sigma_os^2};
% 
%     pg_ml = {@priorClamped};
%     pg_mc = {@priorClamped};
% 
%     prior.cov = {pg_ls, pg_os, pg_ml, pg_mc};
% 
%     mu_lik = log(0.03); sigma_lik = log(2);
%     pg_lik = {@priorGauss, mu_lik, sigma_lik^2};
%     prior.lik = {pg_lik};
end
