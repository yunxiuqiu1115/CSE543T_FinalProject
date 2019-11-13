function [p, fig] = oneCandidatePredict(x, y, parms)
    n = size(x,1);
    nz = 200;
    xs = linspace(x(1,1)-10,0,nz).';
    xs = [xs, zeros(1,nz)',ones(1,nz)', parms.nfirm*ones(1,nz)'];

    % define mean function and kernel
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
    meanfunc = {@meanSum, {mml, mmc}};
    cm = {@covMaterniso, 3};
    cdb = {@covDiag, mb};
    cdf = {@covDiag, md};
    cmf = {@covMask, {firmmask, cdf}};
    cmd = {@covMask, {covmask, cdb}};
    cmm = {@covMask, {meanmask, cm}};
    covfunc = {@covSum, {cmm, cmd}};
    likfunc = @likGauss;
    inffunc = @infExact;

    % setting hyperpriors
    % using linear model as hyperpriors for mean prior
    
    mu_ml = 0; sigma_ml = 0.05/x(1,1);
    mu_mc = mean(y); sigma_mc = 0.05;

%     [b,bint] = regress(y,[x(:,1), ones(n,1)]);
%     mu_ml = b(1); sigma_ml = (bint(1,2) - bint(1,1))/4;
%     mu_mc = b(2); sigma_mc = (bint(2,2) - bint(2,1))/4;
%     if n==2
%       sigma_ml = mu_ml;
%       sigma_mc = 0.05;
%     end       

    pg_ml = {@priorGauss, mu_ml, sigma_ml^2};
    pg_mc = {@priorGauss, mu_mc, sigma_mc^2};
    prior.mean = {pg_ml, pg_mc};
    
%     for i=1:parms.nfirm
%         mu_md = 0; sigma_md = 0.05;
%         pg_md = {@priorGauss, mu_md, sigma_md^2};
%         prior.mean{2+i} = pg_md;
%     end

    % length scale: 30 days in average,
    % no less than 10 days, no more than 90 days
    % output scale: binomial variance
    mu_ls = log(30); sigma_ls = log(90/30)/2;
    mu_os = log(1/40); sigma_os = log(2);
    
%     binom_var = log(sqrt(x(:,2).*(1-x(:,2))./x(:,3)));
%     mu_os = mean(binom_var); sigma_os = std(binom_var);
    
    pg_ls = {@priorGauss, mu_ls, sigma_ls^2};
    pg_os = {@priorGauss, mu_os, sigma_os^2};
    prior.cov = {pg_ls, pg_os};
    
%     for i=1:parms.nfirm
%         mu_vf = log(0.1); sigma_vf = log(2)/2;
%         pg_vf = {@priorGauss, mu_vf, sigma_vf^2};
%         prior.cov{2+i} = pg_vf;
%     end

    % total variance: CI for unknown sigma (chi-squared)
    mu_lik = log(0.05); sigma_lik = log(2)/2;
%     est_var = var(y); 
%     sigma_low = sqrt((n-1)*est_var/chi2inv(0.975,n));
%     sigma_high = sqrt((n-1)*est_var/chi2inv(0.025,n));
%     mu_lik = log(sqrt(est_var)); sigma_lik = (log(sigma_low) - log(sigma_high))/4;
        
    pg_lik = {@priorGauss, mu_lik, sigma_lik^2};
    prior.lik = {pg_lik};

    im = {@infPrior, inffunc, prior};
    par = {meanfunc,covfunc,likfunc, x, y};
    mfun = @minimize_v2;

    % minimize hyperparameters
    if ~parms.mode
        ft = zeros(1,n);
        s2t = zeros(1,n);
        nlZt = zeros(1,n);
        for i = 1:n
            xt = x(1:i,:);
            yt = y(1:i);
            % define hyperparameters and hyperpriors
            hyp.cov = [0, 0]; hyp.mean = mean(yt); hyp.lik = log(1);
            par = {meanfunc,covfunc,likfunc, xt, yt};
            hyp = feval(mfun, hyp,  @gp, -100, im, par{:}); 
            [~, ~, mu, s2] = gp(hyp, inffunc, meanfunc, covfunc, likfunc, xt, yt, [0,0,1]);
            [nlZ, ~] = gp(hyp, inffunc, meanfunc, covfunc, likfunc, x, y);
            ft(i) = mu; s2t(i) = s2; nlZt(i) = nlZ;
        end

        f = [ft+2*sqrt(s2t); flip(ft-2*sqrt(s2t),1)];
        fill([ts.'; flip(ts.',1)], f, [7 7 7]/8);
        hold on; plot(ts, ft);
    else
        % define hyperparameters
        % hyp.cov = [mu_ls; mu_os];
        % hyp.mean = [mu_ml; mu_mc];
        % hyp.lik = [mu_lik];

        bestnlZ = 0;
        for i=1:5
            hyp = sample_prior(prior);
            hyp = feval(mfun, hyp,  @gp, -100, im, par{:});
            [nlZ, ~] = gp(hyp, inffunc, meanfunc, covfunc, likfunc, x, y);
            if nlZ < bestnlZ, bestnlZ = nlZ; besthyp = hyp; end
        end

        % we do not care about ymu, ys2
        [~, ~, fmu, fs2] = gp(besthyp, im, meanfunc, covfunc, likfunc, x, y, xs);

        % plot data and credibal interval of underlying proportion without noise
        fig = plot_posterior(fmu, fs2, x(:,1), y, xs(:,meanmask), parms.i);

        p = fmu(end);
    end

end