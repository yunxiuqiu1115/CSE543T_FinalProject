function fig = plot_posterior(fmu, fs2, x, y, xs, parms)

%
% Plot the posterior mean and credible intervals of underlying voter preference.
% Plot polling data as black dots.
% Plot prior model on intercept.
%

    % define figure object
    fig = figure('visible', 'off');
    tau = parms.tau;
    future = (xs>=-tau);
    pre = (xs<-tau);
    
    
     % prior mean and credible intervals on intercept
     b = parms.prior(1); s = parms.prior(2);
     yp = linspace(b-2*s,b+2*s, 100);
     xp = -normpdf(yp, b,s);
     xp  = xp - max(xp);
     fill(5*xp,yp, [80, 200, 245] / 255,...
     'facealpha', 0.7, ...
     'edgecolor', 'none');
 
  % pre-forcasting posterior mean and credible intervals
   hold on;
    f = [fmu(pre)+2*sqrt(fs2(pre)); flip(fmu(pre)-2*sqrt(fs2(pre)),1)];
    fill([xs(pre); flip(xs(pre),1)], f, [150, 150, 150] / 255, ...
     'facealpha', 0.7, ...
     'edgecolor', 'none');
 
 
    % post-forcasting posterior mean and credible intervals
    f = [fmu(future)+2*sqrt(fs2(future)); flip(fmu(future)-2*sqrt(fs2(future)),1)];
    fill([xs(future); flip(xs(future),1)], f, [200, 200, 200] / 255, ...
     'facealpha', 0.7, ...
     'edgecolor', 'none');
    plot(xs, fmu, "color", [31, 120, 180] / 255); plot(x, y, 'k.'); 
     
    % scale range of polls to [0,1]
    ylim([0,1]);
    BIN = parms.BIN;
    Nx = (abs(min(xs))/BIN);
    XTICK = -BIN*[Nx:-1:0];
    XTICKLABELS = cell(numel(XTICK),1);
    for i=1:numel(XTICK)
        XTICKLABELS{i} = num2str(-XTICK(i));
    end
    
    YT = yticks;
    YTLABELS = cell(numel(YT),1);
    for i=1:numel(YT)
        YTLABELS{i} = num2str(100*YT(i))+"%";
    end
    
    xlim([min(XTICK) ,0]);
    
    set(gca, 'box', 'off', ...
     'tickdir', 'out', ...
     'xtick', XTICK, ...
     'xticklabels', XTICKLABELS, ...
     'ytick', YT, ...
     'yticklabels', YTLABELS);
plot(0, parms.trueVote / 100, 'rd');
legend('Intercept prior', '95% CI pre-forecasting', '95% CI forecasting', 'Posterior mean','Polling data', "Actual Vote Share");
    xlabel("Days to election"); ylabel("Latent voter preference");
end
