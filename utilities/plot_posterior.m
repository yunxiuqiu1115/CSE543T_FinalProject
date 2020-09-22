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
    
    % pre-forcasting posterior mean and credible intervals
    f = [fmu(pre)+2*sqrt(fs2(pre)); flip(fmu(pre)-2*sqrt(fs2(pre)),1)];
    fill([xs(pre); flip(xs(pre),1)], f, [150, 150, 150] / 255, ...
     'facealpha', 0.7, ...
     'edgecolor', 'none');
 
    % post-forcasting posterior mean and credible intervals
    hold on;
    f = [fmu(future)+2*sqrt(fs2(future)); flip(fmu(future)-2*sqrt(fs2(future)),1)];
    fill([xs(future); flip(xs(future),1)], f, [200, 200, 200] / 255, ...
     'facealpha', 0.7, ...
     'edgecolor', 'none');
     plot(xs, fmu, "color", [31, 120, 180] / 255); plot(x, y, 'k.'); 
     
     % prior mean and credible intervals on intercept
     b = parms.prior(1); s = parms.prior(2);
     y = linspace(b-2*s,b+2*s, 100);
     x = -normpdf(y, b,s);
     fill(x,y, [31, 120, 180] / 255,...
     'facealpha', 0.7, ...
     'edgecolor', 'none');
     
    % scale range of polls to [0,1]
    ylim([0,1]);
    legend('95% CI pre-forecasting', '95% CI forecasting', 'mean p*','polling data','intercept prior' ,'Location', 'Best');
    xlabel("Horizon"); ylabel("Voter Preference");
end
