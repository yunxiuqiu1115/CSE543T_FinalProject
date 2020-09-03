function fig = plot_posterior(fmu, fs2, x, y, xs, v, parms)
    fig = figure('visible', 'off');
%     fig = figure(1);
    tau = parms.tau;
    future = (xs>=-tau);
    pre = (xs<-tau);
    f = [fmu(pre)+2*sqrt(fs2(pre)); flip(fmu(pre)-2*sqrt(fs2(pre)),1)];
    fill([xs(pre); flip(xs(pre),1)], f, [150, 150, 150] / 255, ...
     'facealpha', 0.7, ...
     'edgecolor', 'none');
 
    hold on;
    f = [fmu(future)+2*sqrt(fs2(future)); flip(fmu(future)-2*sqrt(fs2(future)),1)];
    fill([xs(future); flip(xs(future),1)], f, [200, 200, 200] / 255, ...
     'facealpha', 0.7, ...
     'edgecolor', 'none');
     plot(xs, fmu, "color", [31, 120, 180] / 255); plot(x, y, 'k.'); 
     
     b = parms.prior(1); s = parms.prior(2);
     y = linspace(b-2*s,b+2*s, 100);
     x = -normpdf(y, b,s);
     fill(x,y, [31, 120, 180] / 255,...
     'facealpha', 0.7, ...
     'edgecolor', 'none');
     
     
    ylim([0,1]);
    % plot(0,v,'b*'); 
%     plot(xs,a*xs+b);
%     scatter(x, y, 16, pb, 'filled');
%     c = jet;
%     colormap(c); colorbar;
%     legend('95% CI','mean p*','data','actual vote','linear trend','pollster biases','Location', 'Best');
    legend('95% CI pre-forecasting', '95% CI forecasting', 'mean p*','polling data','intercept prior' ,'Location', 'Best');
    xlabel("Days Before Election"); ylabel("Voter Preference");
end

% fill([xstar(:,1); flip(xstar(:,1),1)], f, [166, 206, 227] / 255, ...
%      'facealpha', 0.7, ...
%      'edgecolor', 'none');
% hold on; plot(xstar(:,1), fmu, "color", [31, 120, 180] / 255);
% f1 = [fmu1+2*sqrt(fs21); flip(fmu1-2*sqrt(fs21),1)];
%     fill([xstar(:,1); flip(xstar(:,1),1)], f1, [227,166, 206] / 255, ...
%      'facealpha', 0.7, ...
%      'edgecolor', 'none');
%     hold on; plot(xstar(:,1), fmu1, "color", [180,31, 120] / 255);
% legend('95% CI Blunt','mean p* Blunt','95% CI Kander','mean p* Kander','Location', 'Best');
%     xlabel("days before election"); ylabel("voter preference");
% plot_title = year + " " + state;
% title(plot_title);

% scatter(xs{802}(:,1),xs{803}(:,2)-xs{802}(:,2),'k*');
