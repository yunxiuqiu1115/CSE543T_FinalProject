function fig = plot_posterior(fmu, fs2, x, y, xs, v, i)
    fig = figure('visible', 'off');
%     fig = figure(i);
    f = [fmu+2*sqrt(fs2); flip(fmu-2*sqrt(fs2),1)];
    fill([xs; flip(xs,1)], f, [200, 200, 200] / 255, ...
     'facealpha', 0.7, ...
     'edgecolor', 'none');
    hold on; plot(xs, fmu, "color", [31, 120, 180] / 255); plot(x, y, 'k.'); plot(0,v,'b*'); 
%     plot(xs,a*xs+b);
%     scatter(x, y, 16, pb, 'filled');
%     c = jet;
%     colormap(c); colorbar;
%     legend('95% CI','mean p*','data','actual vote','linear trend','pollster biases','Location', 'Best');
    legend('95% CI','mean p*','polling data','vote share','Location', 'Best');
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
