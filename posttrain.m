errors = [];
a = [];
b = [];
for i=1:numel(fn)
    pvs = allRaces.(fn{i});
    ps = pvs(1:2:end);
    vs = pvs(2:2:end);
    [~, p_idx] = max(ps);
    [~, t_idx] = max(vs);
    a = [a, ps];
    b = [b, vs/100];
    errors = [errors, (ps - vs/100)];
end
histogram(abs(errors));

top = 10;
firmsigmas = besthyp.cov(3:end);
[valt, indt] = sort(firmsigmas,'descend');
top_val = valt(1:top);
top_ind = indt(1:top);
worst_firms = unique(CNNdata(ismember(CNNdata.pollsteridx, top_ind),:).pollster);
disp(worst_firms);

bottom = 10;
[valb, indb] = sort(firmsigmas,'ascend');
bottom_val = valb(1:bottom);
bottom_ind = indb(1:bottom);
bottom_ind(1) = 0;
best_firms = unique(CNNdata(ismember(CNNdata.pollsteridx, bottom_ind),:).pollster);
disp(best_firms);

top = 10;
firmbiases = besthyp.mean(end-parms.nfirm:end);
[valt, indt] = sort(firmbiases,'descend');
top_val = valt(1:top);
top_ind = indt(1:top);
R_firms = unique(CNNdata(ismember(CNNdata.pollsteridx, top_ind),:).pollster);
disp(R_firms);

for i=1:top
    disp(unique(CNNdata(CNNdata.pollsteridx == top_ind(i),:).pollster))
end

bottom = 10;
[valb, indb] = sort(firmbiases,'ascend');
bottom_val = valb(1:bottom);
bottom_ind = indb(1:bottom);
D_firms = unique(CNNdata(ismember(CNNdata.pollsteridx, bottom_ind),:).pollster);
disp(D_firms);

for i=1:bottom
    disp(unique(CNNdata(CNNdata.pollsteridx == bottom_ind(i),:).pollster))
end

histogram(exp(firmsigmas));

disp(corr(a',b'));