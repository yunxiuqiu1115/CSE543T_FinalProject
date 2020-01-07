function posttrain(CNNdata, allRaces, besthyp)
    N = 0; nsuc = 0;
    fn = fieldnames(allRaces);
    for i=1:numel(fn)
        pvs = allRaces.(fn{i});
        ps = pvs(1:2:end);
        vs = pvs(2:2:end);
        [~, p_idx] = max(ps);
        [~, t_idx] = max(vs);
        if p_idx == t_idx
           nsuc = nsuc + 1;
        end
        N = N + 1;
    end

    disp(N + " races run.");
    disp(nsuc + " successful predictions.");
    disp("Prediction rate: " + nsuc/N);

    errors = zeros(N,1);
    a = zeros(N,1);
    b = zeros(N,1);
    N = 1;
    for i=1:numel(fn)
        pvs = allRaces.(fn{i});
        ps = pvs(1:2:end);
        vs = pvs(2:2:end)/100;
%         [~, p_idx] = max(ps);
%         [~, t_idx] = max(vs);
        l = size(ps,2);
        a(N:N+l-1) = ps.';
        b(N:N+l-1) = vs.';
        errors(N:N+l-1) = (ps - vs).';
        N = N + l;
    end
    
    fig = histogram(abs(errors));
    xlabel("absolute forecasting error");
    ylabel("count");
    title("Frequency count of absolute forecasting error");
    saveas(fig, "errors.jpg");
    close;

    top = 10;
    firmsigmas = besthyp.cov(3:end);
    [~, indt] = sort(firmsigmas,'descend');
    top_ind = indt(1:top);
    top_ind(1) = [];
    
    for i=1:numel(top_ind)
        disp(unique(CNNdata(CNNdata.pollsteridx == top_ind(i),:).pollster))
        disp((firmsigmas(top_ind(i))))
    end

    bottom = 10;
    [~, indb] = sort(firmsigmas,'ascend');
    bottom_ind = indb(1:bottom);
    bottom_ind(1) = 0;
    
    for i=2:bottom
        disp(unique(CNNdata(CNNdata.pollsteridx == bottom_ind(i),:).pollster))
        disp((firmsigmas(bottom_ind(i))))
    end

    top = 10;
    nfirm = size(besthyp.cov,1) - 2;
    firmbiases = besthyp.mean(end-nfirm+1:end);
    [~, indt] = sort(firmbiases,'descend');
    top_ind = indt(1:top);

    for i=1:top
        disp(unique(CNNdata(CNNdata.pollsteridx == top_ind(i),:).pollster))
    end

    bottom = 10;
    [~, indb] = sort(firmbiases,'ascend');
    bottom_ind = indb(1:bottom);
    bottom_ind(1) = [];

    for i=1:numel(bottom_ind)
        disp(unique(CNNdata(CNNdata.pollsteridx == bottom_ind(i),:).pollster))
    end

    % histogram(exp(firmsigmas));

    disp("Forecasting and actual voting correlation is: " + corr(a,b));
    fig = scatter(a,b, 'k.');
    xlabel("forecasted votes");
    ylabel("actual votes");
    title("Forecast and actual votes correlation");
    saveas(fig, "corr.jpg");
    close;
    
    
%     cycle = cell(760,1);
%     state = cell(760,1);
%     candidate = cell(760,1);
%     posteriormean = cell(760,1);
%     posteriorstd = cell(760,1);
%     for i=1:760
%         cycle{i} = raceinfos{i}{1};
%         state{i} = raceinfos{i}{2};
%         candidate{i} = raceinfos{i}{3};
%         posteriormean{i} = fts(i);
%         posteriorstd{i} = sqrt(s2s(i));
%     end
%     forecast = table(cycle, state, candidate, posteriormean, posteriorstd);
%     writetable(forecast,'forecast1992-2014.csv');

%     validforecast = cell(70, 3+6*2);
%     
%     for i=1:70
%         validforecast{i,1} = validraceinfos{i}{1};
%         validforecast{i,2} = validraceinfos{i}{2};
%         validforecast{i,3} = validraceinfos{i}{3};
%         for j=1:6
%             validforecast{i,3+2*j-1} = validfts{i}(j);
%             validforecast{i,3+2*j} = sqrt(valids2s{i}(j));
%         end
%     end
%     validforecast = array2table(validforecast,'VariableNames',...
%                         {'cycle','state','candidate', 'day90mean', 'day90std',...
%                          'day60mean', 'day60std', 'day30mean', 'day30std',...
%                          'day14mean', 'day14std', 'day7mean', 'day7std',...
%                          'day1mean', 'day1std'});
%     writetable(validforecast,'forecast2016.csv');
end