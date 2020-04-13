function posttrain(raceinfos, fts, s2s, allRaces, besthyp)
    N_train = 0; nsuc_train = 0;
    N_test = 0; nsuc_test = 0;
    fn = fieldnames(allRaces);
    for i=1:numel(fn)
        pvs = allRaces.(fn{i});
        ps = pvs(1:2:end);
        vs = pvs(2:2:end);
        [~, p_idx] = max(ps);
        [~, t_idx] = max(vs);
        year = fn{i}(end-3:end);
        if strcmp(year, '2016')==1
            N_test = N_test + 1;
            if p_idx == t_idx
                nsuc_test = nsuc_test + 1;
            else   
               misclassify(fn{i}, vs, ps, 'test');
            end
        else
            N_train = N_train + 1;
            if p_idx == t_idx
                nsuc_train = nsuc_train + 1;
            else
                % misclassify(fn{i}, vs, ps, 'train');
            end
        end
    end
    
    N = N_train + N_test;
    nsuc = nsuc_train + nsuc_test;
    
    accuracy(N_train, nsuc_train, 'training');
    accuracy(N_test, nsuc_test, 'test');
    accuracy(N, nsuc, 'overall');

    N = numel(raceinfos);
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
    
    disp('correlation of predictive mean and actual vote on training data');
    disp(corr(a(1:760),b(1:760)));
    disp('RMSE of predictive mean and actual vote on training data');
    disp(sqrt(mean((a(1:760)-b(1:760)).^2)));
    
    disp('correlation of predictive mean and actual vote on test data');
    disp(corr(a(761:end),b(761:end)));
    disp('RMSE of predictive mean and actual vote on test data');
    disp(sqrt(mean((a(761:end)-b(761:end)).^2)));
    
%     fig = histogram(abs(errors));
%     xlabel("absolute forecasting error");
%     ylabel("count");
%     title("Frequency count of absolute forecasting error");
%     saveas(fig, "errors.jpg");
%     close;

%     top = 10;
%     firmsigmas = besthyp.cov(3:end);
%     [~, indt] = sort(firmsigmas,'descend');
%     top_ind = indt(1:top);
%     top_ind(1) = [];
%     
%     for i=1:numel(top_ind)
%         disp(unique(CNNdata(CNNdata.pollsteridx == top_ind(i),:).pollster))
%         disp((firmsigmas(top_ind(i))))
%     end
% 
%     bottom = 10;
%     [~, indb] = sort(firmsigmas,'ascend');
%     bottom_ind = indb(1:bottom);
%     bottom_ind(1) = 0;
%     
%     for i=2:bottom
%         disp(unique(CNNdata(CNNdata.pollsteridx == bottom_ind(i),:).pollster))
%         disp((firmsigmas(bottom_ind(i))))
%     end
% 
%     top = 10;
%     nfirm = size(besthyp.cov,1) - 2;
%     firmbiases = besthyp.mean(end-nfirm+1:end);
%     [~, indt] = sort(firmbiases,'descend');
%     top_ind = indt(1:top);
% 
%     for i=1:top
%         disp(unique(CNNdata(CNNdata.pollsteridx == top_ind(i),:).pollster))
%     end
% 
%     bottom = 10;
%     [~, indb] = sort(firmbiases,'ascend');
%     bottom_ind = indb(1:bottom);
%     bottom_ind(1) = [];
% 
%     for i=1:numel(bottom_ind)
%         disp(unique(CNNdata(CNNdata.pollsteridx == bottom_ind(i),:).pollster))
%     end

    % histogram(exp(firmsigmas));

%     disp("Forecasting and actual voting correlation is: " + corr(a,b));
%     fig = scatter(a,b, 'k.');
%     xlabel("forecasted votes");
%     ylabel("actual votes");
%     title("Forecast and actual votes correlation");
%     saveas(fig, "corr.jpg");
%     close;
    
    % get only 1992 - 2014
    N = length(raceinfos);
    cycle = cell(N,1);
    state = cell(N,1);
    candidate = cell(N,1);
    pvi = cell(N,1);
    party = cell(N,1);
    experienced = cell(N,1);
    posteriormean = cell(N,1);
    posteriorstd = cell(N,1);
    vote = cell(N, 1);
    Nout_train = 0;
    Nout_test = 0;
    nlZ = [];
    for i=1:N
        cycle{i} = raceinfos{i}{1};
        state{i} = raceinfos{i}{2};
        candidate{i} = raceinfos{i}{3};
        posteriormean{i} = fts(i);
        posteriorstd{i} = sqrt(s2s(i));
        vote{i} = raceinfos{i}{4};
        pvi{i} = raceinfos{i}{5};
        experienced{i} = raceinfos{i}{6};
        party{i} = raceinfos{i}{7};
        
        u = posteriormean{i} + 1.96*posteriorstd{i};
        l = posteriormean{i} - 1.96*posteriorstd{i};
        v = vote{i}/100;
        
        if cycle{i}==2016
            nlZ = [nlZ, (v-fts(i))^2/2/s2s(i) + log(s2s(i))/2 + log(2*pi)/2];
            if v > u || v < l
%                disp("Posterior for "+cycle{i}+" "+state{i}+" out of 95% CI");
%                disp("95%CI: [" + l + ", " + u + "]" );
%                disp("Actual vote: "+v);
               Nout_test = Nout_test + 1;
            end   
        else
            if v > u || v < l
               Nout_train = Nout_train+1;
            end
        end
    end
    disp(Nout_train + "/760 actual votes out of 95% predictive CI on training data.");
    disp(Nout_test + "/73 actual votes out of 95% predictive CI on test data.");
    disp("Test Average nlZ: " + mean(nlZ));
%     forecast = table(cycle, state, candidate, posteriormean, posteriorstd, vote, pvi, party, experienced);
%     writetable(forecast,'results/forecast1992-2016old14.csv');

    disp("Length Scale: " + exp(besthyp.cov(1)));
    disp("Output Scale: " + exp(besthyp.cov(2)));
    disp("Noise std: " + exp(besthyp.lik));
    
    posteriorstd = cell2mat(posteriorstd);
    disp("Mean of predictive std: " + mean(posteriorstd));
    disp("Median of predictive std: " + median(posteriorstd));
    disp("Std of predictive std: " + std(posteriorstd));
end

function acc = accuracy(N, nsuc, title)
    disp(N + " " + title + " races run.");
    disp(nsuc + " successful " + title + " predictions.");
    acc = nsuc/N;
    disp(title + " accuracy rate: " + acc); 
end

function misclassify(name, v, p, title)
    disp("Inaccurate " + title + " prediction.")
    disp(name);
    disp("Actual voting: " + v);
    disp("Predictive mean: " + p);
%     u = p+1.96*sqrt(s2);
%     l = p-1.96*sqrt(s2);
%     disp("95%CI: [" + l + ", " + u + "]" );
end