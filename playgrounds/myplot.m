% plot samples and credible intervals 

function myplot(x, mu, covariance, num)
    % sample from GP
    R = mvnrnd(mu, covariance, num);
    % plot credible intervals
    n = length(x);
    high = linspace(0,0,n);
    low = linspace(0,0,n);
    for i=1:n
        high(i) = mu(i) + 1.96*covariance(i,i);
        low(i) = mu(i) - 1.96*covariance(i,i);
    end
    fill([x, fliplr(x)], [high, fliplr(low)], [7 7 7]/8);
    hold on;
    % plot samples
    for i=1:num
        plot(x, R(i,:));
        hold on;
    end
end