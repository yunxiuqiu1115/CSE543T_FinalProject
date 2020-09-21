function coefs = priorModel(data, test_year)
%
%   Precompute the coefs of prior linear model for the linear trend intercept
%   The linear model has six coefs with the following order:
%       - 1
%       - pvi
%       - experienced
%       - Democrat
%       - Republican:
%       - pvi*Republican
%
    [~, rows] = unique(data(:, ['Candidateidentifier']), 'rows');
    data = data(rows, :);
    data = data(data.cycle~=test_year, :);
    y = table2array(data(:, ['Percentage_of_Vote_won_x']))/100;
    data = data(:, ["pvi", "experienced", "Democrat", "Republican"]);
    x = ones(size(data,1), 6);
    x(:, 2:5) = table2array(data);
    x(:, 6) = table2array(data(:, "pvi")).*table2array(data(:, "Republican"));
    coefs = regress(y, x);
end