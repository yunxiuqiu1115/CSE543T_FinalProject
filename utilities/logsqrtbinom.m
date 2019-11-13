function A = logsqrtbinom(hyp, x, i)
% Constant mean function. The mean function is parameterized as:
%
% m(x) = p(1-p)/n
% 
% x = [p, n]
%
% The hyperparameter is:
%
% hyp = []
%
    if nargin<2, A = '0'; return; end             % report number of hyperparameters 
    if ~isempty(hyp), error('No hyperparameter needed.'), end
    if nargin==2
      ps = x(:,1); ns = x(:,2);
      A = log(sqrt(ps.*(1-ps)./ns));              % evaluate mean
    else
      error('No derivative')                      % derivative
    end
end