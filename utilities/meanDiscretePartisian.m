function A = meanDiscretePartisian(s, hyp, x, i)

% Mean function for discrete inputs x. Given a function defined on the
% integers 1,2,3,..,s, the mean function is parametrized as:
%
% m(x) = mu_x if x(:,2) == 1
% m(x) = -mu_x if x(:,2) == -1
%
% where mu is a fixed vector of length s.
%
% This implementation assumes that the inputs x are given as integers
% between 1 and s, which simply index the provided vector.
%
% The hyperparameters are:
%
% hyp = [ mu_1
%         mu_2
%         ..
%         mu_s ]
%

if nargin==0, error('s must be specified.'), end           % check for dimension
if nargin<=2, A = num2str(s); return; end     % report number of hyperparameters
mu = hyp(:);
flag = x(:,2);
if nargin==3
  A = mu(x(:,1)).*flag;                                                  % evaluate mean
else
  A = zeros(size(x,1),1);                                            % derivative
  A(x(:,1)==i) = flag(1);
end