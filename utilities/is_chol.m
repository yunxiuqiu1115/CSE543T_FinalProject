function [result] = is_chol(L)
  result = (ismatrix(L) &&  ...                   % is it a matrix?
            (size(L, 1) == size(L, 2)) &&  ...    % is it square?
            isreal(diag(L))  &&       ...         % is the diagonal real?
            all(diag(L) > 0) &&      ...          % is the diagonal positive?
            isequal(L, triu(L)));              % is the matrix upper triangular?
end