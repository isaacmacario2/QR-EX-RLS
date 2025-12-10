function L = Cholesky(P)
    % Check if P is symmetrical
    if ~isequal(P, P')
        error('Matrix P is not symmetric.');
    end

    n = size(P,1);
    L = zeros(n);

    for i = 1:n
        for j = 1:i
            sum = 0;
            for k = 1:j-1
                sum = sum + L(i,k)*L(j,k);
            end

            if i == j
                val = P(i,i) - sum;
                if val <= 0
                    error('The matrix is ​​not positive definite.');
                end
                L(i,j) = sqrt(val);
            else
                L(i,j) = (1 / L(j,j)) * (P(i,j) - sum);
            end
        end
    end
end
