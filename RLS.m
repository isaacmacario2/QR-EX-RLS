function [y, mse, w] = RLS(X, d, lambda, init, type)

    [Nw, N] = size(X);
    y = zeros(N, 1);  % Previsão da saída
    mse = zeros(N,1);
    w = zeros(Nw,N);
    % w = zeros(Nw,1);
    lambda_inv = 1/lambda;

    Sd = init*eye(Nw);

    for n=1:N
        if type == 0
            e = d(n) - X(:,n)'*w(:,n);
        elseif type == 1 
            e = d(n) - w(:,n)'*X(:,n);
        end
        psi = Sd*X(:,n);
        Sd = lambda_inv*(Sd-(psi*psi')/(lambda+psi'*X(:,n)));
        vetk = Sd*X(:,n);
        if type == 0
            w(:,n+1) = w(:,n) + e*vetk;
        elseif type == 1
            w(:,n+1) = w(:,n) + conj(e)*vetk;
        end
        y(n) = w(:,n)'*X(:,n);
        mse(n) = abs(e)^2;

    end 

    w = w(:,end);

end