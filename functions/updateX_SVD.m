function X = updateX_SVD(X, DY, L, R, B, alpha, k, options)
%UPDATEX_SVD_PROJ Update of X in GL_LRD with SVD
%   Default method: GPM. ADMM as an option.
arguments
    X double
    DY double
    L double
    R double
    B double
    alpha double
    k double
    options.solver = 'GPM';
    options.tol = 1e-2;
    options.maxIter = 1000;
end
if strcmp(options.solver, 'GPM')
    isConverge = false;
    isMaxIter = false;
    iter = 1;
    D = @(X) X - R*X*B;
    targetFun = @(X) 1/2*(norm(DY - D(X), 'fro'))^2 + alpha*(trace((D(X))'*L*D(X)));
    gradX = @(X) 2*alpha*(L*D(X)-R*L*D(X)*B') - (DY - D(X) - R*(DY - D(X))*B');
    while ~isConverge && ~isMaxIter
        X_old = X;
        % Linesearch, with armijo
        currentGrad = gradX(X);
        X = lineSearchArminjo(X, currentGrad, targetFun, 1e-4, 1000);
        % k-rank throttling
        X = singularValueThrottling(X, k);
        % Convergence check
        isConverge = norm(X_old - X, 'fro')/norm(X_old, 'fro') < options.tol;
        isMaxIter = iter >= options.maxIter;
        iter = iter + 1;
    end

elseif strcmp(options.solver, 'ADMM')
    isADMMConverge = false;
    isMaxIter = false;

    D = @(X) X - R*X*B;

    iter = 1;
    Phi = X;
    La = X - Phi;
    rho = 1;
    while ~isADMMConverge && ~isMaxIter
        X_old = X;
        % Update of X
        isXConverge = false;
        while ~isXConverge
            X_o = X;
            targetFun = @(X) 1/2*(norm(DY - D(X), 'fro'))^2 + alpha*(trace((D(X))'*L*D(X))) + ...
                trace(La'*(X - Phi)) + rho/2*(norm(X - Phi, 'fro'))^2;
            gradX = @(X) 2*alpha*(L*D(X)-R*L*D(X)*B') - (DY - D(X) - R*(DY - D(X))*B') + ...
                La + rho*(X - Phi);
            X = lineSearchArminjo(X, gradX(X), targetFun, 1e-4, 100);
            isXConverge = norm(X - X_o, 'fro') < 1e-3;
        end
        % Update of Phi
        Phi = singularValueThrottling(X, k);
        % Update of Lambda
        La = La + rho*(X - Phi);
        % Convergence check
        isADMMConverge = norm(X_old - X, 'fro')/norm(X_old, 'fro') < options.tol;
        isMaxIter = iter >= options.maxIter;
        iter = iter + 1;
    end

else
    error('%s is not a vaild solver', options.solver);
end

end

