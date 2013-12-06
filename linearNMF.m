function [D, X, Err, atoms] = linearNMF(Y, atoms0, opts)
[n, T] = size(Y);
opts = initialize_opts(Y, atoms0, opts);

atoms = opts.atoms0;
D = opts.D0;
X = opts.X0;
lambda_factor = repmat(opts.lambda.*ones(opts.m,1), 1, T);

iter = 0;
Err = zeros(opts.max_iter, 1);
while iter < opts.max_iter
    iter = iter + 1;
    
    %update X
    DX = D*X;
    if opts.beta<2
        DX(DX==0) = eps;
    end
    pgradX = D'*(DX.^(opts.beta-1));
    ngradX = D'*((DX.^(opts.beta-2)).*Y);
    X = X .* ngradX ./ (pgradX + lambda_factor);
    
    %update D
    DX = D*X;
    if opts.beta<2
        DX(DX==0) = eps;
    end
    pgradD_m = DX.^(opts.beta-1);
    ngradD_m = (DX.^(opts.beta-2)).*Y;
    
    for j = 1:opts.m
        xj = X(j,:)';
        nom = atoms{j}.Psi' * ngradD_m * xj;
        denom = atoms{j}.Psi' * pgradD_m * xj;
        denom(denom==0) = eps;
        atoms{j}.a = atoms{j}.a .* nom ./ denom;
    end
    for j = 1:opts.m
        D(:, j) = atoms{j}.Psi * atoms{j}.a;
    end
    Err(iter) = beta_divergence(Y, D*X, opts.beta);
    %sparsy = sum(sum(lambda_factor.*X));
    delta = inf;
    if iter>1
        delta = (Err(iter-1)-Err(iter))/Err(iter-1);
    end
    %fprintf('iter = %d, Err = %f, delta = %f, sparsy = %f\n', iter, Err(iter), delta, sparsy);
    if delta < opts.conv_value
        break;
    end
end
Err(iter+1:end) = [];
end

function opts = initialize_opts(Y, atoms0, opts)
[n, T] = size(Y);
opts.m = length(atoms0);
opts.D0 = zeros(n, opts.m);
opts.atoms0 = atoms0;
for j = 1:opts.m
    p = size(opts.atoms0{j}.Psi, 2);
    opts.atoms0{j}.a = rand(p, 1);
    opts.D0(:, j) = opts.atoms0{j}.Psi * opts.atoms0{j}.a; 
end
if ~isfield(opts, 'conv_value') opts.conv_value = 1e-3; end
if ~isfield(opts, 'max_iter') opts.max_iter = 1000; end
if ~isfield(opts, 'lambda') opts.lambda = eps; end
if ~isfield(opts, 'beta') opts.beta = 2; end
if ~isfield(opts, 'X0') opts.X0 = rand(opts.m, T); end
end
