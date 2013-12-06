function [D,X,Err] = simplestNMF(Y, opts)
[n,T] = size(Y);

if nargin<2
    opts = struct();
end
opts = initialize_opts(Y, opts);
D = opts.D0;
X = opts.X0;

iter = 0;
Err = zeros(opts.max_iter, 1);
lambda_factor = repmat(opts.lambda.*ones(opts.m,1), 1, T);
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
    if ~isempty(opts.updateD)
        %DX = D(:, opts.updateD)*X(opts.updateD, :);
        DX = D*X;
        if opts.beta<2
            DX(DX==0) = eps;
        end
        pgradD = (DX.^(opts.beta-1))*X(opts.updateD, :)';
        pgradD(pgradD==0)=eps;
        ngradD = ((DX.^(opts.beta-2)).*Y)*X(opts.updateD, :)';
        D(:, opts.updateD) = D(:, opts.updateD) .* ngradD ./ pgradD;
    end
    
    DX = D*X;
    DX(DX==0) = eps;
    Err(iter) = beta_divergence(Y, DX, opts.beta);
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

function opts = initialize_opts(Y, opts)
[n, T] = size(Y);
if ~isfield(opts, 'm') opts.m = 1; end
if ~isfield(opts, 'conv_value') opts.conv_value = 1e-3; end
if ~isfield(opts, 'max_iter') opts.max_iter = 1000; end
if ~isfield(opts, 'lambda') opts.lambda = eps; end
if ~isfield(opts, 'beta') opts.beta = 2; end
if ~isfield(opts, 'updateD') opts.updateD = 1:opts.m; end
if ~isfield(opts, 'D0') opts.D0 = rand(n, opts.m); end
if ~isfield(opts, 'X0') opts.X0 = rand(opts.m, T); end
end