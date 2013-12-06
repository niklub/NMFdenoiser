function [D, X, Err, atoms] = linearDenseNMF(Y, atoms0, opts)
if nargin<3
    opts = {};
end
[n, T] = size(Y);
opts = initialize_opts(Y, atoms0, opts);
atoms = opts.atoms0;
D = opts.D0;
X = opts.X0;
lambda_factor = repmat(opts.lambda.*ones(opts.m,1), 1, T);
reg_factor = T*opts.alpha.*ones(opts.m, 1);

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
    pgradD = DX.^(opts.beta-1)*X';
    ngradD = (DX.^(opts.beta-2)).*Y*X';
    
    for j = 1:opts.m
        Psi_j = atoms{j}.Psi;
        a_j = atoms{j}.a;
        
        nom = a_j'*Psi_j'*pgradD(:,j) + Psi_j'*ngradD(:,j) + reg_factor(j)*(a_j'*a_j);
        den = Psi_j'*pgradD(:,j) + a_j'*Psi_j'*ngradD(:,j) + reg_factor(j)*a_j;
        
        den(den==0) = eps;
        
        atoms{j}.a = a_j .* nom ./ den;
        atoms{j}.a = atoms{j}.a / sum(atoms{j}.a);
    end
    for j = 1:opts.m
        D(:, j) = atoms{j}.Psi * atoms{j}.a;
    end
    Err(iter) = beta_divergence(Y, D*X, opts.beta);
    %sparsy = sum(sum(X));
    %smooth = .5*smooth_atoms(atoms)*T;
    delta = inf;
    if iter>1
        delta = (Err(iter-1)-Err(iter))/Err(iter-1);
    end
    %imagesc(D);
    %pause(0.01);
    %fprintf('iter = %d, Err = %f, delta = %f, sparsy = %f, smooth = %f\n', iter, Err(iter), delta, sparsy, smooth);
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

if ~isfield(opts, 'conv_value') opts.conv_value = 1e-3; end
if ~isfield(opts, 'max_iter') opts.max_iter = 1000; end
if ~isfield(opts, 'lambda') opts.lambda = 0; end
if ~isfield(opts, 'beta') opts.beta = 1; end
if ~isfield(opts, 'alpha') opts.alpha = 1; end
if ~isfield(opts, 'X0') opts.X0 = rand(opts.m, T); end

for j = 1:opts.m
    p = size(opts.atoms0{j}.Psi, 2);
    opts.atoms0{j}.a = rand(p, 1);
    opts.atoms0{j}.a = opts.atoms0{j}.a / sum(opts.atoms0{j}.a);
    opts.D0(:, j) = opts.atoms0{j}.Psi * opts.atoms0{j}.a; 
end
end
