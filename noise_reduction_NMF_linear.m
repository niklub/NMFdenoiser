function [W] = noise_reduction_NMF_linear(noisy_magnitude, atomsH, atomsN, params)

if nargin < 4
    params = struct();
end
params = initParams(params);
ns = length(atomsH);
atoms0 = {atomsH{:}, atomsN{:}};
params.m = length(atoms0);
params.lambda = [params.speech_sparsity*ones(ns,1); eps*ones(params.m-ns,1)];
switch params.type
    case 'linNMF'
        [Dout, Xout] = linearNMF(noisy_magnitude, atoms0, params);
    case 'denseNMF'
        [Dout, Xout] = linearDenseNMF(noisy_magnitude, atoms0, params);
end

%create Wiener filter
Ys = Dout(:,1:ns)*Xout(1:ns,:);
Yn = Dout(:,(ns+1):end)*Xout((ns+1):end,:);
W = Ys ./ (Ys + Yn);
end

function pout = initParams(p)
pout = p;

if ~isfield(pout, 'alpha') pout.alpha = .2; end
if ~isfield(pout, 'speech_sparsity') pout.speech_sparsity = .2; end
if ~isfield(pout, 'conv_value') pout.conv_value = 0; end
if ~isfield(pout, 'max_iter') pout.max_iter = 25; end
if ~isfield(pout, 'beta') pout.beta = 1; end
if ~isfield(pout, 'type') pout.type = 'denseNMF'; end
end