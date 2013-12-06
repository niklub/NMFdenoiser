function [W] = noise_reduction_NMF(noisy_magnitude, noise_matrix, signal_components)
%"Reduction of Non-stationary Noise for a Robotic Living Assistant using
%Sparse Non-negative Matrix Factorization"

opts.conv_value = 0;
opts.max_iter = 25;
opts.beta = 1;
if nargin<3
    ns = size(noise_matrix, 2);
else
    ns = signal_components;
end
opts.m = ns + size(noise_matrix, 2);
opts.D0 = [rand(size(noisy_magnitude,1), ns), noise_matrix];
opts.updateD = 1:ns;
opts.X0 = rand(size(opts.D0, 2), size(noisy_magnitude, 2));
opts.lambda = [.2*ones(ns,1); eps*ones(opts.m-ns,1)];
[Dout, Xout] = simplestNMF(noisy_magnitude, opts);

%create Wiener filter
Ys = Dout(:,1:ns)*Xout(1:ns,:);
Yn = Dout(:,(ns+1):end)*Xout((ns+1):end,:);
W = Ys ./ (Ys + Yn);
end