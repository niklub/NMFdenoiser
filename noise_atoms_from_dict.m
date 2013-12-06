function [atoms, D] = noise_atoms_from_dict(Dinput, num_atoms, N)

if nargin<3
    N = 2*(size(Dinput,1)-1);
end
p = size(Dinput, 2);
atoms = cell(1,num_atoms);
D = zeros(size(Dinput,1), num_atoms);
for j = 1:num_atoms
    
    atoms{j}.name = 'noise';
    full_spectrum = [D(1:end-1, :); flipdim(D(2:end, :), 1)];
    atoms{j}.signal = real(ifft(full_spectrum.*exp(1i*2*pi*rand(size(full_spectrum)))));
    atoms{j}.signal = atoms{j}.signal(1:N, :);
    atoms{j}.z = rand(p, 1);
    atoms{j}.Psi = Dinput;
    dj = atoms{j}.Psi*atoms{j}.z;
    %subplot(1,2,1), plot(dj);
    %subplot(1,2,2), plot(atoms{j}.z);
    %pause(0.01);
    D(:, j) = dj;
end