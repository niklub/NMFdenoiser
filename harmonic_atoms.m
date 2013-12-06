function [atoms, D] = harmonic_atoms(f0range, sr, window, num_atoms_per_f0, max_harms, fstep)

N = numel(window);
nfft = 2^nextpow2(N);
scale = 2/sum(window);

if numel(f0range)>2
    f0 = f0range;
else
    f0min = f0range(1);
    f0max = f0range(2);
    if nargin<6
        fstep = 5;
    end
    f0 = f0min:fstep:f0max;
end

natom = 0;
atoms = {};
L = numel(f0)*num_atoms_per_f0;
D = zeros(nfft/2+1, L);
for l = 1:numel(f0)
    for j = 1:num_atoms_per_f0
        natom = natom+1;
        p = min(floor(.5*sr/f0(l)), max_harms);
        atoms{natom}.name = 'harmonic';
        atoms{natom}.f0 = f0(l);
        atoms{natom}.z = pulse_model((1:p)*f0(l)/sr);
        atoms{natom}.Psi = zeros(nfft/2+1, p);
        for k = 1:p
            signal = cos(2*pi*k*f0(l)/sr*(0:N-1))';
            spectrum = abs(fft(signal.*window, nfft));
            spectrum = spectrum(1:nfft/2+1)*scale;
            spectrum(1) = spectrum(1)/2;
            spectrum(end) = spectrum(end)/2;
            atoms{natom}.Psi(:,k) = spectrum;
        end
        dj = atoms{natom}.Psi*atoms{natom}.z;
        %subplot(1,2,1), plot(dj);
        %subplot(1,2,2), plot(atoms{natom}.z);
        %pause(0.01);
        D(:, natom) = dj;
    end
end
end

function a = pulse_model(f)
a = sinc(f).^2;
if size(a,1)==1
    a = a';
end
a = a + .01*randn(size(a,1),1);
a(a<0) = eps;
a(a>1) = 1;
end