function SNR = getSNR(clean_file, noisy_file)
clean_audio = wavread(clean_file);
noisy_audio = wavread(noisy_file);
n = min(numel(clean_audio), numel(noisy_audio));
SNR =  20*log10(norm(clean_audio(1:n))/norm(clean_audio(1:n)-noisy_audio(1:n)));
end