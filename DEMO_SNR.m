%================INPUT FILES====================
noisy_speech    = 'test/speech+noise_8khz.wav';
noise           = 'test/noise_8khz.wav';
clean_speech    = 'test/speech_8khz.wav';

%================OUTPUT FILES====================
denoised_ss_nmf     = 'test/denoised_8khz_SS_NMF.wav';
denoised_ss_linnmf  = 'test/denoised_8khz_SS_linNMF.wav';
denoised_ss_densenmf= 'test/denoised_8khz_SS_denseNMF.wav';

denoised_us_nmf     = 'test/denoised_8khz_US_NMF.wav';
denoised_us_linnmf  = 'test/denoised_8khz_US_linNMF.wav';
denoised_us_densenmf= 'test/denoised_8khz_US_denseNMF.wav';

%================INPUT PARAMS==================
params.nwin = 256;
params.show_log = false;

%================PROCESSING====================
fprintf('The input file %s has SNR = %.2f dB\n',...
    noisy_speech, getSNR(clean_speech, noisy_speech));

%-------------SEMI-SUPERVISED---------------
fprintf('\n================\n');
fprintf('First start computing semi-supervised (SS) denoisers...\n');
params.noise    = noise;

%----------NMF-------------
fprintf('Evaluating simple NMF...');
params.type     = 'NMF';
params.output   = denoised_ss_nmf;
NMFdenoiser(noisy_speech, params);
fprintf(' Done! SNR = %.2f dB\n', getSNR(clean_speech, params.output));

%----------linNMF-------------
fprintf('Evaluating linNMF...');
params.type     = 'linNMF';
params.output   = denoised_ss_linnmf;
NMFdenoiser(noisy_speech, params);
fprintf(' Done! SNR = %.2f dB\n', getSNR(clean_speech, params.output));

%----------denseNMF-----------
fprintf('Evaluating denseNMF...');
params.type     = 'denseNMF';
params.output   = denoised_ss_densenmf;
NMFdenoiser(noisy_speech, params);
fprintf(' Done! SNR = %.2f dB\n', getSNR(clean_speech, params.output));


%-------------UNSUPERVISED---------------
fprintf('\n================\n');
fprintf('Then try to solve the unsupervised (US) problem...\n');
params.noise = '';

%----------NMF-------------
fprintf('Evaluating simple NMF...');
params.type     = 'NMF';
params.output   = denoised_us_nmf;
NMFdenoiser(noisy_speech, params);
fprintf(' Done! SNR = %.2f dB\n', getSNR(clean_speech, params.output));

%----------linNMF-------------
fprintf('Evaluating linNMF...');
params.type     = 'linNMF';
params.output   = denoised_us_linnmf;
NMFdenoiser(noisy_speech, params);
fprintf(' Done! SNR = %.2f dB\n', getSNR(clean_speech, params.output));

%----------denseNMF-----------
fprintf('Evaluating denseNMF...');
params.type     = 'denseNMF';
params.output   = denoised_us_densenmf;
NMFdenoiser(noisy_speech, params);
fprintf(' Done! SNR = %.2f dB\n', getSNR(clean_speech, params.output));

fprintf('\n================\n');
fprintf('Please check all output examples:\n%s\n%s\n%s\n%s\n%s\n%s\n',...
    denoised_ss_nmf, denoised_ss_linnmf, denoised_ss_densenmf,...
    denoised_us_nmf, denoised_us_linnmf, denoised_us_densenmf);