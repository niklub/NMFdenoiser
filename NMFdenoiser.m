function [output_file] = NMFdenoiser(input_file, params)
%NMFdenoiser is able to perform audio noise reduction by using several
%Non-negative Matrix Factorization (NMF)-based algorithms.

%For further details please refer to the following paper:
%N. Lyubimov, M. Kotov "Non-negative Matrix Factorization with Linear
%Constraints for Single-Channel Speech Enhancement"
%(http://arxiv.org/abs/1309.6047)

%================USAGE==============================
%output_file = NMFdenoiser(input_file [,params])

%================INPUT PARAMETERS====================
%input_file:                    path to file that should be denoised
%params:                        (optional) input parameters (see below)
%   params.output ('')          output file name (set '' to produce
%                               <input_file>_DENOISED.wav
%   params.nwin (1024)          STFT window size (1024)
%   params.type ('denseNMF')    NMF type to perform:
%                               - 'NMF'     -standart NMF (similar to )
%                               - 'linNMF'  -linear NMF
%                               - 'denseNMF'-dense NMF
%   params.noise ('')           path to file containing noise profile or
%                               empty string if noise profile should 
%                               be estimated from target file
%   params.beta (1)             beta-divergence used in NMF estimations
%   params.m    (32)            number of NMF atoms to estimate
%   params.f0range([80 400])    the range in Hz of expected F0
%   params.fstep (10)           frequency step in Hz
%   params.num_atoms_per_f0 (4) number of harmonic atoms per each F0
%   params.num_noise_atoms (16) number of noise atoms
%   params.max_harms (30)       number of harmonics per each F0
%   params.show_log (true)      show messages
%   params.alpha(.2)            parameter controls envelope smoothness
%   params.energyThr (-40)      logRMS threshold in dB used in energy-based VAD
%   params.speech_sparsity(.2)  parameter controls speech sparsity

if nargin<2
    params = struct();
end

[audio, sr] = wavread(input_file);
params = initParams(params, sr);
if params.show_log
    fprintf('Processing input file %s...\n', input_file);
end
Saudio = m_STFT(audio, sr, params);
Yaudio = abs(Saudio);

if isempty(params.noise)
    if params.show_log
        fprintf('Extracting noise profile from target file...\n');
    end
    if exist('vadsohn')==2
        vad = vadsohn(audio, sr);
    else
        fprintf('WARNING! voicebox is not installed. The results could be inaccurrate\n');
        vad_params.threshold = params.energyThr;
        vad = energyVAD(audio, vad_params);
    end
    noise_profile = audio(vad==0);
    if isempty(noise_profile)
        fprintf('Unable to find noise profile with prespecified params.energyThr. Using first samples instead\n');
        noise_profile = audio(1:round(numel(audio)/10));
    end
        
else
    if params.show_log
        fprintf('Extracting noise profile from file %s...\n', params.noise);
    end
    if ~exist(params.noise, 'file')
        fprintf('ERROR! params.noise==input is specified, but params.noise_profile doesn''t exist\n');
    end
    noise_profile = wavread(params.noise);
end
if params.show_log
    fprintf('Computing noise atoms...\n');
end
Ynoise = abs(m_STFT(noise_profile, sr, params));
if params.m >= size(Ynoise,2)
    fprintf('WARNING! Noise profile is too small\n');
    params.m = size(Ynoise,2)-1;
end
Dnoise = simplestNMF(Ynoise, params);
if strcmp(params.type, 'linNMF') || strcmp(params.type, 'denseNMF')
    if params.show_log
        fprintf('Construct harmonic/noise atoms...\n');
    end
    [atomsH, DH] = harmonic_atoms(params.f0range, sr, params.window, ...
        params.num_atoms_per_f0, params.max_harms, params.fstep);
    [atomsN, DN] = noise_atoms_from_dict(Dnoise, params.num_noise_atoms);
    if params.show_log
        fprintf('Creating signal mask using NMF with linear constraints...\n');
    end
    W = noise_reduction_NMF_linear(Yaudio, atomsH, atomsN, params);
else
    if params.show_log
        fprintf('Creating signal mask using NMF...\n');
    end
    W = noise_reduction_NMF(Yaudio, Dnoise);
end
if params.show_log
    fprintf('Doing denoising transformations...\n');
end
output_audio = m_InverseSTFT(W.*Saudio, params);

if isempty(params.output)
    [pathstr, name] = fileparts(input_file);
    output_file = sprintf('%s/%s_DENOISED.wav', pathstr, name);
else
    output_file = params.output;
end
if params.show_log
    fprintf('Saving output file to %s...\n', output_file);
end
wavwrite(output_audio, sr, output_file);

end

function pout = initParams(p, sr)
pout = p;
pout.sr = sr;
if ~isfield(pout, 'nwin')           pout.nwin = 1024; end
if ~isfield(pout, 'window')         pout.window = hann(pout.nwin); end
if ~isfield(pout, 'type')           pout.type = 'denseNMF'; end
if ~isfield(pout, 'noise')          pout.noise = ''; end
if ~isfield(pout, 'beta')           pout.beta = 1; end
if ~isfield(pout, 'm')              pout.m = 32; end
if ~isfield(pout, 'f0range')        pout.f0range = [80 400]; end
if ~isfield(pout, 'fstep')          pout.fstep = 10; end
if ~isfield(pout, 'num_atoms_per_f0') pout.num_atoms_per_f0 = 4; end
if ~isfield(pout, 'num_noise_atoms') pout.num_noise_atoms = 16; end
if ~isfield(pout, 'max_harms')      pout.max_harms = 30; end
if ~isfield(pout, 'output')         pout.output = ''; end
if ~isfield(pout, 'show_log')       pout.show_log = true; end
if ~isfield(pout, 'energyThr')       pout.energyThr = -40; end
if ~isfield(pout, 'alpha')          pout.alpha = .2; end
if ~isfield(pout, 'speech_sparsity')          pout.speech_sparsity = .2; end
end