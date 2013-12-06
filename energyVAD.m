function [ vad ] = energyVAD( audio, params )
%M_ENERGYVAD calculates voice activity from input audio using simple energy
%threshold-based algorithm
%INPUT:
%   audio:  input audio samples [1xN]
%   vad:    (0,1)-label array [1xN] where 1 indicates voice activity per samples

if nargin<2
    params = struct();
end
params = initParams(params);
frames = buffer(audio, params.nwin);
rms = std(frames);
vad = zeros(size(rms));
vad(20*log10(rms)>params.threshold)=1;
vad = kron(vad, ones(1, params.nwin));
vad = vad(1:numel(audio));
end

function pout = initParams(p)

pout = p;
if ~isfield(pout, 'nwin')       pout.nwin = 64;     end
if ~isfield(pout, 'threshold')  pout.threshold = -40; end
end