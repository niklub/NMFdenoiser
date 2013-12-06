function S = m_STFT(audio, sr, params)
if nargin<3 
    params = struct();
end
params = initParams(params, sr);
frames = buffer(audio, params.nwin, params.nwin-params.nhop, 'nodelay');
frames = frames .* repmat(params.window, 1, size(frames,2));
nfft = 2^nextpow2(params.nwin) * params.npadtimes;
S = fft(frames, nfft, 1);
S = S(1:floor(nfft/2)+1,:);
if params.winscale
    scale = 2/sum(params.window);
    S = S * scale;
    S(1, :) = S(1, :)/2;
    S(end,:) = S(end, :)/2;
end
end


function pout = initParams(p, sr)
pout = p;
pout.sr = sr;
if ~isfield(pout, 'nwin')       pout.nwin = 512;                    end
if ~isfield(pout, 'window')     pout.window = hann(pout.nwin);      end
if ~isfield(pout, 'nhop')       pout.nhop = round(pout.nwin/4);     end
if ~isfield(pout, 'npadtimes')  pout.npadtimes = 1;                 end
if ~isfield(pout, 'winscale')   pout.winscale = false;              end

pout.window = pout.window(:);
end