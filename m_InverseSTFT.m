function audio = m_InverseSTFT(S, params)
if nargin < 2
    params = struct();
end
[nfft, nframes] = size(S);
nfft = 2*(nfft-1);
params = initParams(params, nfft);
if params.winscale
    scale = 2/sum(params.window);
    S = S / scale;
    S(1, :) = S(1,:)*2;
    S(end,:) = S(end,:)*2;
end
S = [S; conj(S(floor((nfft+1)/2):-1:2,:))];
frames = real(ifft(S,[],1));
frames = frames(1:params.nwin, :);
frames = frames .* repmat(params.window, 1, nframes);
audio = 2/3*unbuffer(frames,params.nwin,params.nwin-params.nhop);
end
function pout = initParams(p, nfft)
pout = p;
if ~isfield(pout, 'nwin')       pout.nwin = nfft;                 end
if ~isfield(pout, 'window')     pout.window = hann(pout.nwin);      end
if ~isfield(pout, 'nhop')       pout.nhop = round(pout.nwin/4);     end
if ~isfield(pout, 'winscale')   pout.winscale = false;              end

pout.window = pout.window(:);
end

function y = unbuffer(x,w,o)
y    = [];
skip = w - o;
N    = ceil(w/skip);
L    = (size(x,2) - 1) * skip + size(x,1);
if size(x,1)<skip*N
    x(skip*N,end) = 0; 
end
for i = 1:N
    t = reshape(x(:,i:N:end),1,[]);
    l = length(t);
    y(i,l+(i-1)*skip) = 0;
    y(i,[1:l]+(i-1)*skip) = t;
end;
y = sum(y,1);
y = y(1:L);
end