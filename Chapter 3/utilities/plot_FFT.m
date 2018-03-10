% PLot FFT of signal
% Inputs: 
%   - signal: signal
%   - dt: sampling time
%   - cc: color (optional)
function plot_FFT(signal,dt,cc)

if size(signal,2) > size(signal,1)
    signal  = signal';
end
signal      = signal - mean(signal);
fs          = 1/dt;           
N           = length(signal);
N2          = 2^nextpow2(N);
f           = fs*(0:N2-1)/N2;
s           = [(signal.*hamming(N))',zeros(1,N2-N)];
Y           = fft(s, N2);
specdB      = 20*log10(abs(Y));
specdBn     = specdB;    
if nargin > 2
    plot(f(1:N2/2+1),specdBn(1:N2/2+1),'Color',cc,'Linewidth',2);   
else
    plot(f(1:N2/2+1),specdBn(1:N2/2+1),'Linewidth',2);   
end
axis tight
