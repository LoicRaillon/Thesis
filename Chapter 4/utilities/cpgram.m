% Plot cumulated periodogram of the residuals (res)
%   - lw: LineWidth
%   - col: color
function cpgram(res,lw,col)

if nargin < 2
    lw   = 2;
end

if nargin < 3
    col  = [0 0 1];
end

nobs     = length(res);
rfft     = fft(res);   
n        = length(rfft);
n2       = n/2+1;
power    = (abs(rfft(1:n2)).^2)/n; % calculates intensities
nyquist  = 1/2;
pfreq    = (0:1:n/2);
largofre = length(pfreq)-1;
freq     = pfreq/(largofre)*nyquist;
linea    = freq*2;  
cuper    = cumsum(power); 
cupera   = cuper./max(cuper); 
q        = (nobs/2)+1;
crit     = 1.36/(sqrt(q));

plot(freq,cupera,'Color',col,'linewidth',lw)
axis([0 0.5 0 1])
xlabel('Nyquist Frequency','Fontsize',14)
hold on
plot(freq,linea+crit,'r','linewidth',1) 
plot(freq,linea-crit,'r','linewidth',1)
% axis tight
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',14)
set(gca,'FontName','Times New Roman')
set(gca,'Linewidth',2)
set(gcf,'color','w');

