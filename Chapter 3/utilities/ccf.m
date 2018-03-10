% Compute the correlation between the signals s1 and s2 
%   - Nlags: number of lags
%   - lw: LineWidth in the figure, by default = 2
function [c1,lags1] = ccf(s1,s2,Nlags,lw)

if nargin < 4
    lw = 2;
end

N = length(s1);
if N ~= length(s2)
    disp('Signals with different length')
end
s1m = s1 - mean(s1);
s2m = s2 - mean(s2);
[c,lags] = xcorr(s1m, s2m, Nlags, 'coeff');
stem(lags(Nlags+1:end),c(Nlags+1:end),'markerfacecolor',[0 0 1],'marker','none','linewidth',lw)
line([-2 Nlags+2],[1.96/sqrt(N) 1.96/sqrt(N)],'color',[1 0 0],'linewidth',lw)
line([-2 Nlags+2],[-1.96/sqrt(N) -1.96/sqrt(N)],'color',[1 0 0],'linewidth',lw)
% axis tight
xlim([-2 Nlags+2])
xlabel('lags')
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',16)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
set(gcf,'color','w');

if nargout > 0 
    c1 = c;
    lags1 = lags;
end