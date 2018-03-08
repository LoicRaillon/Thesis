% Integrated autocorrelation time (IACT, section 4.2.2.4
% Inputs: 
%   - pS: Markov Chains
%   - burnin: burnin
%   - Nlags: Number of lags
%   - col: Number of columns in the figure
% Output: 
%   - iact: IACT value for each Markov chains

function iact   = IACT_plot(pS,burnin,Nlags,col)

Nukwn           = size(pS,1);
N               = size(pS,2);
row             = round(Nukwn/col) + mod(Nukwn,col);
bis             = zeros(Nukwn,1);
figure

for i = 1:Nukwn
    subplot(row,col,i)
    chain       = pS(i,burnin:end);
    [acf,~]     = ccf(chain,chain,Nlags);
    acf2        = acf(Nlags+2:end);
    idx         = find(acf2 < 2/sqrt(N-burnin));
    if ~isempty(idx)
        IACT    = 1 + 2*sum(acf2(1:idx(1)));
        bis(i)  = IACT;
        title(['IACT = ', num2str(IACT)]);
    end
end

if any(nargout)
    iact        = bis;
end




