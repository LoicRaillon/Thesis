% Compute and plot all the results of the second chapter
%
% files requested:
%   - var_sens.m : compute sensitivity indices with Quasi Random sampling with uniform distribution [0,1]
%                  to improve the convergence of the estimator
%   - fct.m : List of all equations for each transformations
%   - LPTAU51.m : Quasi Random sampling with uniform distribution [0,1]
%   generator (see inside the function for references)
%
% Functions used: 
%   - subplot_pos: Custom subplot, Pierre Martineau, http://p-martineau.com/perfect-subplot-in-matlab/
%   - export_fig: export figure as it appears on screen, https://github.com/altmany/export_fig
%
% Np = Number of samples for 'propag'
% Ns = Number of iterations for 'sensi'
% mode = propag = error propagation
% mode = sensi = variance based sensitivity analysis
% tfm = choice of the transformation
% tfm = 'rcss' : RC --> SSd
% tfm = 'sstf' : SSd --> TF
% tfm = 'ssrc' : SSd --> RC
% tfm = 'tfss' : TF --> SSd
% tfm = 'rctf' : RC --> TF
% tfm = 'tfrc' : TF --> RC
% ts = sampling time [s]

clear all;
close all;
clc;

Np                              = 1e6;                                      % number of simulation for 'propag'
Ns                              = 2^14;                                     % number of samples for 'sensi'
mode                            = 'sensi';                                  % factor prioritization: 'sensi', error propagation: 'propag'
format shorte                                   
tfm                             = 'tfrc';                                   % choice of transformation 
ts                              = 60;                                       % sampling time
delta                           = 0.1; % 10%                                % delta equation (2.1), scaling of standard deviation        
delta2                          = 0.01; % 1%                                % lower scaling from TF to RC
% Parameter names
pSS                             = {'a11 ','a12 ','a21 ','a22 ','b11 ','b13 ','b21 ','b22 '};
pRC                             = {'Rco ','Rci ','Rv ','Rw ','Cw ','Ca '};
pTF                             = {'n11 ','n12 ','n22 ','n31 ','n32 ','d1 ','d2 '};
% nominal values from R,C
Rco                             = 1/250;
Rci                             = 1/125;
Rv                              = 1/38.3;
Rw                              = 1/1.45;
Cw                              = 4e6;
Ca                              = 82e3;
% SS continuous
a11c                            = -(2*Rci+2*Rv+Rw)/(Ca*Rv*(2*Rci+Rw));
a12c                            = 2/(Ca*(2*Rci+Rw));
a21c                            = 2/(Cw*(2*Rci+Rw));
a22c                            = -(4*(Rci+Rco+Rw))/(Cw*(2*Rco+Rw)*(2*Rci+Rw));
b11c                            = 1/(Rv*Ca);
b13c                            = 1/Ca;
b21c                            = 2/(Cw*(2*Rco+Rw));
b22c                            = (2*Rco)/(Cw*(Rw+2*Rco));
% First order Euler discretization
a11d                            = 1+a11c*ts;
a12d                            = a12c*ts;
a21d                            = a21c*ts;
a22d                            = 1+a22c*ts;
b11d                            = b11c*ts;
b13d                            = b13c*ts;
b21d                            = b21c*ts;
b22d                            = b22c*ts;
Ad                              = [a11d a12d ; a21d a22d];
Bd                              = [b11d 0 b13d ; b21d b22d 0];
C                               = [1 0];
% similarity transformation
T                               = [1 0 ; -a11d/a12d 1/a12d];
Acd                             = T\Ad*T;
Beta                            = T\Bd;
d1                              = -Acd(2,2);
d2                              = -Acd(2,1);
Bcd                             = [1 0 ; d1 1]*Beta;
n11                             = Bcd(1,1);
n12                             = Bcd(2,1);
n22                             = Bcd(2,2);
n31                             = Bcd(1,3);
n32                             = Bcd(2,3);
% All the transformations are coded similarly, so only the first one is
% commented ! 
switch tfm
    case 'rcss'
        Pstr                    = {'a11','a12','a21','a22','b11',...        % Name of output parameters
                                        'b13','b21','b22'};                 % ...
        Pout                    = [a11c a12c a21c a22c b11c b13c b21c b22c];% True values of output parameters
        Pin                     = [Rco Rci Rv Rw Cw Ca];                    % True values of input parameters
        lPin                    = length(Pin);                              % get size
        lPout                   = length(Pout);                             % ....
        num                     = [4 3 3 4 2 1 3 3];                        % Number of input parameters by output parameter, for instance, a11 depends on 4 inputs parameters
        slct                    = [2 3 4 6 ...                              % a11 depends on the parameters n°2,3,4 and 6
                                   2 4 6 ...                                % It could be simplified as shown for the other transformations
                                   2 4 5 ...                                % ...
                                   1 2 4 5 ...                              % ...
                                   3 6 ...                                  % ...
                                   6 ...                                    % ...
                                   1 4 5 ...                                % ...
                                   1 4 5];                                  % ...
        switch mode
            case 'sensi'                                                    % Factor prioritization
                Sf              = zeros(lPout,lPin);                        % Allocation first order sensitivity indices
                STf             = zeros(lPout,lPin);                        % Allocation total sensitivity indices
                k1              = 1;                                        % indexing: start
                for p = 1:lPout     
                    k2          = sum(num(1:p));                            % indexing: stop
                    [S,ST]      = var_sens(tfm,Pstr{p},num(p),ts,Ns,[]);    % Factor prioritization function 
                    idx         = slct(k1:k2);                              % indexing
                    Sf(p,idx)   = S;                                        % save
                    STf(p,idx)  = ST;                                       % ...
                    k1          = k2 + 1;                                   % update indexing
                end
                Tab             = array2table([Sf sum(Sf,2)],'RowNames',... % Table for first order sensitivity indices
                                        Pstr,'VariableNames',[pRC 'sum']);  % ...
                Tabf            = array2table([STf sum(STf,2)],'RowNames'...% Table for total sensitivity indices
                    ,Pstr,'VariableNames',[pRC 'sum']);                     % ...
                disp(['Transformation: ' tfm ' | number of samples: '...    % display
                    num2str(Ns) ' | sampling time: ' num2str(ts) ' s']);    % ...
                fprintf(1,'\n')                                             % ...
                disp('First order sensitivity indices')                     % ...
                fprintf(1,'\n')                                             % ...
                disp(Tab)                                                   % ...
                disp('Total order sensitivity indices')                     % ...
                fprintf(1,'\n')                                             % ...
                disp(Tabf)                                                  % ...
            case 'propag'
                X               = bsxfun(@plus,Pin',bsxfun(@times,...       % Generate random sample (2.1)
                                            delta.*Pin',randn(lPin,Np)))';  % ...
                tab             = zeros(lPout,4);                           % Allocation Table
                y               = zeros(Np,lPout);                          % Allocation model output
                k1              = 1;                                        % count
                for p = 1:lPout
                    k2          = sum(num(1:p));                            % indexing
                    idx         = slct(k1:k2);                              % ...
                    y(:,p)      = fct(tfm,Pstr{p},ts,X(:,idx));             % model evaluation
                    k1          = k2 + 1;                                   % update indexing
                    pd          = fitdist(y(:,p),'normal');                 % fit to Gaussian distribution
                    [~,~,d,~]   = kstest(y(:,p),pd,0.05);                   % Kolmogorov-Smirnov test Figure 2.4
                    [N,bh]      = histcounts(y(:,p));                       % Find expected values
                    [~,idx]     = max(N);                                   % ...
                    E           = (bh(idx) + bh(idx+1))/2;                  % ...
                    MAE         = mean(abs(y(:,p)-Pout(p)));                % Mean absolute error (2.35)
                    tab(p,:)    = [Pout(p) E MAE d];                        % Store in the table
                end
                Tabp            = array2table(tab,'RowNames',Pstr,...       % display table
                                    'VariableNames',{'pout' 'E' 'MAE' 'd'});% ...
                disp(Tabp)                                                  % ...
        end
    case 'sstf'
        Pstr                    = {'n11','n12','n22','n31','n32','d1','d2'};
        Pout                    = [n11 n12 n22 n31 n32 d1 d2];
        Pin                     = [a11c a12c a21c a22c b11c b13c b21c b22c];
        lPin                    = length(Pin);
        lPout                   = length(Pout);
        num                     = 8;
        switch mode
            case 'sensi'
                Sf              = zeros(lPout,lPin);
                STf             = zeros(lPout,lPin);
                for p = 1:lPout
                    [S,ST]      = var_sens(tfm,Pstr{p},num,ts,Ns,[]);
                    Sf(p,:)     = S;
                    STf(p,:)    = ST;
                end
                Tab             = array2table([Sf sum(Sf,2)],'RowNames',Pstr,'VariableNames',[pSS 'sum']);
                Tabf            = array2table([STf sum(STf,2)],'RowNames',Pstr,'VariableNames',[pSS 'sum']);
                disp(['Transformation: ' tfm ' | number of samples: ' num2str(Ns) ' | sampling time: ' num2str(ts) ' s']);
                fprintf(1,'\n')
                disp('First order sensitivity indices')
                fprintf(1,'\n')
                disp(Tab)
                disp('Total order sensitivity indices')
                fprintf(1,'\n')
                disp(Tabf)
            case 'propag'
                X               = bsxfun(@plus,Pin',bsxfun(@times,delta.*Pin',randn(lPin,Np)))';
                tab             = zeros(lPout,4);
                y               = zeros(Np,lPout);
                for p = 1:lPout
                    y(:,p)      = fct(tfm,Pstr{p},ts,X);
                    pd          = fitdist(y(:,p),'normal');
                    [~,~,d,~]   = kstest(y(:,p),pd,0.05);
                    [N,bh]      = histcounts(y(:,p));
                    [~,idx]     = max(N);
                    E           = (bh(idx) + bh(idx+1))/2;
                    MAE         = mean(abs(y(:,p)-Pout(p)));
                    tab(p,:)    = [Pout(p) E MAE d];
                end
                Tabp            = array2table(tab,'RowNames',Pstr,'VariableNames',{'pout' 'E' 'MAE' 'd'});
                disp(Tabp)
        end
    case 'rctf'
        Pstr                    = {'n11','n12','n22','n31','n32','d1','d2'};
        Pout                    = [n11 n12 n22 n31 n32 d1 d2];
        Pin                     = [Rco Rci Rv Rw Cw Ca];
        lPin                    = length(Pin);
        lPout                   = length(Pout);
        num                     = 6;
        switch mode
            case 'sensi'
                Sf              = zeros(lPout,lPin);
                STf             = zeros(lPout,lPin);
                for p = 1:lPout
                    [S,ST]      = var_sens(tfm,Pstr{p},num,ts,Ns,[]);
                    Sf(p,:)     = S;
                    STf(p,:)    = ST;
                end
                Tab             = array2table([Sf sum(Sf,2)],'RowNames',Pstr,'VariableNames',[pRC 'sum']);
                Tabf            = array2table([STf sum(STf,2)],'RowNames',Pstr,'VariableNames',[pRC 'sum']);
                disp(['Transformation: ' tfm ' | number of samples: ' num2str(Ns) ' | sampling time: ' num2str(ts) ' s']);
                fprintf(1,'\n')
                disp('First order sensitivity indices')
                fprintf(1,'\n')
                disp(Tab)
                disp('Total order sensitivity indices')
                fprintf(1,'\n')
                disp(Tabf)
            case 'propag'
                X               = bsxfun(@plus,Pin',bsxfun(@times,delta.*Pin',randn(lPin,Np)))';
                tab             = zeros(lPout,4);
                y               = zeros(Np,lPout);
                for p = 1:lPout
                    y(:,p)      = fct(tfm,Pstr{p},ts,X);
                    pd          = fitdist(y(:,p),'normal');
                    [~,~,d,~]   = kstest(y(:,p),pd,0.05);
                    [N,bh]      = histcounts(y(:,p));
                    [~,idx]     = max(N);
                    E           = (bh(idx) + bh(idx+1))/2;
                    MAE         = mean(abs(y(:,p)-Pout(p)));
                    tab(p,:)    = [Pout(p) E MAE d];
                end
                Tabp            = array2table(tab,'RowNames',Pstr,'VariableNames',{'pout' 'E' 'MAE' 'd'});
                disp(Tabp)
        end
    case 'tfss'
        Pstr                    = {'a21','a22','b11','b13','b21','b22'};
        Pout                    = [a21c a22c b11c b13c b21c b22c];
        Pin                     = [n11 n12 n22 n31 n32 d1 d2];
        lPin                    = length(Pin);
        lPout                   = length(Pout);
        num                     = 7;
        switch mode
            case 'sensi'
                Sf              = zeros(lPout,lPin);
                STf             = zeros(lPout,lPin);
                for p = 1:lPout
                    [S,ST]      = var_sens(tfm,Pstr{p},num,ts,Ns,T);
                    Sf(p,:)     = S;
                    STf(p,:)    = ST;
                end
                Tab             = array2table([Sf sum(Sf,2)],'RowNames',Pstr,'VariableNames',[pTF 'sum']);
                Tabf            = array2table([STf sum(STf,2)],'RowNames',Pstr,'VariableNames',[pTF 'sum']);
                disp(['Transformation: ' tfm ' | number of samples: ' num2str(Ns) ' | sampling time: ' num2str(ts) ' s']);
                fprintf(1,'\n')
                disp('First order sensitivity indices')
                fprintf(1,'\n')
                disp(Tab)
                disp('Total order sensitivity indices')
                fprintf(1,'\n')
                disp(Tabf)
            case 'propag'
                X               = bsxfun(@plus,Pin',bsxfun(@times,abs(delta2.*Pin'),randn(lPin,Np)))';
                tab             = zeros(lPout,4);
                y               = zeros(Np,lPout);
                for p = 1:lPout
                    y(:,p)      = fct(tfm,Pstr{p},ts,X,T);
                    pd          = fitdist(y(:,p),'normal');
                    [~,~,d,~]   = kstest(y(:,p),pd,0.05);
                    [N,bh]      = histcounts(y(:,p));
                    [~,idx]     = max(N);
                    E           = (bh(idx) + bh(idx+1))/2;
                    MAE         = mean(abs(y(:,p)-Pout(p)));
                    tab(p,:)    = [Pout(p) E MAE d];
                end
                Tabp            = array2table(tab,'RowNames',Pstr,'VariableNames',{'pout' 'E' 'MAE' 'd'});
                disp(Tabp)
        end
    case 'ssrc'
        Pstr                    = {'Rco','Rci','Rv','Rw','Cw','Ca'};
        Pout                    = [Rco Rci Rv Rw Cw Ca];
        Pin                     = [a11c a12c a21c a22c b11c b13c b21c b22c];
        lPin                    = length(Pin);
        lPout                   = length(Pout);
        num                     = 8;
        switch mode
            case 'sensi'
                Sf              = zeros(lPout,lPin);
                STf             = zeros(lPout,lPin);
                for p = 1:lPout
                    [S,ST]      = var_sens(tfm,Pstr{p},num,ts,Ns,[]);
                    Sf(p,:)     = S;
                    STf(p,:)    = ST;
                end
                Tab             = array2table([Sf sum(Sf,2)],'RowNames',Pstr,'VariableNames',[pSS 'sum']);
                Tabf            = array2table([STf sum(STf,2)],'RowNames',Pstr,'VariableNames',[pSS 'sum']);
                disp(['Transformation: ' tfm ' | number of samples: ' num2str(Ns) ' | sampling time: ' num2str(ts) ' s']);
                fprintf(1,'\n')
                disp('First order sensitivity indices')
                fprintf(1,'\n')
                disp(Tab)
                disp('Total order sensitivity indices')
                fprintf(1,'\n')
                disp(Tabf)
            case 'propag'
                X               = bsxfun(@plus,Pin',bsxfun(@times,abs(delta2.*Pin'),randn(lPin,Np)))';
                tab             = zeros(lPout,4);
                y               = zeros(Np,lPout);
                for p = 1:lPout
                    y(:,p)      = fct(tfm,Pstr{p},ts,X,[]);
                    pd          = fitdist(y(:,p),'normal');
                    [~,~,d,~]   = kstest(y(:,p),pd,0.05);
                    [N,bh]      = histcounts(y(:,p));
                    [~,idx]     = max(N);
                    E           = (bh(idx) + bh(idx+1))/2;
                    MAE         = mean(abs(y(:,p)-Pout(p)));
                    tab(p,:)    = [Pout(p) E MAE d];
                end
                Tabp            = array2table(tab,'RowNames',Pstr,'VariableNames',{'pout' 'E' 'MAE' 'd'});
                disp(Tabp)
        end
    case 'tfrc'
        Pstr                    = {'Rco','Rci','Rv','Rw','Cw','Ca'};
        Pout                    = [Rco Rci Rv Rw Cw Ca];
        Pin                     = [n11 n12 n22 n31 n32 d1 d2];
        lPin                    = length(Pin);
        lPout                   = length(Pout);
        num                     = 7;
        switch mode
            case 'sensi'
                Sf              = zeros(lPout,lPin);
                STf             = zeros(lPout,lPin);
                for p = 1:lPout
                    [S,ST]      = var_sens(tfm,Pstr{p},num,ts,Ns,T);
                    Sf(p,:)     = S;
                    STf(p,:)    = ST;
                end
                Tab             = array2table([Sf sum(Sf,2)],'RowNames',Pstr,'VariableNames',[pTF 'sum']);
                Tabf            = array2table([STf sum(STf,2)],'RowNames',Pstr,'VariableNames',[pTF 'sum']);
                disp(['Transformation: ' tfm ' | number of samples: ' num2str(Ns) ' | sampling time: ' num2str(ts) ' s']);
                fprintf(1,'\n')
                disp('First order sensitivity indices')
                fprintf(1,'\n')
                disp(Tab)
                disp('Total order sensitivity indices')
                fprintf(1,'\n')
                disp(Tabf)
            case 'propag'
                X               = bsxfun(@plus,Pin',bsxfun(@times,abs(delta2.*Pin'),randn(lPin,Np)))';
                tab             = zeros(lPout,4);
                y               = zeros(Np,lPout);
                for p = 1:lPout
                    y(:,p)      = fct(tfm,Pstr{p},ts,X,T);
                    pd          = fitdist(y(:,p),'normal');
                    [~,~,d,~]   = kstest(y(:,p),pd,0.05);
                    [N,bh]      = histcounts(y(:,p));
                    [~,idx]     = max(N);
                    E           = (bh(idx) + bh(idx+1))/2;
                    MAE         = mean(abs(y(:,p)-Pout(p)));
                    tab(p,:)    = [Pout(p) E MAE d];
                end
                Tabp            = array2table(tab,'RowNames',Pstr,'VariableNames',{'pout' 'E' 'MAE' 'd'});
                disp(Tabp)
        end
end

%% ------------------------------------------------------------------------
% Plot, mode = 'propag'
% -------------------------------------------------------------------------
if strcmp(mode,'propag')
    switch tfm
        %% ------------------------------------------------------------------------
        % Thermal network to state-space
        % -------------------------------------------------------------------------
        case 'rcss'
            Pout                = [a11c a12c a21c a22c b11c b13c b21c b22c];
            plotheight          = 15;
            plotwidth           = 15;
            subplotsx           = 2;
            subplotsy           = 4;
            leftedge            = 1.1;
            rightedge           = 1;
            topedge             = 1.1;
            bottomedge          = 1;
            spacex              = 1.1;
            spacey              = 1.1;
            fontsize            = 16;
            LineWidth           = 2;
            alpha               = 0.3;
            sub_pos             = subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);
            f                   = figure('visible','on','units','normalized','outerposition',[0 0 0.6 1]);
            clf(f);
            set(gcf, 'PaperUnits', 'centimeters');
            set(gcf, 'PaperSize', [plotwidth plotheight]);
            set(gcf, 'PaperPositionMode', 'manual');
            set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
            set(gcf,'color','w');
            % -------------------------- a11 ------------------------------------------
            p                   = 1;
            axes('position',sub_pos{1,4});
            h                   = histogram(y(:,p));
            axis tight
            h.FaceColor         = 'b';
            h.EdgeColor         = 'b';
            h.EdgeAlpha         = 0;
            h.FaceAlpha         = alpha;
            yy                  = max(h.Values);
            hold on
            xh                  = [Pout(p) Pout(p)];
            yh                  = [0 yy];
            plot(xh,yh,'r','LineWidth',LineWidth);
            axis('tight')
            set(gca,'Box','off')
            set(gca,'FontWeight','demi')
            set(gca,'FontSize',fontsize)
            set(gca,'FontName','Cambria')
            set(gca,'Linewidth',LineWidth)
            x1                  = get(gca,'Xtick');
            y1                  = get(gca,'Ytick');
            set(gca,'XTickLabel',x1.*1e4)
            set(gca,'YTickLabel',y1./1e3)
            text(1.01,0.08,'10^{-4}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.95,0.9,'a_1_1','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
            % -------------------------- a12 ------------------------------------------
            p                   = 2;
            axes('position',sub_pos{2,4});
            h                   = histogram(y(:,p));
            axis tight
            h.FaceColor         = 'b';
            h.EdgeColor         = 'b';
            h.EdgeAlpha         = 0;
            h.FaceAlpha         = alpha;
            yy                  = max(h.Values);
            hold on
            xh                  = [Pout(p) Pout(p)];
            yh                  = [0 yy];
            plot(xh,yh,'r','LineWidth',LineWidth);
            axis('tight')
            set(gca,'Box','off')
            set(gca,'FontWeight','demi')
            set(gca,'FontSize',fontsize)
            set(gca,'FontName','Cambria')
            set(gca,'Linewidth',LineWidth)
            x1                  = get(gca,'Xtick');
            y1                  = get(gca,'Ytick');
            set(gca,'XTickLabel',x1.*1e5)
            set(gca,'YTickLabel',y1./1e3)
            text(1.01,0.08,'10^{-5}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.95,0.9,'a_1_2','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
            % -------------------------- a21 ------------------------------------------
            p                   = 3;
            axes('position',sub_pos{1,3});
            h                   = histogram(y(:,p));
            axis tight
            h.FaceColor         = 'b';
            h.EdgeColor         = 'b';
            h.EdgeAlpha         = 0;
            h.FaceAlpha         = alpha;
            yy                  = max(h.Values);
            hold on
            xh                  = [Pout(p) Pout(p)];
            yh                  = [0 yy];
            plot(xh,yh,'r','LineWidth',LineWidth);
            axis('tight')
            set(gca,'Box','off')
            set(gca,'FontWeight','demi')
            set(gca,'FontSize',fontsize)
            set(gca,'FontName','Cambria')
            set(gca,'Linewidth',LineWidth)
            x1                  = get(gca,'Xtick');
            y1                  = get(gca,'Ytick');
            set(gca,'XTickLabel',x1.*1e7)
            set(gca,'YTickLabel',y1./1e3)
            text(1.01,0.08,'10^{-7}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.95,0.9,'a_2_1','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
            % -------------------------- a22 ------------------------------------------
            p                   = 4;
            axes('position',sub_pos{2,3});
            h                   = histogram(y(:,p));
            axis tight
            h.FaceColor         = 'b';
            h.EdgeColor         = 'b';
            h.EdgeAlpha         = 0;
            h.FaceAlpha         = alpha;
            yy                  = max(h.Values);
            hold on
            xh                  = [Pout(p) Pout(p)];
            yh                  = [0 yy];
            plot(xh,yh,'r','LineWidth',LineWidth);
            axis('tight')
            set(gca,'Box','off')
            set(gca,'FontWeight','demi')
            set(gca,'FontSize',fontsize)
            set(gca,'FontName','Cambria')
            set(gca,'Linewidth',LineWidth)
            x1                  = get(gca,'Xtick');
            y1                  = get(gca,'Ytick');
            set(gca,'XTickLabel',x1.*1e6)
            set(gca,'YTickLabel',y1./1e3)
            text(1.01,0.08,'10^{-6}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.95,0.9,'a_2_2','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
            % -------------------------- b11 ------------------------------------------
            p                   = 5;
            axes('position',sub_pos{1,2});
            h                   = histogram(y(:,p));
            axis tight
            h.FaceColor         = 'b';
            h.EdgeColor         = 'b';
            h.EdgeAlpha         = 0;
            h.FaceAlpha         = alpha;
            yy                  = max(h.Values);
            hold on
            xh                  = [Pout(p) Pout(p)];
            yh                  = [0 yy];
            plot(xh,yh,'r','LineWidth',LineWidth);
            axis('tight')
            set(gca,'Box','off')
            set(gca,'FontWeight','demi')
            set(gca,'FontSize',fontsize)
            set(gca,'FontName','Cambria')
            set(gca,'Linewidth',LineWidth)
            x1                  = get(gca,'Xtick');
            y1                  = get(gca,'Ytick');
            set(gca,'XTickLabel',x1.*1e4)
            set(gca,'YTickLabel',y1./1e3)
            text(1.01,0.08,'10^{-4}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.95,0.9,'b_1_1','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
            % -------------------------- b13 ------------------------------------------
            p                   = 6;
            axes('position',sub_pos{2,2});
            h                   = histogram(y(:,p));
            axis tight
            h.FaceColor         = 'b';
            h.EdgeColor         = 'b';
            h.EdgeAlpha         = 0;
            h.FaceAlpha         = alpha;
            yy                  = max(h.Values);
            hold on
            xh                  = [Pout(p) Pout(p)];
            yh                  = [0 yy];
            plot(xh,yh,'r','LineWidth',LineWidth);
            axis('tight')
            set(gca,'Box','off')
            set(gca,'FontWeight','demi')
            set(gca,'FontSize',fontsize)
            set(gca,'FontName','Cambria')
            set(gca,'Linewidth',LineWidth)
            x1                  = get(gca,'Xtick');
            y1                  = get(gca,'Ytick');
            set(gca,'XTickLabel',x1.*1e5)
            set(gca,'YTickLabel',y1./1e3)
            text(1.01,0.08,'10^{-5}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.95,0.9,'b_1_3','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
            % -------------------------- b21 ------------------------------------------
            p                   = 7;
            axes('position',sub_pos{1,1});
            h                   = histogram(y(:,p));
            axis tight
            h.FaceColor         = 'b';
            h.EdgeColor         = 'b';
            h.EdgeAlpha         = 0;
            h.FaceAlpha         = alpha;
            yy                  = max(h.Values);
            hold on
            xh                  = [Pout(p) Pout(p)];
            yh                  = [0 yy];
            plot(xh,yh,'r','LineWidth',LineWidth);
            axis('tight')
            set(gca,'Box','off')
            set(gca,'FontWeight','demi')
            set(gca,'FontSize',fontsize)
            set(gca,'FontName','Cambria')
            set(gca,'Linewidth',LineWidth)
            x1                  = get(gca,'Xtick');
            y1                  = get(gca,'Ytick');
            set(gca,'XTickLabel',x1.*1e7)
            set(gca,'YTickLabel',y1./1e3)
            text(1.01,0.08,'10^{-7}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.95,0.9,'b_2_1','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
            % -------------------------- b22 ------------------------------------------
            p                   = 8;
            axes('position',sub_pos{2,1});
            h                   = histogram(y(:,p));
            axis tight
            h.FaceColor         = 'b';
            h.EdgeColor         = 'b';
            h.EdgeAlpha         = 0;
            h.FaceAlpha         = alpha;
            yy                  = max(h.Values);
            hold on
            xh                  = [Pout(p) Pout(p)];
            yh                  = [0 yy];
            plot(xh,yh,'r','LineWidth',LineWidth);
            axis('tight')
            set(gca,'Box','off')
            set(gca,'FontWeight','demi')
            set(gca,'FontSize',fontsize)
            set(gca,'FontName','Cambria')
            set(gca,'Linewidth',LineWidth)
            x1                  = get(gca,'Xtick');
            y1                  = get(gca,'Ytick');
            set(gca,'XTickLabel',x1.*1e9)
            set(gca,'YTickLabel',y1./1e3)
            text(1.01,0.08,'10^{-9}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.95,0.9,'b_2_2','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
            %% ------------------------------------------------------------------------
            % state-space to transfer function
            % -------------------------------------------------------------------------
        case 'sstf'
            plotheight          = 15;
            plotwidth           = 15;
            subplotsx           = 2;
            subplotsy           = 4;
            leftedge            = 1.1;
            rightedge           = 1;
            topedge             = 1.1;
            bottomedge          = 1;
            spacex              = 1.1;
            spacey              = 1.1;
            fontsize            = 16;
            LineWidth           = 2;
            alpha               = 0.3;
            sub_pos             = subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);
            f                   = figure('visible','on','units','normalized','outerposition',[0 0 0.6 1]);
            clf(f);
            set(gcf, 'PaperUnits', 'centimeters');
            set(gcf, 'PaperSize', [plotwidth plotheight]);
            set(gcf, 'PaperPositionMode', 'manual');
            set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
            set(gcf,'color','w');
            % -------------------------- n11 ------------------------------------------
            p                   = 1;
            axes('position',sub_pos{1,4});
            h                   = histogram(y(:,p));
            axis tight
            h.FaceColor         = 'b';
            h.EdgeColor         = 'b';
            h.EdgeAlpha         = 0;
            h.FaceAlpha         = alpha;
            yy                  = max(h.Values);
            hold on
            xh                  = [Pout(p) Pout(p)];
            yh                  = [0 yy];
            plot(xh,yh,'r','LineWidth',LineWidth);
            axis('tight')
            set(gca,'Box','off')
            set(gca,'FontWeight','demi')
            set(gca,'FontSize',fontsize)
            set(gca,'FontName','Cambria')
            set(gca,'Linewidth',LineWidth)
            x1                  = get(gca,'Xtick');
            y1                  = get(gca,'Ytick');
            set(gca,'XTickLabel',x1.*1e2)
            set(gca,'YTickLabel',y1./1e3)
            text(1.01,0.08,'10^{-2}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.95,0.9,'n_1_1','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
            % -------------------------- n12 ------------------------------------------
            p                   = 2;
            axes('position',sub_pos{2,4});
            h                   = histogram(y(:,p));
            axis tight
            h.FaceColor         = 'b';
            h.EdgeColor         = 'b';
            h.EdgeAlpha         = 0;
            h.FaceAlpha         = alpha;
            yy                  = max(h.Values);
            hold on
            xh                  = [Pout(p) Pout(p)];
            yh                  = [0 yy];
            plot(xh,yh,'r','LineWidth',LineWidth);
            axis('tight')
            set(gca,'Box','off')
            set(gca,'FontWeight','demi')
            set(gca,'FontSize',fontsize)
            set(gca,'FontName','Cambria')
            set(gca,'Linewidth',LineWidth)
            x1                  = get(gca,'Xtick');
            y1                  = get(gca,'Ytick');
            set(gca,'XTickLabel',x1.*1e2)
            set(gca,'YTickLabel',y1./1e3)
            text(1.01,0.08,'10^{-2}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.95,0.9,'n_1_2','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
            % -------------------------- n22 ------------------------------------------
            p                   = 3;
            axes('position',sub_pos{1,3});
            h                   = histogram(y(:,p));
            axis tight
            h.FaceColor         = 'b';
            h.EdgeColor         = 'b';
            h.EdgeAlpha         = 0;
            h.FaceAlpha         = alpha;
            yy                  = max(h.Values);
            hold on
            xh                  = [Pout(p) Pout(p)];
            yh                  = [0 yy];
            plot(xh,yh,'r','LineWidth',LineWidth);
            axis('tight')
            set(gca,'Box','off')
            set(gca,'FontWeight','demi')
            set(gca,'FontSize',fontsize)
            set(gca,'FontName','Cambria')
            set(gca,'Linewidth',LineWidth)
            x1                  = get(gca,'Xtick');
            y1                  = get(gca,'Ytick');
            set(gca,'XTickLabel',x1.*1e10)
            set(gca,'YTickLabel',y1./1e3)
            text(1.01,0.08,'10^{-10}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.95,0.9,'n_2_2','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
            % -------------------------- n31 ------------------------------------------
            p                   = 4;
            axes('position',sub_pos{2,3});
            h                   = histogram(y(:,p));
            axis tight
            h.FaceColor         = 'b';
            h.EdgeColor         = 'b';
            h.EdgeAlpha         = 0;
            h.FaceAlpha         = alpha;
            yy                  = max(h.Values);
            hold on
            xh                  = [Pout(p) Pout(p)];
            yh                  = [0 yy];
            plot(xh,yh,'r','LineWidth',LineWidth);
            axis('tight')
            set(gca,'Box','off')
            set(gca,'FontWeight','demi')
            set(gca,'FontSize',fontsize)
            set(gca,'FontName','Cambria')
            set(gca,'Linewidth',LineWidth)
            x1                  = get(gca,'Xtick');
            y1                  = get(gca,'Ytick');
            set(gca,'XTickLabel',x1.*1e4)
            set(gca,'YTickLabel',y1./1e3)
            text(1.01,0.08,'10^{-4}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.95,0.9,'n_3_1','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
            % -------------------------- n32 ------------------------------------------
            p                   = 5;
            axes('position',sub_pos{1,2});
            h                   = histogram(y(:,p));
            axis tight
            h.FaceColor         = 'b';
            h.EdgeColor         = 'b';
            h.EdgeAlpha         = 0;
            h.FaceAlpha         = alpha;
            yy                  = max(h.Values);
            hold on
            xh                  = [Pout(p) Pout(p)];
            yh                  = [0 yy];
            plot(xh,yh,'r','LineWidth',LineWidth);
            axis('tight')
            set(gca,'Box','off')
            set(gca,'FontWeight','demi')
            set(gca,'FontSize',fontsize)
            set(gca,'FontName','Cambria')
            set(gca,'Linewidth',LineWidth)
            x1                  = get(gca,'Xtick');
            y1                  = get(gca,'Ytick');
            set(gca,'XTickLabel',x1.*1e4)
            set(gca,'YTickLabel',y1./1e3)
            text(1.01,0.08,'10^{-4}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.95,0.9,'n_3_2','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
            % -------------------------- d1 ------------------------------------------
            p                   = 6;
            axes('position',sub_pos{2,2});
            h                   = histogram(y(:,p));
            axis tight
            h.FaceColor         = 'b';
            h.EdgeColor         = 'b';
            h.EdgeAlpha         = 0;
            h.FaceAlpha         = alpha;
            yy                  = max(h.Values);
            hold on
            xh                  = [Pout(p) Pout(p)];
            yh                  = [0 yy];
            plot(xh,yh,'r','LineWidth',LineWidth);
            axis('tight')
            set(gca,'Box','off')
            set(gca,'FontWeight','demi')
            set(gca,'FontSize',fontsize)
            set(gca,'FontName','Cambria')
            set(gca,'Linewidth',LineWidth)
            y1                  = get(gca,'Ytick');
            set(gca,'YTickLabel',y1./1e3)
            text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.95,0.9,'d_1','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
            % -------------------------- d2 ------------------------------------------
            p                   = 7;
            axes('position',sub_pos{1,1});
            h                   = histogram(y(:,p));
            axis tight
            h.FaceColor         = 'b';
            h.EdgeColor         = 'b';
            h.EdgeAlpha         = 0;
            h.FaceAlpha         = alpha;
            yy                  = max(h.Values);
            hold on
            xh                  = [Pout(p) Pout(p)];
            yh                  = [0 yy];
            plot(xh,yh,'r','LineWidth',LineWidth);
            axis('tight')
            set(gca,'Box','off')
            set(gca,'FontWeight','demi')
            set(gca,'FontSize',fontsize)
            set(gca,'FontName','Cambria')
            set(gca,'Linewidth',LineWidth)
            x1                  = get(gca,'Xtick');
            y1                  = get(gca,'Ytick');
            set(gca,'XTickLabel',x1.*1e1)
            set(gca,'YTickLabel',y1./1e3)
            text(1.01,0.08,'10^{-1}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.95,0.9,'d_2','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
            %% --------------------------------------------------------------------
            % Thermal network to transfer function
            % ---------------------------------------------------------------------
        case 'rctf' % plot RC2TF
            plotheight          = 15;
            plotwidth           = 15;
            subplotsx           = 2;
            subplotsy           = 4;
            leftedge            = 1.1;
            rightedge           = 1;
            topedge             = 1.1;
            bottomedge          = 1;
            spacex              = 1.1;
            spacey              = 1.1;
            fontsize            = 16;
            LineWidth           = 2;
            alpha               = 0.3;
            sub_pos             = subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);
            f                   = figure('visible','on','units','normalized','outerposition',[0 0 0.6 1]);
            clf(f);
            set(gcf, 'PaperUnits', 'centimeters');
            set(gcf, 'PaperSize', [plotwidth plotheight]);
            set(gcf, 'PaperPositionMode', 'manual');
            set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
            set(gcf,'color','w');
            % -------------------------- n11 ------------------------------------------
            p                   = 1;
            axes('position',sub_pos{1,4});
            h                   = histogram(y(:,p));
            axis tight
            h.FaceColor         = 'b';
            h.EdgeColor         = 'b';
            h.EdgeAlpha         = 0;
            h.FaceAlpha         = alpha;
            yy                  = max(h.Values);
            hold on
            xh                  = [Pout(p) Pout(p)];
            yh                  = [0 yy];
            plot(xh,yh,'r','LineWidth',LineWidth);
            axis('tight')
            set(gca,'Box','off')
            set(gca,'FontWeight','demi')
            set(gca,'FontSize',fontsize)
            set(gca,'FontName','Cambria')
            set(gca,'Linewidth',LineWidth)
            x1                  = get(gca,'Xtick');
            y1                  = get(gca,'Ytick');
            set(gca,'XTickLabel',x1.*1e2)
            set(gca,'YTickLabel',y1./1e3)
            text(1.01,0.08,'10^{-2}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.95,0.9,'n_1_1','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
            % -------------------------- n12 ------------------------------------------
            p                   = 2;
            axes('position',sub_pos{2,4});
            h                   = histogram(y(:,p));
            axis tight
            h.FaceColor         = 'b';
            h.EdgeColor         = 'b';
            h.EdgeAlpha         = 0;
            h.FaceAlpha         = alpha;
            yy                  = max(h.Values);
            hold on
            xh                  = [Pout(p) Pout(p)];
            yh                  = [0 yy];
            plot(xh,yh,'r','LineWidth',LineWidth);
            axis('tight')
            set(gca,'Box','off')
            set(gca,'FontWeight','demi')
            set(gca,'FontSize',fontsize)
            set(gca,'FontName','Cambria')
            set(gca,'Linewidth',LineWidth)
            x1                  = get(gca,'Xtick');
            y1                  = get(gca,'Ytick');
            set(gca,'XTickLabel',x1.*1e2)
            set(gca,'YTickLabel',y1./1e3)
            text(1.01,0.08,'10^{-2}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.95,0.9,'n_1_2','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
            % -------------------------- n22 ------------------------------------------
            p                   = 3;
            axes('position',sub_pos{1,3});
            h                   = histogram(y(:,p));
            axis tight
            h.FaceColor         = 'b';
            h.EdgeColor         = 'b';
            h.EdgeAlpha         = 0;
            h.FaceAlpha         = alpha;
            yy                  = max(h.Values);
            hold on
            xh                  = [Pout(p) Pout(p)];
            yh                  = [0 yy];
            plot(xh,yh,'r','LineWidth',LineWidth);
            axis('tight')
            set(gca,'Box','off')
            set(gca,'FontWeight','demi')
            set(gca,'FontSize',fontsize)
            set(gca,'FontName','Cambria')
            set(gca,'Linewidth',LineWidth)
            x1                  = get(gca,'Xtick');
            y1                  = get(gca,'Ytick');
            set(gca,'XTickLabel',x1.*1e10)
            set(gca,'YTickLabel',y1./1e3)
            text(1.01,0.08,'10^{-10}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.95,0.9,'n_2_2','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
            % -------------------------- n31 ------------------------------------------
            p                   = 4;
            axes('position',sub_pos{2,3});
            h                   = histogram(y(:,p));
            axis tight
            h.FaceColor         = 'b';
            h.EdgeColor         = 'b';
            h.EdgeAlpha         = 0;
            h.FaceAlpha         = alpha;
            yy                  = max(h.Values);
            hold on
            xh                  = [Pout(p) Pout(p)];
            yh                  = [0 yy];
            plot(xh,yh,'r','LineWidth',LineWidth);
            axis('tight')
            set(gca,'Box','off')
            set(gca,'FontWeight','demi')
            set(gca,'FontSize',fontsize)
            set(gca,'FontName','Cambria')
            set(gca,'Linewidth',LineWidth)
            x1                  = get(gca,'Xtick');
            y1                  = get(gca,'Ytick');
            set(gca,'XTickLabel',x1.*1e4)
            set(gca,'YTickLabel',y1./1e3)
            text(1.01,0.08,'10^{-4}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.95,0.9,'n_3_1','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
            % -------------------------- n32 ------------------------------------------
            p                   = 5;
            axes('position',sub_pos{1,2});
            h                   = histogram(y(:,p));
            axis tight
            h.FaceColor         = 'b';
            h.EdgeColor         = 'b';
            h.EdgeAlpha         = 0;
            h.FaceAlpha         = alpha;
            yy                  = max(h.Values);
            hold on
            xh                  = [Pout(p) Pout(p)];
            yh                  = [0 yy];
            plot(xh,yh,'r','LineWidth',LineWidth);
            axis('tight')
            set(gca,'Box','off')
            set(gca,'FontWeight','demi')
            set(gca,'FontSize',fontsize)
            set(gca,'FontName','Cambria')
            set(gca,'Linewidth',LineWidth)
            x1                  = get(gca,'Xtick');
            y1                  = get(gca,'Ytick');
            set(gca,'XTickLabel',x1.*1e4)
            set(gca,'YTickLabel',y1./1e3)
            text(1.01,0.08,'10^{-4}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.95,0.9,'n_3_2','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
            % -------------------------- d1 ------------------------------------------
            p                   = 6;
            axes('position',sub_pos{2,2});
            h                   = histogram(y(:,p));
            axis tight
            h.FaceColor         = 'b';
            h.EdgeColor         = 'b';
            h.EdgeAlpha         = 0;
            h.FaceAlpha         = alpha;
            yy                  = max(h.Values);
            hold on
            xh                  = [Pout(p) Pout(p)];
            yh                  = [0 yy];
            plot(xh,yh,'r','LineWidth',LineWidth);
            axis('tight')
            set(gca,'Box','off')
            set(gca,'FontWeight','demi')
            set(gca,'FontSize',fontsize)
            set(gca,'FontName','Cambria')
            set(gca,'Linewidth',LineWidth)
            y1                  = get(gca,'Ytick');
            set(gca,'YTickLabel',y1./1e3)
            text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.95,0.9,'d_1','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
            % -------------------------- d2 ------------------------------------------
            p                   = 7;
            axes('position',sub_pos{1,1});
            h                   = histogram(y(:,p));
            axis tight
            h.FaceColor         = 'b';
            h.EdgeColor         = 'b';
            h.EdgeAlpha         = 0;
            h.FaceAlpha         = alpha;
            yy                  = max(h.Values);
            hold on
            xh                  = [Pout(p) Pout(p)];
            yh                  = [0 yy];
            plot(xh,yh,'r','LineWidth',LineWidth);
            axis('tight')
            set(gca,'Box','off')
            set(gca,'FontWeight','demi')
            set(gca,'FontSize',fontsize)
            set(gca,'FontName','Cambria')
            set(gca,'Linewidth',LineWidth)
            x1                  = get(gca,'Xtick');
            y1                  = get(gca,'Ytick');
            set(gca,'XTickLabel',x1.*1e1)
            set(gca,'YTickLabel',y1./1e3)
            text(1.01,0.08,'10^{-1}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.95,0.9,'d_2','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
    %% --------------------------------------------------------------------
    % Transfer function to state-space
    % ---------------------------------------------------------------------
        case 'tfss'
            plotheight          = 15;
            plotwidth           = 15;
            subplotsx           = 2;
            subplotsy           = 3;
            leftedge            = 1.1;
            rightedge           = 1;
            topedge             = 1.1;
            bottomedge          = 1;
            spacex              = 1.1;
            spacey              = 1.1;
            fontsize            = 16;
            LineWidth           = 2;
            alpha               = 0.3;
            sub_pos             = subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);
            f                   = figure('visible','on','units','normalized','outerposition',[0 0 0.6 1]);
            clf(f);
            set(gcf, 'PaperUnits', 'centimeters');
            set(gcf, 'PaperSize', [plotwidth plotheight]);
            set(gcf, 'PaperPositionMode', 'manual');
            set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
            set(gcf,'color','w');
            % -------------------------- a21 ------------------------------------------
            p                   = 1;
            axes('position',sub_pos{1,3});
            h                   = histogram(y(:,p));
            axis tight
            h.FaceColor         = 'b';
            h.EdgeColor         = 'b';
            h.EdgeAlpha         = 0;
            h.FaceAlpha         = alpha;
            yy                  = max(h.Values);
            hold on
            xh                  = [Pout(p) Pout(p)];
            yh                  = [0 yy];
            plot(xh,yh,'r','LineWidth',LineWidth);
            axis('tight')
            set(gca,'Box','off')
            set(gca,'FontWeight','demi')
            set(gca,'FontSize',fontsize)
            set(gca,'FontName','Cambria')
            set(gca,'Linewidth',LineWidth)
            x1                  = get(gca,'Xtick');
            y1                  = get(gca,'Ytick');
            set(gca,'XTickLabel',x1.*1e1)
            set(gca,'YTickLabel',y1./1e3)
            text(1.01,0.08,'10^{-1}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.95,0.9,'a_2_1','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
            % -------------------------- a22 ------------------------------------------
            p                   = 2;
            axes('position',sub_pos{2,3});
            h                   = histogram(y(:,p));
            axis tight
            h.FaceColor         = 'b';
            h.EdgeColor         = 'b';
            h.EdgeAlpha         = 0;
            h.FaceAlpha         = alpha;
            yy                  = max(h.Values);
            hold on
            xh                  = [Pout(p) Pout(p)];
            yh                  = [0 yy];
            plot(xh,yh,'r','LineWidth',LineWidth);
            axis('tight')
            set(gca,'Box','off')
            set(gca,'FontWeight','demi')
            set(gca,'FontSize',fontsize)
            set(gca,'FontName','Cambria')
            set(gca,'Linewidth',LineWidth)
            x1                  = get(gca,'Xtick');
            y1                  = get(gca,'Ytick');
            set(gca,'XTickLabel',x1.*1e3)
            set(gca,'YTickLabel',y1./1e3)
            text(1.01,0.08,'10^{-3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.95,0.9,'a_2_2','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
            % -------------------------- b11 ------------------------------------------
            p                   = 3;
            axes('position',sub_pos{1,2});
            h                   = histogram(y(:,p));
            axis tight
            h.FaceColor         = 'b';
            h.EdgeColor         = 'b';
            h.EdgeAlpha         = 0;
            h.FaceAlpha         = alpha;
            yy                  = max(h.Values);
            hold on
            xh                  = [Pout(p) Pout(p)];
            yh                  = [0 yy];
            plot(xh,yh,'r','LineWidth',LineWidth);
            axis('tight')
            set(gca,'Box','off')
            set(gca,'FontWeight','demi')
            set(gca,'FontSize',fontsize)
            set(gca,'FontName','Cambria')
            set(gca,'Linewidth',LineWidth)
            x1                  = get(gca,'Xtick');
            y1                  = get(gca,'Ytick');
            set(gca,'XTickLabel',x1.*1e4)
            set(gca,'YTickLabel',y1./1e3)
            text(1.01,0.08,'10^{-4}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.95,0.9,'b_1_1','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
            % -------------------------- b13 ------------------------------------------
            p                   = 4;
            axes('position',sub_pos{2,2});
            h                   = histogram(y(:,p));
            axis tight
            h.FaceColor         = 'b';
            h.EdgeColor         = 'b';
            h.EdgeAlpha         = 0;
            h.FaceAlpha         = alpha;
            yy                  = max(h.Values);
            hold on
            xh                  = [Pout(p) Pout(p)];
            yh                  = [0 yy];
            plot(xh,yh,'r','LineWidth',LineWidth);
            axis('tight')
            set(gca,'Box','off')
            set(gca,'FontWeight','demi')
            set(gca,'FontSize',fontsize)
            set(gca,'FontName','Cambria')
            set(gca,'Linewidth',LineWidth)
            x1                  = get(gca,'Xtick');
            y1                  = get(gca,'Ytick');
            set(gca,'XTickLabel',x1.*1e5)
            set(gca,'YTickLabel',y1./1e3)
            text(1.01,0.08,'10^{-5}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.95,0.9,'b_1_3','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
            % -------------------------- b21 ------------------------------------------
            p                   = 5;
            axes('position',sub_pos{1,1});
            h                   = histogram(y(:,p));
            axis tight
            h.FaceColor         = 'b';
            h.EdgeColor         = 'b';
            h.EdgeAlpha         = 0;
            h.FaceAlpha         = alpha;
            yy                  = max(h.Values);
            hold on
            xh                  = [Pout(p) Pout(p)];
            yh                  = [0 yy];
            plot(xh,yh,'r','LineWidth',LineWidth);
            axis('tight')
            set(gca,'Box','off')
            set(gca,'FontWeight','demi')
            set(gca,'FontSize',fontsize)
            set(gca,'FontName','Cambria')
            set(gca,'Linewidth',LineWidth)
            x1                  = get(gca,'Xtick');
            y1                  = get(gca,'Ytick');
            set(gca,'XTickLabel',x1.*1e2)
            set(gca,'YTickLabel',y1./1e3)
            text(1.01,0.08,'10^{-2}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.95,0.9,'b_2_1','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
            % -------------------------- b22 ------------------------------------------
            p                   = 6;
            axes('position',sub_pos{2,1});
            h                   = histogram(y(:,p));
            axis tight
            h.FaceColor         = 'b';
            h.EdgeColor         = 'b';
            h.EdgeAlpha         = 0;
            h.FaceAlpha         = alpha;
            yy                  = max(h.Values);
            hold on
            xh                  = [Pout(p) Pout(p)];
            yh                  = [0 yy];
            plot(xh,yh,'r','LineWidth',LineWidth);
            axis('tight')
            set(gca,'Box','off')
            set(gca,'FontWeight','demi')
            set(gca,'FontSize',fontsize)
            set(gca,'FontName','Cambria')
            set(gca,'Linewidth',LineWidth)
            x1                  = get(gca,'Xtick');
            y1                  = get(gca,'Ytick');
            set(gca,'XTickLabel',x1.*1e9)
            set(gca,'YTickLabel',y1./1e3)
            text(1.01,0.08,'10^{-9}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.95,0.9,'b_2_2','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
    %% --------------------------------------------------------------------
    % State-space to thermal network
    % ---------------------------------------------------------------------
        case 'ssrc' 
            plotheight          = 15;
            plotwidth           = 15;
            subplotsx           = 2;
            subplotsy           = 3;
            leftedge            = 1.1;
            rightedge           = 1;
            topedge             = 1.1;
            bottomedge          = 1;
            spacex              = 1.1;
            spacey              = 1.1;
            fontsize            = 16;
            LineWidth           = 2;
            alpha               = 0.3;
            sub_pos             = subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);
            f                   = figure('visible','on','units','normalized','outerposition',[0 0 0.6 1]);
            clf(f);
            set(gcf, 'PaperUnits', 'centimeters');
            set(gcf, 'PaperSize', [plotwidth plotheight]);
            set(gcf, 'PaperPositionMode', 'manual');
            set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
            set(gcf,'color','w');
            % -------------------------- Rco ------------------------------------------
            p                   = 1;
            axes('position',sub_pos{1,3});
            h                   = histogram(y(:,p));
            axis tight
            h.FaceColor         = 'b';
            h.EdgeColor         = 'b';
            h.EdgeAlpha         = 0;
            h.FaceAlpha         = alpha;
            yy                  = max(h.Values);
            hold on
            xh                  = [Pout(p) Pout(p)];
            yh                  = [0 yy];
            plot(xh,yh,'r','LineWidth',LineWidth);
            axis('tight')
            xlim([2e-3 7e-3])
            set(gca,'Box','off')
            set(gca,'FontWeight','demi')
            set(gca,'FontSize',fontsize)
            set(gca,'FontName','Cambria')
            set(gca,'Linewidth',LineWidth)
            x1                  = get(gca,'Xtick');
            y1                  = get(gca,'Ytick');
            set(gca,'XTickLabel',x1.*1e3)
            set(gca,'YTickLabel',y1./1e3)
            text(1.01,0.08,'10^{-3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.95,0.9,'R_c_o','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
            % -------------------------- Rci ------------------------------------------
            p                   = 2;
            axes('position',sub_pos{2,3});
            h                   = histogram(y(:,p));
            axis tight
            h.FaceColor         = 'b';
            h.EdgeColor         = 'b';
            h.EdgeAlpha         = 0;
            h.FaceAlpha         = alpha;
            yy                  = max(h.Values);
            hold on
            xh                  = [Pout(p) Pout(p)];
            yh                  = [0 yy];
            plot(xh,yh,'r','LineWidth',LineWidth);
            axis('tight')
            xlim([-0.3 0.3])
            set(gca,'Box','off')
            set(gca,'FontWeight','demi')
            set(gca,'FontSize',fontsize)
            set(gca,'FontName','Cambria')
            set(gca,'Linewidth',LineWidth)
            x1                  = get(gca,'Xtick');
            y1                  = get(gca,'Ytick');
            set(gca,'XTickLabel',x1.*1e1)
            set(gca,'YTickLabel',y1./1e3)
            text(1.01,0.08,'10^{-1}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.95,0.9,'R_c_i','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
            % -------------------------- Rv ------------------------------------------
            p                   = 3;
            axes('position',sub_pos{1,2});
            h                   = histogram(y(:,p));
            axis tight
            h.FaceColor         = 'b';
            h.EdgeColor         = 'b';
            h.EdgeAlpha         = 0;
            h.FaceAlpha         = alpha;
            yy                  = max(h.Values);
            hold on
            xh                  = [Pout(p) Pout(p)];
            yh                  = [0 yy];
            plot(xh,yh,'r','LineWidth',LineWidth);
            axis('tight')
            set(gca,'Box','off')
            set(gca,'FontWeight','demi')
            set(gca,'FontSize',fontsize)
            set(gca,'FontName','Cambria')
            set(gca,'Linewidth',LineWidth)
            x1                  = get(gca,'Xtick');
            y1                  = get(gca,'Ytick');
            set(gca,'XTickLabel',x1.*1e2)
            set(gca,'YTickLabel',y1./1e3)
            text(1.01,0.08,'10^{-2}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.95,0.9,'R_v','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
            % -------------------------- Rw ------------------------------------------
            p                   = 4;
            axes('position',sub_pos{2,2});
            h                   = histogram(y(:,p));
            axis tight
            h.FaceColor         = 'b';
            h.EdgeColor         = 'b';
            h.EdgeAlpha         = 0;
            h.FaceAlpha         = alpha;
            yy                  = max(h.Values);
            hold on
            xh                  = [Pout(p) Pout(p)];
            yh                  = [0 yy];
            plot(xh,yh,'r','LineWidth',LineWidth);
            axis('tight')
            xlim([-2 3])
            set(gca,'Box','off')
            set(gca,'FontWeight','demi')
            set(gca,'FontSize',fontsize)
            set(gca,'FontName','Cambria')
            set(gca,'Linewidth',LineWidth)
            y1                  = get(gca,'Ytick');
            set(gca,'YTickLabel',y1./1e3)
            text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.95,0.9,'R_w','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
            % -------------------------- Cw ------------------------------------------
            p                   = 5;
            axes('position',sub_pos{1,1});
            h                   = histogram(y(:,p));
            axis tight
            h.FaceColor         = 'b';
            h.EdgeColor         = 'b';
            h.EdgeAlpha         = 0;
            h.FaceAlpha         = alpha;
            yy                  = max(h.Values);
            hold on
            xh                  = [Pout(p) Pout(p)];
            yh                  = [0 yy];
            plot(xh,yh,'r','LineWidth',LineWidth);
            axis('tight')
            set(gca,'Box','off')
            set(gca,'FontWeight','demi')
            set(gca,'FontSize',fontsize)
            set(gca,'FontName','Cambria')
            set(gca,'Linewidth',LineWidth)
            x1                  = get(gca,'Xtick');
            y1                  = get(gca,'Ytick');
            set(gca,'XTickLabel',x1./1e6)
            set(gca,'YTickLabel',y1./1e3)
            text(1.01,0.08,'10^{6}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.95,0.9,'C_w','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
            % -------------------------- Ca ------------------------------------------
            p                   = 6;
            axes('position',sub_pos{2,1});
            h                   = histogram(y(:,p));
            axis tight
            h.FaceColor         = 'b';
            h.EdgeColor         = 'b';
            h.EdgeAlpha         = 0;
            h.FaceAlpha         = alpha;
            yy                  = max(h.Values);
            hold on
            xh                  = [Pout(p) Pout(p)];
            yh                  = [0 yy];
            plot(xh,yh,'r','LineWidth',LineWidth);
            axis('tight')
            set(gca,'Box','off')
            set(gca,'FontWeight','demi')
            set(gca,'FontSize',fontsize)
            set(gca,'FontName','Cambria')
            set(gca,'Linewidth',LineWidth)
            x1                  = get(gca,'Xtick');
            y1                  = get(gca,'Ytick');
            set(gca,'XTickLabel',x1./1e4)
            set(gca,'YTickLabel',y1./1e3)
            text(1.01,0.08,'10^{4}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
            text(0.95,0.9,'C_a','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
    %% --------------------------------------------------------------------
    % Transfer function to thermal network
    % ---------------------------------------------------------------------
        case 'tfrc'
        plotheight              = 15;
        plotwidth               = 15;
        subplotsx               = 2;
        subplotsy               = 3;
        leftedge                = 1.1;
        rightedge               = 1;
        topedge                 = 1.1;
        bottomedge              = 1;
        spacex                  = 1.1;
        spacey                  = 1.1;
        fontsize                = 16;
        LineWidth               = 2;
        alpha                   = 0.3;
        sub_pos                 = subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);
        f                       = figure('visible','on','units','normalized','outerposition',[0 0 0.6 1]);
        clf(f);
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperSize', [plotwidth plotheight]);
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
        set(gcf,'color','w');
        % -------------------------- Rco ------------------------------------------
        p                       = 1;
        axes('position',sub_pos{1,3});
        h                       = histogram(y(:,p));
        axis tight
        h.FaceColor             = 'b';
        h.EdgeColor             = 'b';
        h.EdgeAlpha             = 0;
        h.FaceAlpha             = alpha;
        yy                      = max(h.Values);
        hold on
        xh                      = [Pout(p) Pout(p)];
        yh                      = [0 yy];
        plot(xh,yh,'r','LineWidth',LineWidth);
        axis('tight')
        xlim([-20e-6 20e-6])
        set(gca,'Box','off')
        set(gca,'FontWeight','demi')
        set(gca,'FontSize',fontsize)
        set(gca,'FontName','Cambria')
        set(gca,'Linewidth',LineWidth)
        x1                      = get(gca,'Xtick');
        y1                      = get(gca,'Ytick');
        set(gca,'XTickLabel',x1.*1e5)
        set(gca,'YTickLabel',y1./1e5)
        text(1.01,0.08,'10^{-5}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
        text(0.01,1.07,'10^{5}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
        text(0.95,0.9,'R_c_o','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
        % -------------------------- Rci ------------------------------------------
        p                       = 2;
        axes('position',sub_pos{2,3});
        h                       = histogram(y(:,p));
        axis tight
        h.FaceColor             = 'b';
        h.EdgeColor             = 'b';
        h.EdgeAlpha             = 0;
        h.FaceAlpha             = alpha;
        yy                      = max(h.Values);
        hold on
        xh                      = [Pout(p) Pout(p)];
        yh                      = [0 yy];
        plot(xh,yh,'r','LineWidth',LineWidth);
        axis('tight')
        xlim([-6 11])
        set(gca,'Box','off')
        set(gca,'FontWeight','demi')
        set(gca,'FontSize',fontsize)
        set(gca,'FontName','Cambria')
        set(gca,'Linewidth',LineWidth)
        y1                      = get(gca,'Ytick');
        set(gca,'YTickLabel',y1./1e4)
        text(0.01,1.07,'10^{4}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
        text(0.95,0.9,'R_c_i','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
        % -------------------------- Rv ------------------------------------------
        p                       = 3;
        axes('position',sub_pos{1,2});
        h                       = histogram(y(:,p));
        axis tight
        h.FaceColor             = 'b';
        h.EdgeColor             = 'b';
        h.EdgeAlpha             = 0;
        h.FaceAlpha             = alpha;
        yy                      = max(h.Values);
        hold on
        xh                      = [Pout(p) Pout(p)];
        yh                      = [0 yy];
        plot(xh,yh,'r','LineWidth',LineWidth);
        axis('tight')
        set(gca,'Box','off')
        set(gca,'FontWeight','demi')
        set(gca,'FontSize',fontsize)
        set(gca,'FontName','Cambria')
        set(gca,'Linewidth',LineWidth)
        x1                      = get(gca,'Xtick');
        y1                      = get(gca,'Ytick');
        set(gca,'XTickLabel',x1.*1e2)
        set(gca,'YTickLabel',y1./1e3)
        text(1.01,0.08,'10^{-2}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
        text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
        text(0.95,0.9,'R_v','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
        % -------------------------- Rw ------------------------------------------
        p                       = 4;
        axes('position',sub_pos{2,2});
        h                       = histogram(y(:,p));
        axis tight
        h.FaceColor             = 'b';
        h.EdgeColor             = 'b';
        h.EdgeAlpha             = 0;
        h.FaceAlpha             = alpha;
        yy                      = max(h.Values);
        hold on
        xh                      = [Pout(p) Pout(p)];
        yh                      = [0 yy];
        plot(xh,yh,'r','LineWidth',LineWidth);
        axis('tight')
        xlim([-7 5])
        set(gca,'Box','off')
        set(gca,'FontWeight','demi')
        set(gca,'FontSize',fontsize)
        set(gca,'FontName','Cambria')
        set(gca,'Linewidth',LineWidth)
        y1                      = get(gca,'Ytick');
        set(gca,'YTickLabel',y1./1e4)
        text(0.01,1.07,'10^{4}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
        text(0.95,0.9,'R_w','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
        % -------------------------- Cw ------------------------------------------
        p                       = 5;
        axes('position',sub_pos{1,1});
        h                       = histogram(y(:,p));
        axis tight
        h.FaceColor             = 'b';
        h.EdgeColor             = 'b';
        h.EdgeAlpha             = 0;
        h.FaceAlpha             = alpha;
        yy                      = max(h.Values);
        hold on
        xh                      = [Pout(p) Pout(p)];
        yh                      = [0 yy];
        plot(xh,yh,'r','LineWidth',LineWidth);
        axis('tight')
        xlim([-5e3 5e3])
        set(gca,'Box','off')
        set(gca,'FontWeight','demi')
        set(gca,'FontSize',fontsize)
        set(gca,'FontName','Cambria')
        set(gca,'Linewidth',LineWidth)
        x1                      = get(gca,'Xtick');
        y1                      = get(gca,'Ytick');
        set(gca,'XTickLabel',x1./1e3)
        set(gca,'YTickLabel',y1./1e5)
        text(1.01,0.08,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
        text(0.01,1.07,'10^{5}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
        text(0.95,0.9,'C_w','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
        % -------------------------- Ca ------------------------------------------
        p                       = 6;
        axes('position',sub_pos{2,1});
        h                       = histogram(y(:,p));
        axis tight
        h.FaceColor             = 'b';
        h.EdgeColor             = 'b';
        h.EdgeAlpha             = 0;
        h.FaceAlpha             = alpha;
        yy                      = max(h.Values);
        hold on
        xh                      = [Pout(p) Pout(p)];
        yh                      = [0 yy];
        plot(xh,yh,'r','LineWidth',LineWidth);
        axis('tight')
        set(gca,'Box','off')
        set(gca,'FontWeight','demi')
        set(gca,'FontSize',fontsize)
        set(gca,'FontName','Cambria')
        set(gca,'Linewidth',LineWidth)
        x1                      = get(gca,'Xtick');
        y1                      = get(gca,'Ytick');
        set(gca,'XTickLabel',x1./1e4)
        set(gca,'YTickLabel',y1./1e3)
        text(1.01,0.08,'10^{4}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
        text(0.01,1.07,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
        text(0.95,0.9,'C_a','FontSize',18,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria','Color','red');
    end
end
