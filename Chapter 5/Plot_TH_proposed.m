% Functions: 
%   - linspecer: Beautiful and distinguishable line colors + colormap, Jonathan C. Lansey: https://fr.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-+-colormap?focused=5372538&tab=function
%   - subplot_pos: Custom subplot, Pierre Martineau, http://p-martineau.com/perfect-subplot-in-matlab/
%   - export_fig: export figure as it appears on screen, https://github.com/altmany/export_fig
load('TH_proposed_data.mat')
%% ------------------------------------------------------------------------
% Prior and posterior of parameters, Figure 5.10
% -------------------------------------------------------------------------
plotheight  = 15;
plotwidth   = 15;
subplotsx   = 2;
subplotsy   = 5;
leftedge    = 1.1;
rightedge   = 1;
topedge     = 1.1;
bottomedge  = 1;
spacex      = 1.1;
spacey      = 1;
fontsize    = 16;
sub_pos     = subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);
f           = figure('visible','on','units','normalized','outerposition',[0 0 1 1]);
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'color','w');
% --------------------------------------- % 
axes('position',sub_pos{1,5});
p           = 1;
[nb,xb]     = hist(sys.m(p,:),sqrt(Nprt));
bh          = bar(xb,nb,'hist');
set(bh,'facecolor','w','edgecolor','k');
hold on
[nb,xb]     = hist(m(p,:),sqrt(Nprt));
bh          = bar(xb,nb,'hist');
set(bh,'facecolor','b','edgecolor','b');
ylim([0 270])
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
x1          = get(gca,'Xtick');
set(gca,'XTickLabel',x1.*100)
text(0.05,0.9,'R_o','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
text(1.01,0.08,'10^{-2}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
% --------------------------------------- % 
axes('position',sub_pos{2,5});
p           = 2;
[nb,xb]     = hist(sys.m(p,:),sqrt(Nprt));
bh          = bar(xb,nb,'hist');
set(bh,'facecolor','w','edgecolor','k');
hold on
[nb,xb]     = hist(m(p,:),sqrt(Nprt));
bh          = bar(xb,nb,'hist');
set(bh,'facecolor','b','edgecolor','b');
ylim([0 270])
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
x1          = get(gca,'Xtick');
set(gca,'XTickLabel',x1.*100)
text(0.05,0.9,'R_i','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
text(1.01,0.08,'10^{-2}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
% --------------------------------------- % 
axes('position',sub_pos{1,4});
p           = 3;
[nb,xb]     = hist(sys.m(p,:),sqrt(Nprt));
bh          = bar(xb,nb,'hist');
set(bh,'facecolor','w','edgecolor','k');
hold on
[nb,xb]     = hist(m(p,:),sqrt(Nprt));
bh          = bar(xb,nb,'hist');
set(bh,'facecolor','b','edgecolor','b');
ylim([0 270])
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
x1          = get(gca,'Xtick');
set(gca,'XTickLabel',x1.*1000)
text(0.05,0.9,'R_m','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
text(1.01,0.08,'10^{-3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
% --------------------------------------- % 
axes('position',sub_pos{2,4});
p           = 4;
[nb,xb]     = hist(sys.m(p,:),sqrt(Nprt));
bh          = bar(xb,nb,'hist');
set(bh,'facecolor','w','edgecolor','k');
hold on
[nb,xb]     = hist(m(p,:),sqrt(Nprt));
bh          = bar(xb,nb,'hist');
set(bh,'facecolor','b','edgecolor','b');
ylim([0 270])
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
x1          = get(gca,'Xtick');
set(gca,'XTickLabel',x1.*1000)
text(0.05,0.9,'R_z','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
text(1.01,0.08,'10^{-3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
% --------------------------------------- % 
axes('position',sub_pos{1,3});
p           = 9;
[nb,xb]     = hist(sys.m(p,:),sqrt(Nprt));
bh          = bar(xb,nb,'hist');
set(bh,'facecolor','w','edgecolor','k');
hold on
[nb,xb]     = hist(m(p,:),sqrt(Nprt));
bh          = bar(xb,nb,'hist');
set(bh,'facecolor','b','edgecolor','b');
ylim([0 270])
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
x1          = get(gca,'Xtick');
set(gca,'XTickLabel',x1.*1000)
text(0.05,0.9,'R_s_o','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
text(1.01,0.08,'10^{-3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
% --------------------------------------- % 
axes('position',sub_pos{2,3});
p           = 5;
[nb,xb]     = hist(sys.m(p,:),sqrt(Nprt));
bh          = bar(xb,nb,'hist');
set(bh,'facecolor','w','edgecolor','k');
hold on
[nb,xb]     = hist(m(p,:),sqrt(Nprt));
bh          = bar(xb,nb,'hist');
set(bh,'facecolor','b','edgecolor','b');
ylim([0 270])
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
x1          = get(gca,'Xtick');
set(gca,'XTickLabel',x1.*10)
text(0.05,0.9,'C_w','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
text(1.01,0.08,'10^{7}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
% --------------------------------------- % 
axes('position',sub_pos{1,2});
p           = 6;
[nb,xb]     = hist(sys.m(p,:),sqrt(Nprt));
bh          = bar(xb,nb,'hist');
set(bh,'facecolor','w','edgecolor','k');
hold on
[nb,xb]     = hist(m(p,:),sqrt(Nprt));
bh          = bar(xb,nb,'hist');
set(bh,'facecolor','b','edgecolor','b');
ylim([0 270])
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
x1          = get(gca,'Xtick');
set(gca,'XTickLabel',x1.*1000)
text(0.05,0.9,'C_i','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
text(1.01,0.08,'10^{5}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
% --------------------------------------- % 
axes('position',sub_pos{2,2});
p           = 7;
[nb,xb]     = hist(sys.m(p,:),sqrt(Nprt));
bh          = bar(xb,nb,'hist');
set(bh,'facecolor','w','edgecolor','k');
hold on
[nb,xb]     = hist(m(p,:),sqrt(Nprt));
bh          = bar(xb,nb,'hist');
set(bh,'facecolor','b','edgecolor','b');
ylim([0 270])
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
x1          = get(gca,'Xtick');
set(gca,'XTickLabel',x1.*10)
text(0.05,0.9,'C_m','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
text(1.01,0.08,'10^{7}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
% --------------------------------------- % 
axes('position',sub_pos{1,1});
p           = 8;
[nb,xb]     = hist(sys.m(p,:),sqrt(Nprt));
bh=bar(xb,nb,'hist');
set(bh,'facecolor','w','edgecolor','k');
hold on
[nb,xb]     = hist(m(p,:),sqrt(Nprt));
bh          = bar(xb,nb,'hist');
set(bh,'facecolor','b','edgecolor','b');
ylim([0 270])
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
x1          = get(gca,'Xtick');
set(gca,'XTickLabel',x1.*10)
text(0.05,0.9,'a_I','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
text(1.01,0.08,'10^{-1}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
%% ------------------------------------------------------------------------
% plot states, Figure 5.9
% -------------------------------------------------------------------------
startDate   = datenum('28-04-2014 09', 'dd-mm-yy HH');
endDate     = datenum('14-05-2014 00', 'dd-mm-yy HH');
xData       = linspace(startDate,endDate,T);

plotheight  = 15;
plotwidth   = 15;
subplotsx   = 1;
subplotsy   = 3;
leftedge    = 1.1;
rightedge   = 1;
topedge     = 1.1;
bottomedge  = 1;
spacex      = 1.1;
spacey      = 1;
fontsize    = 16;
sub_pos     = subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);
f=figure('visible','on','units','normalized','outerposition',[0 0 0.5 1]);
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'color','w');
% --------------------------------------- % 
axes('position',sub_pos{1,3});
p           = 1;
plot(xData,xS(p,:),'b','linewidth',1.5)
datetick('x','dd mmm','keepticks')
axis tight
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
ylabel('°C','fontsize',fontsize)
text(0.025,0.9,'\theta_w','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
% --------------------------------------- % 
axes('position',sub_pos{1,2});
p           = 2;
plot(xData,xS(p,:),'b',xData,sys.y,'r','linewidth',1.5)
datetick('x','dd mmm','keepticks')
axis tight
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
ylabel('°C','fontsize',fontsize)
text(0.025,0.9,'\theta_i','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
% --------------------------------------- % 
axes('position',sub_pos{1,1});
p           = 3;
plot(xData,xS(p,:),'b','linewidth',1.5)
datetick('x','dd mmm','keepticks')
axis tight
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
ylabel('°C','fontsize',fontsize)
text(0.025,0.9,'\theta_m','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
%% ------------------------------------------------------------------------
% plot simulation, Figure 5.13
% -------------------------------------------------------------------------
Np                          = 9;                                            
Nx                          = 3;                                            
Ny                          = 1;                                            
Nu                          = 5;                                            
sys.dim                     = [Np Nx Ny Nu];                                
cutp                        = cut(end)+1:cut(end)+234;
up                          = [Text(cutp) Tzone2_all(cutp) Srw(cutp) Sri(cutp) Qhvac1(cutp)]';
Tp                          = length(cutp);
xp                          = zeros(Nx,Tp);
sig                         = sys.sig;
xp(:,1)                     = xS(1:Nx,end);

mdl                         = feval(sys.fun,xS(Nx+1:end,end),sys);          
F                           = expm(mdl.A*dt);                               
bis                         = F - eye(Nx);                                  
G(:,1:Nu)                   = mdl.A\bis*mdl.B;                              
G(:,Nu+1:2*Nu)              = mdl.A\(-mdl.A\bis + F*dt)*mdl.B;              
Cs                          = sig*sig';                                     
V                           = Cs - F*Cs*F';                                 
Q                           = diag(diag(lyap(mdl.A,V)));                    
Qstd                        = sqrt(diag(Q));
for k = 2:Tp
    afoh                    = (up(:,k) - up(:,k-1))./sys.dt;
    xp(:,k)                 = F*xp(:,k-1) + G(:,1:Nu)*(up(:,k-1) + ...
                    afoh*sys.dt) - G(:,Nu+1:2*Nu)*afoh + Qstd.*randn(Nx,1);
end

startDate                   = datenum('14-05-2014 01', 'dd-mm-yy HH');
endDate                     = datenum('23-05-2014 18', 'dd-mm-yy HH');
xData                       = linspace(startDate,endDate,Tp);

plotheight  = 15;
plotwidth   = 15;
subplotsx   = 1;
subplotsy   = 1;
leftedge    = 1.1;
rightedge   = 1;
topedge     = 1.1;
bottomedge  = 1.2;
spacex      = 1.1;
spacey      = 1;
fontsize    = 16;
sub_pos     = subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);
f           = figure('visible','on','units','normalized','outerposition',[0 0 0.5 0.5]);
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'color','w');
% --------------------------------------- % 
axes('position',sub_pos{1,1});
plot(xData,xp(2,:),'b',xData,Tzone1_all(cutp),'r','linewidth',1.5)
datetick('x','dd mmm','keepticks')
axis tight
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
ylabel('°C','fontsize',16)
text(0.3,0.5,'\theta_i measured','color','r','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
text(0.3,0.4,'\theta_i simulated','color','b','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');