% Functions: 
%   - linspecer: Beautiful and distinguishable line colors + colormap, Jonathan C. Lansey: https://fr.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-+-colormap?focused=5372538&tab=function
%   - subplot_pos: Custom subplot, Pierre Martineau, http://p-martineau.com/perfect-subplot-in-matlab/
%   - export_fig: export figure as it appears on screen, https://github.com/altmany/export_fig
load('MC_simu.mat')
MC          = size(xS1,3);
cut         = 465:840; 
T           = length(cut);
dt          = 60*60;
t_cut       = 0:dt:(T-1)*dt;
tdays       = t_cut/3600/24;
%% ------------------------------------------------------------------------
% Plot Figure 5.7
% -------------------------------------------------------------------------
plotheight  = 15;
plotwidth   = 15;
subplotsx   = 2;
subplotsy   = 5;
leftedge    = 1;
rightedge   = 1;
topedge     = 1;
bottomedge  = 1;
spacex      = 1.2;
spacey      = 1.1;
fontsize    = 16;
sub_pos     = subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);
f           = figure('visible','on','units','normalized','outerposition',[0 0 0.5 1]);
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'color','w');
% ------------ 1/5 -------------------- % 
axes('position',sub_pos{1,5});
mmx         = 0;
for mc = 1:MC
   [f,x]    = ksdensity(mS1(1,:,mc));
   plot(x,f,'b','Linewidth',2)
   hold on
   if max(f) > mmx
       mmx  = max(f);
   end
end
axis tight
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
yy1         = get(gca,'Ytick');
x1          = [par(1) par(1)];
y1          = [0 mmx];
hold on
plot(x1,y1,'r','LineWidth',2)
xx          = get(gca,'Xtick');
set(gca,'Ytick',yy1);
set(gca,'XTickLabel',xx.*100)
set(gca,'YTickLabel',yy1./10)
text(0.05,0.7,'R_o','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
text(1.01,0.05,'10^{-2}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
text(-0.025,1.15,'10^{1}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
% ------------ 2/5 -------------------- % 
axes('position',sub_pos{2,5});
mmx         = 0;
for mc = 1:MC
   [f,x]    = ksdensity(mS1(2,:,mc));
   plot(x,f,'b','Linewidth',2)
   hold on
   if max(f) > mmx
       mmx  = max(f);
   end
end
axis tight
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
yy1         = get(gca,'Ytick');
x1          = [par(2) par(2)];
y1          = [0 mmx];
hold on
plot(x1,y1,'r','LineWidth',2)
xx          = get(gca,'Xtick');
set(gca,'Ytick',yy1);
set(gca,'XTickLabel',xx.*1000)
set(gca,'YTickLabel',yy1./1000)
text(0.05,0.7,'R_i','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
text(1.01,0.05,'10^{-3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
text(-0.025,1.15,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
% ------------ 1/4 -------------------- % 
axes('position',sub_pos{1,4});
mmx         = 0;
for mc = 1:MC
   [f,x]    = ksdensity(mS1(3,:,mc));
   plot(x,f,'b','Linewidth',2)
   hold on
   if max(f) > mmx
       mmx  = max(f);
   end
end
axis tight
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
yy1         = get(gca,'Ytick');
x1          = [par(3) par(3)];
y1          = [0 mmx];
hold on
plot(x1,y1,'r','LineWidth',2)
xx          = get(gca,'Xtick');
set(gca,'Ytick',yy1);
set(gca,'XTickLabel',xx.*1000)
set(gca,'YTickLabel',yy1./1000)
text(0.05,0.7,'R_m','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
text(1.01,0.05,'10^{-3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
text(-0.025,1.15,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
% ------------ 2/4 -------------------- % 
axes('position',sub_pos{2,4});
mmx         = 0;
for mc = 1:MC
   [f,x]    = ksdensity(mS1(4,:,mc));
   plot(x,f,'b','Linewidth',2)
   hold on
   if max(f) > mmx
       mmx  = max(f);
   end
end
axis tight
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
yy1         = get(gca,'Ytick');
x1          = [par(4) par(4)];
y1          = [0 mmx];
hold on
plot(x1,y1,'r','LineWidth',2)
xx          = get(gca,'Xtick');
set(gca,'XTickLabel',xx.*1000)
set(gca,'Ytick',yy1);
set(gca,'YTickLabel',yy1./1000)
text(0.05,0.7,'R_z','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
text(1.01,0.05,'10^{-3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
text(-0.025,1.15,'10^{3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
% ------------ 1/3 -------------------- % 
axes('position',sub_pos{1,3});
mmx         = 0;
for mc = 1:MC
   [f,x]    = ksdensity(mS1(9,:,mc));
   plot(x,f,'b','Linewidth',2)
   hold on
   if max(f) > mmx
       mmx  = max(f);
   end
end
axis tight
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
yy1         = get(gca,'Ytick');
x1          = [par(9) par(9)];
y1          = [0 mmx];
hold on
plot(x1,y1,'r','LineWidth',2)
xx          = get(gca,'Xtick');
set(gca,'XTickLabel',xx.*1000)
set(gca,'Ytick',yy1);
set(gca,'YTickLabel',yy1./100)
text(0.05,0.7,'R_s_o','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
text(1.01,0.05,'10^{-3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
text(-0.025,1.15,'10^{2}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
% ------------ 2/3 -------------------- % 
axes('position',sub_pos{2,3});
mmx         = 0;
for mc = 1:MC
   [f,x]    = ksdensity(mS1(5,:,mc));
   plot(x,f,'b','Linewidth',2)
   hold on
   if max(f) > mmx
       mmx  = max(f);
   end
end
axis tight
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
yy1         = get(gca,'Ytick');
x1          = [par(5) par(5)];
y1          = [0 mmx];
hold on
plot(x1,y1,'r','LineWidth',2)
xx          = get(gca,'Xtick');
set(gca,'XTickLabel',xx.*100)
set(gca,'Ytick',yy1);
set(gca,'YTickLabel',yy1./10)
text(0.05,0.7,'C_w','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
text(1.01,0.05,'10^{6}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
text(-0.025,1.15,'10^{1}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
% ------------ 1/2 -------------------- % 
axes('position',sub_pos{1,2});
mmx         = 0;
for mc = 1:MC
   [f,x]    = ksdensity(mS1(6,:,mc));
   plot(x,f,'b','Linewidth',2)
   hold on
   if max(f) > mmx
       mmx  = max(f);
   end
end
axis tight
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
yy1         = get(gca,'Ytick');
x1          = [par(6) par(6)];
y1          = [0 mmx];
hold on
plot(x1,y1,'r','LineWidth',2)
xx          = get(gca,'Xtick');
set(gca,'XTickLabel',xx.*1000)
set(gca,'Ytick',yy1);
set(gca,'YTickLabel',yy1./100)
text(0.05,0.7,'C_i','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
text(1.01,0.05,'10^{5}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
text(-0.025,1.15,'10^{2}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
% ------------ 2/2 -------------------- % 
axes('position',sub_pos{2,2});
mmx         = 0;
for mc = 1:MC
   [f,x]    = ksdensity(mS1(7,:,mc));
   plot(x,f,'b','Linewidth',2)
   hold on
   if max(f) > mmx
       mmx  = max(f);
   end
end
axis tight
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
yy1         = get(gca,'Ytick');
x1          = [par(7) par(7)];
y1          = [0 mmx];
hold on
plot(x1,y1,'r','LineWidth',2)
xx          = get(gca,'Xtick');
set(gca,'XTickLabel',xx.*10)
set(gca,'Ytick',yy1);
set(gca,'YTickLabel',yy1./10)
text(0.05,0.7,'C_m','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
text(1.01,0.05,'10^{7}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
text(-0.025,1.15,'10^{1}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
% ------------ 1/1 -------------------- % 
axes('position',sub_pos{1,1});
mmx         = 0;
for mc = 1:MC
   [f,x]    = ksdensity(mS1(8,:,mc));
   plot(x,f,'b','Linewidth',2)
   hold on
   if max(f) > mmx
       mmx  = max(f);
   end
end
axis tight
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
yy1         = get(gca,'Ytick');
x1          = [par(8) par(8)];
y1          = [0 mmx];
hold on
plot(x1,y1,'r','LineWidth',2)
xx          = get(gca,'Xtick');
set(gca,'XTickLabel',xx.*10)
set(gca,'Ytick',yy1);
set(gca,'YTickLabel',yy1./10)
text(0.05,0.7,'a_I','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
text(1.01,0.05,'10^{-1}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
text(-0.025,1.15,'10^{1}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');

%% ------------------------------------------------------------------------
% Plot Figure 5.4
% -------------------------------------------------------------------------
aa          = 0.1;
plotheight  = 15;
plotwidth   = 15;
subplotsx   = 2;
subplotsy   = 3;
leftedge    = 1;
rightedge   = 1;
topedge     = 1;
bottomedge  = 1;
spacex      = 1.1;
spacey      = 1.3;
fontsize    = 16;
sub_pos     = subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);
f=figure('visible','on','units','normalized','outerposition',[0 0 0.5 1]);
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'color','w');
% -------------- 2/3 --------------------- % 
axes('position',sub_pos{2,3});
for i = 1:MC
    z       = zeros(size(xS1(1,:,i)));
    S       = surface([tdays;tdays],[xS1(1,:,i);xS1(1,:,i)],[z;z],...
                'facecol','no',...
                'edgecol','interp',...
                'linew',2,...
                'edgealpha',aa,...
                'edgecolor','b');
    hold on
end
axis tight
plot(tdays,xsim(1,:),'r','Linewidth',2);
ylabel('°C','FontName','Cambria','FontWeight','demi','FontSize',fontsize)
xlabel('days','FontName','Cambria','FontWeight','demi','FontSize',fontsize)
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
text(0.03,0.95,'\theta_w','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
% -------------- 2/2 --------------------- % 
axes('position',sub_pos{2,2});
for i = 1:MC
    z       = zeros(size(xS1(2,:,i)));
    S       = surface([tdays;tdays],[xS1(2,:,i);xS1(2,:,i)],[z;z],...
                'facecol','no',...
                'edgecol','interp',...
                'linew',2,...
                'edgealpha',aa,...
                'edgecolor','b');
    hold on
end
axis tight
plot(tdays,xsim(2,:),'r','Linewidth',2);
ylabel('°C','FontName','Cambria','FontWeight','demi','FontSize',fontsize)
xlabel('days','FontName','Cambria','FontWeight','demi','FontSize',fontsize)
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
text(0.03,0.95,'\theta_i','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
% -------------- 2/1 --------------------- % 
axes('position',sub_pos{2,1});
for i = 1:MC
    z       = zeros(size(xS1(3,:,i)));
    S       = surface([tdays;tdays],[xS1(3,:,i);xS1(3,:,i)],[z;z],...
                'facecol','no',...
                'edgecol','interp',...
                'linew',2,...
                'edgealpha',aa,...
                'edgecolor','b');
   hold on
end
axis tight
plot(tdays,xsim(3,:),'r','Linewidth',2);
ylabel('°C','FontName','Cambria','FontWeight','demi','FontSize',fontsize)
xlabel('days','FontName','Cambria','FontWeight','demi','FontSize',fontsize)
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
text(0.03,0.95,'\theta_m','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
% -------------- 1/3 --------------------- % 
axes('position',sub_pos{1,3});
for i = 1:MC
    z       = zeros(size(xS2(1,:,i)));
    S       = surface([tdays;tdays],[xS2(1,:,i);xS2(1,:,i)],[z;z],...
                'facecol','no',...
                'edgecol','interp',...
                'linew',2,...
                'edgealpha',aa,...
                'edgecolor','b');
   hold on
end
axis tight
plot(tdays,xsim(1,:),'r','Linewidth',2);
ylabel('°C','FontName','Cambria','FontWeight','demi','FontSize',fontsize)
xlabel('days','FontName','Cambria','FontWeight','demi','FontSize',fontsize)
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
text(0.03,0.95,'\theta_w','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
% -------------- 1/2 --------------------- % 
axes('position',sub_pos{1,2});
for i = 1:MC
    z       = zeros(size(xS2(2,:,i)));
    S       = surface([tdays;tdays],[xS2(2,:,i);xS2(2,:,i)],[z;z],...
                'facecol','no',...
                'edgecol','interp',...
                'linew',2,...
                'edgealpha',aa,...
                'edgecolor','b');
    hold on
end
axis tight
plot(tdays,xsim(2,:),'r','Linewidth',2);
ylabel('°C','FontName','Cambria','FontWeight','demi','FontSize',fontsize)
xlabel('days','FontName','Cambria','FontWeight','demi','FontSize',fontsize)
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
text(0.03,0.95,'\theta_i','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
% -------------- 1/1 --------------------- % 
axes('position',sub_pos{1,1});
for i = 1:MC
    z       = zeros(size(xS2(3,:,i)));
    S       = surface([tdays;tdays],[xS2(3,:,i);xS2(3,:,i)],[z;z],...
                'facecol','no',...
                'edgecol','interp',...
                'linew',2,...
                'edgealpha',aa,...
                'edgecolor','b');
   hold on
end
axis tight
plot(tdays,xsim(3,:),'r','Linewidth',2);
ylabel('°C','FontName','Cambria','FontWeight','demi','FontSize',fontsize)
xlabel('days','FontName','Cambria','FontWeight','demi','FontSize',fontsize)
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
text(0.03,0.95,'\theta_m','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');

%% ------------------------------------------------------------------------
% Plot FIgure 5.6
% -------------------------------------------------------------------------
plotheight  = 15;
plotwidth   = 15;
subplotsx   = 2;
subplotsy   = 5;
leftedge    = 1;
rightedge   = 1;
topedge     = 1;
bottomedge  = 1;
spacex      = 1.2;
spacey      = 1.1;
fontsize    = 16;
sub_pos     = subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);
f           = figure('visible','on','units','normalized','outerposition',[0 0 0.5 1]);
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
set(gcf,'color','w');
% ------------ 1/5 -------------------- % 
axes('position',sub_pos{1,5});
mmx         = 0;
for mc = 1:MC
   [f,x]    = ksdensity(mS2(1,:,mc));
   plot(x,f,'b','Linewidth',2)
   hold on
   if max(f) > mmx
       mmx  = max(f);
   end
end
axis tight
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
yy1         = get(gca,'Ytick');
x1          = [par(1) par(1)];
y1          = [0 mmx];
hold on
plot(x1,y1,'r','LineWidth',2)
xx          = get(gca,'Xtick');
set(gca,'Ytick',yy1);
set(gca,'XTickLabel',xx.*100)
set(gca,'YTickLabel',yy1./1e4)
text(0.9,0.7,'R_o','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
text(1.01,0.05,'10^{-2}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
text(-0.025,1.15,'10^{4}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
% ------------ 2/5 -------------------- % 
axes('position',sub_pos{2,5});
mmx         = 0;
for mc = 1:MC
   [f,x]    = ksdensity(mS2(2,:,mc));
   plot(x,f,'b','Linewidth',2)
   hold on
   if max(f) > mmx
       mmx  = max(f);
   end
end
axis tight
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
yy1         = get(gca,'Ytick');
x1          = [par(2) par(2)];
y1          = [0 mmx];
hold on
plot(x1,y1,'r','LineWidth',2)
xx          = get(gca,'Xtick');
set(gca,'Ytick',yy1);
set(gca,'XTickLabel',xx.*1000)
set(gca,'YTickLabel',yy1./1e6)
text(0.9,0.7,'R_i','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
text(1.01,0.05,'10^{-3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
text(-0.025,1.15,'10^{6}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
% ------------ 1/4 -------------------- % 
axes('position',sub_pos{1,4});
mmx         = 0;
for mc = 1:MC
   [f,x]    = ksdensity(mS2(3,:,mc));
   plot(x,f,'b','Linewidth',2)
   hold on
   if max(f) > mmx
       mmx  = max(f);
   end
end
axis tight
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
yy1         = get(gca,'Ytick');
x1          = [par(3) par(3)];
y1          = [0 mmx];
hold on
plot(x1,y1,'r','LineWidth',2)
xx          = get(gca,'Xtick');
set(gca,'Ytick',yy1);
set(gca,'XTickLabel',xx.*1000)
set(gca,'YTickLabel',yy1./1e5)
text(0.9,0.7,'R_m','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
text(1.01,0.05,'10^{-3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
text(-0.025,1.15,'10^{5}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
% ------------ 2/4 -------------------- % 
axes('position',sub_pos{2,4});
mmx         = 0;
for mc = 1:MC
   [f,x]    = ksdensity(mS2(4,:,mc));
   plot(x,f,'b','Linewidth',2)
   hold on
   if max(f) > mmx
       mmx  = max(f);
   end
end
axis tight
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
yy1         = get(gca,'Ytick');
x1          = [par(4) par(4)];
y1          = [0 mmx];
hold on
plot(x1,y1,'r','LineWidth',2)
xx          = get(gca,'Xtick');
set(gca,'XTickLabel',xx.*1000)
set(gca,'Ytick',yy1);
set(gca,'YTickLabel',yy1./1e5)
text(0.9,0.7,'R_z','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
text(1.01,0.05,'10^{-3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
text(-0.025,1.15,'10^{5}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
% ------------ 1/3 -------------------- % 
axes('position',sub_pos{1,3});
mmx         = 0;
for mc = 1:MC
   [f,x]    = ksdensity(mS2(9,:,mc));
   plot(x,f,'b','Linewidth',2)
   hold on
   if max(f) > mmx
       mmx  = max(f);
   end
end
axis tight
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
yy1         = get(gca,'Ytick');
x1          = [par(9) par(9)];
y1          = [0 mmx];
hold on
plot(x1,y1,'r','LineWidth',2)
xx          = get(gca,'Xtick');
set(gca,'XTickLabel',xx.*1000)
set(gca,'Ytick',yy1);
set(gca,'YTickLabel',yy1./1e5)
text(0.9,0.7,'R_s_o','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
text(1.01,0.05,'10^{-3}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
text(-0.025,1.15,'10^{5}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
% ------------ 2/3 -------------------- % 
axes('position',sub_pos{2,3});
mmx         = 0;
for mc = 1:MC
   [f,x]    = ksdensity(mS2(5,:,mc));
   plot(x,f,'b','Linewidth',2)
   hold on
   if max(f) > mmx
       mmx  = max(f);
   end
end
axis tight
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
yy1         = get(gca,'Ytick');
x1          = [par(5) par(5)];
y1          = [0 mmx];
hold on
plot(x1,y1,'r','LineWidth',2)
xx          = get(gca,'Xtick');
set(gca,'XTickLabel',xx.*100)
set(gca,'Ytick',yy1);
set(gca,'YTickLabel',yy1./1e4)
text(0.9,0.7,'C_w','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
text(1.01,0.05,'10^{6}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
text(-0.025,1.15,'10^{4}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
% ------------ 1/2 -------------------- % 
axes('position',sub_pos{1,2});
mmx         = 0;
for mc = 1:MC
   [f,x]    = ksdensity(mS2(6,:,mc));
   plot(x,f,'b','Linewidth',2)
   hold on
   if max(f) > mmx
       mmx  = max(f);
   end
end
axis tight
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
yy1         = get(gca,'Ytick');
x1          = [par(6) par(6)];
y1          = [0 mmx];
hold on
plot(x1,y1,'r','LineWidth',2)
xx          = get(gca,'Xtick');
set(gca,'XTickLabel',xx.*1000)
set(gca,'Ytick',yy1);
set(gca,'YTickLabel',yy1./1e5)
text(0.9,0.7,'C_i','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
text(1.01,0.05,'10^{5}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
text(-0.025,1.15,'10^{5}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
% ------------ 2/2 -------------------- % 
axes('position',sub_pos{2,2});
mmx         = 0;
for mc = 1:MC
   [f,x]    = ksdensity(mS2(7,:,mc));
   plot(x,f,'b','Linewidth',2)
   hold on
   if max(f) > mmx
       mmx  = max(f);
   end
end
axis tight
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
yy1         = get(gca,'Ytick');
x1          = [par(7) par(7)];
y1          = [0 mmx];
hold on
plot(x1,y1,'r','LineWidth',2)
xx          = get(gca,'Xtick');
set(gca,'XTickLabel',xx.*10)
set(gca,'Ytick',yy1);
set(gca,'YTickLabel',yy1./1e4)
text(0.9,0.7,'C_m','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
text(1.01,0.05,'10^{7}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
text(-0.025,1.15,'10^{4}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
% ------------ 1/1 -------------------- % 
axes('position',sub_pos{1,1});
mmx         = 0;
for mc = 1:MC
   [f,x]    = ksdensity(mS2(8,:,mc));
   plot(x,f,'b','Linewidth',2)
   hold on
   if max(f) > mmx
       mmx  = max(f);
   end
end
axis tight
set(gca,'Box','off')
set(gca,'FontWeight','demi')
set(gca,'FontSize',fontsize)
set(gca,'FontName','Cambria')
set(gca,'Linewidth',2)
yy1         = get(gca,'Ytick');
x1          = [par(8) par(8)];
y1          = [0 mmx];
hold on
plot(x1,y1,'r','LineWidth',2)
xx          = get(gca,'Xtick');
set(gca,'XTickLabel',xx.*10)
set(gca,'Ytick',yy1);
set(gca,'YTickLabel',yy1./1e4)
text(0.9,0.7,'a_I','FontSize',20,'FontWeight','demi','Units','Normalized','FontAngle','italic','FontName','Cambria');
text(1.01,0.05,'10^{-1}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');
text(-0.025,1.15,'10^{4}','FontSize',fontsize,'FontWeight','demi','Units','Normalized','FontName','Cambria');