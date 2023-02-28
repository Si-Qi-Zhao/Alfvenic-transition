clear; close all;

date   = '20031203200000';
trange = '300';
numk   = 800;
B0     = 106;

angle_threshod = 10; 
[P_perp_10,P_para_10,kperp_10,kpara_10,kpara_new_10,Ekperp_10,Ekperp_pre2_10,XiB_new_10,vk_new_10,XiB_mat_sc_10,vk2_sc_10] = nonlinearity_V(date,trange,numk,angle_threshod,B0);

angle_threshod = 15; 
[P_perp_15,P_para_15,kperp_15,kpara_15,kpara_new_15,Ekperp_15,Ekperp_pre2_15,XiB_new_15,vk_new_15,XiB_mat_sc_15,vk2_sc_15] = nonlinearity_V(date,trange,numk,angle_threshod,B0);

angle_threshod = 20; 
[P_perp_20,P_para_20,kperp_20,kpara_20,kpara_new_20,Ekperp_20,Ekperp_pre2_20,XiB_new_20,vk_new_20,XiB_mat_sc_20,vk2_sc_20] = nonlinearity_V(date,trange,numk,angle_threshod,B0);

angle_threshod = 25; 
[P_perp_25,P_para_25,kperp_25,kpara_25,kpara_new_25,Ekperp_25,Ekperp_pre2_25,XiB_new_25,vk_new_25,XiB_mat_sc_25,vk2_sc_25] = nonlinearity_V(date,trange,numk,angle_threshod,B0);

angle_threshod = 30; 
[P_perp_30,P_para_30,kperp_30,kpara_30,kpara_new_30,Ekperp_30,Ekperp_pre2_30,XiB_new_30,vk_new_30,XiB_mat_sc_30,vk2_sc_30] = nonlinearity_V(date,trange,numk,angle_threshod,B0);

%%
figure ; 
loglog(kperp_10,P_perp_10,'o-','linewidth',2,'Markersize',3); hold on;
loglog(kperp_15,P_perp_15,'o-','linewidth',2,'Markersize',3); 
loglog(kperp_20,P_perp_20,'o-','linewidth',2,'Markersize',3); 
loglog(kperp_25,P_perp_25,'o-','linewidth',2,'Markersize',3); 
loglog(kperp_30,P_perp_30,'o-','linewidth',2,'Markersize',3); 
legend('\theta_{kk''}=10^\circ','\theta_{kk''}=15^\circ',...
       '\theta_{kk''}=20^\circ','\theta_{kk''}=25^\circ','\theta_{kk''}=30^\circ'); legend('boxoff','Location','best');
xlabel('k_{\perp} (1/km)');ylabel('P(k_{\perp}) (km/s)^2/Hz');
set(gca,'Fontsize',25,'tickDir','both','linewidth',2);
xlim([5e-5,1e-3]);set(gcf,'Position',[10 10 740 550]);

figure;
loglog(kperp_10,Ekperp_pre2_10,'d-','linewidth',2,'Markersize',3); hold on;
loglog(kperp_15,Ekperp_pre2_15,'d-','linewidth',2,'Markersize',3); 
loglog(kperp_20,Ekperp_pre2_20,'d-','linewidth',2,'Markersize',3); 
loglog(kperp_25,Ekperp_pre2_25,'d-','linewidth',2,'Markersize',3); 
loglog(kperp_30,Ekperp_pre2_30,'d-','linewidth',2,'Markersize',3); 
yyy2 = 0.18*kperp_10.^(-1/3);
loglog(kperp_10,yyy2,':k','linewidth',3);
legend('\theta_{kk''}=10^\circ','\theta_{kk''}=15^\circ',...
       '\theta_{kk''}=20^\circ','\theta_{kk''}=25^\circ','\theta_{kk''}=30^\circ','\propto k_{\perp}^{5/3-2}','NumColumns',2,'Location','best'); legend('boxoff');
set(gca,'Fontsize',22,'tickDir','both','linewidth',2);
xlim([5e-5,1e-3]); ylim([1.5,5]);set(gcf,'Position',[10 10 840 350]);

figure ; 
loglog(kperp_10,XiB_new_10,'o-','linewidth',2,'Markersize',3); hold on;
loglog(kperp_15,XiB_new_15,'o-','linewidth',2,'Markersize',3);
loglog(kperp_20,XiB_new_20,'o-','linewidth',2,'Markersize',3);
loglog(kperp_25,XiB_new_25,'o-','linewidth',2,'Markersize',3);
loglog(kperp_30,XiB_new_30,'o-','linewidth',2,'Markersize',3);
yline(1,':','linewidth',3);
legend('\theta_{kk''}=10^\circ','\theta_{kk''}=15^\circ',...
       '\theta_{kk''}=20^\circ','\theta_{kk''}=25^\circ','\theta_{kk''}=30^\circ','y=1','NumColumns',2,'Location','best'); legend('boxoff');
ylabel('\chi_B'); xlabel('k_{\perp}');
set(gca,'Fontsize',20,'tickDir','both','linewidth',2,'xscale','log','yscale','lin')
xlim([5e-5,1e-3 ]);set(gcf,'Position',[10 10 740 450]);

figure ; 
loglog(kperp_10,kpara_new_10,'+-','linewidth',2,'Markersize',3); hold on;
loglog(kperp_15,kpara_new_15,'+-','linewidth',2,'Markersize',3); 
loglog(kperp_20,kpara_new_20,'+-','linewidth',2,'Markersize',3); 
loglog(kperp_25,kpara_new_25,'+-','linewidth',2,'Markersize',3); 
loglog(kperp_30,kpara_new_30,'+-','linewidth',2,'Markersize',3); 
yyy = 0.023*kperp_10.^(2/3);
yyy2 = 0.3*kperp_10;
loglog(kperp_10,yyy,':k','linewidth',3);
loglog(kperp_10,yyy2,':r','linewidth',3);
legend('\theta_{kk''}=10^\circ','\theta_{kk''}=15^\circ',...
       '\theta_{kk''}=20^\circ','\theta_{kk''}=25^\circ','\theta_{kk''}=30^\circ','\propto k_{\perp}^{2/3}','\propto k_{\perp}','NumColumns',3,'Location','southeast'); legend('boxoff');
ylabel('k_{||}'); xlabel('k_{\perp}');
set(gca,'Fontsize',20,'tickDir','both','linewidth',2);
xlim([5e-5,1e-3]);ylim([2.5e-5,3e-4]);set(gcf,'Position',[10 10 740 450]);


%% %%%%%%%%%%%%%%%%%%%%%%%%
figure;
h1 = subplot(3,1,1);
h1.Position = [0.1700    0.673    0.70    0.25];
loglog(kperp_10,Ekperp_pre2_10,'d-','linewidth',2,'Markersize',3); hold on;
loglog(kperp_15,Ekperp_pre2_15,'d-','linewidth',2,'Markersize',3); 
loglog(kperp_20,Ekperp_pre2_20,'d-','linewidth',2,'Markersize',3); 
loglog(kperp_25,Ekperp_pre2_25,'d-','linewidth',2,'Markersize',3); 
loglog(kperp_30,Ekperp_pre2_30,'d-','linewidth',2,'Markersize',3); 
yyy2 = 0.08*kperp_10.^(-1/3);
loglog(kperp_10,yyy2,'--k','linewidth',3);
xline(2e-4,':','linewidth',1.5);
xline(3e-4,':','linewidth',1.5);
xline(7e-4,':','linewidth',1.5);
legend('','','','','','\propto k_{\perp}^{5/3-2}','','','NumColumns',3,'Location','best'); legend('boxoff');
set(gca,'Fontsize',23,'tickDir','both','linewidth',2,'XTickLabel','');
ylabel({'k_{\perp}^{5/3}E_V_{A}(k_{\perp})','(km s^{-1})^2 km^{-2/3}'});
xlim([5e-5,1.0e-3]); ylim([0.5,3])

h2 = subplot(3,1,2);
h2.Position = [0.1700    0.40    0.70     0.25];
loglog(kperp_10,kpara_new_10,'+-','linewidth',2,'Markersize',3); hold on;
loglog(kperp_15,kpara_new_15,'+-','linewidth',2,'Markersize',3); 
loglog(kperp_20,kpara_new_20,'+-','linewidth',2,'Markersize',3); 
loglog(kperp_25,kpara_new_25,'+-','linewidth',2,'Markersize',3); 
loglog(kperp_30,kpara_new_30,'+-','linewidth',2,'Markersize',3); 
xx1 = 5e-6:1e-5:1e-3;
yyy = 0.037*xx1.^(2/3);
loglog(xx1,yyy,':m','linewidth',3);
xline(2e-4,':','linewidth',1.5);
xline(3e-4,':','linewidth',1.5);
xline(7e-4,':','linewidth',1.5);
yline(7e-5,':','linewidth',2);
yline(1e-4,':','linewidth',2);
legend('','','','','','\propto k_{\perp}^{2/3}','','','','','','NumColumns',3,'Location','northwest'); legend('boxoff');
ylabel({'k_{||}','(km s^{-1})'}); 
set(gca,'Fontsize',23,'tickDir','both','linewidth',2,'XTickLabel','');
xlim([5e-5,1e-3]);ylim([2.5e-5,5e-4]);

h3 = subplot(3,1,3);
h3.Position = [0.1700    0.125     0.70     0.25];
h3 = pcolor(kperp_30,kpara_30,XiB_mat_sc_30); hold on;
shading(gca,'flat'); colormap('jet'); c=colorbar;caxis(gca,[0 1]);
c.Location = 'eastoutside';
c.Label.String = '\chi_V_{A}';
xlim([5e-5,1.0e-3]); ylim([2.5e-5,5e-4]);
xline(2e-4,':','linewidth',1.5);
xline(3e-4,':','linewidth',1.5);
xline(7e-4,':','linewidth',1.5);
xlabel('k_{\perp} (km^{-1})');ylabel({'k_{||}','(km s^{-1})'}); 
set(gca,'Fontsize',23,'tickDir','both','linewidth',2,'xscale','log','yscale','log');
set(gcf,'Position',[10 10 800 820]);

xx1 = 5e-6:1e-5:1e-3;
yyy = 0.26*xx1.^(1);
yyy2 = 0.017*xx1.^(2/3);
loglog(xx1,yyy,':w','linewidth',3);
loglog(xx1,yyy2,':m','linewidth',3);
text(8e-5,2e-4,'k_{||}\propto k_{\perp}^{2/3}','FontSize',25,'Color','m');
text(8e-5,1e-4,'k_{||}\propto k_{\perp}','FontSize',25,'Color','w');

% loglog(kperp_10,XiB_new_10,'o-','linewidth',2,'Markersize',3); hold on;
% loglog(kperp_15,XiB_new_15,'o-','linewidth',2,'Markersize',3);
% loglog(kperp_20,XiB_new_20,'o-','linewidth',2,'Markersize',3);
% loglog(kperp_25,XiB_new_25,'o-','linewidth',2,'Markersize',3);
% loglog(kperp_30,XiB_new_30,'o-','linewidth',2,'Markersize',3);
% yline(1,':k','linewidth',2);
% yline(0.5,':k','linewidth',2);
% xline(1.6e-4,':','linewidth',1.5);
% xline(3e-4,':','linewidth',1.5);
% xline(6e-4,':','linewidth',1.5);
% legend('\theta_{kk''}=10^\circ','\theta_{kk''}=15^\circ',...
%        '\theta_{kk''}=20^\circ','\theta_{kk''}=25^\circ','\theta_{kk''}=30^\circ','','','','','NumColumns',3,'Location','best'); legend('boxoff');
% ylabel('\chi_V_{A}'); xlabel('k_{\perp} (1/km)');
% set(gca,'Fontsize',18,'tickDir','both','linewidth',2,'xscale','log','yscale','lin')
% xlim([5e-5,1.0e-3]);  ylim([0,2])
% set(gcf,'Position',[10 10 800 820]);



%
figure;
for i=26:10:536%length(kpara_20)
    kpara_line = log10(kpara_20(i)).*ones(1,length(kperp_20));
    scatter(kperp_20,XiB_mat_sc_20(i,:),5,kpara_line,'filled'); hold on;
end
yline(1,':','linewidth',1.5);
ylim([1e-3,5e1]);colormap jet;
xlabel('k_{\perp} (km^{-1})');ylabel('\chi_V');
set(gca,'Fontsize',18,'tickDir','both','linewidth',2,'xscale','log','yscale','log');
ccc = colorbar;
ccc.Location = 'south';
ccc.Position = [0.5433    0.2308    0.25    0.0724];


% 
% 
% figure;
% h1 = subplot(4,1,1);
% h1.Position = [0.1700    0.7373    0.70    0.2];
% loglog(kperp_10,Ekperp_10,'d-','linewidth',2,'MarkerSize',2);hold on;
% loglog(kperp_15,Ekperp_15,'d-','linewidth',2,'MarkerSize',2);
% loglog(kperp_20,Ekperp_20,'d-','linewidth',2,'MarkerSize',2);
% loglog(kperp_25,Ekperp_25,'d-','linewidth',2,'MarkerSize',2);
% loglog(kperp_30,Ekperp_30,'d-','linewidth',2,'MarkerSize',2);
% % idx_p1 = 71:135;  idx_p2 = 136:355;
% % [p1,S1] = polyfit(log(kperp(idx_p1)),log(Ekperp(idx_p1)),1);
% % yy1 = exp(polyval(p1,log(kperp(idx_p1)),S1));
% % loglog(kperp(idx_p1),yy1,'r--','linewidth',3);  
% % [p2,S2] = polyfit(log(kperp(idx_p2)),log(Ekperp(idx_p2)),1);
% % yy2 = exp(polyval(p2,log(kperp(idx_p2)),S2));
% % loglog(kperp(idx_p2),yy2,'--','linewidth',3,'color','m'); 
% xline(2e-4,':','linewidth',1.5);
% xline(3e-4,':','linewidth',1.5);
% xline(6e-4,':','linewidth',1.5);
% ylabel({'E_V_{A}(k_{\perp})','(km s^{-1})^2 km'});
% set(gca,'Fontsize',18,'linewidth',2,'XTickLabel','');
% xlim([5e-5,1.0e-3]); ylim([5e4,2e7]);
% 
% h2 = subplot(4,1,2);
% h2.Position = [0.1700    0.5232    0.70     0.2];
% loglog(kperp_10,Ekperp_pre2_10,'d-','linewidth',2,'Markersize',3); hold on;
% loglog(kperp_15,Ekperp_pre2_15,'d-','linewidth',2,'Markersize',3); 
% loglog(kperp_20,Ekperp_pre2_20,'d-','linewidth',2,'Markersize',3); 
% loglog(kperp_25,Ekperp_pre2_25,'d-','linewidth',2,'Markersize',3); 
% loglog(kperp_30,Ekperp_pre2_30,'d-','linewidth',2,'Markersize',3); 
% yyy2 = 0.08*kperp_10.^(-1/3);
% loglog(kperp_10,yyy2,'--k','linewidth',3);
% xline(2e-4,':','linewidth',1.5);
% xline(3e-4,':','linewidth',1.5);
% xline(6e-4,':','linewidth',1.5);
% legend('','','','','','\propto k_{\perp}^{5/3-2}','','','NumColumns',3,'Location','best'); legend('boxoff');
% set(gca,'Fontsize',18,'tickDir','both','linewidth',2,'XTickLabel','');
% ylabel({'k_{\perp}^{5/3}E_V_{A}(k_{\perp})','(km s^{-1})^2 km^{-2/3}'});
% xlim([5e-5,1.0e-3]); ylim([0.5,2.2])
% 
% h3 = subplot(4,1,3);
% h3.Position = [0.1700    0.3091     0.70     0.2];
% loglog(kperp_10,kpara_new_10,'+-','linewidth',2,'Markersize',3); hold on;
% loglog(kperp_15,kpara_new_15,'+-','linewidth',2,'Markersize',3); 
% loglog(kperp_20,kpara_new_20,'+-','linewidth',2,'Markersize',3); 
% loglog(kperp_25,kpara_new_25,'+-','linewidth',2,'Markersize',3); 
% loglog(kperp_30,kpara_new_30,'+-','linewidth',2,'Markersize',3); 
% xx1 = 3e-4:1e-5:6e-4;
% yyy = 0.037*xx1.^(2/3);
% loglog(xx1,yyy,':k','linewidth',3);
% xline(2e-4,':','linewidth',1.5);
% xline(3e-4,':','linewidth',1.5);
% xline(6e-4,':','linewidth',1.5);
% yline(7e-5,':','linewidth',2);
% yline(1e-4,':','linewidth',2);
% legend('','','','','','\propto k_{\perp}^{2/3}','','','','','','NumColumns',3,'Location','northwest'); legend('boxoff');
% ylabel({'k_{||}','(km s^{-1})'}); 
% set(gca,'Fontsize',18,'tickDir','both','linewidth',2,'XTickLabel','');
% xlim([5e-5,1e-3]);ylim([2.5e-5,5e-4]);
% 
% h4 = subplot(4,1,4);
% h4.Position = [0.1700    0.0950     0.70     0.2];
% h4 = pcolor(kperp_30,kpara_30,XiB_mat_sc_30); hold on;
% shading(gca,'flat'); colormap('jet'); c=colorbar;caxis(gca,[0 1]);
% c.Location = 'eastoutside';
% c.Label.String = '\chi_V_{A}';
% xlim([5e-5,1.0e-3]); ylim([2.5e-5,5e-4]);
% xline(2e-4,':','linewidth',1.5);
% xline(3e-4,':','linewidth',1.5);
% xline(6e-4,':','linewidth',1.5);
% xlabel('k_{\perp} (km^{-1})');ylabel({'k_{||}','(km s^{-1})'}); 
% set(gca,'Fontsize',18,'tickDir','both','linewidth',2,'xscale','log','yscale','log');
% set(gcf,'Position',[10 10 800 820]);
% 
% xx1 = 3e-4:1e-5:6e-4;
% yyy = 0.26*xx1.^(1);
% yyy2 = 0.017*xx1.^(2/3);
% loglog(xx1,yyy,':w','linewidth',3);
% loglog(xx1,yyy2,':k','linewidth',3);