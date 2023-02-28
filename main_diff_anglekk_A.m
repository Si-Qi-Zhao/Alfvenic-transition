clear; close all;

date   = '20031203200000';%'20021216050000';%'20011219000000';%'20011216120000';%'20171126214800';
trange = '300';
numk   = 800;
B0     = 106;

angle_threshod = 10; 
[P_perp_10,P_para_10,kperp_10,kpara_10,kpara_new_10,Ekperp_10,Ekperp_pre2_10,XiB_new_10,bk_new_10,XiB_mat_sc_10,bk2_sc_10] = nonlinearity_A(date,trange,numk,angle_threshod,B0);

angle_threshod = 15; 
[P_perp_15,P_para_15,kperp_15,kpara_15,kpara_new_15,Ekperp_15,Ekperp_pre2_15,XiB_new_15,bk_new_15,XiB_mat_sc_15,bk2_sc_15] = nonlinearity_A(date,trange,numk,angle_threshod,B0);

angle_threshod = 20; 
[P_perp_20,P_para_20,kperp_20,kpara_20,kpara_new_20,Ekperp_20,Ekperp_pre2_20,XiB_new_20,bk_new_20,XiB_mat_sc_20,bk2_sc_20] = nonlinearity_A(date,trange,numk,angle_threshod,B0);

angle_threshod = 25; 
[P_perp_25,P_para_25,kperp_25,kpara_25,kpara_new_25,Ekperp_25,Ekperp_pre2_25,XiB_new_25,bk_new_25,XiB_mat_sc_25,bk2_sc_25] = nonlinearity_A(date,trange,numk,angle_threshod,B0);

angle_threshod = 30; 
[P_perp_30,P_para_30,kperp_30,kpara_30,kpara_new_30,Ekperp_30,Ekperp_pre2_30,XiB_new_30,bk_new_30,XiB_mat_sc_30,bk2_sc_30] = nonlinearity_A(date,trange,numk,angle_threshod,B0);

%%
figure ; 
loglog(kperp_10,P_perp_10,'o-','Linewidth',2,'Markersize',2,'Color','#0072BD'); hold on;
loglog(kperp_15,P_perp_15,'o-','Linewidth',2,'Markersize',2,'Color','#D95319'); 
loglog(kperp_20,P_perp_20,'o-','Linewidth',2,'Markersize',2,'Color','#EDB120'); 
loglog(kperp_25,P_perp_25,'o-','Linewidth',2,'Markersize',2,'Color','#7E2F8E'); 
loglog(kperp_30,P_perp_30,'o-','Linewidth',2,'Markersize',2,'Color','#77AC30'); 
legend('\theta_{kk''}=10^\circ','\theta_{kk''}=15^\circ',...
       '\theta_{kk''}=20^\circ','\theta_{kk''}=25^\circ','\theta_{kk''}=30^\circ','Location','best'); legend('boxoff');
xlabel({'k_{\perp}','(km^{-1})'}); ylabel({'P(k_{\perp})','(km s^{-1})^2/Hz'});
set(gca,'Fontsize',25,'tickDir','both','Linewidth',2);
xlim([5e-5,1e-3]);set(gcf,'Position',[10 10 640 550]);

figure ; 
loglog(kpara_10,P_para_10,':','Linewidth',4,'Markersize',2,'Color','#0072BD'); hold on;
loglog(kpara_15,P_para_15,':','Linewidth',4,'Markersize',2,'Color','#D95319'); 
loglog(kpara_20,P_para_20,':','Linewidth',4,'Markersize',2,'Color','#EDB120'); 
loglog(kpara_25,P_para_25,':','Linewidth',4,'Markersize',2,'Color','#7E2F8E'); 
loglog(kpara_30,P_para_30,':','Linewidth',4,'Markersize',2,'Color','#77AC30'); 
legend('\theta_{kk''}=10^\circ','\theta_{kk''}=15^\circ',...
       '\theta_{kk''}=20^\circ','\theta_{kk''}=25^\circ','\theta_{kk''}=30^\circ','Location','best'); legend('boxoff');
xlabel({'k_{||}','(km^{-1})'});ylabel({'P(k_{||})','(km s^{-1})^2/Hz'});
set(gca,'Fontsize',25,'tickDir','both','Linewidth',2);
xlim([8e-6,1e-3]);  ylim([1e3,1e10]);  set(gcf,'Position',[10 10 640 550]);

figure ; 
loglog(kperp_10,P_perp_10,'o-','Linewidth',2,'Markersize',2,'Color','#0072BD'); hold on;
loglog(kperp_15,P_perp_15,'o-','Linewidth',2,'Markersize',2,'Color','#D95319'); 
loglog(kperp_20,P_perp_20,'o-','Linewidth',2,'Markersize',2,'Color','#EDB120'); 
loglog(kperp_25,P_perp_25,'o-','Linewidth',2,'Markersize',2,'Color','#7E2F8E'); 
loglog(kperp_30,P_perp_30,'o-','Linewidth',2,'Markersize',2,'Color','#77AC30'); 
loglog(kpara_10,P_para_10,':','Linewidth',4,'Markersize',2,'Color','#0072BD'); hold on;
loglog(kpara_15,P_para_15,':','Linewidth',4,'Markersize',2,'Color','#D95319'); 
loglog(kpara_20,P_para_20,':','Linewidth',4,'Markersize',2,'Color','#EDB120'); 
loglog(kpara_25,P_para_25,':','Linewidth',4,'Markersize',2,'Color','#7E2F8E'); 
loglog(kpara_30,P_para_30,':','Linewidth',4,'Markersize',2,'Color','#77AC30'); 
% legend('\theta_{kk''}=10^\circ','\theta_{kk''}=15^\circ',...
%        '\theta_{kk''}=20^\circ','\theta_{kk''}=25^\circ','\theta_{kk''}=30^\circ'); legend('boxoff','Location','best');
xlabel('k (km^{-1})');   ylabel({'P(k_{\perp}) & P(k_{||})','(km s^{-1})^2'});
set(gca,'Fontsize',25,'tickDir','both','Linewidth',2);
xlim([6e-5,1e-3]);   ylim([1e-2,5e1]);
set(gcf,'Position',[10 10 640 500]);

%% normarlize bu psd at k ~ 6e-3 1/km
P_perp_10_norm = P_perp_10./P_perp_10(17);
P_perp_15_norm = P_perp_15./P_perp_15(20);
P_perp_20_norm = P_perp_20./P_perp_20(17);
P_perp_25_norm = P_perp_25./P_perp_25(18);
P_perp_30_norm = P_perp_30./P_perp_30(17);
P_para_10_norm = P_para_10./P_para_10(36);
P_para_15_norm = P_para_15./P_para_15(39);
P_para_20_norm = P_para_20./P_para_20(36);
P_para_25_norm = P_para_25./P_para_25(37);
P_para_30_norm = P_para_30./P_para_30(36);
figure ; 
loglog(kperp_10,P_perp_10_norm,'o-','Linewidth',2,'Markersize',2,'Color','#0072BD'); hold on;
loglog(kperp_15,P_perp_15_norm,'o-','Linewidth',2,'Markersize',2,'Color','#D95319'); 
loglog(kperp_20,P_perp_20_norm,'o-','Linewidth',2,'Markersize',2,'Color','#EDB120'); 
loglog(kperp_25,P_perp_25_norm,'o-','Linewidth',2,'Markersize',2,'Color','#7E2F8E'); 
loglog(kperp_30,P_perp_30_norm,'o-','Linewidth',2,'Markersize',2,'Color','#77AC30'); 
loglog(kpara_10,P_para_10_norm,':','Linewidth',4,'Markersize',2,'Color','#0072BD'); hold on;
loglog(kpara_15,P_para_15_norm,':','Linewidth',4,'Markersize',2,'Color','#D95319'); 
loglog(kpara_20,P_para_20_norm,':','Linewidth',4,'Markersize',2,'Color','#EDB120'); 
loglog(kpara_25,P_para_25_norm,':','Linewidth',4,'Markersize',2,'Color','#7E2F8E'); 
loglog(kpara_30,P_para_30_norm,':','Linewidth',4,'Markersize',2,'Color','#77AC30'); 
legend('\theta_{kk''}=10^\circ','\theta_{kk''}=15^\circ',...
       '\theta_{kk''}=20^\circ','\theta_{kk''}=25^\circ','\theta_{kk''}=30^\circ','Location','best'); legend('boxoff');
xlabel('k (km^{-1})');   ylabel('P(k_{\perp}) & P(k_{||})');
set(gca,'Fontsize',25,'tickDir','both','Linewidth',2);
xlim([6e-5,1e-3]);   ylim([5e-4,1e1]);
set(gcf,'Position',[10 10 640 500]);

figure;
loglog(kperp_10,Ekperp_pre2_10,'d-','Linewidth',2,'Markersize',2); hold on;
loglog(kperp_15,Ekperp_pre2_15,'d-','Linewidth',2,'Markersize',2); 
loglog(kperp_20,Ekperp_pre2_20,'d-','Linewidth',2,'Markersize',2); 
loglog(kperp_25,Ekperp_pre2_25,'d-','Linewidth',2,'Markersize',2); 
loglog(kperp_30,Ekperp_pre2_30,'d-','Linewidth',2,'Markersize',2); 
yyy2 = 0.18*kperp_10.^(-1/3);
loglog(kperp_10,yyy2,':k','Linewidth',3);
legend('\theta_{kk''}=10^\circ','\theta_{kk''}=15^\circ',...
       '\theta_{kk''}=20^\circ','\theta_{kk''}=25^\circ','\theta_{kk''}=30^\circ','\propto k_{\perp}^{5/3-2}','NumColumns',2,'Location','best'); legend('boxoff');
xlim([5e-5,1e-3]);   ylim([1.5,5]);   set(gcf,'Position',[10 10 840 350]);
set(gca,'Fontsize',22,'tickDir','both','Linewidth',2);

figure ; 
loglog(kperp_10,XiB_new_10,'o-','Linewidth',2,'Markersize',2); hold on;
loglog(kperp_15,XiB_new_15,'o-','Linewidth',2,'Markersize',2);
loglog(kperp_20,XiB_new_20,'o-','Linewidth',2,'Markersize',2);
loglog(kperp_25,XiB_new_25,'o-','Linewidth',2,'Markersize',2);
loglog(kperp_30,XiB_new_30,'o-','Linewidth',2,'Markersize',2);
yline(1,':','Linewidth',3);
legend('\theta_{kk''}=10^\circ','\theta_{kk''}=15^\circ',...
       '\theta_{kk''}=20^\circ','\theta_{kk''}=25^\circ','\theta_{kk''}=30^\circ','y=1','NumColumns',2,'Location','best'); legend('boxoff');
ylabel('\chi_B');    xlabel('k_{\perp}');
set(gca,'Fontsize',20,'tickDir','both','Linewidth',2,'xscale','log','yscale','lin')
xlim([5e-5,1e-3 ]);  set(gcf,'Position',[10 10 740 450]);

figure ; 
loglog(kperp_10,kpara_new_10,'+-','Linewidth',2,'Markersize',2);   hold on;
loglog(kperp_15,kpara_new_15,'+-','Linewidth',2,'Markersize',2); 
loglog(kperp_20,kpara_new_20,'+-','Linewidth',2,'Markersize',2); 
loglog(kperp_25,kpara_new_25,'+-','Linewidth',2,'Markersize',2); 
loglog(kperp_30,kpara_new_30,'+-','Linewidth',2,'Markersize',2); 
yyy = 0.023*kperp_10.^(2/3);
yyy2 = 0.3*kperp_10;
loglog(kperp_10,yyy, ':k','Linewidth',3);
loglog(kperp_10,yyy2,':r','Linewidth',3);
legend('\theta_{kk''}=10^\circ','\theta_{kk''}=15^\circ',...
       '\theta_{kk''}=20^\circ','\theta_{kk''}=25^\circ','\theta_{kk''}=30^\circ','\propto k_{\perp}^{2/3}','\propto k_{\perp}','NumColumns',3,'Location','southeast'); legend('boxoff');
ylabel('k_{||}'); xlabel('k_{\perp}');
set(gca,'Fontsize',20,'tickDir','both','Linewidth',2);
xlim([5e-5,1e-3]);ylim([2.5e-5,3e-4]);
set(gcf,'Position',[10 10 740 450]);

%% %%%%%%%%%%%%%%%%%%%%%%%%
figure;
h1 = subplot(3,1,1);
h1.Position = [0.1700    0.673    0.70    0.25];
loglog(kperp_10,Ekperp_pre2_10,'d-','Linewidth',2,'Markersize',2); hold on;
loglog(kperp_15,Ekperp_pre2_15,'d-','Linewidth',2,'Markersize',2); 
loglog(kperp_20,Ekperp_pre2_20,'d-','Linewidth',2,'Markersize',2); 
loglog(kperp_25,Ekperp_pre2_25,'d-','Linewidth',2,'Markersize',2); 
loglog(kperp_30,Ekperp_pre2_30,'d-','Linewidth',2,'Markersize',2); 
yyy2 = 0.18*kperp_10.^(-1/3);
loglog(kperp_10,yyy2,'--k','Linewidth',3);
xline(1.6e-4,':','Linewidth',1.5);
xline(3.0e-4,':','Linewidth',1.5);
xline(7.0e-4,':','Linewidth',1.5);
legend('','','','','','\propto k_{\perp}^{5/3-2}','','','NumColumns',3,'Location','best'); legend('boxoff');
set(gca,'Fontsize',23,'tickDir','both','Linewidth',2,'XTickLabel','');
ylabel({'k_{\perp}^{5/3}E_B_{A}(k_{\perp})','(km s^{-1})^2 km^{-2/3}'});
xlim([5e-5,1.0e-3]);

h2 = subplot(3,1,2);
h2.Position = [0.1700    0.40    0.70     0.25];
loglog(kperp_10,kpara_new_10,'+-','Linewidth',2,'Markersize',2); hold on;
loglog(kperp_15,kpara_new_15,'+-','Linewidth',2,'Markersize',2); 
loglog(kperp_20,kpara_new_20,'+-','Linewidth',2,'Markersize',2); 
loglog(kperp_25,kpara_new_25,'+-','Linewidth',2,'Markersize',2); 
loglog(kperp_30,kpara_new_30,'+-','Linewidth',2,'Markersize',2); 
xx1 = 5e-5:1e-5:1e-3;
yyy = 0.03*xx1.^(2/3);
loglog(xx1,yyy,':m','Linewidth',3);
xline(1.6e-4,':','Linewidth',1.5);
xline(3.0e-4,':','Linewidth',1.5);
xline(7.0e-4,':','Linewidth',1.5);
yline(7.0e-5,':','Linewidth',2);
legend('','','','','','\propto k_{\perp}^{2/3}','','','','NumColumns',3,'Location','northwest'); legend('boxoff');
ylabel({'k_{||}','(km^{-1})'}); 
set(gca,'Fontsize',23,'tickDir','both','Linewidth',2,'XTickLabel','');
xlim([5e-5,1.0e-3]);ylim([2.5e-5,4e-4]);

h3 = subplot(3,1,3);
h3.Position = [0.1700    0.125     0.70     0.25];
h3 = pcolor(kperp_30,kpara_30,XiB_mat_sc_30); hold on;
shading(gca,'flat'); colormap('jet'); c=colorbar; caxis(gca,[0 1]);
c.Location = 'eastoutside';
c.Label.String = '\chi_B_A';
xlim([5e-5,1.0e-3]); ylim([2.5e-5,4e-4]);
xline(1.6e-4,':','Linewidth',1.5);
xline(3.0e-4,':','Linewidth',1.5);
xline(7.0e-4,':','Linewidth',1.5);
xlabel('k_{\perp} (km^{-1})');ylabel({'k_{||}','(km^{-1})'}); 
set(gca,'Fontsize',23,'tickDir','both','Linewidth',2,'xscale','log','yscale','log');
set(gcf,'Position',[10 10 800 820]);
xx1 = 5e-5:1e-5:1e-3;
yyy = 0.36*xx1.^(1);
yyy2 = 0.023*xx1.^(2/3);
loglog(xx1,yyy, ':w','Linewidth',3);
loglog(xx1,yyy2,':m','Linewidth',3);
text(8e-5,2e-4,'k_{||}\propto k_{\perp}^{2/3}','FontSize',25,'Color','m');
text(8e-5,1e-4,'k_{||}\propto k_{\perp}','FontSize',25,'Color','w');

%% all xi pdf
h_xi = histogram(XiB_mat_sc_30(:));
yyy= h_xi.Values./length(kpara_30)./length(kperp_30);
xxx= h_xi.BinEdges;
plot(xxx(2:end),yyy,'Linewidth',4);
set(gca,'Fontsize',30,'tickDir','both','linewidth',2,'xscale','lin','yscale','log');
xlabel('\chi_B');ylabel('PDF(\chi_B)');xlim([4e-2,40]);

%% Zone 1
xi_1 = XiB_mat_sc_30(4:160,1:65);
xi_2 = XiB_mat_sc_30(4:160,66:140);
xi_3 = XiB_mat_sc_30(4:160,141:356);
xi_4 = XiB_mat_sc_30(4:160,141:517);

% zone 1
figure;
h_xi_1 =  histogram(xi_1(:));
yyy1= h_xi_1.Values./size(xi_1,1)./size(xi_1,2);
xxx1= h_xi_1.BinEdges;
plot(xxx1(2:end),yyy1,'Linewidth',4); hold on;
set(gca,'Fontsize',30,'tickDir','both','linewidth',2,'xscale','lin','yscale','log');
xlabel('\chi_B');ylabel('PDF(\chi_B)');xlim([4e-2,10]);
title('Zone 1');set(gcf,'Position',[10 10 800 820]);

% zone 2
figure;
h_xi_2 =  histogram(xi_2(:));
yyy2= h_xi_2.Values./size(xi_2,1)./size(xi_2,2);
xxx2= h_xi_2.BinEdges;
plot(xxx2(2:end),yyy2,'Linewidth',4); hold on;
set(gca,'Fontsize',30,'tickDir','both','linewidth',2,'xscale','lin','yscale','log');
xlabel('\chi_B');ylabel('PDF(\chi_B)');xlim([4e-2,20]);
title('Zone 2');set(gcf,'Position',[10 10 800 820]);

% zone 3
figure;
h_xi_3 =  histogram(xi_3(:));
yyy3= h_xi_3.Values./size(xi_3,1)./size(xi_3,2);
xxx3= h_xi_3.BinEdges;
plot(xxx3(2:end),yyy3,'Linewidth',4); hold on;
set(gca,'Fontsize',30,'tickDir','both','linewidth',2,'xscale','lin','yscale','log');
xlabel('\chi_B');ylabel('PDF(\chi_B)');xlim([4e-2,30]);
title('Zone 3');set(gcf,'Position',[10 10 800 820]);

% zone 4
figure;
h_xi_4 =  histogram(xi_4(:));
yyy4= h_xi_4.Values./size(xi_4,1)./size(xi_4,2);
xxx4= h_xi_4.BinEdges;
plot(xxx4(2:end),yyy4,'Linewidth',4); hold on;
set(gca,'Fontsize',30,'tickDir','both','linewidth',2,'xscale','lin','yscale','log');
xlabel('\chi_B');ylabel('PDF(\chi_B)');xlim([4e-2,30]);
title('Zone 4');set(gcf,'Position',[10 10 800 820]);

xx1 = logspace(-2,1);
yy1 = xx1.^(-2)/20;
yy2 = xx1.^(-5/3)/20;
yy3 = xx1.^(-3/2)/20;
hold on
loglog(xx1,yy1,'--','LineWidth',3);
loglog(xx1,yy2,':','LineWidth',3);
loglog(xx1,yy3,'.-','LineWidth',3);
legend('','-2','-5/3','-3/2'); legend('boxoff')
%% output data 
%      outputarray_new = strings(size(kperp_10,2),11);
%      outputarray_new(:,1)  = kperp_10;
%      outputarray_new(:,2)  = Ekperp_pre2_10;  % without 0.5
%      outputarray_new(:,3)  = Ekperp_pre2_15;
%      outputarray_new(:,4)  = Ekperp_pre2_20;  % without 0.5
%      outputarray_new(:,5)  = Ekperp_pre2_25;  % without 0.5
%      outputarray_new(:,6)  = Ekperp_pre2_30;
%      outputarray_new(:,7)  = kpara_new_10;  % without 0.5
%      outputarray_new(:,8)  = kpara_new_15;
%      outputarray_new(:,9)  = kpara_new_20;  % without 0.5
%      outputarray_new(:,10) = kpara_new_25;  % without 0.5
%      outputarray_new(:,11) = kpara_new_30;
% 
%      file_path = ['/Users/zhaosiqi/matlab_Figure/',date,'/',trange,'min-window/'];
% 
%      filename =[file_path,'data.txt'];
%      writematrix(outputarray_new,filename);
% 

% 
% XiB_mat_s = zeros(length(kpara_20),length(kperp_20));
% for i=1:1:length(kpara_20)-1
%     XiB_mat_s(i,:) = kperp_20.*bk_new_20./kpara_20(i)./106;
% end
% 
% figure; 
% pcolor(kperp_20,kpara_20,XiB_mat_s); hold on;
% shading(gca,'flat'); colormap('jet'); c=colorbar;
% c.Location = 'eastoutside';
% c.Label.String = '\chi_B_{\perp1}';
% set(gca,'Linewidth',2,'Fontsize',30,'YScale','log','XScale','log','tickDir','both');
% caxis(gca,[0 1]);
% xlim([4.5e-5,1.0e-3]);   ylim([2.5e-5,4e-4]);
% xline(1.6e-4,':','Linewidth',1.5);
% xline(3e-4,':','Linewidth',1.5);
% xline(6e-4,':','Linewidth',1.5);
% xlabel('k_{\perp} (km^{-1})');ylabel('k_{||} (km^{-1})');
% xx1 = 3e-4:1e-5:6e-4;
% yyy = 0.36*xx1.^(1);
 
% loglog(kperp_10,XiB_new_10,'o-','Linewidth',2,'Markersize',2); hold on;
% loglog(kperp_15,XiB_new_15,'o-','Linewidth',2,'Markersize',2);
% loglog(kperp_20,XiB_new_20,'o-','Linewidth',2,'Markersize',2);
% loglog(kperp_25,XiB_new_25,'o-','Linewidth',2,'Markersize',2);
% loglog(kperp_30,XiB_new_30,'o-','Linewidth',2,'Markersize',2);
% yline(1,':k','Linewidth',2);
% xline(1.6e-4,':','Linewidth',1.5);
% xline(3e-4,':','Linewidth',1.5);
% xline(6e-4,':','Linewidth',1.5);
% legend('\theta_{kk''}=10^\circ','\theta_{kk''}=15^\circ',...
%        '\theta_{kk''}=20^\circ','\theta_{kk''}=25^\circ','\theta_{kk''}=30^\circ','','','','NumColumns',3,'Location','best'); legend('boxoff');
% ylabel('\chi_B_{\perp1}');

%%
figure;
for i=26:10:536%length(kpara_20)    
    kpara_line = log10(kpara_20(i)).*ones(1,length(kperp_20));
    scatter(kperp_20,XiB_mat_20(i,:),5,kpara_line,'filled'); hold on;
end
yline(1,':','Linewidth',1.5);
ylim([1e-3,5e1]);   colormap jet;
xlabel('k_{\perp} (km^{-1})');   ylabel('\chi_B');
set(gca,'Fontsize',18,'tickDir','both','Linewidth',2,'xscale','log','yscale','log');
ccc = colorbar;
ccc.Location = 'south';
ccc.Position = [0.5433    0.2308    0.25    0.0724];
%%
figure;
for i=26:10:536%length(kpara_20)
    kpara_line = log10(kpara_20(i)).*ones(1,length(kperp_20));
    scatter(kperp_20,XiB_mat_sc_20(i,:),5,kpara_line,'filled'); hold on;
end
yline(1,':','Linewidth',1.5);
ylim([1e-3,1e1]);  colormap jet;   
xlabel('k_{\perp} (km^{-1})');   ylabel('\chi_B');
set(gca,'Fontsize',18,'tickDir','both','Linewidth',2,'xscale','log','yscale','log');
ccc = colorbar;
ccc.Location = 'south';
ccc.Position = [0.5433    0.2308    0.25    0.0724];

%% test Yan 2008
L0_para = 7.2e4;
L0_perp = 1.1e5;

MA = 0.4;

kpara_mat = repmat(kpara_20',1,length(kperp_20)); 
kperp_mat = repmat(kperp_20, length(kpara_20),1);  
E_mat = 0.5*bk2_sc_20./kpara_mat./kperp_mat;

I = zeros(length(kpara_20),length(kperp_20)); 
for i=1:1:length(kpara_20)
    kperp_norm = kperp_20*L0_perp;
    kpara_norm = kpara_20(i)*L0_para;
    A = exp(-1*L0_para^(1/3)*abs(kpara_20(i))./(MA^(4/3).*kperp_20.^(2/3)));
%     A = exp(-1*abs(kpara_norm)./MA^(4/3)./kperp_norm.^(2/3));
%     I(i,:) = kperp_norm.^(-10/3);%.*L0_perp^(-1/3).*MA^(4/3)/6/pi.*A;
    I(i,:) = kperp_20.^(-7/3).*A;
end

%  I(I<1e2) = nan;
% figure;
% pcolor(kperp_20,kpara_20,log10(I));       hold on;
% shading(gca,'flat');   colormap('jet');   c=colorbar; 
% caxis(gca,[1 15]);     c.Label.String = {'log10(I)';'(Yan et al. 2008)'};
% set(gca,'Fontsize',30,'tickDir','both','Linewidth',2,'xscale','log','yscale','log');
% xlabel('k_{\perp} (km^{-1})');   ylabel('k_{||} (km^{-1})'); 
% xlim([4.5e-5,1.0e-3]);           ylim([1e-5,1.0e-3]);
% xx1  = 2e-4:1e-5:1e-3;
% yyy  = 0.36*xx1.^(1);
% yyy2 = 0.015*xx1.^(2/3);
% loglog(xx1,yyy, ':w','Linewidth',3);
% loglog(xx1,yyy2,':k','Linewidth',3);
% set(gcf,'Position',[10 10 800 820]);
% 
% [M,color] = contour(kperp_20,kpara_20,log10(I),10);
% color.LineWidth = 2;   % color.LineColor = 'k'; 
% circles(0,0,5e-4,'facecolor','none','Linewidth',3,'linestyle',':','edgecolor','w');  

%% 
file_path = ['/Users/zhaosiqi/matlab_Figure/',date,'/',trange,'min-window/'];
angle_threshod = 20; 
%% loadin perp
load([file_path,'Pkk_2D_numk',num2str(numk),'_kklt',num2str(angle_threshod),'_VB_nonlinear'],...
      'PkkB_A','xvecA');

PkkB_A_norm = PkkB_A./max(max(PkkB_A));
PkkB_A_norm(PkkB_A_norm<1e-6) = nan;
figure;
pcolor(xvecA,xvecA,log10(PkkB_A_norm)); hold on;
shading(gca,'flat'); colormap('jet'); c=colorbar; caxis(gca,[-6 0]);  
c.Label.String = {'log_{10}P_{A}(k_{\perp},k_{||}) & log_{10}I_{A}(k_{\perp},k_{||})'};
set(gca,'Fontsize',30,'tickDir','both','linewidth',2,'xscale','log','yscale','log');
xlim([5e-5,1.0e-3]);ylim([1e-5,1.0e-3]);
% caxis(gca,[2 8]);
xlabel('k_{\perp} (km^{-1})');    ylabel('k_{||} (km^{-1})'); 

I_norm = I./max(max(I));
I_norm(I_norm<1e-6) = nan;
[M,color] = contour(kperp_20,kpara_20,log10(I_norm),10);
color.LineWidth = 5;   % color.LineColor = 'k'; 
[M,color2] = contour(kperp_20,kpara_20,log10(I_norm),10);
color2.LineStyle = '--';   color2.LineColor = 'k'; color2.LineWidth = 1.5;
circles(0,0,1.4e-04,'facecolor','none','Linewidth',2,'linestyle',':','edgecolor','k');  
circles(0,0,3*1.4e-04,'facecolor','none','Linewidth',2,'linestyle',':','edgecolor','k');  
xx1  = 1.6e-4:1e-5:1e-3;
yyy2 = 0.02*xx1.^(2/3);
loglog(xx1,yyy2,':k','Linewidth',3);
set(gcf,'Position',[10 10 800 650]);

%% %%%%%%%%%%%%%%%%%%%%%%%%%
temp = XiB_mat_sc_20;
temp(temp<1) = 0;
d_xi = abs(temp-1);
[aaa, ind_l] = min(d_xi,[],1);

kpara_line = zeros(length(kperp_20),1);
for i=1:1:length(kperp_20)
    kpara_line(i) = kpara_20(ind_l(i));
end

loglog(kperp_20,kpara_line,':k','Linewidth',3);

idx_perp = xvecA >= 4e-5;
idx_para = xvecA >= 5e-6;
kpara = xvecA(idx_para);
kperp = xvecA(idx_perp);
PkkB_A = PkkB_A(idx_para,idx_perp);

E_mat(PkkB_A<1e1) = nan;


figure; 
pcolor(kperp_20,kpara_20,log10(E_mat));   hold on;
shading(gca,'flat');   colormap('jet');   c=colorbar; caxis(gca,[8 12]);
   c.Label.String = {'log_{10}(E(k_{\perp},k_{||}))'};
set(gca,'Fontsize',30,'tickDir','both','Linewidth',2,'xscale','log','yscale','log');
xlabel('k_{\perp} (km^{-1})');   ylabel({'k_{||}','(km^{-1})'}); 
xlim([5e-5,1.0e-3]);   ylim([1e-5,1.0e-3]);
set(gcf,'Position',[10 10 800 820]);
xx1  = 3e-4:1e-5:6e-4;
yyy  = 0.36*xx1.^(1);
yyy2 = 0.023*xx1.^(2/3);
loglog(xx1,yyy, ':w','Linewidth',3);
loglog(xx1,yyy2,':k','Linewidth',3);

[M,color]=contour(kperp_20,kpara_20,log10(E_mat),10);
color.LineWidth = 2;    color.LineColor = 'k'; 
circles(0,0,5e-4,'facecolor','none','Linewidth',3,'linestyle',':','edgecolor','w');  

%% %%%%%%%%%%%%%%%%%%%%%%%%
% figure;
% h1 = subplot(4,1,1);
% h1.Position = [0.1700    0.7373    0.70    0.2];
% loglog(kperp_10,Ekperp_10,'d-','Linewidth',2,'MarkerSize',2);hold on;
% loglog(kperp_15,Ekperp_15,'d-','Linewidth',2,'MarkerSize',2);
% loglog(kperp_20,Ekperp_20,'d-','Linewidth',2,'MarkerSize',2);
% loglog(kperp_25,Ekperp_25,'d-','Linewidth',2,'MarkerSize',2);
% loglog(kperp_30,Ekperp_30,'d-','Linewidth',2,'MarkerSize',2);
% idx_p1 = 65:140;    idx_p2 = 141:301;
% [p1_10,S1] = polyfit(log(kperp_10(idx_p1)),log(Ekperp_10(idx_p1)),1); 
% [p2_10,S2] = polyfit(log(kperp_10(idx_p2)),log(Ekperp_10(idx_p2)),1);
% [p1_15,S1] = polyfit(log(kperp_15(idx_p1)),log(Ekperp_15(idx_p1)),1); 
% [p2_15,S2] = polyfit(log(kperp_15(idx_p2)),log(Ekperp_15(idx_p2)),1);
% [p1_20,S1] = polyfit(log(kperp_20(idx_p1)),log(Ekperp_20(idx_p1)),1); 
% [p2_20,S2] = polyfit(log(kperp_20(idx_p2)),log(Ekperp_20(idx_p2)),1);
% [p1_25,S1] = polyfit(log(kperp_25(idx_p1)),log(Ekperp_25(idx_p1)),1); 
% [p2_25,S2] = polyfit(log(kperp_25(idx_p2)),log(Ekperp_25(idx_p2)),1);
% [p1_30,S1] = polyfit(log(kperp_30(idx_p1)),log(Ekperp_30(idx_p1)),1); 
% [p2_30,S2] = polyfit(log(kperp_30(idx_p2)),log(Ekperp_30(idx_p2)),1);
% text(2.0e-4,5e7,['\alpha_1=',num2str(p1_10(1),3),',\alpha_2=',num2str(p2_10(1),3)],'FontSize',18,'Color','#0072BD');
% text(4.5e-4,5e7,['\alpha_1=',num2str(p1_15(1),3),',\alpha_2=',num2str(p2_15(1),3)],'FontSize',18,'Color','#D95319');
% text(2.0e-4,2e7,['\alpha_1=',num2str(p1_20(1),3),',\alpha_2=',num2str(p2_20(1),3)],'FontSize',18,'Color','#EDB120');
% text(4.5e-4,2e7,['\alpha_1=',num2str(p1_25(1),3),',\alpha_2=',num2str(p2_25(1),3)],'FontSize',18,'Color','#7E2F8E');
% text(2.0e-4,8e6,['\alpha_1=',num2str(p1_30(1),3),',\alpha_2=',num2str(p2_30(1),3)],'FontSize',18,'Color','#77AC30');
% xline(1.6e-4,':','Linewidth',1.5);
% xline(3.0e-4,':','Linewidth',1.5);
% xline(6.0e-4,':','Linewidth',1.5);
% ylabel({'E_B_{A}(k_{\perp})','(km s^{-1})^2 km'});
% set(gca,'Fontsize',23,'Linewidth',2,'XTickLabel','');
% xlim([5e-5,1.0e-3]); ylim([1e5,1e8]);
% 
% h2 = subplot(4,1,2);
% h2.Position = [0.1700    0.5232    0.70     0.2];
% loglog(kperp_10,Ekperp_pre2_10,'d-','Linewidth',2,'Markersize',2); hold on;
% loglog(kperp_15,Ekperp_pre2_15,'d-','Linewidth',2,'Markersize',2); 
% loglog(kperp_20,Ekperp_pre2_20,'d-','Linewidth',2,'Markersize',2); 
% loglog(kperp_25,Ekperp_pre2_25,'d-','Linewidth',2,'Markersize',2); 
% loglog(kperp_30,Ekperp_pre2_30,'d-','Linewidth',2,'Markersize',2); 
% yyy2 = 0.18*kperp_10.^(-1/3);
% loglog(kperp_10,yyy2,'--k','Linewidth',3);
% xline(1.6e-4,':','Linewidth',1.5);
% xline(3.0e-4,':','Linewidth',1.5);
% xline(6.0e-4,':','Linewidth',1.5);
% legend('','','','','','\propto k_{\perp}^{5/3-2}','','','NumColumns',3,'Location','best'); legend('boxoff');
% set(gca,'Fontsize',23,'tickDir','both','Linewidth',2,'XTickLabel','');
% ylabel({'k_{\perp}^{5/3}E_B_{A}(k_{\perp})','(km s^{-1})^2 km^{-2/3}'});
% xlim([5e-5,1.0e-3]);
% 
% h3 = subplot(4,1,3);
% h3.Position = [0.1700    0.3091     0.70     0.2];
% loglog(kperp_10,kpara_new_10,'+-','Linewidth',2,'Markersize',2); hold on;
% loglog(kperp_15,kpara_new_15,'+-','Linewidth',2,'Markersize',2); 
% loglog(kperp_20,kpara_new_20,'+-','Linewidth',2,'Markersize',2); 
% loglog(kperp_25,kpara_new_25,'+-','Linewidth',2,'Markersize',2); 
% loglog(kperp_30,kpara_new_30,'+-','Linewidth',2,'Markersize',2); 
% xx1 = 3e-4:1e-5:6e-4;
% yyy = 0.03*xx1.^(2/3);
% loglog(xx1,yyy,':k','Linewidth',3);
% xline(1.6e-4,':','Linewidth',1.5);
% xline(3.0e-4,':','Linewidth',1.5);
% xline(6.0e-4,':','Linewidth',1.5);
% yline(7.0e-5,':','Linewidth',2);
% legend('','','','','','\propto k_{\perp}^{2/3}','','','','NumColumns',3,'Location','northwest'); legend('boxoff');
% ylabel({'k_{||}','(km^{-1})'}); 
% set(gca,'Fontsize',23,'tickDir','both','Linewidth',2,'XTickLabel','');
% xlim([5e-5,1.0e-3]);ylim([2.5e-5,4e-4]);
% 
% h4 = subplot(4,1,4);
% h4.Position = [0.1700    0.0950     0.70     0.2];
% h4 = pcolor(kperp_30,kpara_30,XiB_mat_sc_30); hold on;
% shading(gca,'flat'); colormap('jet'); c=colorbar; caxis(gca,[0 1]);
% c.Location = 'eastoutside';
% c.Label.String = '\chi_B_A';
% xlim([5e-5,1.0e-3]);
% ylim([2.5e-5,4e-4]);
% xline(1.6e-4,':','Linewidth',1.5);
% xline(3.0e-4,':','Linewidth',1.5);
% xline(6.0e-4,':','Linewidth',1.5);
% xlabel('k_{\perp} (km^{-1})');ylabel({'k_{||}','(km^{-1})'}); 
% set(gca,'Fontsize',23,'tickDir','both','Linewidth',2,'xscale','log','yscale','log');
% set(gcf,'Position',[10 10 800 820]);
% xx1 = 3e-4:1e-5:6e-4;
% yyy = 0.36*xx1.^(1);
% yyy2 = 0.023*xx1.^(2/3);
% loglog(xx1,yyy, ':w','Linewidth',3);
% loglog(xx1,yyy2,':k','Linewidth',3);
