%% 3d sc frame
%% theta difference from k.B=0  and timing 
clear; close all;

%% output path
date = '20031203200000';%'20011219000000';%'20011216120000';
trange = '300';
file_path = ['/Users/zhaosiqi/matlab_Figure/',date,'/',trange,'min-window/'];
data_path = ['/Users/zhaosiqi/matlab_Figure/',date,'/data/nonlinear/',trange,'min-data/'];
%% parameters
numk = 400;   rci  = 74;
fci  = 0.24;  angle_threshod = 30;  d_sc = 2000;
flim = 0.5*fci;  
fmin = 0.00025;

%% load in sub window
filename  = 'cluster_list.txt';
list      = readtable(filename); 
tstart    = char(list.Var1);    tstop = char(list.Var2);
lines     = size(tstart,1);
data_name = strings(lines,1);

Power_MA =[]; Power_KA  =[];
kx_A =[];  ky_A =[];  kz_A =[];
kxav =[];  kyav =[];  kzav =[];
Bav = [];  a =[];  vxyz = [];  Va = []; 
%% Load data
for i=1:1:lines
  tinput = [tstart(i,1:10),'T',tstart(i,12:19),'.00Z/',tstop(i,1:10),'T',tstop(i,12:19),'.00Z'];
  Tints  = irf.tint(tinput);           % timing interval
  tstart0_str    = Tints(1).utc;       tend0_str = Tints(2).utc;
  data_name(i,:) = [tstart0_str(1:10),'-',tstart0_str(12:13),tstart0_str(15:16),tstart0_str(18:19),'-',tend0_str(12:13),tend0_str(15:16),tend0_str(18:19),'VB'];
  S = load(strcat(data_path,data_name(i)));
  kx_A = [kx_A;S.kx_A];
  ky_A = [ky_A;S.ky_A];
  kz_A = [kz_A;S.kz_A];
  kxav = [kxav; S.kxav]; 
  kyav = [kyav; S.kyav]; 
  kzav = [kzav; S.kzav]; 
  Power_MA = [Power_MA; S.Power_MA];
  Power_KA = [Power_KA; S.Power_KA ];
  Bav  = [Bav; S.Bav];
  a    = [a;   S.a];
  Va   = [Va;  S.Va];
  vxyz = [vxyz;S.vxyz];
end
numf = S.numf;  fre1 = S.fre;  t_length = str2num(trange)*60;     

%%
Bavxmat = Bav(:,1)*ones(1,numf);
Bavymat = Bav(:,2)*ones(1,numf);
Bavzmat = Bav(:,3)*ones(1,numf);
Bavabsmat = sqrt(Bav(:,1).^2+Bav(:,2).^2+Bav(:,3).^2)*ones(1,numf);

%% Doppler shift
[k0_A_pf,kx_A_pf,ky_A_pf,kz_A_pf,freq_mean_A,freq_A_A,freq_f_A,freq_s_A,vp_mean_A,theta_A,f_doppler_A] = doppler_shift(kx_A,ky_A,kz_A,fre1,Bav,vxyz,Va,a); 

kmag_A = sqrt(kx_A.*kx_A + ky_A.*ky_A + kz_A.*kz_A); 
%%  MHD threthod
f_2_t = round(2./t_length,3);
Power_MA(freq_mean_A>flim)   = 0; 
Power_MA(freq_mean_A<f_2_t)  = 0; 
Power_MA(isnan(freq_mean_A)) = 0; 
Power_KA(freq_mean_A>flim)   = 0; 
Power_KA(freq_mean_A<f_2_t)  = 0; 
Power_KA(isnan(freq_mean_A)) = 0; 

%% wave vectors: spacecraft frame
kpar_A     = (kx_A.*Bavxmat + ky_A.*Bavymat + kz_A.*Bavzmat)./Bavabsmat;
kperp_A    = sqrt(kmag_A.^2 - kpar_A.^2);
kpar_A_abs = abs(kpar_A);
  
%% theta from timing
kmag_av  = sqrt(kxav.*kxav + kyav.*kyav + kzav.*kzav);
dtheta_A = acosd((kx_A.*kxav + ky_A.*kyav + kz_A.*kzav)./kmag_av./kmag_A);

obtuse_theta_A = dtheta_A > 90;      % theta < 90
dtheta_A(obtuse_theta_A)  = 180 - dtheta_A(obtuse_theta_A);

%% f threshod
dtheta_A(freq_mean_A>flim)   = 0; 
dtheta_A(freq_mean_A<f_2_t)  = 0; 
dtheta_A(isnan(freq_mean_A)) = 90;   % freq_mean is nan delete 
%% deNAN 
vp_mean_A(freq_mean_A>flim)     = 0; 
freq_mean_A(freq_mean_A>flim)   = 0; 
freq_mean_A(isnan(freq_mean_A)) = 0;

Power_MA(freq_mean_A==0) = 0; 
Power_KA(freq_mean_A==0) = 0; 
%% k threshod: k rho_ci < 0.1
klim = 0.1/rci;
Power_MA(kmag_A   > klim) = 0; % threshod |k|>klim, Power = 0
Power_KA(kmag_A   > klim) = 0; 
kpar_A_abs(kmag_A > klim) = 0; % all k(|k|>klim) = 0 
kperp_A(kmag_A    > klim) = 0;
dtheta_A(kmag_A   > klim) = 0; % |k|>klim, theta = 0  
kmin = 1/10/d_sc; % /2 sqzhao
Power_MA(kmag_A   < kmin) = 0;
Power_KA(kmag_A   < kmin) = 0; 
kmag_A(kmag_A     > klim) = 0; % add 20220710
% use same length scale
kmax_A = klim*1.1; 
kmin_A = -kmax_A;
% kmagvec_A = linspace(0,kmax_A,numk);
% dkmag_A = kmax_A/numk;

%% theta kk' threshod
% angle_threshod = 360; % if without angle threshod
% angle_threshod = angle_threshod;
ind_dtheta_gt  = dtheta_A > angle_threshod;   % phi_kk' > 10
Power_MA(ind_dtheta_gt) = 0;   % angle threshod, Power = 0
Power_KA(ind_dtheta_gt) = 0; 

%% normalize to time points
numt = size(Power_MA,1);
f_mat0 = repmat(fre1,1,numt);  f_mat0 = f_mat0';
% fpower_MA = Power_MA.*f_mat0;
% fpower_KA = Power_KA.*f_mat0;

idx_M = Power_MA~=0;
numt_data_M = nansum(idx_M,1);  numt_data_M(numt_data_M==0)=1;
numt_mat_M  = repmat(numt_data_M',1,numt); numt_mat_M = numt_mat_M';
mean_Power_MA = Power_MA./numt_mat_M;
% mean_fpower_MA = fpower_MA./numt_mat_M;

idx_K = Power_KA~=0;
numt_data_K = nansum(idx_K,1);  numt_data_K(numt_data_K==0) = 1;
numt_mat_K  = repmat(numt_data_K',1,numt); numt_mat_K = numt_mat_K';
mean_Power_KA = Power_KA./numt_mat_K;
% mean_fpower_KA = fpower_KA./numt_mat_K;

% df*P 
df_sc = abs(diff(fre1));    df_sc = [df_sc;0];
df_sc_mat = repmat(df_sc,1,numt); df_sc_mat = df_sc_mat';
df_mean_Power_MA = df_sc_mat.*mean_Power_MA;
df_mean_Power_KA = df_sc_mat.*mean_Power_KA;
df_Power_MA = df_sc_mat.*Power_MA;
df_Power_KA = df_sc_mat.*Power_KA;
%%
clear Bavabsmat Bavxmat Bavymat Bavzmat dtheta_A f_doppler_A ... 
      freq_f_A freq_s_A freq_A_A idx_K idx_M ind_dtheta_gt ...
      kmag_av kx_A kx_A_pf ky_A ky_A_pf kz_A kz_A_pf kxav kyav kzav
% frequency map
numf0 = 200; % 800=>1600
afmin = log10(fmin);
afmax = log10(flim);
af = logspace(afmin,afmax,numf0);

diff_af = abs(diff(af)); diff_af = [diff_af,0];
% k 
amax=log10(kmax_A); 
amin=log10(1e-6);   % 1e-6
a = logspace(amin,amax,numk);
xvecA = a;
diff_a = abs(diff(a)); diff_a = [diff_a,0];
%% %%%%%%%%%%%%%%%      3D      %%%%%%%%%%%%%%%%
%% Plasma frame:  fplasma  |kpar| |kperp|
df_mean_Power_MA_fpkk = zeros(numf0,numk,numk);
df_mean_Power_KA_fpkk = zeros(numf0,numk,numk); 
fpf_mat = [];  kpara_mat = [];  kperp_mat = [];  PfpB_mat = [];
for mm = 1:numt
	for nn = 1:numf
      [kpara_min,idx_kpara] = min(abs(kpar_A_abs(mm,nn) - a)); 
      [kperp_min,idx_kperp] = min(abs(kperp_A(mm,nn)    - a));  
      [af_min,idx_af]       = min(abs(freq_mean_A(mm,nn)- af));  
      df_mean_Power_MA_fpkk(idx_af,idx_kpara,idx_kperp) = df_mean_Power_MA_fpkk(idx_af,idx_kpara,idx_kperp) + df_mean_Power_MA(mm,nn)./diff_a(idx_kpara)./diff_a(idx_kperp)./diff_af(idx_af);
      df_mean_Power_KA_fpkk(idx_af,idx_kpara,idx_kperp) = df_mean_Power_KA_fpkk(idx_af,idx_kpara,idx_kperp) + df_mean_Power_KA(mm,nn)./diff_a(idx_kpara)./diff_a(idx_kperp)./diff_af(idx_af);
    end  
end

% Pkf = zeros(length(af),length(xvecA));
% for  i=1:1:length(af)-1
%     for j = 1:1:length(xvecA)-1
%        Pkf(i,j) = trapz(flipud(xvecA),flipud(df_mean_Power_MA_fpkk(i,:,j)));
%     end
% end

% test 
Pkf = zeros(length(af),length(xvecA));
for  i=1:1:length(af)-1
    for j = 1:1:length(xvecA)-1
       Pkf(i,j) = nansum(df_mean_Power_MA_fpkk(i,:,j).*diff_a);
    end
end

Pkf = abs(Pkf);
Pkf(isinf(Pkf)) = nan;
Pkf_norm = Pkf./max(max(Pkf));
Pkf_norm(Pkf_norm<1e-6) = nan;
figure;
pcolor(xvecA,af,log10(Pkf_norm)); hold on;
shading(gca,'flat'); colormap('jet'); c=colorbar; 
c.Label.String = 'log_{10}P_{B_A}(k_{\perp},f_{rest})' ;
set(gca,'Fontsize',35,'tickDir','both','linewidth',2,'xscale','log','yscale','log');
xx1 = logspace(-4,-2.8);
yyy = 2*xx1.^(2/3);
yyy2 = 42*xx1;
loglog(xx1,yyy,'--','linewidth',4,'color','k');
loglog(xx1,yyy2,':','linewidth',4,'color','k');
xlabel('k_{\perp} (km^{-1})');ylabel('f_{rest} (Hz)');
xlim([1e-5,1.2e-3]);ylim([5e-4,3e-2]);

k_para = 7e-5;  Va = 106;  dVa = 11;dk = 0;%3e-6;
fA   = k_para*Va/2/pi;
xxA  = logspace(-6,-4,8);
yyyA = fA*ones(length(xxA),1);
f_err_2 = dk^2*Va^2/(2*pi)^2 + k_para^2/(2*pi)^2*dVa^2;
f_err = sqrt(f_err_2);
y_err = f_err*ones(length(xxA),1);caxis(gca,[-4 -0.5]);
e1 = errorbar(xxA,yyyA,y_err,':','linewidth',2,'color','k');hold on;

k_para = 1e-4;  Va = 106;  dVa = 11; dk = 0;%3e-6;
fA   = k_para*Va/2/pi;
xxA  = logspace(-6,-4,8);
yyyA = fA*ones(length(xxA),1);
f_err_2 = dk^2*Va^2/(2*pi)^2 + k_para^2/(2*pi)^2*dVa^2;
f_err = sqrt(f_err_2);
y_err = f_err*ones(length(xxA),1);caxis(gca,[-4 -0.5]);
e2 = errorbar(xxA,yyyA,y_err,':','linewidth',2,'color','k');hold on;
caxis(gca,[-6 -2]);
set(gcf,'Position',[10 10 800 650]);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% %% fix k_kerp  X zone
% [m, idx_s] = min(abs(xvecA - 3.5e-4));
% [m, idx_e] = min(abs(xvecA - 4e-4));
% 
% figure;
% dnumk = 100;
% kpara_color = zeros(1,length(f_mat));
% for i = 1:1:(length(xvecA(7:163))-dnumk) % 270
%     sum_kpara_tmp       = squeeze(nansum(df_mean_Power_MA_fpkk(:,i:(i+dnumk),:),2));
%     sum_kpara_kperp_tmp = squeeze(nansum(sum_kpara_tmp(:,idx_s:idx_e),2));
%     sum_kpara_kperp_tmp = sum_kpara_kperp_tmp./sum(sum_kpara_kperp_tmp);
%     kpara_color(1,:) = 106*xvecA(i+dnumk/2-1)/2/pi;
%     scatter(f_mat,sum_kpara_kperp_tmp,15,kpara_color,'filled'); hold on;
%     f_mat_part_idx = f_mat<106*xvecA(i+dnumk/2-1)/2/pi;
% %     scatter(f_mat(f_mat_part_idx),sum_kpara_kperp_tmp(f_mat_part_idx),50,'k','filled'); 
%     scatter(f_mat(f_mat_part_idx),sum_kpara_kperp_tmp(f_mat_part_idx),40,kpara_color(f_mat_part_idx)); 
%     scatter(f_mat(f_mat_part_idx),sum_kpara_kperp_tmp(f_mat_part_idx),15,'k','filled'); 
% end
% % text(1e-3,1e-4,'\propto f_{rest}^{-2}','FontSize',35);
% text(1e-2,1e-2,'\propto f_{rest}^{-2}','FontSize',35);
% xlabel('f_{rest} (Hz)');  zlabel('f_A (Hz)');%ylabel('P_B_A');
% xlim([4e-4,3.0e-2]);ylim([5e-6,1.0e-1]);
% c=colorbar; colormap jet; c.Label.String = 'f_A (Hz)';
% set(gca,'Fontsize',30,'tickDir','both','linewidth',2,'xscale','log','yscale','log','zscale','lin');
% 
% xxx1 = 2.5e-3:1e-4:1e-2;
% P_t1 = 1e-7*xxx1.^(-2);
% loglog(xxx1,P_t1,'--k','LineWidth',3); hold on;
% 
% %% %% fix k_kpara  X zone
% [m, idx_s] = min(abs(xvecA - 6e-5));
% [m, idx_e] = min(abs(xvecA - 1e-4));
% 
% figure;
% dnumk = 10;
% kperp_color = zeros(1,length(f_mat));
% for i = 1:1:(length(xvecA(1:55))-dnumk) %270
%     sum_kperp_tmp       = squeeze(nansum(df_mean_Power_MA_fpkk(:,:,i:(i+dnumk)),3));
%     sum_kpara_kperp_tmp = squeeze(nansum(sum_kperp_tmp(:,idx_s:idx_e),2));
%     sum_kpara_kperp_tmp_norm = sum_kpara_kperp_tmp./sum(sum_kpara_kperp_tmp);
%     kperp_color(1,:) = xvecA(i+dnumk/2-1);
%     scatter(f_mat,sum_kpara_kperp_tmp_norm,15,log10(kperp_color),'filled'); hold on;
% end
% set(gca,'Fontsize',30,'tickDir','both','linewidth',2,'xscale','log','yscale','log');
% c=colorbar; colormap jet; 
% xlabel('f_{rest} (Hz)'); ylabel('P_B_A');  
% xlim([2e-4,3.0e-2]);ylim([1e-4,1.0e1]);
% 
% fA1 = 106*xvecA(idx_s)/2/pi;
% fA2 = 106*xvecA(idx_e)/2/pi;
% xline(fA1,':','linewidth',3,'color','#EDB120');
% xline(fA2,':','linewidth',3,'color','#0072BD');
% title('normarlized P; k_{||}~[6*10^{-5},1*10^{-5}] km^{-1}')
% 
% xxx1 = 2e-4:1e-4:1e-2;
% P_t1  = 1e-5*xxx1.^(-5/3);
% P_t2  = 1e-5*xxx1.^(-3/2);
% P_t3  = 1e-5*xxx1.^(-2);
% loglog(xxx1,P_t1,'LineWidth',2);
% loglog(xxx1,P_t2,'LineWidth',2);
% loglog(xxx1,P_t3,'LineWidth',2);
% 
% %% fix k_kerp  X zone
% sum_kpara_1 = squeeze(nansum(df_mean_Power_MA_fpkk(:,28:39,:),2));
% sum_kpara_2 = squeeze(nansum(df_mean_Power_MA_fpkk(:,39:55,:),2));
% sum_kpara_3 = squeeze(nansum(df_mean_Power_MA_fpkk(:,55:71,:),2));
% sum_kpara_4 = squeeze(nansum(df_mean_Power_MA_fpkk(:,71:109,:),2));
% sum_kpara_5 = squeeze(nansum(df_mean_Power_MA_fpkk(:,109:135,:),2));
% sum_kpara_6 = squeeze(nansum(df_mean_Power_MA_fpkk(:,135:163,:),2));
% sum_kpara_kperp_1 = squeeze(nansum(sum_kpara_1(:,163:270),2));
% sum_kpara_kperp_2 = squeeze(nansum(sum_kpara_2(:,163:270),2));
% sum_kpara_kperp_3 = squeeze(nansum(sum_kpara_3(:,163:270),2));
% sum_kpara_kperp_4 = squeeze(nansum(sum_kpara_4(:,163:270),2));
% sum_kpara_kperp_5 = squeeze(nansum(sum_kpara_5(:,163:270),2));
% sum_kpara_kperp_6 = squeeze(nansum(sum_kpara_6(:,163:270),2));
% figure;
% loglog(f_mat,sum_kpara_kperp_1,'*','linewidth',2,'MarkerSize',8); hold on;
% loglog(f_mat,sum_kpara_kperp_2,'*','linewidth',2,'MarkerSize',8); hold on;
% loglog(f_mat,sum_kpara_kperp_3,'*','linewidth',2,'MarkerSize',8); hold on;
% loglog(f_mat,sum_kpara_kperp_4,'*','linewidth',2,'MarkerSize',8); hold on;
% loglog(f_mat,sum_kpara_kperp_5,'*','linewidth',2,'MarkerSize',8); hold on;
% loglog(f_mat,sum_kpara_kperp_6,'*','linewidth',2,'MarkerSize',8); hold on;
% 
% xlabel('f_{rest} (Hz)'); ylabel('E_B');  legend('boxoff','Location','bestoutside')
% set(gca,'linewidth',2,'Fontsize',30,'tickDir','both');
% legend('k_{||}~[5*10^{-5},7*10^{-5}]; k_{\perp}~[3*10^{-4},5*10^{-4}]',...
%        'k_{||}~[7*10^{-5},1*10^{-4}]; k_{\perp}~[3*10^{-4},5*10^{-4}]',...
%        'k_{||}~[1*10^{-4},1.3*10^{-4}]; k_{\perp}~[3*10^{-4},5*10^{-4}]',...
%        'k_{||}~[1.3*10^{-4},2*10^{-4}]; k_{\perp}~[3*10^{-4},5*10^{-4}]',...
%        'k_{||}~[2*10^{-4},2.5*10^{-4}]; k_{\perp}~[3*10^{-4},5*10^{-4}]',...
%        'k_{||}~[2.5*10^{-4},3*10^{-4}]; k_{\perp}~[3*10^{-4},5*10^{-4}]');
% % 
% sum_kpara_1 = squeeze(nansum(df_mean_Power_MA_fpkk(:,34:55,:),2));
% sum_kpara_2 = squeeze(nansum(df_mean_Power_MA_fpkk(:,55:109,:),2));
% sum_kpara_3 = squeeze(nansum(df_mean_Power_MA_fpkk(:,109:163,:),2));
% sum_kpara_kperp_1 = squeeze(nansum(sum_kpara_1(:,163:270),2));
% sum_kpara_kperp_2 = squeeze(nansum(sum_kpara_2(:,163:270),2));
% sum_kpara_kperp_3 = squeeze(nansum(sum_kpara_3(:,163:270),2));
% figure;
% loglog(f_mat,sum_kpara_kperp_1,'*','linewidth',2,'MarkerSize',8,'color','#0072BD'); hold on;
% loglog(f_mat,sum_kpara_kperp_2,'*','linewidth',2,'MarkerSize',8,'color','#D95319'); hold on;
% loglog(f_mat,sum_kpara_kperp_3,'*','linewidth',2,'MarkerSize',8,'color','#EDB120'); hold on;
% xlabel('f_{rest} (Hz)'); ylabel('E_B');  legend('boxoff','Location','bestoutside')
% set(gca,'linewidth',2,'Fontsize',30,'tickDir','both');
% legend('k_{||}~[6*10^{-5},1*10^{-4}]; k_{\perp}~[3*10^{-4},5*10^{-4}]',...
%        'k_{||}~[1*10^{-4},2*10^{-4}]; k_{\perp}~[3*10^{-4},5*10^{-4}]',...
%        'k_{||}~[2*10^{-4},3*10^{-4}]; k_{\perp}~[3*10^{-4},5*10^{-4}]');
% 
% dVa = 11; Va = 106; 
% fA1 = xvecA(34)*Va/2/pi;
% fA2 = xvecA(55)*Va/2/pi;
% fA3 = xvecA(109)*Va/2/pi;
% fA4 = xvecA(163)*Va/2/pi;
% xline(fA1,':','linewidth',3,'color','#0072BD');
% xline(fA2,':','linewidth',3,'color','#D95319');
% xline(fA3,':','linewidth',3,'color','#EDB120');
% xline(fA4,':','linewidth',3,'color','#EDB120');
% %% fix kpara  4 zone
% sum_kpara_1 = squeeze(nansum(df_mean_Power_MA_fpkk(:,34:109,:),2));
% sum_kpara_kperp1 = squeeze(nansum(sum_kpara_1(:,28:55),2)); %28:87
% sum_kpara_kperp2 = squeeze(nansum(sum_kpara_1(:,55:163),2));
% sum_kpara_kperp3 = squeeze(nansum(sum_kpara_1(:,163:270),2));
% sum_kpara_kperp4 = squeeze(nansum(sum_kpara_1(:,270:378),2));
% 
% figure;
% loglog(f_mat,sum_kpara_kperp1,'*','linewidth',2,'MarkerSize',8); hold on;
% loglog(f_mat,sum_kpara_kperp2,'*','linewidth',2,'MarkerSize',8); hold on;
% loglog(f_mat,sum_kpara_kperp3,'*','linewidth',2,'MarkerSize',8); hold on;
% loglog(f_mat,sum_kpara_kperp4,'*','linewidth',2,'MarkerSize',8); hold on;
% xlabel('f_{rest} (Hz)'); ylabel('E_B');  legend('boxoff','Location','bestoutside')
% set(gca,'linewidth',2,'Fontsize',30,'tickDir','both');
% legend('k_{||}~[6*10^{-5},2*10^{-4}]; k_{\perp}~[5*10^{-5},1*10^{-4}]',...
%        'k_{||}~[6*10^{-5},2*10^{-4}]; k_{\perp}~[1*10^{-4},3*10^{-4}]',...
%        'k_{||}~[6*10^{-5},2*10^{-4}]; k_{\perp}~[3*10^{-4},5*10^{-4}]',...
%        'k_{||}~[6*10^{-5},2*10^{-4}]; k_{\perp}~[5*10^{-4},7*10^{-4}]');
% 
% %
% sum_kpara_1 = squeeze(nansum(df_mean_Power_MA_fpkk(:,34:163,:),2));
% sum_kpara_kperp1 = squeeze(nansum(sum_kpara_1(:,28:55),2)); %28:87
% sum_kpara_kperp2 = squeeze(nansum(sum_kpara_1(:,55:163),2));
% sum_kpara_kperp3 = squeeze(nansum(sum_kpara_1(:,163:270),2));
% sum_kpara_kperp4 = squeeze(nansum(sum_kpara_1(:,270:378),2));
% 
% figure;
% loglog(f_mat,sum_kpara_kperp1,'*','linewidth',2,'MarkerSize',8); hold on;
% loglog(f_mat,sum_kpara_kperp2,'*','linewidth',2,'MarkerSize',8); hold on;
% loglog(f_mat,sum_kpara_kperp3,'*','linewidth',2,'MarkerSize',8); hold on;
% loglog(f_mat,sum_kpara_kperp4,'*','linewidth',2,'MarkerSize',8); hold on;
% xlabel('f_{rest} (Hz)'); ylabel('E_B');  legend('boxoff','Location','bestoutside')
% set(gca,'linewidth',2,'Fontsize',30,'tickDir','both');
% legend('k_{||}~[6*10^{-5},3*10^{-4}]; k_{\perp}~[5*10^{-5},1*10^{-4}]',...
%        'k_{||}~[6*10^{-5},3*10^{-4}]; k_{\perp}~[1*10^{-4},3*10^{-4}]',...
%        'k_{||}~[6*10^{-5},3*10^{-4}]; k_{\perp}~[3*10^{-4},5*10^{-4}]',...
%        'k_{||}~[6*10^{-5},3*10^{-4}]; k_{\perp}~[5*10^{-4},7*10^{-4}]');
% 
% 
% % fix kpara  X zone
% sum_kpara_1 = squeeze(nansum(df_mean_Power_MA_fpkk(:,34:109,:),2));
% sum_kpara_kperp1 = squeeze(nansum(sum_kpara_1(:,28:39),2)); %28:87
% sum_kpara_kperp2 = squeeze(nansum(sum_kpara_1(:,39:55),2));
% sum_kpara_kperp3 = squeeze(nansum(sum_kpara_1(:,55:71),2));
% sum_kpara_kperp4 = squeeze(nansum(sum_kpara_1(:,71:109),2));
% sum_kpara_kperp5 = squeeze(nansum(sum_kpara_1(:,109:135),2)); %28:87
% sum_kpara_kperp6 = squeeze(nansum(sum_kpara_1(:,135:163),2));
% sum_kpara_kperp7 = squeeze(nansum(sum_kpara_1(:,163:270),2));
% sum_kpara_kperp8 = squeeze(nansum(sum_kpara_1(:,270:378),2));
% 
% 
% figure;
% loglog(f_mat,sum_kpara_kperp1,'*','linewidth',2,'MarkerSize',8); hold on;
% loglog(f_mat,sum_kpara_kperp2,'*','linewidth',2,'MarkerSize',8); hold on;
% loglog(f_mat,sum_kpara_kperp3,'*','linewidth',2,'MarkerSize',8); hold on;
% loglog(f_mat,sum_kpara_kperp4,'*','linewidth',2,'MarkerSize',8); hold on;
% loglog(f_mat,sum_kpara_kperp5,'*','linewidth',2,'MarkerSize',8); hold on;
% loglog(f_mat,sum_kpara_kperp6,'*','linewidth',2,'MarkerSize',8); hold on;
% loglog(f_mat,sum_kpara_kperp7,'*','linewidth',2,'MarkerSize',8); hold on;
% loglog(f_mat,sum_kpara_kperp8,'*','linewidth',2,'MarkerSize',8); hold on;
% xlabel('f_{rest} (Hz)'); ylabel('E_B');  legend('boxoff','Location','bestoutside')
% set(gca,'linewidth',2,'Fontsize',30,'tickDir','both');
% legend('k_{||}~[6*10^{-5},2*10^{-4}]; k_{\perp}~[5*10^{-5},1*10^{-4}]',...
%        'k_{||}~[6*10^{-5},2*10^{-4}]; k_{\perp}~[1*10^{-4},3*10^{-4}]',...
%        'k_{||}~[6*10^{-5},2*10^{-4}]; k_{\perp}~[3*10^{-4},5*10^{-4}]',...
%        'k_{||}~[6*10^{-5},2*10^{-4}]; k_{\perp}~[5*10^{-4},7*10^{-4}]');
% 
% %
% sum_kpara_1 = squeeze(nansum(df_mean_Power_MA_fpkk(:,34:163,:),2));
% sum_kpara_kperp1 = squeeze(nansum(sum_kpara_1(:,28:55),2)); %28:87
% sum_kpara_kperp2 = squeeze(nansum(sum_kpara_1(:,55:163),2));
% sum_kpara_kperp3 = squeeze(nansum(sum_kpara_1(:,163:270),2));
% sum_kpara_kperp4 = squeeze(nansum(sum_kpara_1(:,270:378),2));
% 
% figure;
% loglog(f_mat,sum_kpara_kperp1,'*','linewidth',2,'MarkerSize',8); hold on;
% loglog(f_mat,sum_kpara_kperp2,'*','linewidth',2,'MarkerSize',8); hold on;
% loglog(f_mat,sum_kpara_kperp3,'*','linewidth',2,'MarkerSize',8); hold on;
% loglog(f_mat,sum_kpara_kperp4,'*','linewidth',2,'MarkerSize',8); hold on;
% xlabel('f_{rest} (Hz)'); ylabel('E_B');  legend('boxoff','Location','bestoutside')
% set(gca,'linewidth',2,'Fontsize',30,'tickDir','both');
% legend('k_{||}~[6*10^{-5},3*10^{-4}]; k_{\perp}~[5*10^{-5},1*10^{-4}]',...
%        'k_{||}~[6*10^{-5},3*10^{-4}]; k_{\perp}~[1*10^{-4},3*10^{-4}]',...
%        'k_{||}~[6*10^{-5},3*10^{-4}]; k_{\perp}~[3*10^{-4},5*10^{-4}]',...
%        'k_{||}~[6*10^{-5},3*10^{-4}]; k_{\perp}~[5*10^{-4},7*10^{-4}]');
% 
% 
% %% test E-f
% % diff Zone   small range k_para
% sum_kpara_PB_Z1 = squeeze(nansum(df_mean_Power_MA_fpkk(:,34:55,:),2));
% sum_kpara_kperp_PB_Z1 = squeeze(nansum(sum_kpara_PB_Z1(:,28:87),2)); %28:87
% sum_kpara_kperp_PB_Z2 = squeeze(nansum(sum_kpara_PB_Z1(:,87:163),2));
% sum_kpara_kperp_PB_Z3 = squeeze(nansum(sum_kpara_PB_Z1(:,163:324),2));
% % diff Zone   all range k_para
% sum_kpara_PB_Z1_a = squeeze(nansum(df_mean_Power_MA_fpkk(:,34:163,:),2));
% sum_kpara_kperp_PB_Z1_a = squeeze(nansum(sum_kpara_PB_Z1_a(:,28:87),2));
% sum_kpara_kperp_PB_Z2_a = squeeze(nansum(sum_kpara_PB_Z1_a(:,87:163),2));
% sum_kpara_kperp_PB_Z3_a = squeeze(nansum(sum_kpara_PB_Z1_a(:,163:324),2));
% 
% % diff Zone   right range k_para
% sum_kpara_PB_Z1_b = squeeze(nansum(df_mean_Power_MA_fpkk(:,55:163,:),2));
% sum_kpara_kperp_PB_Z1_b = squeeze(nansum(sum_kpara_PB_Z1_b(:,28:87),2));
% sum_kpara_kperp_PB_Z2_b = squeeze(nansum(sum_kpara_PB_Z1_b(:,87:163),2));
% sum_kpara_kperp_PB_Z3_b = squeeze(nansum(sum_kpara_PB_Z1_b(:,163:324),2));
% 
% figure;
% loglog(f_mat,sum_kpara_kperp_PB_Z1,'*','linewidth',2,'MarkerSize',10); hold on;
% loglog(f_mat,sum_kpara_kperp_PB_Z2,'*','linewidth',2,'MarkerSize',10);
% loglog(f_mat,sum_kpara_kperp_PB_Z3,'*','linewidth',2,'MarkerSize',10);
% loglog(f_mat,sum_kpara_kperp_PB_Z1_a,'s','linewidth',2,'MarkerSize',10); hold on;
% loglog(f_mat,sum_kpara_kperp_PB_Z2_a,'s','linewidth',2,'MarkerSize',10);
% loglog(f_mat,sum_kpara_kperp_PB_Z3_a,'s','linewidth',2,'MarkerSize',10);
% loglog(f_mat,sum_kpara_kperp_PB_Z1_b,'+','linewidth',2,'MarkerSize',10); hold on;
% loglog(f_mat,sum_kpara_kperp_PB_Z2_b,'+','linewidth',2,'MarkerSize',10);
% loglog(f_mat,sum_kpara_kperp_PB_Z3_b,'+','linewidth',2,'MarkerSize',10);
% 
% legend('k_{||}~[6*10^{-5},1*10^{-4}]; k_{\perp}~[5*10^{-5},1.6*10^{-4}]',...
%        'k_{||}~[6*10^{-5},1*10^{-4}]; k_{\perp}~[1.6*10^{-4},3*10^{-4}]',...
%        'k_{||}~[6*10^{-5},1*10^{-4}]; k_{\perp}~[3*10^{-4},6*10^{-4}]',...
%        'k_{||}~[6*10^{-5},3*10^{-4}]; k_{\perp}~[5*10^{-5},1.6*10^{-4}]',...
%        'k_{||}~[6*10^{-5},3*10^{-4}]; k_{\perp}~[1.6*10^{-4},3*10^{-4}]',...
%        'k_{||}~[6*10^{-5},3*10^{-4}]; k_{\perp}~[3*10^{-4},6*10^{-4}]',...
%        'k_{||}~[1*10^{-4},3*10^{-4}]; k_{\perp}~[5*10^{-5},1.6*10^{-4}]',...
%        'k_{||}~[1*10^{-4},3*10^{-4}]; k_{\perp}~[1.6*10^{-4},3*10^{-4}]',...
%        'k_{||}~[1*10^{-4},3*10^{-4}]; k_{\perp}~[3*10^{-4},6*10^{-4}]');
% xlabel('f_{rest} (Hz)'); ylabel('P_B');  legend('boxoff','Location','bestoutside')
% set(gca,'linewidth',2,'Fontsize',30,'tickDir','both');
% 
% %
% sum_kpara_PB_Z1 = squeeze(nansum(df_mean_Power_MA_fpkk(:,28:44,:),2));
% sum_kpara_kperp_PB_Z1 = squeeze(nansum(sum_kpara_PB_Z1(:,28:44),2)); %28:87
% figure;
% loglog(f_mat,sum_kpara_kperp_PB_Z1,'*','linewidth',2,'MarkerSize',8); hold on;
% 
% legend('k_{||}~[5*10^{-5},8*10^{-5}]; k_{\perp}~[5*10^{-5},8*10^{-5}]');
% xlabel('f_{rest} (Hz)'); ylabel('P_B');  legend('boxoff','Location','bestoutside')
% set(gca,'linewidth',2,'Fontsize',30,'tickDir','both');
% 
% 
% Fi2 = fit(f_mat(:), sum_kpara_kperp_PB_Z1(:),'smoothingspline');
% 
% plot(Fi2,'--',f_mat(:), sum_kpara_kperp_PB_Z1(:));
% set(gca,'linewidth',2,'Fontsize',30,'tickDir','both','xscale','log','yscale','log');
% ylim([1e-4,5e0])
% 
% 
% sum_kpara_PB_Z1_b = squeeze(nansum(df_mean_Power_MA_fpkk(:,85:135,:),2));
% sum_kpara_kperp_PB_Z3_b = squeeze(nansum(sum_kpara_PB_Z1_b(:,163:217),2));
% 
% figure;
% loglog(f_mat,sum_kpara_kperp_PB_Z3_b,'*','linewidth',2,'MarkerSize',8); hold on;
% xlabel('f_{rest} (Hz)'); ylabel('P_B');  legend('boxoff','Location','bestoutside')
% set(gca,'linewidth',2,'Fontsize',30,'tickDir','both');
% legend('k_{||}~[1.5*10^{-4},2.5*10^{-4}]; k_{\perp}~[3*10^{-4},4*10^{-4}]');
% 
% %%
% f_m     = zeros(numf,numk,numk); 
% kpara_m = zeros(numf,numk,numk);
% kperp_m = zeros(numf,numk,numk); 
% 
% fre1_mat  = repmat(f_mat',1,numk); % per kperp same
% xvecA_mat = repmat(xvecA,numf,1);  % per kperp same
% 
% for nnkperp = 1:numk
%    f_m(:,:,nnkperp) = fre1_mat;
%    kpara_m(:,:,nnkperp) = xvecA_mat;
%    kperp_m(:,:,nnkperp) = xvecA(nnkperp);
% end
% 
% f_temp     = reshape(f_m,[],1);
% kpara_temp = reshape(kpara_m,[],1);
% kperp_temp = reshape(kperp_m,[],1);
% % PB_temp    = reshape(df_Power_MA_fpkk,[],1);
% % PV_temp    = reshape(df_Power_KA_fpkk,[],1);
% PB_temp    = reshape(df_mean_Power_MA_fpkk,[],1);
% PV_temp    = reshape(df_mean_Power_KA_fpkk,[],1);
% 
% PB_temp_norm = PB_temp./max(PB_temp);
% PB_temp_norm(PB_temp_norm<1e-6) = nan;
% 
% figure;
% scatter3(kperp_temp,kpara_temp,f_temp,9,log10(PB_temp_norm),'filled','s'); hold on;
% c = colorbar;  colormap jet; %caxis(gca,[4 7]); 
% c.Label.String = 'log_{10} P_{B_{A}}(k_{\perp},k_{||},f_{rest})';
% % c.Position = [ 0.8051    0.326    0.0170    0.45];
% xlabel('k_{\perp} (km^{-1})'); ylabel('k_{||} (km^{-1})');  zlabel('f_{rest} (Hz)');
% set(gca,'linewidth',2,'Fontsize',30,'xdir','normal','YScale','log','XScale','log','ZScale','log');
% % set(gca,'linewidth',2,'Fontsize',20,'xdir','reverse','YScale','lin','XScale','lin','ZScale','lin');
% xlim([1e-5,1e-3]); ylim([1e-5,1e-3]); zlim([2e-4,0.02]);
% view(-75,10);% view(10,10);view(-30,30);
% set(gca,'xtick',[1e-5,1e-4,1e-3],'ytick',[1e-5,1e-4,1e-3])
% set(gcf,'Position',[0 0 740 580]);
% 
% [X,Y] = meshgrid(1e-6:5e-4:5e-3,1e-6:5e-4:1e-2);
% Z = Y.*106/2/pi;    
% surf(X,Y,Z,'Facecolor','k','FaceAlpha',0.15,'EdgeColor','none');
% %%   V
% figure;
% PV_temp(PV_temp<1e1) = nan;
% scatter3(kperp_temp,kpara_temp,f_temp,6,log10(PV_temp),'filled'); hold on;
% c = colorbar;caxis(gca,[3 8]);colormap jet;
% c.Label.String = 'log_{10} P_{V_{A}}(k_{\perp},f_{rest})';
% c.Position = [ 0.8051    0.2426    0.0170    0.55];
% xlabel('k_{\perp} (km^{-1})'); ylabel('k_{||} (km^{-1})');  zlabel('f (Hz)');
% set(gca,'linewidth',2,'Fontsize',20,'xdir','reverse','YScale','log','XScale','log','ZScale','log');
% xlim([5e-5,1e-3]); ylim([1e-6,1e-3]); zlim([1e-4,0.04]);
% 
% [X,Y] = meshgrid(1e-6:5e-4:5e-3,1e-6:5e-4:1e-2);
% Z = Y.*109/2/pi;    
% surf(X,Y,Z,'Facecolor','k','FaceAlpha',0.3,'EdgeColor','none');
% 
% %% B
% sum_PB = squeeze(nansum(df_mean_Power_MA_fpkk(:,:,:),2));
% sum_PB_norm = sum_PB./max(max(sum_PB));
% sum_PB_norm(sum_PB_norm<1e-4) = nan;
% figure;
% pcolor(xvecA,f_mat',log10(sum_PB_norm)); hold on;
% shading(gca,'flat'); c=colorbar; colormap('jet');ylabel(c,'log_{10} P_{B_{A}}(k_{\perp},f_{rest})');  
% % caxis(gca,[3 7.5]);
% set(gca,'linewidth',2,'Fontsize',30,'YScale','log','XScale','log','tickDir','both');
% xx1 = logspace(-4,-3);
% yyy = 1.8*xx1.^(2/3);
% yyy2 = 38*xx1;
% loglog(xx1,yyy,'--','linewidth',3,'color','k');
% loglog(xx1,yyy2,':','linewidth',3,'color','k');
% xlabel('k_{\perp} (km^{-1})');ylabel('f_{rest} (Hz)');
% xlim([1e-5,1e-3]);ylim([5e-4,3e-2]);
% set(gcf,'Position',[10 10 800 400]);
% k_para = 7e-5;  Va = 106;  dVa = 11;dk = 0;%3e-6;
% fA   = k_para*Va/2/pi;
% xxA  = logspace(-6,-4,15);
% yyyA = fA*ones(length(xxA),1);
% f_err_2 = dk^2*Va^2/(2*pi)^2 + k_para^2/(2*pi)^2*dVa^2;
% f_err = sqrt(f_err_2);
% y_err = f_err*ones(length(xxA),1);caxis(gca,[-4 -0.5]);
% errorbar(xxA,yyyA,y_err,'linewidth',3,'color','k');hold on;
% 
% k_para = 1e-4;  Va = 106;  dVa = 11;dk = 0;%3e-6;
% fA   = k_para*Va/2/pi;
% xxA  = logspace(-6,-4,15);
% yyyA = fA*ones(length(xxA),1);
% f_err_2 = dk^2*Va^2/(2*pi)^2 + k_para^2/(2*pi)^2*dVa^2;
% f_err = sqrt(f_err_2);
% y_err = f_err*ones(length(xxA),1);caxis(gca,[-4 -0.5]);
% errorbar(xxA,yyyA,y_err,'linewidth',3,'color','m');hold on;
% 
% 
% 
% colormap jet
% 
% % B part
% sum_PB = squeeze(nansum(df_mean_Power_MA_fpkk(:,34:44,:),2));
% sum_PB_norm = sum_PB./max(max(sum_PB));
% sum_PB_norm(sum_PB_norm<5e-5) = nan;
% figure;
% pcolor(xvecA,f_mat',log10(sum_PB_norm)); hold on;
% shading(gca,'flat'); c=colorbar; colormap('jet');ylabel(c,'log_{10} P_{B_{A}}(k_{\perp},f_{rest})');  
% set(gca,'linewidth',2,'Fontsize',30,'YScale','log','XScale','log','tickDir','both');
% xx1 = logspace(-4,-3);
% yyy = 1.1*xx1.^(2/3);
% yyy2 = 23*xx1;
% loglog(xx1,yyy,'--','linewidth',3,'color','k');
% loglog(xx1,yyy2,':','linewidth',3,'color','k');
% xlabel('k_{\perp} (km^{-1})');ylabel('f_{rest} (Hz)');
% xlim([1e-5,1e-3]);ylim([5e-4,3e-2]);
% set(gcf,'Position',[10 10 800 400]);caxis(gca,[-4.5 -1]);
% k_para = 7e-5;  Va = 106; dk = 3e-6; dVa = 11;
% fA   = k_para*Va/2/pi;
% xxA  = logspace(-6,-4,15);
% yyyA = fA*ones(length(xxA),1);
% f_err_2 = dk^2*Va^2/(2*pi)^2 + k_para^2/(2*pi)^2*dVa^2;
% f_err = sqrt(f_err_2);
% y_err = f_err*ones(length(xxA),1);caxis(gca,[-3 0]);
% errorbar(xxA,yyyA,y_err,'linewidth',3,'color','k');hold on;
% 
% 
% 
% % V
% sum_PV = squeeze(sum(PfpkkV_A(:,:,:),2));
% sum_PV(sum_PV<1e1) = nan;
% figure;
% pcolor(xvecA,f_mat',log10(sum_PV)); hold on;
% shading(gca,'flat'); c=colorbar; colormap('jet');caxis(gca,[4 9]);ylabel(c,'log_{10} P_{B_{A}}(k_{\perp},f_{sc})');  
% set(gca,'linewidth',2,'Fontsize',20,'YScale','log','XScale','log','tickDir','both');
% xx1 = logspace(-4,-3);
% yyy = 1.8*xx1.^(2/3);
% yyy2 = 38*xx1;
% loglog(xx1,yyy,'--','linewidth',3,'color','k');
% loglog(xx1,yyy2,':','linewidth',3,'color','k');
% xlabel('k_{\perp} (km^{-1})');ylabel('f_{rest} (Hz)');
% xlim([5e-6,1.2e-3]);ylim([6e-4,3e-2]);
% title(['V','-',trange,'min']);
% colormap hsv
% 
% % V part
% sum_PV = squeeze(sum(PfpkkV_A(:,28:39,:),2));
% sum_PV(sum_PV<1e1) = nan;
% figure;
% pcolor(xvecA,f_mat',log10(sum_PV)); hold on;
% shading(gca,'flat'); c=colorbar; colormap('jet');caxis(gca,[3 7]);ylabel(c,'log_{10} P_{B_{A}}(k_{\perp},f_{sc})');  
% set(gca,'linewidth',2,'Fontsize',20,'YScale','log','XScale','log','tickDir','both');
% xx1 = logspace(-4,-3);
% yyy = 1.2*xx1.^(2/3);
% yyy2 = 27*xx1;
% loglog(xx1,yyy,'--','linewidth',3,'color','k');
% loglog(xx1,yyy2,':','linewidth',3,'color','k');
% xlabel('k_{\perp} (km^{-1})');ylabel('f_{rest} (Hz)');
% xlim([5e-6,1.2e-3]);ylim([6e-4,1.2e-2]);
% title(['V','-',trange,'min']);
% colormap jet;
% 
% %% 2d scatter
% 
% f_2d     = repmat(f_mat',1,numk); % per kperp same
% xvecA_2d = repmat(xvecA,numf,1); % per kperp same
% 
% f_2d_tmp = reshape(f_2d,[],1);
% xvecA_2d_tmp = reshape(xvecA_2d,[],1);
% sum_PB_tmp = reshape(sum_PB,[],1);
% 
% figure; 
% scatter(xvecA_2d_tmp,f_2d_tmp,15,sum_PB_tmp,'filled');
% 
% xlim([5e-6,1e-3]); ylim([8e-4,3e-2]);
% 
% 
% 
% %% %%%  B
% % sum_PB = squeeze(sum(PfpkkB_A(:,:,109:150),3));
% % sum_PB = squeeze(sum(PfpkkB_A(:,:,87:163),3));
% sum_PB = squeeze(sum(df_mean_Power_MA_fpkk(:,:,87:164),3));
% sum_PB_norm = sum_PB./max(max(sum_PB));
% sum_PB_norm(sum_PB_norm<1e-6) = nan;
% figure;
% pcolor(xvecA,f_mat',log10(sum_PB_norm)); hold on;%caxis(gca,[3 6]);
% shading(gca,'flat'); c=colorbar; colormap('jet');ylabel(c,'log_{10} P_{B_{A}}(k_{\perp},f_{rest})');  
% set(gca,'linewidth',2,'Fontsize',30,'YScale','lin','XScale','lin','tickDir','both');
% xx1 = logspace(-5,-3);
% yyy2 = 106*xx1/2/pi;
% y_err = 11*xx1/2/pi;
% errorbar(xx1,yyy2,y_err,'linewidth',3,'color','k');hold on;
% xlabel('k_{||} (km^{-1})');ylabel('f_{rest} (Hz)');
% xlim([5e-6,5e-4]);  ylim([6e-4,1e-2]);
% set(gcf,'Position',[10 10 900 400]);
% colormap hsv;
% 
% %%
% sum_PB = squeeze(sum(PfpkkB_A(:,:,195:297),3));
% sum_PB(sum_PB<1e1) = nan;
% figure;
% pcolor(xvecA,f_mat',log10(sum_PB)); hold on;
% shading(gca,'flat'); c=colorbar; colormap('jet');caxis(gca,[3 6]);ylabel(c,'log_{10} P_{B_{A}}(k_{\perp},f_{rest})');  
% set(gca,'linewidth',2,'Fontsize',30,'YScale','lin','XScale','lin','tickDir','both');
% xx1 = logspace(-5,-3);
% yyy2 = 106*xx1/2/pi;
% y_err = 11*xx1/2/pi;
% errorbar(xx1,yyy2,y_err,'linewidth',3,'color','k');hold on;
% xlabel('k_{||} (km^{-1})');ylabel('f_{rest} (Hz)');
% xlim([5e-6,5e-4]);ylim([6e-4,1.2e-2]);
% set(gcf,'Position',[10 10 800 400]);
% 
% %% %%%
% sum_PV = squeeze(sum(PfpkkV_A(:,:,109:151),3));
% sum_PV(sum_PV<1e1) = nan;
% figure;
% pcolor(xvecA,f_mat',log10(sum_PV)); hold on;
% shading(gca,'flat'); c=colorbar; colormap('jet');caxis(gca,[4 7]);ylabel(c,'log_{10} P_{V_{A}}(k_{\perp},f_{sc})');  
% set(gca,'linewidth',2,'Fontsize',20,'YScale','lin','XScale','lin','tickDir','both');
% xx1 = logspace(-5,-3);
% yyy2 = 109*xx1/2/pi;
% loglog(xx1,yyy2,':','linewidth',3,'color','k');
% xlabel('k_{||} (km^{-1})');ylabel('f_{rest} (Hz)');
% xlim([5e-6,6e-4]);ylim([6e-4,1.2e-2]);
% title(['V','-',trange,'min']);
% colormap hsv;
% 
% pause
% % 
% % 
% % kpara_m = zeros(numf,numk,numk);
% %        kperp_m = zeros(numf,numk,numk); 
% %     for nnkperp = 1:numk
% %        kpara_m(:,:,nnkperp) = xvecA_mat;
% %        kperp_m(:,:,nnkperp) = xvecA(nnkperp);
% %     end
% %        kpara_temp = reshape(kpara_m,[],1);
% %        kperp_temp = reshape(kperp_m,[],1);
% %        fpf_temp   = reshape(fre_pf_A,[],1);
% %        PfpB_temp  = reshape(PfpkkB_A,[],1);
% %        idx = find(PfpB_temp ~= 0);
% %        kpara_temp = kpara_temp(idx);
% %        kperp_temp = kperp_temp(idx); 
% %        PfpB_temp  = PfpB_temp(idx);              
% %        fpf_temp   = fpf_temp(idx);        
% %        kpara_mat  = [kpara_mat; kpara_temp];
% %        kperp_mat  = [kperp_mat; kperp_temp];
% %        PfpB_mat   = [PfpB_mat;  PfpB_temp];
% %        fpf_mat    = [fpf_mat;   fpf_temp]; 
% clear kpara_m kperp_m fre_pf_A PfpkkB_A Power_MA freq_mean_A 
% 
% save([file_path,'3D_fpfkk_kklt',num2str(angle_threshod),'3D'],...
%     'kperp_mat','kpara_mat','fpf_mat','Va','PfpB_mat','xvecA','f_mat');
% 
% figure;
% PfpB_mat(PfpB_mat<1e1) = nan;
% scatter3(kperp_mat,kpara_mat,fpf_mat,15,log10(PfpB_mat),'filled'); hold on;
% colorbar;caxis(gca,[3 6]);colormap jet;
% xlabel('k_{\perp} (km^{-1})'); ylabel('k_{||} (km^{-1})');  zlabel('f (Hz)');
% set(gca,'linewidth',2,'Fontsize',20,'xdir','reverse');%,'YScale','log','XScale','log','ZScale','log');
% xlim([1e-5,1e-3]); ylim([1e-5,1e-3]);zlim([0.0001,0.05]);
% 
% %% Three-view Picture
% P_fpkpara = sum(PfpkkB_A,3);
% P_fpkperp = sum(PfpkkB_A,2); P_fpkperp = squeeze(P_fpkperp);
% P_kk      = sum(PfpkkB_A,1); P_kk = squeeze(P_kk);
% 
% P_fpkpara(P_fpkpara<1e1) = nan;
% P_fpkperp(P_fpkperp<1e1) = nan;
% P_kk(P_kk<1e1) = nan;
% 
% figure;
% pcolor(xvecA,xvecA,log10(P_kk)); hold on;
% shading(gca,'flat'); c=colorbar; colormap('jet'); ylabel(c,'log_{10} P_{B_{A}}(k_{\perp},k_{||})');  
% % [M,color]=contour(kperp,kpara,log10(PkkB_A_norm),1);
% % color.LineWidth = 2;  color.LineColor = 'k';  
% % caxis(gca,[-5 0]); 
% circles(0,0,0.1*klim,'facecolor','none','linewidth',3,'linestyle',':','edgecolor','k');  
% circles(0,0,0.3*klim,'facecolor','none','linewidth',3,'linestyle',':','edgecolor','k');  
% circles(0,0,0.5*klim,'facecolor','none','linewidth',3,'linestyle',':','edgecolor','k');  
% xlabel(gca,'k_{\perp} (km^{-1})');   ylabel(gca,'k_{||} (km^{-1})');
% set(gcf,'Position',[0 0 740 580]); axis(gca,[5e-6 1e-3 5e-6 1e-3]); 
% set(gca,'xscale','log','yscale','log','Fontsize',30,'tickDir','both','linewidth',2,'XTick',[1e-4,1e-3],'YTick',[1e-4,1e-3])   
% % filename_fig = [file_path,'kk_power_2D_AnormA_kklt',num2str(angle_threshod),'BA'];
% % irf_print_fig(filename_fig,'png'); saveas(gca,[filename_fig,'.fig'])
% 
% % caxis(gca,[-3 -1]);
% xx1 = logspace(-4,-3);
% yyy = xx1.^(2/3);
% plot(xx1,yyy*0.03,'--k','linewidth',2)
% 
% 
% 
% figure;
% pcolor(xvecA,f_mat',log10(P_fpkpara)); hold on;
% shading(gca,'flat'); c=colorbar; colormap('jet'); ylabel(c,'log_{10} P_{B_{A}}(k_{\perp},k_{||})');  
% xlabel(gca,'k_{||} (km^{-1})');   ylabel(gca,'f_{rest} (Hz)');
% set(gcf,'Position',[0 0 740 580]); axis(gca,[5e-6 1e-3 5e-6 1e-3]); 
% set(gca,'xscale','log','yscale','log','Fontsize',30,'tickDir','both','linewidth',2,'XTick',[1e-4,1e-3],'YTick',[1e-4,1e-3])   
% % filename_fig = [file_path,'kk_power_2D_AnormA_kklt',num2str(angle_threshod),'BA'];
% % irf_print_fig(filename_fig,'png'); saveas(gca,[filename_fig,'.fig'])
% 
% xlim([1e-5,1e-3]); ylim([7e-4,0.02]);
% 
% xx1 = logspace(-5,-3);
% yyy2 = 106*xx1/2/pi;
% y_err = 11*xx1/2/pi;
% errorbar(xx1,yyy2,y_err,'linewidth',3,'color','k');hold on;
% 
% 
% 
% 
% 
% figure;
% pcolor(xvecA,f_mat',log10(P_fpkperp)); hold on;
% shading(gca,'flat'); c=colorbar; colormap('jet'); ylabel(c,'log_{10} P_{B_{A}}(k_{\perp},k_{||})');  
% xlabel(gca,'k_{\perp} (km^{-1})');   ylabel(gca,'f_{rest} (Hz)');
% set(gcf,'Position',[0 0 740 580]); axis(gca,[5e-6 1e-3 5e-6 1e-3]); 
% set(gca,'xscale','log','yscale','log','Fontsize',30,'tickDir','both','linewidth',2,'XTick',[1e-4,1e-3],'YTick',[1e-4,1e-3])   
% % filename_fig = [file_path,'kk_power_2D_AnormA_kklt',num2str(angle_threshod),'BA'];
% % irf_print_fig(filename_fig,'png'); saveas(gca,[filename_fig,'.fig'])
% 
% xlim([5e-6,1e-3]); ylim([7e-4,0.02]);
% 
% xx1 = logspace(-4,-3);
% yyy = 1.8*xx1.^(2/3);
% yyy2 = 38*xx1;
% loglog(xx1,yyy,'--','linewidth',3,'color','k');
% loglog(xx1,yyy2,':','linewidth',3,'color','k');
% 
% %% 1D
% P_kk_perp = nansum(P_kk,1);
% P_kk_para = nansum(P_kk,2);
% 
% figure;
% loglog(xvecA,P_kk_perp,'o-','linewidth',2,'Markersize',2);
% hold on 
% loglog(xvecA,P_kk_para,':', 'linewidth',4,'Markersize',2);
% 
% xlabel('k (km^{-1})');
% ylabel({'P_{B_{A}}(k_{\perp}) & P_{B_{A}}(k_{||})','(km s^{-1})^2 Hz^{-1}'});
% set(gca,'Fontsize',30,'tickDir','both','linewidth',2);
% xlim([4.5e-5,1e-3]);set(gcf,'Position',[0 0 740 580]); 
% ylim([1e4,5e9]);
% 
% [a_B1, b_B, MSE_B, R2_B, S_B] = logfit(xvecA(163:325),P_kk_perp(163:325),'loglog','color','m','linewidth',2);
% % %% test I and P
% % figure;
% % for i=50:50:539
% %      I = 1e7^(-1/3).*10.^((-(1e7)^(1/3).*xvecA./xvecA(i).^(2/3)));
% %      kperp_line = log10(xvecA(i)).*ones(1,length(xvecA));
% %      yyaxis left
% %      scatter(xvecA,I,10,kperp_line,'filled'); hold on;     
% %      set(gca,'Fontsize',18,'tickDir','both','linewidth',2,'xscale','log','yscale','log');
% %      ylim([1e-30,1]);
% %      yyaxis right
% %      scatter(xvecA,P_kk(:,i)./max(P_kk(:,i)),25,kperp_line,'+','LineWidth',2); hold on;
% %      set(gca,'Fontsize',18,'tickDir','both','linewidth',2,'xscale','log','yscale','log');
% % %      ylim([1e-50,1]);
% % end
% % colorbar
% % colormap jet
% % 
% % 
% % 
% % figure;
% % for i=50:50:539
% %      I = 1e4^(-1/3).*exp((-(1e4)^(1/3).*xvecA(i)./xvecA.^(2/3)));
% %      kpara_line = log(xvecA(i)).*ones(1,length(xvecA));
% %      yyaxis left
% %      scatter(xvecA,I,10,kpara_line,'filled'); hold on;     
% %      set(gca,'Fontsize',18,'tickDir','both','linewidth',2,'xscale','log','yscale','log');
% %      ylim([1e-5,1]);
% %      yyaxis right
% %      scatter(xvecA,P_kk(i,:),25,kpara_line,'+','LineWidth',2); hold on;
% %      set(gca,'Fontsize',18,'tickDir','both','linewidth',2,'xscale','log','yscale','log');
% % %      ylim([1e-5,1]);
% % end
% % colorbar
% % colormap jet
% % xlim([1e-5,1e-3])