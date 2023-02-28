%% theta difference from k.B=0  and timing 
%% control angle<k,k'>
clear; close all;

%% output path
date = '20031203200000';%'20021216050000';%'20011219000000';%'20011216120000';
trange = '300';
file_path = ['/Users/zhaosiqi/matlab_Figure/',date,'/',trange,'min-window/'];
data_path = ['/Users/zhaosiqi/matlab_Figure/',date,'/data/nonlinear/',trange,'min-data/'];
%% parameters
numk = 400;     rci = 74;
fci  = 0.24;    angle_threshod = 25;    d_sc = 2000;
flim = 0.5*fci;  
fmin = 0.00025;
%% load in sub window
filename  = 'cluster_list.txt';
list      = readtable(filename); 
tstart    = char(list.Var1);     tstop = char(list.Var2);
lines     = size(tstart,1);
data_name = strings(lines,1);

Power_MA = [];  Power_KA = [];
kx_A =[]; ky_A =[]; kz_A =[];  kxav = []; kyav = []; kzav = [];
Bav = []; a = [];  vxyz = []; Va = []; 
%% Load data
for i=1:1:lines
  tinput = [tstart(i,1:10),'T',tstart(i,12:19),'.00Z/',tstop(i,1:10),'T',tstop(i,12:19),'.00Z'];
  Tints  = irf.tint(tinput);          % timing interval
  tstart0_str    = Tints(1).utc;      tend0_str = Tints(2).utc;
  data_name(i,:) = [tstart0_str(1:10),'-',tstart0_str(12:13),tstart0_str(15:16),tstart0_str(18:19),'-',tend0_str(12:13),tend0_str(15:16),tend0_str(18:19),'VB'];
  S = load(strcat(data_path,data_name(i)));
  kx_A = [kx_A; S.kx_A];
  ky_A = [ky_A; S.ky_A];
  kz_A = [kz_A; S.kz_A];
  kxav = [kxav; S.kxav]; 
  kyav = [kyav; S.kyav]; 
  kzav = [kzav; S.kzav]; 
  Power_MA = [Power_MA; S.Power_MA];
  Power_KA = [Power_KA; S.Power_KA];
  Bav  = [Bav; S.Bav];
  a    = [a;   S.a];
  Va   = [Va;  S.Va];
  vxyz = [vxyz;S.vxyz];
end
numf0 = S.numf;  fre1 = S.fre;  t_length = str2num(trange)*60;     

%%
Bavxmat = Bav(:,1)*ones(1,numf0);
Bavymat = Bav(:,2)*ones(1,numf0);
Bavzmat = Bav(:,3)*ones(1,numf0);
Bavabsmat = sqrt(Bav(:,1).^2+Bav(:,2).^2+Bav(:,3).^2)*ones(1,numf0);

%% Doppler shift
[k0_A_pf,kx_A_pf,ky_A_pf,kz_A_pf,freq_mean_A,freq_A_A,freq_f_A,freq_s_A,vp_mean_A,theta_A,f_doppler_A] = doppler_shift(kx_A,ky_A,kz_A,fre1,Bav,vxyz,Va,a); 

kmag_A = sqrt(kx_A.*kx_A + ky_A.*ky_A + kz_A.*kz_A); 

%%  MHD threthod
f_2_t = round(4./t_length,5);
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

obtuse_theta_A = dtheta_A > 90;   % theta < 90
dtheta_A(obtuse_theta_A)  = 180 - dtheta_A(obtuse_theta_A);

%% f threshod
dtheta_A(freq_mean_A>flim)   = 0; 
dtheta_A(freq_mean_A<f_2_t)  = 0; 
dtheta_A(isnan(freq_mean_A)) = 90; % freq_mean is nan delete 

%% deNAN 
vp_mean_A(freq_mean_A>flim)     = 0; 
freq_mean_A(freq_mean_A>flim)   = 0; 
freq_mean_A(isnan(freq_mean_A)) = 0;

%% k threshod: k rho_ci < 0.1
klim = 0.1/rci;
Power_MA(kmag_A   > klim) = 0; % threshod |k|>klim, Power = 0
Power_KA(kmag_A   > klim) = 0; 
kpar_A_abs(kmag_A > klim) = 0; % all k(|k|>klim) = 0 
kperp_A(kmag_A    > klim) = 0;
dtheta_A(kmag_A   > klim) = 0; % |k|>klim, theta = 0  
kmin = 1/10/d_sc; % /2 sqzhao  100/d_sc
Power_MA(kmag_A   < kmin) = 0;
Power_KA(kmag_A   < kmin) = 0; 
kmag_A(kmag_A     > klim) = 0; % add 20220710

% use same length scale
% kmax_A = klim*1.1; 
% kmin_A = -kmax_A;
% kmagvec_A = linspace(kmin,kmax_A,numk);
% dkmag_A = kmax_A/numk;

%% theta kk' threshod
%   angle_threshod = 360;      % if without angle threshod
% angle_threshod = angle_threshod;
ind_dtheta_gt  = dtheta_A > angle_threshod;   % phi_kk' 
Power_MA(ind_dtheta_gt) = 0;   % angle threshod, Power = 0
Power_KA(ind_dtheta_gt) = 0; 

%% normalize to time points
numt = size(Power_MA,1);
fmat = repmat(fre1,1,numt);  fmat = fmat';

% for M
idx_M = Power_MA~=0;
numt_data_M = nansum(idx_M,1);  numt_data_M(numt_data_M==0) = 1;
numt_mat_M  = repmat(numt_data_M',1,numt); numt_mat_M = numt_mat_M';
mean_Power_MA = Power_MA./numt_mat_M;

% for K
idx_K = Power_KA~=0;
numt_data_K = nansum(idx_K,1);  numt_data_K(numt_data_K==0) = 1;
numt_mat_K  = repmat(numt_data_K',1,numt); numt_mat_K = numt_mat_K';
mean_Power_KA = Power_KA./numt_mat_K;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kmax_A = klim*1.1; 
amin=log10(1e-6); 
amax=log10(kmax_A); 
a = logspace(amin,amax,numk);
xvecA = a;

diff_a = abs(diff(a)); diff_a = [diff_a,0];
%% %%%%%%%%%%%%%%%%%%%    3D    %%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fsc,|kperp|,|kpar|  sqzhao  spacecraft frame 
irf.log('notice','Computing power versus |kperp|,|kpar|')
PfkkB_A      = zeros(numf0,numk,numk);         
PfkkV_A      = zeros(numf0,numk,numk);
PfkkB_mean_A = zeros(numf0,numk,numk);    
PfkkV_mean_A = zeros(numf0,numk,numk);

for mm = 1:numt
	for nn = 1:numf0
      [kpara_min,idx_kpara] = min(abs(kpar_A_abs(mm,nn) - a)); 
      [kperp_min,idx_kperp] = min(abs(kperp_A(mm,nn)    - a));  
      PfkkB_A(nn,idx_kpara,idx_kperp) = PfkkB_A(nn,idx_kpara,idx_kperp) + Power_MA(mm,nn)./diff_a(idx_kpara)./diff_a(idx_kperp);
      PfkkV_A(nn,idx_kpara,idx_kperp) = PfkkV_A(nn,idx_kpara,idx_kperp) + Power_KA(mm,nn)./diff_a(idx_kpara)./diff_a(idx_kperp);
      PfkkB_mean_A(nn,idx_kpara,idx_kperp) = PfkkB_mean_A(nn,idx_kpara,idx_kperp) + mean_Power_MA(mm,nn)./diff_a(idx_kpara)./diff_a(idx_kperp);
      PfkkV_mean_A(nn,idx_kpara,idx_kperp) = PfkkV_mean_A(nn,idx_kpara,idx_kperp) + mean_Power_KA(mm,nn)./diff_a(idx_kpara)./diff_a(idx_kperp);
    end
end

%%
idx_perp = xvecA >= 5e-5;
idx_para = xvecA >= 9e-6;
kpara = xvecA(idx_para);
kperp = xvecA(idx_perp);

PfkkB_mean_A = PfkkB_mean_A(:,idx_para,idx_perp);

% integrate
% Pkk = zeros(length(kpara),length(kperp));
% for  i=1:1:length(kpara)-1
%     for j = 1:1:length(kperp)-1
%        Pkk(i,j) = trapz(flipud(fre1),flipud(PfkkB_mean_A(:,i,j)));
%     end
% end

% sum
diff_f = abs(diff(fre1)); diff_f = [diff_f;0];
Pkk = zeros(length(kpara),length(kperp));
for  i=1:1:length(kpara)-1
    for j = 1:1:length(kperp)-1
       Pkk(i,j) = nansum(PfkkB_mean_A(:,i,j).*diff_f);
    end
end


%%
L0_para = 4.6e4;%7.6e4; % dvperp * Tcor(B)
L0_perp = 1.1e5;
MA = 0.34;

kpara_mat = repmat(kpara',1,length(kperp)); 
kperp_mat = repmat(kperp, length(kpara),1);  

I = zeros(length(kpara),length(kperp)); 
for i=1:1:length(kpara)
    kperp_norm = kperp*L0_perp;
    kpara_norm = kpara(i)*L0_para;
    A = exp(-1*L0_para^(1/3)*abs(kpara(i))./(MA^(4/3).*kperp.^(2/3)));
%     I(i,:) = A;
    I(i,:) = kperp.^(-7/3).*A;
%       I(i,:) = kperp.^(-7/3).*exp(-1*abs(kpara(i))./(MA^(4/3)));
end

%% Fig. 2a
Pkk_norm = Pkk./max(max(Pkk));
Pkk_norm(Pkk_norm<1e-5) = nan;
figure;
pcolor(kperp,kpara,log10(Pkk_norm)); hold on;
shading(gca,'flat'); colormap('jet'); c=colorbar;% caxis(gca,[0 1]);
c.Label.String = {'log_{10}P_{A}(k_{\perp},k_{||}) & log_{10}I_{A}(k_{\perp},k_{||})'};
set(gca,'Fontsize',35,'tickDir','both','linewidth',2,'xscale','log','yscale','log');
xlim([5e-5,1.2e-3]); ylim([1e-5,1.2e-3]);
xlabel('k_{\perp} (km^{-1})');    ylabel('k_{||} (km^{-1})'); 

I_norm = I./max(max(I));
I_norm(I_norm<1e-5) = nan;
[M,color] = contour(kperp,kpara,log10(I_norm),10);
color.LineWidth = 5;   % color.LineColor = 'k'; 
[M,color2] = contour(kperp,kpara,log10(I_norm),10);
color2.LineStyle = '--';   color2.LineColor = 'k'; color2.LineWidth = 3;
circles(0,0,1.4e-04,'facecolor','none','Linewidth',3,'linestyle',':','edgecolor','k');  
circles(0,0,3*1.4e-04,'facecolor','none','Linewidth',3,'linestyle',':','edgecolor','k');  
xx1  = 2e-4:1e-5:1.2e-3;
yyy2 = 0.03*xx1.^(2/3);
yyy1 = 0.52*xx1;
loglog(xx1,yyy2,':m','Linewidth',5);
loglog(xx1,yyy1,':k','Linewidth',5);
set(gcf,'Position',[10 10 800 650]);
caxis(gca,[-6 0]);

%% contour Fig. 2b
figure;
[MP,colorP] = contourf(kperp,kpara,log10(Pkk_norm),10); hold on;
colorP.LineWidth = 0.1;    %color.LineColor = ' '; 
set(gca,'Fontsize',35,'tickDir','both','linewidth',2,'xscale','log','yscale','log');
xlim([5e-5,1.2e-3]); ylim([1e-5,1.2e-3]);
xlabel('k_{\perp} (km^{-1})');    ylabel('k_{||} (km^{-1})'); 

I_norm = I./max(max(I));
I_norm(I_norm<1e-6) = nan;
[M,color] = contour(kperp,kpara,log10(I_norm),10);
color.LineWidth = 5;   % color.LineColor = 'k'; 
[M,color2] = contour(kperp,kpara,log10(I_norm),10);
color2.LineStyle = '--';   color2.LineColor = 'k'; color2.LineWidth = 3;
circles(0,0,1.4e-04,'facecolor','none','Linewidth',3,'linestyle',':','edgecolor','k');  
circles(0,0,3*1.4e-04,'facecolor','none','Linewidth',3,'linestyle',':','edgecolor','k');  
xx1  = 2e-4:1e-5:1.2e-3;
yyy2 = 0.03*xx1.^(2/3);
yyy1 = 0.52*xx1;
loglog(xx1,yyy2,':m','Linewidth',5);
loglog(xx1,yyy1,':b','Linewidth',5);
set(gcf,'Position',[10 10 800 650]);
c=colorbar;colormap hsv; caxis(gca,[-6 -1]);
c.Label.String = {'log_{10}P_{A}(k_{\perp},k_{||}) & log_{10}I_{A}(k_{\perp},k_{||})'};


% smooth
% Pkk_norm_filt = imgaussfilt(Pkk_norm,0.01);
Pkk_norm_filt = Pkk_norm;
figure;
[MP,colorP_filt] = contourf(kperp,kpara,log10(Pkk_norm_filt),10); hold on;
% colorP_filt.LineWidth = 0.1;   % color.LineColor = 'k'; 
set(gca,'Fontsize',35,'tickDir','both','linewidth',2,'xscale','log','yscale','log');
xlim([5e-5,1.2e-3]); ylim([1e-5,1.2e-3]);
xlabel('k_{\perp} (km^{-1})');    ylabel('k_{||} (km^{-1})'); 

I_norm = I./max(max(I));
I_norm(I_norm<1e-5) = nan;
[M,color] = contour(kperp,kpara,log10(I_norm.*3),10);
color.LineWidth = 5;   % color.LineColor = 'k'; 
[M,color2] = contour(kperp,kpara,log10(I_norm),10);
color2.LineStyle = '--';   color2.LineColor = 'k'; color2.LineWidth = 3;
circles(0,0,1.4e-04,'facecolor','none','Linewidth',3,'linestyle',':','edgecolor','k');  
circles(0,0,3*1.4e-04,'facecolor','none','Linewidth',3,'linestyle',':','edgecolor','k');  
xx1  = 2e-4:1e-5:1.2e-3;
yyy2 = 0.03*xx1.^(2/3);
yyy1 = 0.52*xx1;
loglog(xx1,yyy2,':m','Linewidth',5);
loglog(xx1,yyy1,':b','Linewidth',5);
set(gcf,'Position',[10 10 800 650]);
c=colorbar;colormap jet; caxis(gca,[-6 0]);
c.Label.String = {'log_{10}P_{A}(k_{\perp},k_{||}) & log_{10}I_{A}(k_{\perp},k_{||})'};

