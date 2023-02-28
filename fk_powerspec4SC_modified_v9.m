% all sc have the same coordinate at each time and fre
% including V and B projection   % nonlinear f
% we decrease the svd resolution to speed up the program. And average in
% time
%% corresponding to global phase list code
function [] = fk_powerspec4SC_modified_v9(varargin)

fmin = 0.00025;  fmax = 0.12;% f_ci;  fmax = 0.5fci
date = ['20031203','200000'];
trange = '300'; % min
numf = 400;
dt_svd = 4;
dt   = 256;%;64; % average over dt? seconds
res  = 8;%1/8;%1;%1./4;  %1/dt_B ; % Hz

file_path = ['/Users/zhaosiqi/matlab_Figure/',date,'/data/nonlinear/',trange,'min-data/'];

% first time clip then project & timing function sqzhao
ic = 1:4;
c_eval('E?=evalin(''base'',irf_ssub(varargin{1},?));',ic);  % fluctuating manetic field
c_eval('R?=evalin(''base'',irf_ssub(varargin{2},?));',ic);
c_eval('B?=evalin(''base'',irf_ssub(varargin{3},?));',ic);  % mean magnetic field
Tint      = varargin{4};
data_name = varargin{5};
Va_ts     = varargin{6};
a_ts      = varargin{7};
vxyz_ts   = varargin{8};
Vp1       = varargin{9};
Vp3       = varargin{10};
norm_cc   = varargin{11};

Bav = irf.ts_vec_xyz(B1.time,(B1.data+B2.data+B3.data+B4.data)/4);
%% wavelet fluctuating magnetic field  %% nonlinear map
W1 = irf_wavelet(E1,'returnpower',0,'cutedge',0,'nf',numf,'wavelet_width',5.36*1,'f',[fmin fmax]);
W2 = irf_wavelet(E2,'returnpower',0,'cutedge',0,'nf',numf,'wavelet_width',5.36*1,'f',[fmin fmax]);
W3 = irf_wavelet(E3,'returnpower',0,'cutedge',0,'nf',numf,'wavelet_width',5.36*1,'f',[fmin fmax]);
W4 = irf_wavelet(E4,'returnpower',0,'cutedge',0,'nf',numf,'wavelet_width',5.36*1,'f',[fmin fmax]);
W1_mat_x = transpose(cell2mat(W1.p(1))); W2_mat_x = transpose(cell2mat(W2.p(1)));
W1_mat_y = transpose(cell2mat(W1.p(2))); W2_mat_y = transpose(cell2mat(W2.p(2)));
W1_mat_z = transpose(cell2mat(W1.p(3))); W2_mat_z = transpose(cell2mat(W2.p(3)));
W3_mat_x = transpose(cell2mat(W3.p(1))); W4_mat_x = transpose(cell2mat(W4.p(1)));
W3_mat_y = transpose(cell2mat(W3.p(2))); W4_mat_y = transpose(cell2mat(W4.p(2)));
W3_mat_z = transpose(cell2mat(W3.p(3))); W4_mat_z = transpose(cell2mat(W4.p(3)));

fre = W1.f;
clear W1 W2 W3 W4 B1 B2 B3 B4
%% wavelet fluctuating velocity  
Wv1 = irf_wavelet(Vp1,'returnpower',0,'cutedge',0,'nf',numf,'wavelet_width',5.36*1,'f',[fmin fmax]);
Wv3 = irf_wavelet(Vp3,'returnpower',0,'cutedge',0,'nf',numf,'wavelet_width',5.36*1,'f',[fmin fmax]);

Wv1_mat_x = transpose(cell2mat(Wv1.p(1)));
Wv1_mat_y = transpose(cell2mat(Wv1.p(2)));
Wv1_mat_z = transpose(cell2mat(Wv1.p(3)));
Wv3_mat_x = transpose(cell2mat(Wv3.p(1)));
Wv3_mat_y = transpose(cell2mat(Wv3.p(2)));
Wv3_mat_z = transpose(cell2mat(Wv3.p(3)));

clear  Wv1 Wv3
%% interp to lower resolution dt
Etime = E1.time.epochUnix;
c_eval('Edata? = double(E?.data);',ic);

t_start = Etime(1);   t_end = Etime(end);    
t_mat   = t_start:dt_svd:t_end;    %  low res to perform svd
t_mat_epoch = irf_time(t_mat,'epoch>utc');  t_mat = t_mat';   
c_eval('E?_data = interp1(Etime,Edata?,t_mat);',ic);
E1_ts = irf.ts_vec_xyz(t_mat_epoch,E1_data);

clear Edata1 Edata2 Edata3 Edata4
%% get k - spectra
[fre1,WaveVector_B1] = SVD_B_2022([t_mat,E1_data],fmin,fmax,numf);
[fre2,WaveVector_B2] = SVD_B_2022([t_mat,E2_data],fmin,fmax,numf);
[fre3,WaveVector_B3] = SVD_B_2022([t_mat,E3_data],fmin,fmax,numf);
[fre4,WaveVector_B4] = SVD_B_2022([t_mat,E4_data],fmin,fmax,numf);

c_eval('kx_? = squeeze(WaveVector_B?(:,:,1));',ic);
c_eval('ky_? = squeeze(WaveVector_B?(:,:,2));',ic);
c_eval('kz_? = squeeze(WaveVector_B?(:,:,3));',ic);

clear WaveVector_B1 WaveVector_B2 WaveVector_B3 WaveVector_B4 ...
      fre1 fre2 fre3 fre4 E1_data E2_data E3_data E4_data

%% time clip; cutoff edge effects
time_h = E1.time;      idx_h = tlim(E1.time,Tint);  % tclip 
idx_svd = tlim(E1_ts.time,Tint);

% If odd, remove last data point (as is done in irf_wavelet)
if mod(length(idx_h),2)
    idx_h(end)=[];
end
L_h = length(idx_h);    %times_h = time_h(idx_h);

if mod(length(idx_svd),2)
    idx_svd(end)=[];
end
L_s = length(idx_svd);   %times_h = time_h(idx_h);

%%
c_eval('W?_mat_x = W?_mat_x(:,idx_h);',ic);  % dB
c_eval('W?_mat_y = W?_mat_y(:,idx_h);',ic);
c_eval('W?_mat_z = W?_mat_z(:,idx_h);',ic);
Wv1_mat_x = Wv1_mat_x(:,idx_h);   % dV1
Wv1_mat_y = Wv1_mat_y(:,idx_h);
Wv1_mat_z = Wv1_mat_z(:,idx_h);
Wv3_mat_x = Wv3_mat_x(:,idx_h);   % dV3
Wv3_mat_y = Wv3_mat_y(:,idx_h);
Wv3_mat_z = Wv3_mat_z(:,idx_h);
c_eval('kx_? = kx_?(:,idx_svd);',ic);
c_eval('ky_? = ky_?(:,idx_svd);',ic);
c_eval('kz_? = kz_?(:,idx_svd);',ic);

kxav = (kx_1 + kx_2 + kx_3 + kx_4)/4;
kyav = (ky_1 + ky_2 + ky_3 + ky_4)/4;
kzav = (kz_1 + kz_2 + kz_3 + kz_4)/4;

clear kx_1 kx_2 kx_3 kx_4 ky_1 ky_2 ky_3 ky_4 kz_1 kz_2 kz_3 kz_4 ...
      E1 E2 E3 E4

%% projection k b0
kxav  = kxav';    kyav = kyav';   kzav = kzav';
cav_s = dt/dt_svd;   % average in dt s; points = dt / dt_svd;
N_s   = floor(L_s/cav_s)-1;

kxav  = Power_average(kxav,numf,cav_s,N_s);   kxav = kxav';
kyav  = Power_average(kyav,numf,cav_s,N_s);   kyav = kyav';
kzav  = Power_average(kzav,numf,cav_s,N_s);   kzav = kzav';
 
cav_h = res*dt;      % average in dt s;
%% project
c_eval('dB_A_?    = zeros(numf,cav_h);',ic);
c_eval('dB_para_? = zeros(numf,cav_h);',ic);
c_eval('dB_perp_? = zeros(numf,cav_h);',ic);
dV1_A    = zeros(numf,cav_h);       dV3_A    = zeros(numf,cav_h); 
dV1_para = zeros(numf,cav_h);       dV3_para = zeros(numf,cav_h); 
dV1_perp = zeros(numf,cav_h);       dV3_perp = zeros(numf,cav_h);  

c_eval('W_A?    = [];',ic);
c_eval('W_para? = [];',ic);
c_eval('W_perp? = [];',ic);
Wv_A1 = [];   Wv_para1 = [];   Wv_perp1 = [];   
Wv_A3 = [];   Wv_para3 = [];   Wv_perp3 = [];
B0 = mean((Bav.data),1); 
%% 
for j = 1:N_s+1
    for i = 1:numf
        kxyz = [kxav(i,j),kyav(i,j),kzav(i,j)];
        kxyz_unit = kxyz./norm(kxyz);
        B_bg_unit = B0./norm(B0);
    
        %% axis 
        kpara = B_bg_unit; 
        kperp = cross(kpara,cross(kxyz_unit,kpara))./norm(cross(kpara,cross(kxyz_unit,kpara)));%bx(kxb)
        kesai_A = cross(kperp,kpara);
       
        % project to kB axes
        %% ms: mat start point      me: mat end point
        ms = (j-1)*res*dt + 1;      me = j*res*dt;
        % dB
        c_eval('dB_k_? = [W?_mat_x(i,ms:me);W?_mat_y(i,ms:me);W?_mat_z(i,ms:me)];',ic);
        c_eval('dB_A_?(i,:)    = dB_k_?(1,:)*kesai_A(1) + dB_k_?(2,:)*kesai_A(2) + dB_k_?(3,:)*kesai_A(3);',ic);
        c_eval('dB_para_?(i,:) = dB_k_?(1,:)*kpara(1)   + dB_k_?(2,:)*kpara(2)   + dB_k_?(3,:)*kpara(3);'  ,ic);
        c_eval('dB_perp_?(i,:) = dB_k_?(1,:)*kperp(1)   + dB_k_?(2,:)*kperp(2)   + dB_k_?(3,:)*kperp(3);'  ,ic);             
        % dV
        dV1_k = [Wv1_mat_x(i,ms:me);Wv1_mat_y(i,ms:me);Wv1_mat_z(i,ms:me)];
        dV3_k = [Wv3_mat_x(i,ms:me);Wv3_mat_y(i,ms:me);Wv3_mat_z(i,ms:me)];
        dV1_A(i,:)    = dV1_k(1,:)*kesai_A(1) + dV1_k(2,:)*kesai_A(2) + dV1_k(3,:)*kesai_A(3);
        dV1_para(i,:) = dV1_k(1,:)*kpara(1)   + dV1_k(2,:)*kpara(2)   + dV1_k(3,:)*kpara(3);
        dV1_perp(i,:) = dV1_k(1,:)*kperp(1)   + dV1_k(2,:)*kperp(2)   + dV1_k(3,:)*kperp(3);         
        dV3_A(i,:)    = dV3_k(1,:)*kesai_A(1) + dV3_k(2,:)*kesai_A(2) + dV3_k(3,:)*kesai_A(3);
        dV3_para(i,:) = dV3_k(1,:)*kpara(1)   + dV3_k(2,:)*kpara(2)   + dV3_k(3,:)*kpara(3);
        dV3_perp(i,:) = dV3_k(1,:)*kperp(1)   + dV3_k(2,:)*kperp(2)   + dV3_k(3,:)*kperp(3);         
    end
        c_eval('W_A?    = [W_A?,   dB_A_?];'   ,ic);
        c_eval('W_para? = [W_para?,dB_para_?];',ic);
        c_eval('W_perp? = [W_perp?,dB_perp_?];',ic);
        Wv_A1    = [Wv_A1,   dV1_A];
        Wv_A3    = [Wv_A3,   dV3_A];
        Wv_para1 = [Wv_para1,dV1_para];
        Wv_para3 = [Wv_para3,dV3_para];
        Wv_perp1 = [Wv_perp1,dV1_perp];
        Wv_perp3 = [Wv_perp3,dV3_perp];
end
L_h  = size(Wv_A1,2);

clear dB_A_1 dB_A_2 dB_A_3 dB_A_4 dB_para_1 dB_para_2 dB_para_3 dB_para_4 dB_perp_1 dB_perp_2 dB_perp_3 dB_perp_4 ...
      dV1_A  dV3_A  dV1_para dV3_para dV1_perp dV3_perp
%% %%%%%%%%%%%%%%%    timing   %%%%%%%%%%%%%%%
c_eval('W_A?    = transpose(W_A?);'   ,ic);
c_eval('W_para? = transpose(W_para?);',ic);
c_eval('W_perp? = transpose(W_perp?);',ic);
Wv_A1    = transpose(Wv_A1);
Wv_para1 = transpose(Wv_para1);
Wv_perp1 = transpose(Wv_perp1);
Wv_A3    = transpose(Wv_A3);
Wv_para3 = transpose(Wv_para3);
Wv_perp3 = transpose(Wv_perp3);

clear W1_mat_x W1_mat_y W1_mat_z W2_mat_x W2_mat_y W2_mat_z ...
      W3_mat_x W3_mat_y W3_mat_z W4_mat_x W4_mat_y W4_mat_z ...
      Wv1_mat_x Wv1_mat_y Wv1_mat_z Wv3_mat_x Wv3_mat_y Wv3_mat_z

%% Power
fkPower_A    = 0.25*(W_A1.*conj(W_A1) + W_A2.*conj(W_A2) ...
                   + W_A3.*conj(W_A3) + W_A4.*conj(W_A4));
fkPower_para = 0.25*(W_para1.*conj(W_para1) + W_para2.*conj(W_para2) ...
                   + W_para3.*conj(W_para3) + W_para4.*conj(W_para4));
fkPower_perp = 0.25*(W_perp1.*conj(W_perp1) + W_perp2.*conj(W_perp2) ...
                   + W_perp3.*conj(W_perp3) + W_perp4.*conj(W_perp4));

%% tclip data
R1 = R1.data;      R2 = R2.data;      R3 = R3.data;      R4 = R4.data;
R1 = R1(idx_h,:);  R2 = R2(idx_h,:);  R3 = R3(idx_h,:);  R4 = R4(idx_h,:);
R1 = R1(1:L_h,:);  R2 = R2(1:L_h,:);  R3 = R3(1:L_h,:);  R4 = R4(1:L_h,:);

%% timing 
N_h = floor(L_h/cav_h)-1; % L_h is new

[kx_A,ky_A,kz_A,Powerav_A] = timing_2022_v2(fkPower_A,numf,W_A1,W_A2,W_A3,W_A4,R1,R2,R3,R4,fre,cav_h,N_h);
clear W_A2 W_A3 W_A4 
[kx_para,ky_para,kz_para,Powerav_para] = timing_2022_v2(fkPower_para,numf,W_para1,W_para2,W_para3,W_para4,R1,R2,R3,R4,fre,cav_h,N_h);
clear W_para1 W_para2 W_para3 W_para4 
[kx_perp,ky_perp,kz_perp,Powerav_perp] = timing_2022_v2(fkPower_perp,numf,W_perp1,W_perp2,W_perp3,W_perp4,R1,R2,R3,R4,fre,cav_h,N_h);
clear W_perp1 W_perp2 W_perp3 W_perp4 ...
      R1 R2 R3 R4

N_k  = size(kx_A,1);
kxav = kxav(:,1:N_k);    kxav = kxav'; 
kyav = kyav(:,1:N_k);    kyav = kyav'; 
kzav = kzav(:,1:N_k);    kzav = kzav';
    
% PM: unit (km/s)^2/Hz
Power_MA    = Powerav_A./norm_cc*1e-18/1e6;     %fkPower_A;  %(km/s)^2/Hz  
Power_Mpara = Powerav_para./norm_cc*1e-18/1e6;  %fkPower_para;
Power_Mperp = Powerav_perp./norm_cc*1e-18/1e6;  %fkPower_perp;%

%%
Power_KA    = 0.5*(Wv_A1.*conj(Wv_A1) + Wv_A3.*conj(Wv_A3));
Power_Kpara = 0.5*(Wv_para1.*conj(Wv_para1) + Wv_para3.*conj(Wv_para3));
Power_Kperp = 0.5*(Wv_perp1.*conj(Wv_perp1) + Wv_perp3.*conj(Wv_perp3));

Power_KA    = Power_average(Power_KA,   numf,cav_h,N_h);
Power_Kpara = Power_average(Power_Kpara,numf,cav_h,N_h);
Power_Kperp = Power_average(Power_Kperp,numf,cav_h,N_h);

%% normalized cross-helicity
% W_A1_new = W_A1./sqrt(norm_cc)*1e-12; 
% PB1 = W_A1_new.*conj(W_A1_new);
% PV1 = Wv_A1.*conj(Wv_A1);
% % sigma_c = 2*W_A1_new.*conj(Wv_A1)./(PB1 + PV1);
% sigma_r = (PV1 - PB1)./(PB1 + PV1);
% r_A = PV1./PB1;

clear Wv_A1 Wv_para1 Wv_perp1 Wv_A3 Wv_para3 Wv_perp3
%%
len = size(kx_A,1);
Va = Va_ts.data;   a = a_ts.data;   vxyz = vxyz_ts.data;   Bav = Bav.data;
Va = Va(1:len,:);  a = a(1:len,:);  vxyz = vxyz(1:len,:);  Bav = Bav(1:len,:);

%% output data
save([file_path,data_name],...
    'kx_A','ky_A','kz_A','kx_para','ky_para','kz_para','kx_perp','ky_perp','kz_perp',...
    'kxav','kyav','kzav',...
    'Powerav_A','Powerav_para','Powerav_perp',...
    'Power_MA','Power_Mpara','Power_Mperp',...
    'Power_KA','Power_Kpara','Power_Kperp',...
    'numf','fre','Bav','Va','a','vxyz','res');

clear 'kx_A' 'ky_A' 'kz_A' 'kx_para' 'ky_para' 'kz_para' 'kx_perp' 'ky_perp' 'kz_perp' ...
      'kxav' 'kyav' 'kzav'...
      'Powerav_A' 'Powerav_para' 'Powerav_perp' ...
      'Power_MA'  'Power_Mpara'  'Power_Mperp' ...
      'Power_KA'  'Power_Kpara'  'Power_Kperp' ...
      'numf' 'fre' 'Bav' 'Va' 'a' 'vxyz' 'res'
end
