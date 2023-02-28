%% modified by sqzhao 2021-12-15
close all;  clear;

%% time list
filename  = 'cluster_list.txt';
list   = readtable(filename); 
tstart = char(list.Var1);    tstop = char(list.Var2);
lines  = size(tstart,1);

%% Load data
for i=1:1:lines
tinput   = [tstart(i,1:10),'T',tstart(i,12:19),'.00Z/',tstop(i,1:10),'T',tstop(i,12:19),'.00Z'];
Tints    = irf.tint(tinput);                      % timing interval
t_length = Tints(2) - Tints(1);
Tint1    = Tints(1).epochUnix - t_length/2;  
Tint2    = Tints(2).epochUnix + t_length/2;
Tint_t   = irf_time(Tint1,'epoch>utc');
Tint_e   = irf_time(Tint2,'epoch>utc');
Tint     = irf.tint([Tint_t,'/',Tint_e]);        % longer interval for devoid of edge effects

tstart0_str = Tints(1).utc;    tend0_str = Tints(2).utc;
data_name   = [tstart0_str(1:10),'-',tstart0_str(12:13),tstart0_str(15:16),tstart0_str(18:19),'-',tend0_str(12:13),tend0_str(15:16),tend0_str(18:19),'VB'];

%% Load data
ic = 1:4; 
%% load fgm
mag_data1 = load_mag(Tint,'c1'); mag_data2 = load_mag(Tint,'c2'); mag_data3 = load_mag(Tint,'c3'); mag_data4 = load_mag(Tint,'c4');
Bxyz1 = mag_data1.Bxyz;          Bxyz2 = mag_data2.Bxyz;          Bxyz3 = mag_data3.Bxyz;          Bxyz4 = mag_data4.Bxyz;
r1    = mag_data1.r;             r2    = mag_data2.r;             r3    = mag_data3.r;             r4    = mag_data4.r;

%% denan
Bxyz1 = denan_sq(Bxyz1);  r1 = denan_sq(r1);
Bxyz2 = denan_sq(Bxyz2);  r2 = denan_sq(r2); 
Bxyz3 = denan_sq(Bxyz3);  r3 = denan_sq(r3); 
Bxyz4 = denan_sq(Bxyz4);  r4 = denan_sq(r4);

Btime1 = Bxyz1(:,1);   Bdata1 = Bxyz1(:,2:4);   rtime1 = r1(:,1);   rdata1 = r1(:,2:4);
Btime2 = Bxyz2(:,1);   Bdata2 = Bxyz2(:,2:4);   rtime2 = r2(:,1);   rdata2 = r2(:,2:4);
Btime3 = Bxyz3(:,1);   Bdata3 = Bxyz3(:,2:4);   rtime3 = r3(:,1);   rdata3 = r3(:,2:4);   
Btime4 = Bxyz4(:,1);   Bdata4 = Bxyz4(:,2:4);   rtime4 = r4(:,1);   rdata4 = r4(:,2:4);

%% resample:  B => interp higher resolution 4 sample/s 
t_start = Btime1(1);   t_end = Btime1(end);   dt = 1/8;%8;%1/4;% 8;%4;
t_mat   = t_start:dt:t_end; 
t_mat_epoch = irf_time(t_mat,'epoch>utc');

% interp 
c_eval('Bxyz?_data_temp = interp1(Btime?,Bdata?,t_mat);',ic);
c_eval('r?_data_temp    = interp1(rtime?,rdata?,t_mat);',ic);

c_eval('Bxyz? = irf.ts_vec_xyz(t_mat_epoch,Bxyz?_data_temp);',ic);
c_eval('Rxyz? = irf.ts_vec_xyz(t_mat_epoch,r?_data_temp);'   ,ic);

%% B mean timeclip in Tint; whole window
c_eval('Bxyz?_temp = Bxyz?.tlim(Tints);',ic);
c_eval('Bxyz?_mean = nanmean(double(Bxyz?_temp.data));',ic);

clear Bxyz1_temp   Bxyz2_temp   Bxyz3_temp   Bxyz4_temp ...
      r1_data_temp r2_data_temp r3_data_temp r4_data_temp

%% plasma data    !!!!  only ues 'c1' !!!!!!!!!!!!!!!!!!!!!!! 
cis_data1   = load_cis(Tint,'c1');     cis_data3 = load_cis(Tint,'c1'); % here is 'c3'
Np_HIA1     = cis_data1.Np_HIA;          Vp_HIA3 = cis_data3.Vp_HIA;
Vp_HIA1     = cis_data1.Vp_HIA;
Tp_HIA_par  = cis_data1.Tp_HIA_par;
Tp_HIA_perp = cis_data1.Tp_HIA_perp;
T = (Tp_HIA_par(:,2)+2.*Tp_HIA_perp(:,2))/3;   % MK

Vp_HIA1 = denan_sq(Vp_HIA1);
Vp_HIA3 = denan_sq(Vp_HIA3);

Vp_temp1 = interp1(Vp_HIA1(:,1),Vp_HIA1(:,2:4),t_mat);
Vp_temp3 = interp1(Vp_HIA3(:,1),Vp_HIA3(:,2:4),t_mat);
Vp1 = irf.ts_vec_xyz(t_mat_epoch,Vp_temp1);
% Vp3 = irf.ts_vec_xyz(t_mat_epoch,Vp_temp3);
Vp3 = Vp1; 

%% mean; tclip
tts = Tints(1).epochUnix; tte = Tints(2).epochUnix;
Vp_m = irf_tlim(Vp_HIA1,tts,tte);
Np_m = irf_tlim(Np_HIA1,tts,tte);
T_m  = irf_tlim([Tp_HIA_par(:,1),T],tts,tte);
vxyz_mean = nanmean(Vp_m(:,2:4));
Np = nanmean(Np_m(:,2));
T_mean = nanmean(T_m(:,2));

clear Bdata1 Bdata2 Bdata3 Bdata4 Btime1 Btime2 Btime3 Btime4 t_mat t_mat_epoch...
      Bxyz1_data_temp Bxyz2_data_temp Bxyz3_data_temp Bxyz4_data_temp Vp_temp1 Vp_temp3...
      mag_data1 mag_data2 mag_data3 mag_data4 Vp_HIA1 Vp_HIA3 ...
      r1 r2 r3 r4 rdata1 rdata2 rdata3 rdata4 rtime1 rtime2 rtime3 rtime4 r1_data_temp r2_data_temp r3_data_temp r4_data_temp

%% parameters
gamma = 5./3;    miu0 = 4*pi*1e-7;    e_charge = 1.6*1e-19;    mp = 1.67*1e-27; % kg  
Tp   = T_mean./11604.505.*1e6;          % eV  thermal velocity Mk => Tp
Bt   = norm(Bxyz1_mean);
Pmag = (Bt.^2/(2*miu0))*1e-9; 
Pth  = Np.*e_charge.*Tp*1e6*1e9;            % Pth = Np*T   nPa
beta = Pth./Pmag;
rho  = Np*mp*1e6;                           % kg/m^3
a    = sqrt(gamma.*Pth.*1e-9./rho)*1e-3;    % km/s
V_A  = Bt./sqrt(miu0*mp*1e6.*Np)*1e-9*1e-3; % km/s
f_ci_Hz = 0.01525.*mean(Bt);
rci  = 1.02*100*sqrt(Tp)./Bt;               
di   = 228./sqrt(Np);
alpha= beta.*(gamma/2); 
norm_cc = miu0*mp*Np*1e6;

%% resample
btime = Bxyz1.time; 
c_eval('B?_m_ts = irf.ts_vec_xyz(btime(1),Bxyz?_mean);',ic);
c_eval('B?_m_ts = B?_m_ts.resample(Bxyz?);'            ,ic);

Va_ts   = irf.ts_scalar(btime(1),V_A);            Va_ts =   Va_ts.resample(Bxyz1);
a_ts    = irf.ts_scalar(btime(1),a);               a_ts =    a_ts.resample(Bxyz1);
vxyz_ts = irf.ts_vec_xyz(btime(1),vxyz_mean);   vxyz_ts = vxyz_ts.resample(Bxyz1);

%% Compute dispersion relation     % 4 sc with the same coordinate
% fk_powerspec4SC_modified_v8('Bxyz?','Rxyz?','B?_m_ts',Tints,data_name,Va_ts,a_ts,vxyz_ts,Vp1,Vp3,norm_cc);
fk_powerspec4SC_modified_v9('Bxyz?','Rxyz?','B?_m_ts',Tints,data_name,Va_ts,a_ts,vxyz_ts,Vp1,Vp3,norm_cc);
disp(['i=',num2str(i)]);
end

disp('done!'); 
load chirp  
sound(y,Fs)
% 
% % test
% ttt = Bxyz1.time.epochUnix;
% ddd = Bxyz1.data;
% E = [ttt,ddd(:,2)];
% [fre_fft,Ekx,n,Fs] = fft_zsq_2022(E);
% 
% Power = abs(Ekx.*conj(Ekx))./(n*Fs);
% 
% Power_sth = smoothdata(Power,7);
% loglog(fre_fft,Power_sth,'o-');
% 
% Wv1 = irf_wavelet(E,'returnpower',0,'nf',600,'wavelet_width',5.36*1,'f',[0.0001 0.05]);
% Power_wv = nanmean(abs(cell2mat(Wv1.p).*conj(cell2mat(Wv1.p))));
% fre_wv = Wv1.f;
% t = Wv1.t; dt = t(2)-t(1);
% 
% lll = zeros(1,600-1);
% kkk = zeros(1,600-1);
% for i =1:1:600-1
%     lll(i) = trapz(flipud(fre_wv(i:end)),flipud(Power_wv(i:end)));
%     kkk(i) = sum(flipud(transpose(fre_wv(i:end))).*flipud(Power_wv(i:end)));
% end