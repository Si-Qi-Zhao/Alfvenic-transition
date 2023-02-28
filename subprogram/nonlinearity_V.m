function [P_perp,P_para,kperp,kpara,kpara_new,Ekperp,Ekperp_pre2,XiB_new,vk_new,XiB_mat_sc,vk2_sc] = nonlinearity_V(date,trange,numk,angle_threshod,B0)

file_path = ['/Users/zhaosiqi/matlab_Figure/',date,'/',trange,'min-window/'];

%% loadin perp
load([file_path,'Pkk_2D_numk',num2str(numk),'_kklt',num2str(angle_threshod),'_VB_nonlinear'],...
      'P_mean_fkperp_VA','PfkkV_mean_A','xvecA','fre1','P_mean_fkpara_VA');

idx_perp = xvecA >= 4e-5;
idx_para = xvecA >= 5e-6;
kpara = xvecA(idx_para);
kperp = xvecA(idx_perp);

%% test
diff_f = abs(diff(fre1)); diff_f = [diff_f;0];
dfmat  = repmat(diff_f,1,length(xvecA));  
Pperp  = nansum(P_mean_fkperp_VA.*dfmat,1);
Ppara  = nansum(P_mean_fkpara_VA.*dfmat,1);

P_mean_fkpara_VA = P_mean_fkpara_VA(:,idx_perp);        % sc frame
PfkkV_mean_A     = PfkkV_mean_A(:,idx_para,idx_perp);   % sc frame
Pperp = Pperp(idx_perp);
Ppara = Ppara(idx_para);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx_para = 1:400;
% P_para = sum(Pkk,2);
% P_perp = sum(Pkk,1);
P_para = Ppara;
P_perp = Pperp;

kpara_new = interp1(P_para(idx_para),kpara(idx_para),P_perp);
%% sc frame
P_sumk_sc = zeros(length(fre1),length(kpara),length(kperp));
for i = 1:1:length(fre1)
  tmp = squeeze(PfkkV_mean_A(i,:,:));
  tmp_para = cumsum(tmp,2,'reverse');
  tmp_all  = cumsum(tmp_para,1,'reverse');
  P_sumk_sc(i,:,:) = tmp_all;   
end

vk2_sc = zeros(length(kpara),length(kperp));
for  i=1:1:length(kpara)-1
    for j = 1:1:length(kperp)-1
       vk2_sc(i,j) = 2*trapz(flipud(fre1),flipud(P_sumk_sc(:,i,j)));
    end
end
vk_sc = sqrt(vk2_sc); 

XiB_mat_sc = zeros(length(kpara),length(kperp));
for i=1:1:length(kpara)-1
    XiB_mat_sc(i,:) = kperp.*vk_sc(i,:)./kpara(i)./B0;
end

%%   %%%%%%%%%%%%%%%%
P_sumk_perp = zeros(length(fre1),length(kperp));
for i=1:1:length(kperp)-1
     P_sumk_perp(:,i) = sum(P_mean_fkperp_VA(:,i:end),2);
end

% int
% vk_new2 = 2*trapz(flipud(fre1),flipud(P_sumk_perp));  % during ourputing P, we do not time 0.5

% test sum
dfmat2  = repmat(diff_f,1,length(kperp));  
vk_new2 = 2*nansum(dfmat2.*P_sumk_perp,1);


vk_new  = sqrt(vk_new2);

Ekperp = 0.5*vk_new2./kperp;
Ekperp_pre2 = kperp.^(5/3).*Ekperp;

XiB_new = kperp.*vk_new./kpara_new./B0;
end

% %% test
% Pkk = zeros(length(kpara),length(kperp));
% for  i=1:1:length(kpara)-1
%     for j = 1:1:length(kperp)-1
%        Pkk(i,j) = trapz(flipud(fre1),flipud(PfkkV_A(:,i,j)));
%     end
% end

% load([file_path,'Pkk_2D_numk',num2str(numk),'_kklt',num2str(angle_threshod),'_VB_nonlinear'],...
%     'countk_A',...
%     'PkkB_A','PkkV_A','PkkB_mean_A','PkkV_mean_A','fPkkB_mean_A','fPkkV_mean_A',...
%     'xvecA','yvecA','klim','fre1',...
%     'PkperpB_A','PkperpV_A','Pkperp_mean_B_A','Pkperp_mean_V_A','fPkperp_mean_B_A','fPkperp_mean_V_A',...
%     'PkparaB_A','PkparaV_A','Pkpara_mean_B_A','Pkpara_mean_V_A','fPkpara_mean_B_A','fPkpara_mean_V_A',...
%     'PkabsB_A' ,'PkabsV_A' ,'Pkabs_mean_B_A' ,'Pkabs_mean_V_A' ,'fPkabs_mean_B_A' ,'fPkabs_mean_V_A',...
%     'P_fkperp_BA','P_fkperp_VA','P_mean_fkperp_BA','P_mean_fkperp_VA','fP_fkperp_BA','fP_fkperp_VA',...
%     'PfpkkB_A','PfpkkV_A','PfpkkB_mean_A','PfpkkV_mean_A','f_mat',...
%     'PfkkB_A','PfkkV_A','PfkkB_mean_A','PfkkV_mean_A');