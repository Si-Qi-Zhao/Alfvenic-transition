function [P_perp,P_para,kperp,kpara,kpara_new,Ekperp,Ekperp_pre2,XiB_new,bk_new,XiB_mat_sc,bk2_sc] = nonlinearity_A(date,trange,numk,angle_threshod,B0)

file_path = ['/Users/zhaosiqi/matlab_Figure/',date,'/',trange,'min-window/'];

%% loadin 
load([file_path,'Pkk_2D_numk',num2str(numk),'_kklt',num2str(angle_threshod),'_VB_nonlinear'],...
      'P_mean_fkperp_BA','PfkkB_mean_A','xvecA','fre1','P_mean_fkpara_BA');

idx_perp = xvecA >= 4e-5;
idx_para = xvecA >= 5e-6;
kpara = xvecA(idx_para);
kperp = xvecA(idx_perp);

%% test
diff_f = abs(diff(fre1)); diff_f = [diff_f;0];
dfmat  = repmat(diff_f,1,length(xvecA));  
Pperp  = nansum(P_mean_fkperp_BA.*dfmat,1);
Ppara  = nansum(P_mean_fkpara_BA.*dfmat,1);

P_mean_fkperp_BA = P_mean_fkperp_BA(:,idx_perp);        % sc frame
PfkkB_mean_A     = PfkkB_mean_A(:,idx_para,idx_perp);   % sc frame
Pperp = Pperp(idx_perp);
Ppara = Ppara(idx_para);
% 
% Pkk = zeros(length(kpara),length(kperp));
% for  i=1:1:length(kpara)
%     for j = 1:1:length(kperp)
%        Pkk(i,j) = trapz(flipud(fre1),flipud(PfkkB_mean_A(:,i,j)));
%     end
% end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
idx_para = 1:400;
% P_para = sum(Pkk,2);
% P_perp = sum(Pkk,1);
P_para = Ppara;
P_perp = Pperp;

kpara_new = interp1(P_para(idx_para),kpara(idx_para),P_perp);
%% sc frame
P_sumk_sc = zeros(length(fre1),length(kpara),length(kperp));
for i = 1:1:length(fre1)
    tmp = squeeze(PfkkB_mean_A(i,:,:));
    tmp_para = cumsum(tmp,2,'reverse');
    tmp_all  = cumsum(tmp_para,1,'reverse');
    P_sumk_sc(i,:,:) = tmp_all;   
end

bk2_sc = zeros(length(kpara),length(kperp));
for i=1:1:length(kpara)-1
    for j = 1:1:length(kperp)-1
       bk2_sc(i,j) = 2*trapz(flipud(fre1),flipud(P_sumk_sc(:,i,j)));
    end
end
bk_sc = sqrt(bk2_sc); 

XiB_mat_sc = zeros(length(kpara),length(kperp));
for i=1:1:length(kpara)-1
    XiB_mat_sc(i,:) = kperp.*bk_sc(i,:)./kpara(i)./B0;
end

%% test
% bk2_sc_para = nansum(bk2_sc,2);
% Ekpara = 0.5*bk2_sc_para'./kpara;
% Ekpara_pre = kpara.^(2).*Ekpara;
% 
% kpara_mat = repmat(kpara',1,length(kperp)); 
% kperp_mat = repmat(kperp, length(kpara),1);  
% % E_mat = 0.5*bk2_sc./kpara_mat./kperp_mat; 
% E_mat = 0.5*bk2_sc./kperp_mat;
% E_mat_norm = E_mat./max(max(E_mat));
%%
P_sumk_perp = zeros(length(fre1),length(kperp));
for i=1:1:length(kperp)-1
    P_sumk_perp(:,i) = sum(P_mean_fkperp_BA(:,i:end),2);
end

% test sum
dfmat2  = repmat(diff_f,1,length(kperp));  
bk_new2 = 2*nansum(dfmat2.*P_sumk_perp,1);

% integrate
% bk_new2 = 2*trapz(flipud(fre1),flipud(P_sumk_perp));

bk_new  = sqrt(bk_new2);

Ekperp = 0.5*bk_new2./kperp;
Ekperp_pre2 = kperp.^(5/3).*Ekperp;

XiB_new = kperp.*bk_new./kpara_new./B0;
end

%% plasma frame
% P_sumk = zeros(length(f_mat),length(kpara),length(kperp));
% for i = 1:1:length(f_mat)
%   tmp = squeeze(PfpkkB_mean_A(i,:,:));
%   tmp_para = cumsum(tmp,2,'reverse');
%   tmp_all  = cumsum(tmp_para,1,'reverse');
%   P_sumk(i,:,:) = tmp_all;   
% end
% 
% bk2 = zeros(length(kpara),length(kperp));
% for  i=1:1:length(kpara)-1
%     for j = 1:1:length(kperp)-1
%        bk2(i,j) = 2*trapz(flipud(f_mat),flipud(P_sumk(:,i,j)));
%     end
% end
% bk = sqrt(bk2); 
% 
% XiB_mat = zeros(length(kpara),length(kperp));
% for i=1:1:length(kpara)-1
%     XiB_mat (i,:) = kperp.*bk(i,:)./kpara(i)./B0;
% end


%% sc frame
% P_sumk_sc = zeros(length(fre1),length(kpara),length(kperp));
% for i = 1:1:length(fre1)
%     tmp = squeeze(PfkkB_mean_A_lin(i,:,:));
%     tmp_para = cumsum(tmp,2,'reverse');
%     tmp_all  = cumsum(tmp_para,1,'reverse');
%     P_sumk_sc(i,:,:) = tmp_all;   
% end
% 
% bk2_sc = zeros(length(kpara),length(kperp));
% for i=1:1:length(kpara)-1
%     for j = 1:1:length(kperp)-1
% %        bk2_sc(i,j) = 2.48e-04*2*sum((P_sumk_sc(:,i,j)),1);
%        bk2_sc(i,j) = 2*trapz(flipud(f_mat),flipud(P_sumk_sc(:,i,j)));
%     end
% end
% bk_sc = sqrt(bk2_sc); 
% 
% XiB_mat_sc = zeros(length(kpara),length(kperp));
% for i=1:1:length(kpara)-1
%     XiB_mat_sc(i,:) = kperp.*bk_sc(i,:)./kpara(i)./B0;
% end

%% sum
% df_sc = abs(diff(fre1));    df_sc = [df_sc;0];
% Pkk = zeros(length(kpara),length(kperp));
% for  i=1:1:length(kpara)-1
%     for j = 1:1:length(kperp)-1
%        Pkk(i,j) = nansum(df_sc.*PfkkB_A(:,i,j));
%     end
% end


% PfkkB_A_lin = PfkkB_A_lin(:,idx_para,idx_perp);      % sc frame
% P_mean_fkperp_BA_lin = P_mean_fkperp_BA_lin(:,idx_perp);
% PfpkkB_mean_A = PfpkkB_mean_A(:,idx_para,idx_perp);  % pl frame
% PfkkB_A = PfkkB_A(:,idx_para,idx_perp);              % sc frame
% PfkkB_mean_A_lin = PfkkB_mean_A_lin(:,idx_para,idx_perp);