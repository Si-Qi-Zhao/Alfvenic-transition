function [fre_SVD,WaveVector_B] = SVD_B_2022(dB_ts,fmin,fmax,numf)
  %% attention: xyz coordinates: theta is not kB0 angle
    ts = dB_ts(:,1);
    dcBx = [ts dB_ts(:,2)];
    dcBy = [ts dB_ts(:,3)];
    dcBz = [ts dB_ts(:,4)];
 
 [Q_11,f] = xwt_2022(dcBx,dcBx,fmin,fmax,numf);
 [Q_12,f] = xwt_2022(dcBx,dcBy,fmin,fmax,numf);
 [Q_13,f] = xwt_2022(dcBx,dcBz,fmin,fmax,numf);

%  [Q_21,f] = xwt_2022(dcBy,dcBx,fmin,fmax,numf);   % no use in A matrix
 [Q_22,f] = xwt_2022(dcBy,dcBy,fmin,fmax,numf);
 [Q_23,f] = xwt_2022(dcBy,dcBz,fmin,fmax,numf);

%  [Q_31,f] = xwt_2022(dcBz,dcBx,fmin,fmax,numf);   % no use in A matrix
%  [Q_32,f] = xwt_2022(dcBz,dcBy,fmin,fmax,numf);   % no use in A matrix
 [Q_33,f] = xwt_2022(dcBz,dcBz,fmin,fmax,numf);
 
fre_SVD = f;
length_f=size(Q_11,1);
length_t=size(Q_11,2); 
WaveVector_B   = zeros(length_f,length_t,3);
% Lp_only_B      = zeros(length_f,length_t);
% theta_only_B   = zeros(length_f,length_t);

for j=1:length_t
    for i=1:length_f
          A_only_B=[real(Q_11(i,j)), real(Q_12(i,j)), real(Q_13(i,j));...
                    real(Q_12(i,j)), real(Q_22(i,j)), real(Q_23(i,j));...
                    real(Q_13(i,j)), real(Q_23(i,j)), real(Q_33(i,j));...
                         0,        -imag(Q_12(i,j)), -imag(Q_13(i,j));...
                    imag(Q_12(i,j)),       0,        -imag(Q_23(i,j));...
                    imag(Q_13(i,j)), imag(Q_23(i,j)),      0       ];
              
         [U_only_B, W_only_B, V_only_B] = svd(A_only_B,0); 
         W_B_reciprocal = [1/W_only_B(1,1) 0 0;0 1/W_only_B(2,2) 0;0 0 1/W_only_B(3,3)];  % 1/n  reciprocal
         ifnan = W_B_reciprocal;
         if isinf(ifnan(1))
            WaveVector_B(i,j,:) = [nan,nan,nan];
         else
            tt = V_only_B';
            k_only_B = tt(3,:);
            % cannot identifly two antiparallel wave vectors
            if k_only_B(1) > 0
               k_only_B = -1.*k_only_B;
            end
            WaveVector_B(i,j,:) = k_only_B;
         end
    end
end         
end

%             B0 = mean(Bav.data,1);
%             theta_temp = acosd(dot(k_only_B,B0)/norm(B0)/norm(k_only_B));
%             if theta_temp < 90
%               theta_only_B(i,j) = theta_temp;
%             else
%               theta_only_B(i,j) = 180 - theta_temp;
%             end
         
       % phi_only_B(i,j)=atan2(k_only_B(1),k_only_B(2)).*180./3.14159;      
%          W_only_B_sort = sort([W_only_B(1,1),W_only_B(2,2),W_only_B(3,3)]);
%        F_only_B(i,j)=1-sqrt(W_only_B_sort(1)/W_only_B_sort(3));
%          Lp_only_B(i,j) = W_only_B_sort(2)/W_only_B_sort(3);
 
