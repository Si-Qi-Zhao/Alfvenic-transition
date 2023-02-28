%% timing to get wave vectors
function [kx,ky,kz,Powerav] = timing_2022_v2(fkPower,numf,W1,W2,W3,W4,R1,R2,R3,R4,fre1,cav,N)

posav = cav/2 + (0:1:N)*cav;

cx12 = zeros(N+1,numf);
cx13 = zeros(N+1,numf);
cx14 = zeros(N+1,numf);
cx23 = zeros(N+1,numf);
cx24 = zeros(N+1,numf);
cx34 = zeros(N+1,numf);
Powerav = zeros(N+1,numf);
R1av = zeros(N+1,3);
R2av = zeros(N+1,3);
R3av = zeros(N+1,3);
R4av = zeros(N+1,3);
for m = 1:N+1
    cx12(m,:) = irf.nanmean(W1([posav(m)-cav/2+1:posav(m)+cav/2],:).*conj(W2([posav(m)-cav/2+1:posav(m)+cav/2],:)));
    cx13(m,:) = irf.nanmean(W1([posav(m)-cav/2+1:posav(m)+cav/2],:).*conj(W3([posav(m)-cav/2+1:posav(m)+cav/2],:)));
    cx14(m,:) = irf.nanmean(W1([posav(m)-cav/2+1:posav(m)+cav/2],:).*conj(W4([posav(m)-cav/2+1:posav(m)+cav/2],:)));
    cx23(m,:) = irf.nanmean(W2([posav(m)-cav/2+1:posav(m)+cav/2],:).*conj(W3([posav(m)-cav/2+1:posav(m)+cav/2],:)));
    cx24(m,:) = irf.nanmean(W2([posav(m)-cav/2+1:posav(m)+cav/2],:).*conj(W4([posav(m)-cav/2+1:posav(m)+cav/2],:)));
    cx34(m,:) = irf.nanmean(W3([posav(m)-cav/2+1:posav(m)+cav/2],:).*conj(W4([posav(m)-cav/2+1:posav(m)+cav/2],:)));
    Powerav(m,:) = irf.nanmean(fkPower([posav(m)-cav/2+1:posav(m)+cav/2],:));
    R1av(m,:)  = irf.nanmean(R1([posav(m)-cav/2+1:posav(m)+cav/2],:));
    R2av(m,:)  = irf.nanmean(R2([posav(m)-cav/2+1:posav(m)+cav/2],:));
    R3av(m,:)  = irf.nanmean(R3([posav(m)-cav/2+1:posav(m)+cav/2],:));
    R4av(m,:)  = irf.nanmean(R4([posav(m)-cav/2+1:posav(m)+cav/2],:)); 
end

% Compute phase differences between each spacecraft pair
th12 = atan2(imag(cx12),real(cx12));
th13 = atan2(imag(cx13),real(cx13));
th14 = atan2(imag(cx14),real(cx14));
th23 = atan2(imag(cx23),real(cx23));
th24 = atan2(imag(cx24),real(cx24));
th34 = atan2(imag(cx34),real(cx34));

wmat = 2*pi*ones(N+1,1)*(fre1)';

% Convert phase difference to time delay
dt12 = th12./wmat;
dt13 = th13./wmat;
dt14 = th14./wmat;
dt23 = th23./wmat;
dt24 = th24./wmat;
dt34 = th34./wmat;

% Weighted averaged time delay using all spacecraft pairs
dt2 = 0.5*dt12 + 0.2*(dt13 - dt23) + 0.2*(dt14 - dt24) + 0.1*(dt14 - dt34 - dt23);
dt3 = 0.5*dt13 + 0.2*(dt12 + dt23) + 0.2*(dt14 - dt34) + 0.1*(dt12 + dt24 - dt34);
dt4 = 0.5*dt14 + 0.2*(dt12 + dt24) + 0.2*(dt13 + dt34) + 0.1*(dt12 + dt23 + dt34);

% Compute phase speeds
kx = zeros(N+1,numf);
ky = zeros(N+1,numf);
kz = zeros(N+1,numf);

for ii = 1:N+1
	dR = [R2av(ii,:);R3av(ii,:);R4av(ii,:)]-[R1av(ii,:);R1av(ii,:);R1av(ii,:)];
  for jj = 1:numf
    % Volumetric tensor with SC1 as center.
    m = dR\[dt2(ii,jj);dt3(ii,jj);dt4(ii,jj)]; % "1/v vector"
    kx(ii,jj) = 2*pi*fre1(jj)*m(1);
    ky(ii,jj) = 2*pi*fre1(jj)*m(2);
    kz(ii,jj) = 2*pi*fre1(jj)*m(3);  
  end
end

% kx = kx/1e3; % m
% ky = ky/1e3;
% kz = kz/1e3;
end