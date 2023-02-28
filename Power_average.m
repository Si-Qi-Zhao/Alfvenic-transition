function Powerav = Power_average(fkPower,numf,cav,N)

posav = cav/2 + (0:1:N)*cav;
Powerav = zeros(N+1,numf);
for m = 1:N+1
     Powerav(m,:) = irf.nanmean(fkPower([posav(m)-cav/2+1:posav(m)+cav/2],:));
end
end