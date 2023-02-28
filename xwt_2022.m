function varargout=xwt_2022(x,y,fmin,fmax,numf)
%% Cross wavelet transform

% USAGE: [Wxy,period,scale,coi,sig95]=xwt(x,y,[,settings])
% x & y: two time series

%% 
x_wavelet = irf_wavelet(x,'returnpower',0,'f',[fmin fmax],'nf',numf); 
y_wavelet = irf_wavelet(y,'returnpower',0,'f',[fmin fmax],'nf',numf);
%-----------:::::::::::::--------- ANALYZE ----------::::::::::::------------

X = cell2mat(x_wavelet.p); X = X.';
Y = cell2mat(y_wavelet.p); Y = Y.';

% -------- Cross
Wxy = X.*conj(Y); 
f   = y_wavelet.f;

ind_NaN = find(isnan(Wxy));  Wxy(ind_NaN) = 0;

varargout={Wxy,f};




