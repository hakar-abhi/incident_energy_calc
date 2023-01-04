% Script for plotting spectrum irradiance and solar insolation  

%  solar measured data is in data.m
% downloaded data from https://www.nrel.gov/grid/solar-resource/spectra-am1.5.html


clear
data %loaded data from data.m
Lambda_nm = sol(:,1);
irrad= sol(:,2);

plot(Lambda_nm,irrad),grid,xlabel('wavelength-nm'),ylabel('irradiance'),title('Solar Irradiance Spectrum')

%  Integrate Irradiance to obtain insolation, L 
% trapezoidal integration is used
%yint = yint+ dL*(x1+x2)/2

yint =0;
insol(1)=0;
for i = 2:length(Lambda_nm)
    yint = yint+(Lambda_nm(i)-Lambda_nm(i-1))*(irrad(i)+irrad(i-1))/2;
    insol(i)=yint;
end
insol=insol;

plot(Lambda_nm,irrad,Lambda_nm,insol/1000),grid,xlabel('wavelength-nm'),ylabel('irradiance and insolation'),title('Solar Irradiance Spectrum and Solar Insolation')






