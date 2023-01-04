%   Description: This file calculates the incident energy on a horizontal/tilted surface given the
%   latitude and day number with atmospheric effects taken into consideration.
 


%   INPUTS: N = Day Number, N=1 for Jan 1st and N=365 for December 31st, Feb 29th is not included.

%           Q = Latitude of the place in degrees (need to convert to radians, places on earth are on
%               certain arc lengths and arc length are measured in radians as per standards.

%   Author: Abhisek Adhikari
%   Created on: Jun 2022

clc
clear

%Inputs
%Locality = 27 degrees, 39 minutes; (Just used my office location, you can
%           use your required site location)

Q = 12.97;  % expressed in degrees, took different location to adjust with data gathered for modeling

x = Q-35; %Fourier Curve Fit variable
Q = Q * pi/180; %now converted latitude into radians

B = 10*pi/180; %tilt angle in radians as found optimal without atmosphere

% constants

Lsc = 1.37; % kW/m2 mean solar constant (does not depend on earth's position relative to sun)

rho = 0.2; %reflection coeffecient, small value since we assume the site is in plains

G =[ 1.6225   -0.2288   -0.0068;     % water vapor estimate equation coeffeients
    1.6881   -0.0069   -0.0020;
   -1.2463   -0.0816   -0.0010;
    0.5667    0.0714    0.0019;
   -1.2174    0.0059    0.0014;
   -0.0975   -0.0131   -0.0004;
    0.4983    0.0607    0.0015 ];
% Clearness index estimate quation coeffecients

 A = [ 0.5564  0.0089    0.0002    0.0743   -0.0089;
   -0.2351    0.0119    0.0004    0.1475   -0.0238;
   -0.1007   -0.0090   -0.0004    0.1027   -0.0200;
    0.0134    0.0040    0.0002   -0.0070    0.0010;
    0.1300   -0.0133   -0.0003   -0.0848    0.0098;
   -0.0601    0.0048    0.0002    0.0734   -0.0133;
    0.0968    0.0058    0.0002   -0.0282    0.0010];


%Calculation of insolation and energy on all days of the year
    
for N = 1:365
    %Claculation of Declination, Declination accounts for days in a year/ time of a year.
    t = 2*pi*(N-80)/365;       %equinox is at N=80  so N-80 is done.
    d = 23.45*sin(t)*(pi/180); % declination expressed in radians, this formula gives best approximation for declination angle, found in literature. 
    
    %calculation of extraterrestrial insolation scale factor
    k = 1 + 0.033*cos(2*pi*N/365); % for scaling Lsc depending on number of day in year, derived
                                   % empiricaly, I took from certain literatures
    
    wsr = acos(-1 * tan(Q)*tan(d)); %solar hours calculation
          
    wsrB = acos(-1 * tan(Q-B)*tan(d)); %hour angle based on tilt angle B
    
    wsrt = min(wsr, wsrB);
    
    Ho(N) = (24*k*Lsc/pi)*(cos(Q)*cos(d)*sin(wsr)+wsr*sin(Q)*sin(d));
    
    Hot(N) = (24*k*Lsc/pi)*(cos(Q-B)*cos(d)*sin(wsrt)+wsrt*sin(Q-B)*sin(d)); %energy on tilted surface
    
    days(N) = N;
    
    %Introduction of atmospheric effects
    %estimate water vapor content
    
    XX=[1;x;x*x];
    G1= G(1,1:3)*XX;
    G2= G(2,1:3)*XX;
    G3= G(3,1:3)*XX;
    G4= G(4,1:3)*XX;
    G5= G(5,1:3)*XX;
    G6= G(6,1:3)*XX;
    G7= G(7,1:3)*XX;
   W(N) = G1+G2*sin(t)+G3*sin(2*t)+G4*sin(4*t)+G5*cos(t)+G6*cos(2*t)+G7*cos(3*t); 
    
    % estimated clearness index
    
    YY=[1;x;x*x;W(N);W(N)*W(N)];
    A1= A(1,1:5)*YY;
    A2= A(2,1:5)*YY;
    A3= A(3,1:5)*YY;
    A4= A(4,1:5)*YY;
    A5= A(5,1:5)*YY;
    A6= A(6,1:5)*YY;
    A7= A(7,1:5)*YY;
   Kte(N) = A1+A2*sin(t)+A3*sin(2*t)+A4*sin(4*t)+A5*cos(t)+A6*cos(2*t)+A7*cos(3*t); 
    
    Rd(N) = Hot(N)/Ho(N); %tilt factor
    
    Kd(N) = 1-1.13*Kte(N); %Kd(N)=Hd/Ha diffuse radiation factor, approximated it using the equation.
    
    rt(N) = ((1-Kd(N))*Rd(N)) + (Kd(N)*(1+cos(B))/2) + (rho*(1-cos(B))/2); %overall titlt factor
    
end

Hat_direct = Rd.*Kte.*Ho; % diffused and reflected factor not included

Hat = rt.*Kte.*Ho; % tilted surface with atmospheric effects  

%show results

Hat_trans = (Hat'/6).*1000


figure(1),plot(days,Ho, days, Hot); 
h=plot(days,Ho, days, Hot);
legend(h,'Ho','Hot');
grid, xlabel('Day Number, N'), ylabel('kWh/m2/day'), title('Ho and Hat versus Day of the year');

figure(2),plot(days,Rd,'b',days,rt,'g');

g= plot(days,Rd,'b',days,rt,'g');
legend(g,'Rd','rt');
grid, xlabel('Day Number, N'), ylabel('tilt factor'), title('Rd vs rt');


figure(3),plot(days,Kte);
i = plot(days,Kte);
legend(i,'Kte');
grid, xlabel('Day Number, N'), ylabel('Kte'), title('Clearness Index vs Day of Year');    
