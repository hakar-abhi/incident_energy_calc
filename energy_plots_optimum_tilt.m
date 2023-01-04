%   Description: This file calculates the incident energy on a horizontal/tilted surface given the
%   latitude and day number
%     
%   INPUTS: N = Day Number, N=1 for Jan 1st and N=365 for December 31st, Feb 29th is not included.

%           Q = Latitude of the place in degrees (need to convert to radians, places on earth are on
%               certain arc lengths and arc length are measured in radians as per standards.
%   Author: Abhisek Adhikari
%   Created on: Jun 2022

clc
clear

%Inputs
%Locality = 27 degrees, 39 minutes; (Just used my office location, you can
%           use your own)

Q = 27.65;  % expressed in degrees
Q = Q * pi/180;

% B = Q;      %tilt angle in radians

% constants

Lsc = 1.37; % kW/m2 mean solar constant (does not depend on earth's position relative to sun)

fprintf('Hot-MIN \t Hot-MAX \t H-RIPPLE \t LATITUDE,deg \t TILTANGLE,deg\n');


Bmax = 30 * pi/180;  %just making an assumption for maximum tilt angle, you can have your own

for B = 0:Bmax/30:Bmax
    
    for N = 1:365,
    %Claculation of Declination
    t = 2*pi*(N-80)/365;
    d = 23.45*sin(t)*(pi/180); % declination expressed in radians
    
    %calculation of extraterrestrial insolation scale factor
    k = 1 + 0.033*cos(2*pi*N/365); % for scaling Lsc depending on number of day in year, derived
                                   % empiricaly, I took from certain literatures
    
    wsr = acos(-1 * tan(Q)*tan(d)); %solar hours calculation
    
    wsrB = acos(-1 * tan(Q-B)*tan(d)); %hour angle based on tilt angle B
    
    wsrt = min(wsr, wsrB);
    
    Ho(N) = (24*k*Lsc/pi)*(cos(Q)*cos(d)*sin(wsr)+wsr*sin(Q)*sin(d));
    
    Hot(N) = (24*k*Lsc/pi)*(cos(Q-B)*cos(d)*sin(wsrt)+wsrt*sin(Q-B)*sin(d)); %energy on tilted surface
    
    days(N) = N;
    end % of day number loop


    
fprintf('%f \t %f \t %f \t %f \t \t %f\n', min(Hot), max(Hot), max(Hot)-min(Hot), Q*180/pi, B*180/pi);
    
%show results

plot(days,Ho, days, Hot);
grid, xlabel('Day Number, N'), ylabel('kWh/m2/day'), title('Ho and Hot');

pause

end 
 
    
    
    
