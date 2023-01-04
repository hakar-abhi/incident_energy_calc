%   Description: This file calculates the incident energy on a horizontal/tilted surface given the
%   latitude and day number and also finds the best tilt angle for each day
%   such that it tracks the sun each day for best incident energy.
%    
%   Disclaimer : No Atmospheric effect is taken into consideration, will do
%              later.

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

Q = 27.65;  % expressed in degrees
Q = Q * pi/180;



% constants

Lsc = 1.37; % kW/m2 mean solar constant (does not depend on earth's position relative to sun)






    
for N = 1:365
    %Claculation of Declination, Declination accounts for days in a year/ time of a year.
    t = 2*pi*(N-80)/365;       %equinox is at N=80 and at 365-80 so N-80 is done.
    d = 23.45*sin(t)*(pi/180); % declination expressed in radians, this formula gives best approximation for declination angle, found in literature. 
    
    %calculation of extraterrestrial insolation scale factor
    k = 1 + 0.033*cos(2*pi*N/365); % for scaling Lsc depending on number of day in year, derived
                                   % empiricaly, I took from certain literatures
    
    wsr = acos(-1 * tan(Q)*tan(d)); %solar hours calculation
    
    B = Q-d;    %tilt angle required is zero when Sun is in our latitude location
    
    wsrB = acos(-1 * tan(Q-B)*tan(d)); %hour angle based on tilt angle B
    
    wsrt = min(wsr, wsrB);
    
    Ho(N) = (24*k*Lsc/pi)*(cos(Q)*cos(d)*sin(wsr)+wsr*sin(Q)*sin(d));
    
    Hot(N) = (24*k*Lsc/pi)*(cos(Q-B)*cos(d)*sin(wsrt)+wsrt*sin(Q-B)*sin(d)); %energy on tilted surface
    
    days(N) = N;
    end % of day number loop


    

%show results

plot(days,Ho, days, Hot);
grid, xlabel('Day Number, N'), ylabel('kWh/m2/day'), title('Ho and Hot');

 
    
    
    
