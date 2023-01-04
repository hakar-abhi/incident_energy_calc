% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Steps to find incident energy taking atmospheric effects into
%   consideration
%   Determine Q,N,B (we have script to do B calculation now, Q,N are inputs)
%   Estimate Ho (without tilt angle, scripts is already there to do it)
%   Estimate Hot (With B calculated, we have script to calculate Hot)
%   Estimate Kt=f(Q,N,watervapor,etc), Kt - Clearness Index
%   Find Rds - Direct radiation factor = Hot/Ho
%   Find rt (diffused, reflected, direct radiation factor)
%   Calculate Hat = rt * Kt * Ho (Hatmin is used to size PV panels which is
%   our company's goal)

% Kt estimation (Kte) for certain place for certain time also requires modeling
% as it is a multivariate function; where Q, N are known parameters but
% watervapor is seasonal/periodic quantity so we need to estimate it,
% it is periodic so we use sine waves to estimate it / we use Fourier
% Series from Signal/Systems to estimate a periodic function.

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Fourier Series Curve Fit for Clearness Index for Nepal

% Clearness Index Estimate Model:

% we use model Kt = Kte + error where Kte is estimated Kt

% Kte/f(x) = A1 + A2.sint + A3.sin2t + A4.sin3t + A5.cost + A6.cos2t + A7.cos3t

% where 
%   Ai = ai1 + ai2.x + ai3.x^2 + ai4.w + ai5.w^2
%   x = (Q -35)    in degrees, Q-35 is used as this gives better accuracy
%   for the model.
%   t = 2*pi*(N-80)/365
%   w = water vapor content in estimated from another Fourier Curve fit
%       model
%   w = G1 + G2.sint + G3.sin2t + G4.sin3t + G5.cost + G6.cos2t + G7.cos3t
%   Gi = gi1 + gi2.x + gi3.x^2, I have another script for that.

%   The program gives the coeffecient ai1,..... ai5 for Clearness Index
%   The data taken into consideration are Kt (measured clearness Index), N
%   (day number of year), W (Water content in atmosphere) and latitude Q.
%   The data I took are for locations in India, because I found
%   Kt(measured) for India(not Nepal) but data if found or calculated for
%   Nepal/any place can be used here too.

% Kt can be calculated using formula Kt = Ha/Ho
%  where Ho is incident energy on horizontal flat surface, we have script
%  to calculate it, Hoa is incident energy on horizontal flat plate but
%  taking atmospheric consideration into account, we use pyranometer to
%  measure insolation (W/m2). Each day of a year we measure it, or certain
%  number of days we measure insolation then integrate particular day
%  insolation for particular day 'N' which gives Wh/m2/day i.e,
%  Ha, then find Kt for that particular day using Ha/Ho.

clc 
clear

% measured Kt data of India for specific locations

Kt = [0.7070;0.7164;0.7198;0.7073;0.6958;0.5998;0.4554;0.4393;0.6100;0.6952;0.7157;0.7051;
    0.7231;0.7297;0.7242;0.7007;0.6991;0.5693;0.4330;0.4167;0.5969;0.7041;0.7239;0.7095;
    0.6694;0.6680;0.6766;0.6745;0.6811;0.5305;0.3836;0.3851;0.5229;0.6376;0.6730;0.6661;
    0.6012;0.6410;0.6441;0.6268;0.6366;0.4617;0.4270;0.4101;0.4740;0.5464;0.5924;0.6083;
    0.7394;0.7349;0.7289;0.7087;0.6887;0.6467;0.5592;0.5537;0.6660;0.7269;0.7501;0.7334;
    0.7707;0.7634;0.7009;0.6259;0.5755;0.5086;0.4208;0.4325;0.4565;0.4381;0.5251;0.6662;
    0.6427;0.6882;0.6955;0.6720;0.6467;0.5882;0.5385;0.5574;0.5953;0.5546;0.5453;0.5590;
    0.6906;0.6887;0.6716;0.6611;0.6493;0.5438;0.4021;0.3969;0.5470;0.6770;0.7087;0.6997;
    0.6792;0.6998;0.7108;0.6902;0.6751;0.5937;0.5044;0.5159;0.6504;0.7072;0.7244;0.7034;
    0.6904;0.7098;0.7034;0.6902;0.6972;0.5488;0.4216;0.4315;0.5434;0.6668;0.6888;0.6932;
    0.6809;0.6748;0.6679;0.6165;0.5771;0.5508;0.5108;0.5549;0.6043;0.5683;0.5719;0.6377;
    0.7068;0.6999;0.6770;0.6453;0.6375;0.5026;0.4483;0.4827;0.5643;0.6685;0.7029;0.7030];

% Day Number when measurement was done

N = [15;46;74;105;135;166;196;227;258;288;319;349;
    15;46;74;105;135;166;196;227;258;288;319;349;
    15;46;74;105;135;166;196;227;258;288;319;349;
    15;46;74;105;135;166;196;227;258;288;319;349;
    15;46;74;105;135;166;196;227;258;288;319;349;
    15;46;74;105;135;166;196;227;258;288;319;349;
    15;46;74;105;135;166;196;227;258;288;319;349;
    15;46;74;105;135;166;196;227;258;288;319;349;
    15;46;74;105;135;166;196;227;258;288;319;349;
    15;46;74;105;135;166;196;227;258;288;319;349;
    15;46;74;105;135;166;196;227;258;288;319;349;
    15;46;74;105;135;166;196;227;258;288;319;349];

% Latitude in Degrees

Q = [23.07;23.07;23.07;23.07;23.07;23.07;23.07;23.07;23.07;23.07;23.07;23.07;
    21.75;21.75;21.75;21.75;21.75;21.75;21.75;21.75;21.75;21.75;21.75;21.75;
    19.12;19.12;19.12;19.12;19.12;19.12;19.12;19.12;19.12;19.12;19.12;19.12;
    22.65;22.65;22.65;22.65;22.65;22.65;22.65;22.65;22.65;22.65;22.65;22.65;
    26.30;26.30;26.30;26.30;26.30;26.30;26.30;26.30;26.30;26.30;26.30;26.30;
    10.23;10.23;10.23;10.23;10.23;10.23;10.23;10.23;10.23;10.23;10.23;10.23;
    13.00;13.00;13.00;13.00;13.00;13.00;13.00;13.00;13.00;13.00;13.00;13.00;
    21.15;21.15;21.15;21.15;21.15;21.15;21.15;21.15;21.15;21.15;21.15;21.15;
    28.58;28.58;28.58;28.58;28.58;28.58;28.58;28.58;28.58;28.58;28.58;28.58;
    18.53;18.53;18.53;18.53;18.53;18.53;18.53;18.53;18.53;18.53;18.53;18.53;
    08.48;08.48;08.48;08.48;08.48;08.48;08.48;08.48;08.48;08.48;08.48;08.48;
    17.72;17.72;17.72;17.72;17.72;17.72;17.72;17.72;17.72;17.72;17.72;17.72];


% water vapor content


w = [1.67;1.82;2.16;2.73;3.33;4.66;5.37;5.32;4.40;3.27;2.00;1.91;
    1.84;1.85;2.30;2.88;3.42;4.81;5.32;5.16;4.59;3.31;2.29;2.03;
    2.70;2.55;2.90;3.30;3.95;5.17;5.36;5.05;4.78;4.01;2.94;2.93;
    2.10;2.22;2.63;3.49;4.30;5.62;6.28;6.17;5.69;4.55;2.86;2.06;
    1.59;1.60;1.83;2.29;2.86;4.22;5.23;5.25;3.91;2.53;1.71;1.73;
    0.74;0.77;0.77;1.12;1.58;1.86;1.88;1.77;1.82;1.63;1.04;0.91;
    2.73;2.66;2.51;3.25;4.36;5.07;5.09;4.94;4.82;4.64;3.45;3.12;
    1.75;2.09;2.08;2.58;2.99;4.66;5.32;5.41;4.59;3.53;1.99;1.89;
    1.40;1.38;1.66;1.93;2.55;4.01;5.56;5.74;4.19;2.67;1.50;1.42;
    1.77;1.75;2.06;2.45;2.87;3.94;4.26;4.12;3.86;3.15;2.29;2.09;
    3.02;3.16;3.49;4.30;4.47;4.47;4.46;4.29;4.39;4.37;4.02;3.67;
    2.73;3.21;3.08;3.98;4.59;5.50;5.55;5.48;5.39;4.91;3.18;2.67];


x = (Q-35);
t = ((2*pi)/365).*(N-80);

n = length(x);

c(1:n,1)=1;
c(1:n,2)=x;
c(1:n,3)=x.*x;
c(1:n,4)=w;
c(1:n,5)=w.*w;

c(1:n,6)=sin(t);
c(1:n,7)=x.*sin(t);
c(1:n,8)=x.*x.*sin(t);
c(1:n,9)=w.*sin(t);
c(1:n,10)=w.*w.*sin(t);

c(1:n,11)=sin(2*t);
c(1:n,12)=x.*sin(2*t);
c(1:n,13)=x.*x.*sin(2*t);
c(1:n,14)=w.*sin(2*t);
c(1:n,15)=w.*w.*sin(2*t);

c(1:n,16)=sin(3*t);
c(1:n,17)=x.*sin(3*t);
c(1:n,18)=x.*x.*sin(3*t);
c(1:n,19)=w.*sin(3*t);
c(1:n,20)=w.*w.*sin(3*t);

c(1:n,21)=cos(t);
c(1:n,22)=x.*cos(t);
c(1:n,23)=x.*x.*cos(t);
c(1:n,24)=w.*cos(t);
c(1:n,25)=w.*w.*cos(t);


c(1:n,26)=cos(2*t);
c(1:n,27)=x.*cos(2*t);
c(1:n,28)=x.*x.*cos(2*t);
c(1:n,29)=w.*cos(2*t);
c(1:n,30)=w.*w.*cos(2*t);

c(1:n,31)=cos(3*t);
c(1:n,32)=x.*cos(3*t);
c(1:n,33)=x.*x.*cos(3*t);
c(1:n,34)=w.*cos(3*t);
c(1:n,35)=w.*w.*cos(3*t);


[row, col] = size(c);

for i = 1:col
    for j = 1:col
        A(i,j) = sum(c(1:n,j).*c(1:n,i));
        b(i) = sum(Kt.*c(1:n,i));
    end
end

cof = A\b'
k = 1;
for i = 1:7
    for j =1:5
        coeff(i,j) = cof(k);
        k = k+1;
end
end

coeff













