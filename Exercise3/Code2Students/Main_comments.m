% Positioning and Location Based Services
% A.A. 2023/2024
% 3rd Exercise: Ionospheric delay computation
%
% @author: Marianna Alghisi


%load parameters
line1 = 'GPSA 7.4506D-09 1.4901D-08 -5.9605D-08 -1.1921D-07 IONOSPHERIC CORR';
line2 = 'GPSB 9.2160D+04 1.3107D+05 -6.5536D+04 -5.2429D+05 IONOSPHERIC CORR';
ionoparams = [cell2mat(textscan(line1, '%*s %f %f %f %f %*s')) ...
cell2mat(textscan(line2, '%*s %f %f %f %f %*s'))];
% zenithal maps of Iono corrections

%initialize values for the zenital cycle

%initialize matrix

%time, phi and lambda cycle

% plots

%%polar map in Milano
% Milano position in degrees
phi2 = 45 + 28 / 60 + 38.28 / 60^2; %degrees
lambda2 = 9 + 10 / 60 + 53.40 / 60^2; %degrees

%inizialize values for the cycle

% matrix inizialization

%time, elevation and azimuth cycle 

%plots
