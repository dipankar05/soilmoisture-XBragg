%% Inversion X-bragg compact-pol alpha method
% Ref: Ponnurangam et al. (2016) "Soil Moisture Estimation Using Hybrid Polarimetric SAR Data of RISAT-1", IEEE TGRS. 
% @author: Dr. Dipankar Mandal
%%  ---------------------------------------------------------------------------------------
%   ---------------------------------------------------------------------------------------
%   Copyright (C) 2021 by Microwave Remote Sensing Lab, IITBombay http://www.mrslab.in
%
%   This program is free software; you can redistribute it and/or modify it
%   under the terms of the GNU General Public License as published by the Free
%   Software Foundation; either version 3 of the License, or (at your option)
%   any later version.
%   This program is distributed in the hope that it will be useful, but WITHOUT
%   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
%   FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
%   more details.
%
%   You should have received a copy of the GNU General Public License along
%   with this program; if not, see http://www.gnu.org/licenses/
%   ---------------------------------------------------------------------------------------
%%
clear;
load('LUT_CP.mat'); %Load mv-H-alpha LUT generated from forward X-bragg model

%% Observed covariance matrix C2 of RH-RV
c11c = 0.0032;
c12r = -0.000065;
c12img = 0.0013;
c22c = 0.0039;
c12c = c12r + c12img*1i;
c21c = conj(c12c);
C2 = [c11c c12c; c21c c22c];
% transmit chi
% chi_transmit = -45; % degree -ve for R and +ve for L
%%
% Observed Stokes LH-LV/RH-RV
g0 = c11c + c22c;
g1 = c11c - c22c;
g2 = c12c + c21c;
g3 = -(1i.*(c12c - c21c)); %sign for RC case

%% conformity coefficient
cc = (2*c12img)/(c11c+c22c);

%% Observed alpha_c
% Ref: Cloude et al. 2012, "Compact Decomposition Theory", IEEE GRSL, p. 31
alphaC = 0.5 * atan2d(g3,sqrt(g1.^2 + g2.^2)); %%tan2d(Y,X)

if cc >0.5
    mv_percent = NaN; %% mask regions with cc>0.5
else
%% Soil moisture
[mindiff, row] = min(0.5*sqrt((LUT_CP(:,2)-alphaC).^2));
mv = LUT_CP(row,1);
mv_percent = mv*100;
end


%% Print values
fprintf('Volumetric soil moisture (percent) = %0.2f \n',mv_percent);