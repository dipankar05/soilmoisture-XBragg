%% Inversion X-bragg H-A-alpha method
% Ref: Irena Hajnsek, "INVERSION OF SURFACE PARAMETERS USING POLARIMETRIC SAR", Dissertation, pp.165-175
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
c11c = 0.0025;
c12r = -0.00465;
c12img = 0.00163;
c12c = c12r + c12img*1i;
c21c = conj(c12c);
c22c = 0.00970;
C2 = [c11c c12c; c21c c22c];
% transmit chi
chi_transmit = -45; % degree -ve for R and +ve for L

% Observed Stokes LH-LV/RH-RV
g0 = c11c + c22c;
g1 = c11c - c22c;
g2 = c12c + c21c;
if (chi_transmit >= 0)
    g3 = (1i.*(c12c - c21c)); %The sign is according to RC or LC sign !!
end

if (chi_transmit < 0)
    g3 = -(1i.*(c12c - c21c)); %The sign is according to RC or LC sign !!
end

%% Observed alpha_c
alphaC = rad2deg(0.5 * atan(real(g1./g3)));



%% Soil moisture
[mindiff, row] = min(0.5*sqrt((LUT(:,2)-alphaC).^2));
mv = LUT(row,1);
mv_percent = mv*100;


%% Print values
fprintf('Soil roughness rms (cm)= %0.2f \n',rms_cm);
fprintf('Volumetric soil moisture (percent) = %0.2f \n',mv_percent);