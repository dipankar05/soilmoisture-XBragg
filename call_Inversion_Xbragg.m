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
load('LUT.mat'); %Load mv-H-alpha LUT generated from forward X-bragg model

%% Observed Coherency matrix
T3 = [1.8236 -0.1310-0.0053i 0; -0.1310+0.0053i 0.0239 0; 0 0 0.0165];

%% H-A-Alpha decomposition
[alpha, H, A] = HAalphadecomp(T3);

if (H>0.5 && alpha>45)
    rms = 0;
    mvpercent = 0;
    % masked out infeasible region
else
%% Soil roughness
f = 1.5; %radar frequency GHz
k = 2*pi*f/0.3; %calculate the wavenumber 
ks = 1-A; 
% Ref: CLOUDE, S. R., ‘Eigenvalue Parameters for Surface Roughness Studies’, Proceedings of SPIE Conference on Polarization: Measurement, Analysis and Remote Sensing II, vol. 3754,Denver, Colorado, USA, July 1999.
rms_meter = ks./k;
rms_cm = rms_meter*100; %%cm
clear rms_meter k ks f;

%% Soil moisture
[mindiff, row] = min(0.5.*sqrt(((LUT(:,2)-H)./max(LUT(:,2))).^2+((LUT(:,3)-alpha)./max(LUT(:,3))).^2));
mv = LUT(row,1);
mv_percent = mv*100;
end

%% Print values
fprintf('Soil roughness rms (cm)= %0.2f \n',rms_cm);
fprintf('Volumetric soil moisture (percent) = %0.2f \n',mv_percent);