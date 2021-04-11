%% Forward X-bragg model to generate H-alpha LUT for given volumetric soil moisture mv
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

%% Incidence angle (SAR geometry)
inc = 35; %in degree
%% Volumetric soil moisture given real dielectric constant of soil
epr = 2:0.1:40;
eps = complex(epr,0.00);
[xx,yy] = size(epr);
for i = 1:yy
mv(i) = 0.01*(-5.3+(2.92.*epr(i))-(0.055.*epr(i).*epr(i))+(0.0004.*epr(i).*epr(i).*epr(i)));  %% Topp et al. 1980 model
end

%% Plotting with different width beta (0<b<90) distribution
b10 = 10; % in degree
for i = 1:yy
alphaC10(i) = xbragg_compactpol(b10,eps(i),inc);
end
%% Plotting with varying eps
figure
line(mv,alphaC10,'Color','red','LineStyle','-','LineWidth',1.2)
hold on
text(max(mv),max(alphaC10),num2str('\beta=10'))

%%
b20 = 20;
for i = 1:yy
alphaC20(i) = xbragg_compactpol(b20,eps(i),inc);
end
%% Plotting with varying eps
line(mv,alphaC20,'Color','red','LineStyle','-','LineWidth',1.2)
text(max(mv),max(alphaC20),num2str('\beta=20'))

%%
b40 = 40;
for i = 1:yy
alphaC40(i) = xbragg_compactpol(b40,eps(i),inc);
end
%% Plotting with varying eps
line(mv,alphaC40,'Color','red','LineStyle','-','LineWidth',1.2)
text(max(mv),max(alphaC40),num2str('\beta=40'))

%%
b60 = 60;
for i = 1:yy
alphaC60(i) = xbragg_compactpol(b60,eps(i),inc);
end
%% Plotting with varying eps
line(mv,alphaC60,'Color','red','LineStyle','-','LineWidth',1.2)
text(max(mv),max(alphaC60),num2str('\beta=60'))

%%
b70 = 70;
for i = 1:yy
alphaC70(i) = xbragg_compactpol(b70,eps(i),inc);
end
%% Plotting with varying eps
line(mv,alphaC70,'Color','red','LineStyle','-','LineWidth',1.2)
text(max(mv),max(alphaC70),num2str('\beta=70'))

%%
b80 = 80;
for i = 1:yy
alphaC80(i) = xbragg_compactpol(b80,eps(i),inc);
end
%% Plotting with varying eps
line(mv,alphaC80,'Color','red','LineStyle','-','LineWidth',1.2)
text(max(mv),max(alphaC80),num2str('\beta=80'))

%%
b90 = 90;
for i = 1:yy
alphaC90(i) = xbragg_compactpol(b90,eps(i),inc);
end
%% Plotting with varying eps
line(mv,alphaC90,'Color','red','LineStyle','-','LineWidth',1.2)
text(max(mv),max(alphaC90),num2str('\beta=90'))

%%-----------------------------------------------------------------------------------
%% Build LUT mv, H, alpha --with stack of different beta
LUT_CP = [mv',alphaC10'; mv',alphaC20'; mv',alphaC40'; mv',alphaC60'; mv',alphaC70'; mv',alphaC80'; mv',alphaC90'];

%% Saving LUT 
save('F:\Github_Dipankar\SAR-processing-matlab\XBragg_surfacemodel\LUT_CP.mat','LUT_CP')


%% Loading LUT from file
% clear
% load('test.mat')