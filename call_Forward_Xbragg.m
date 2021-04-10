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
%% H-Alpha curve
m = 0:0.01:1;
[x,y] = size(m);
%% Curve -1 --Lower curve
for i = 1:y
    H1(i) = ((-1)/(1+2*m(i)))*(log10(power(m(i),2*m(i))/power((1+2*m(i)),(1+2*m(i))))/log10(3));
    a1(i) = (m(i)*180)/(2*m(i)+1);
end
figure
line(H1,a1)
hold on
%% Curve -II
%%limit-1 0<m<0.5
m = 0:0.01:0.5;
[x,y] = size(m);
for i = 1:y
    H21(i) = ((-1)/(1+2*m(i)))*(log10(power(2*m(i),2*m(i))/power((1+2*m(i)),(1+2*m(i))))/log10(3));
    a21(i) = 180/2;
end
line(H21,a21)

%%limit-1 0.5<m<1
m = 0.5:0.01:1;
[x,y] = size(m);
for i = 1:y
    H22(i) = ((-1)/(1+2*m(i)))*(log10(power((2*m(i)-1),(2*m(i)-1))/power((1+2*m(i)),(1+2*m(i))))/log10 (3));
    a22(i) = 180/(1+2*m(i));
end
line(H22,a22)

%% Zone boundaries
axis([0 1 0 90])
% line([x1 x2],[y1 y2])
hold on
plot([0.5 0.5],[0 90])
hold on
line([0.9 0.9],[0 90])
hold on
line([0 0.5],[42.5 42.5])
hold on
line([0 0.5],[47.5 47.5])
hold on
line([0.5 0.9],[40 40])
hold on
line([0.5 0.9],[50 50])
hold on
line([0.9 1.0],[40 40])
hold on
line([0.9 1.0],[60 60])
%-------------------------------------------------------------------------
% zone tagging
text(0.95,75,{'Z1'});
text(0.95,50,{'Z2'});
text(0.95,20,{'Z3'});
text(0.70,75,{'Z4'});
text(0.70,45,{'Z5'});
text(0.70,20,{'Z6'});
text(0.25,75,{'Z7'});
text(0.25,44,{'Z8'});
text(0.25,20,{'Z9'});
%-----------------------------------
grid on;
xlabel('Entropy');
ylabel('Alpha');

%% Volumetric soil moisture given real dielectric constant of soil
epr = 2:0.1:40;
[xx,yy] = size(epr);
for i = 1:yy
mv(i) = 0.01*(-5.3+(2.92.*epr(i))-(0.055.*epr(i).*epr(i))+(0.0004.*epr(i).*epr(i).*epr(i)));  %% Topp et al. 1980 model
end

%% Plotting with different width beta (0<b<90) distribution
b10 = 10;
[alphaC10,HC10] = xbragg(b10);
%% Plotting with varying eps
line(HC10,alphaC10,'Color','red','LineStyle','-','LineWidth',1.2)
text(max(HC10),max(alphaC10),num2str('\beta=10'))

%%
b20 = 20;
[alphaC20,HC20] = xbragg(b20);
%% Plotting with varying eps
line(HC20,alphaC20,'Color','red','LineStyle','-','LineWidth',1.2)
text(max(HC20),max(alphaC20),num2str('\beta=20'))

%%
b40 = 40;
[alphaC40,HC40] = xbragg(b40);
%% Plotting with varying eps
line(HC40,alphaC40,'Color','red','LineStyle','-','LineWidth',1.2)
text(max(HC40),max(alphaC40),num2str('\beta=40'))

%%
b60 = 60;
[alphaC60,HC60] = xbragg(b60);
%% Plotting with varying eps
line(HC60,alphaC60,'Color','red','LineStyle','-','LineWidth',1.2)
text(max(HC60),max(alphaC60),num2str('\beta=60'))

%%
b70 = 70;
[alphaC70,HC70] = xbragg(b70);
%% Plotting with varying eps
line(HC70,alphaC70,'Color','red','LineStyle','-','LineWidth',1.2)
text(max(HC70),max(alphaC70),num2str('\beta=70'))

%%
b80 = 80;
[alphaC80,HC80] = xbragg(b80);
%% Plotting with varying eps
line(HC80,alphaC80,'Color','red','LineStyle','-','LineWidth',1.2)
text(max(HC80),max(alphaC80),num2str('\beta=80'))

%%
b90 = 90;
[alphaC90,HC90] = xbragg(b90);
%% Plotting with varying eps
line(HC90,alphaC90,'Color','red','LineStyle','-','LineWidth',1.2)
text(max(HC90),max(alphaC90),num2str('\beta=90'))

%%-----------------------------------------------------------------------------------
%% Build LUT mv, H, alpha --with stack of different beta
LUT = [mv',HC10',alphaC10'; mv',HC20',alphaC20'; mv',HC40',alphaC40'; mv',HC60',alphaC60'; mv',HC70',alphaC70'; mv',HC80',alphaC80'; mv',HC90',alphaC90'];

%% Saving LUT 
save('F:\Github_Dipankar\SAR-processing-matlab\XBragg_surfacemodel\LUT.mat','LUT')


%% Loading LUT from file
% clear
% load('test.mat')