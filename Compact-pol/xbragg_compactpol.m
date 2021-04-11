%% X-bragg compact-pol model
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
function [alphaC] = xbragg_compactpol(b,epr,inc)
%% Orientation angle -distribution
%% Considering uniform distribution in 2Beta
%b = 35; %degree

%% Relative dielectric constant
% epr = 2:0.1:40; %% real part
eps = complex(epr,0.00);

%% Incidence angle (SAR geometry)
% inc = 35; %in degree
incr = deg2rad(inc);  % in radian



%% Bragg coefficients for horizontal (Rh) and vertical (Rv) polarizations
Rh = (cos(incr) - sqrt(eps - (sin(incr)*sin(incr))))/(cos(incr) + sqrt(eps - (sin(incr)*sin(incr))));
Rv = ((eps*cos(incr)) - sqrt(eps - (sin(incr)*sin(incr))))/((eps*cos(incr)) + sqrt(eps - (sin(incr)*sin(incr))));


%% X-Bragg coefficients
C1 = abs((Rh + Rv))^2;    %% |Rh + Rv|^2;
C2 = (Rh + Rv)*(conj(Rh) - conj(Rv));
C3 = 0.5 * (abs((Rh - Rv))^2);    %% |Rh - Rv|^2;


%% Stokes vector of X-bragg
g0 = 0.5* (C1 + 2*C3);
g1 = C2*(sin(2*b*pi/180)/(2*b*pi/180));
g2 = 0;
g3 = 0.5* (C1 - 2*C3);


%% H-A-Alpha decomposition
alphaC = rad2deg(0.5 * atan(real(g1./g3)));


end




