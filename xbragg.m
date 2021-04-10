%% X-bragg model
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
function [alphaC,HC] = xbragg(b)
%% Orientation angle -distribution
%% Considering uniform distribution in 2Beta
%b = 35; %degree

%% Relative dielectric constant
epr = 2:0.1:40; %% real part
eps = complex(epr,0.00);
[x,y] = size(epr);

%% Incidence angle (SAR geometry)
inc = 35; %in degree
incr = deg2rad(inc);  % in radian
    
for i = 1:y
   
    %% Bragg coefficients for horizontal (Rh) and vertical (Rv) polarizations
    Rh = (cos(incr) - sqrt(eps(i) - (sin(incr)*sin(incr))))/(cos(incr) + sqrt(eps(i) - (sin(incr)*sin(incr))));
    Rv = ((eps(i)*cos(incr)) - sqrt(eps(i) - (sin(incr)*sin(incr))))/((eps(i)*cos(incr)) + sqrt(eps(i) - (sin(incr)*sin(incr))));
    
    
    %% X-Bragg coefficients
    C1 = abs((Rh - Rv))^2;    %% |Rh + Rv|^2;
    C2 = (Rh + Rv)*(conj(Rh) - conj(Rv));
    C3 = abs((Rh + Rv))^2;    %% |Rh - Rv|^2;
    
    
    %% T3 elements of X-bragg
    t11 = 0.5*C1;
    t12 = 0.5*C2*(sin(2*b*pi/180)/(2*b*pi/180));
    t13 = 0;
    t21 = conj(t12);
    t22 = 0.25*C3*(1+(sin(4*b*pi/180)/(4*b*pi/180)));
    t23 = 0;
    t31 = 0;
    t32 = 0;
    t33 = 0.25*C3*(1-(sin(4*b*pi/180)/(4*b*pi/180)));
    
    T_XB = [t11 t12 t13;t21 t22 t23;t31 t32 t33];
    
    %% H-A-Alpha decomposition
    [alphaC(i), HC(i), AC(i)] = HAalphadecomp(T_XB);
    
    
end

end


