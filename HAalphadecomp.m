%% H-A-alpha decomposition 
% Ref: CLOUDE, S. R. & POTTIER, E.,‘An Entropy Based Classification Scheme for Land Applications of Polarimetric SAR’, IEEE Transactions on Geoscience and Remote Sensing, vol. 35, no. 1, pp. 68-78, 1997.
% @author: Dr. Dipankar Mandal and Dr. Avik Bhattacharya
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
function [alpha,H, A] = HAalphadecomp(T3)
%% H-A-Alpha decomposition
[evec_v, eval] = eig(T3);

%% Eigenvalues

eval_diag = (sort(diag(eval)))';

if (eval_diag(1) < 0)
    eval_diag(1) = 0;
end

if (eval_diag(2) < 0)
    eval_diag(2) = 0;
end

if (eval_diag(3) < 0)
    eval_diag(3) = 0;
end

%Lambda 1
eval_norm1 = (eval_diag(3))./(eval_diag(1) + eval_diag(2) + eval_diag(3));

eval_norm1(eval_norm1 < 0) = 0;
eval_norm1(eval_norm1 > 1) = 1;

%Lambda 2
eval_norm2 = (eval_diag(2))./(eval_diag(1) + eval_diag(2) + eval_diag(3));

eval_norm2(eval_norm2 < 0) = 0;
eval_norm2(eval_norm2 > 1) = 1;

%Lambda 3
eval_norm3 = (eval_diag(1))./(eval_diag(1) + eval_diag(2) + eval_diag(3));

eval_norm3(eval_norm3 < 0) = 0;
eval_norm3(eval_norm3 > 1) = 1;


%% Eigenvectors

%Alpha 1
eig_vec_r1 = real(evec_v(1,3));
eig_vec_c1 = imag(evec_v(1,3));
alpha1 = acos(sqrt(eig_vec_r1*eig_vec_r1 + eig_vec_c1*eig_vec_c1));

%Alpha 2
eig_vec_r2 = real(evec_v(1,2));
eig_vec_c2 = imag(evec_v(1,2));
alpha2 = acos(sqrt(eig_vec_r2*eig_vec_r2 + eig_vec_c2*eig_vec_c2));

%Alpha 3
eig_vec_r3 = real(evec_v(1,1));
eig_vec_c3 = imag(evec_v(1,1));
alpha3 = acos(sqrt(eig_vec_r3*eig_vec_r3 + eig_vec_c3*eig_vec_c3));


%Cloude Alpha
alpha = (eval_norm1*alpha1*180./pi + eval_norm2*alpha2*180./pi + ...
    eval_norm3*alpha3*180./pi);

%Entropy
H = -eval_norm1*log10(eval_norm1)./log10(3) - ...
    eval_norm2*log10(eval_norm2)./log10(3) - ...
    eval_norm3*log10(eval_norm3)./log10(3);

A = (eval_norm2 - eval_norm3)./(eval_norm2 + eval_norm3);

end