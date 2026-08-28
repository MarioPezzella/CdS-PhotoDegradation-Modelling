function [C,t] = GH_Cal_CdS_PCTrap_Adim_Denser(xi,mu,nu)
% =========================================================================
% GH_CAL_CDS_PCTRAP_ADIM_Denser
%
% Software and code developed by: MARIO PEZZELLA
%
% This software is associated with the manuscript
%       "Early detection and long-term prediction of pigment degradation in
%        paintings through data-informed mathematical modelling"
%
% Authors:
%       Sara Mattana* and Mario Pezzella*,
%       Francesca Rosi, Marine Cotte, Aldo Romani, Costanza Miliani,
%       Roberto Natalini† and Letizia Monico†
%                         * These authors contributed equally to this work.
%                         † Corresponding authors.
%
% =========================================================================
% This function solves the adimensional CdS degradation model by means of
% a predictor-corrector scheme with the trapezoidal rule.
%
% The experimental data and model-related quantities required by the
% simulation are loaded from Smoothed_Data.mat. If the smoothed data are
% not available, the GH_Data_Smoothing script is executed to generate them.
%
% INPUTS:
%   xi  - Adimensional reaction parameter
%   mu  - Adimensional optical absorption parameter
%   nu  - Ratio between CdS and CdSO4 reference absorptivities
%
% OUTPUT:
%   C - CdS relative concentration profiles at the selected times t
%         corresponding to the experimental samples.
% =========================================================================

% =========================================================================
% Load smoothed data
% =========================================================================
dataPath = fileparts(mfilename('fullpath'));
dataFile = fullfile(dataPath,'Smoothed_Data.mat');
smoothingFile = fullfile(dataPath,'GH_Data_Smoothing.m');
if isfile(dataFile)
    fprintf('Loading smoothed data.\n');
    load(dataFile);
else
    fprintf('Smoothed data not found. Running GH_Data_Smoothing.m.\n');
    run(smoothingFile);
    close all
    load(dataFile);
end
% =========================================================================
% Model parameters
% =========================================================================
L = x_crop_all(end);
c_c_ref = 4.82/144.46;      % CdS reference value [mol/cm^3]
T = (3600*24)*365*30+1;     % Reference time [s] (30 years)
lambda_g = 512;             % CdS band-gap wavelength [nm]
lambda_0 = 380;             % Lower wavelength threshold [nm]
lambda_g0 = lambda_g - lambda_0;
% =========================================================================
% Illumination function
% =========================================================================
I_func_adim = @(x) spline(lambda_adim,I_mean_adim,x);
I = @(x) subplus(I_func_adim(x));
% =========================================================================
% Water profile
% =========================================================================
w = @(x) (1 - wbm/w_ref).*(1 - npsi*x);
% =========================================================================
% CdS and CdSO4 absorptivity functions
% =========================================================================
eps_c_values = -log(Refl_un)/(c_c_ref*L);
eps_ref_c = max(eps_c_values);

eps_c = @(l) subplus(spline(WL,eps_c_values, ...
             lambda_0 + l*lambda_g0))/eps_ref_c;

eps_g_values = eps_c_values/nu;
eps_ref_g = max(eps_g_values);

eps_g = @(l) subplus(spline(WL,eps_g_values, ...
             lambda_0 + l*lambda_g0))/eps_ref_g;

% =========================================================================
% Numerical discretization
% =========================================================================
% The problem is adimensionalized and solved in [0,1]^2,
% with lambda in [0,1].
%
% D_z = Spatial stepsize
% D_t = Temporal stepsize
% D_L = Wavelength stepsize
%
% xi  = Adimensional reaction parameter
% mu  = Adimensional optical absorption parameter
% nu  = Ratio between CdS and CdSO4 reference absorptivities
%
% w   = Adimensionalized humidity function
% e_c = Adimensionalized molar absorptivity function for CdS
% e_g = Adimensionalized molar absorptivity function for CdSO4
% I   = Adimensionalized illumination function
% =========================================================================
D_z=0.006622516556291/100;  D_t=24*3600/T;  D_L=1e-2;
z = 0:D_z:1;    t = 0:D_t:1;    L = 0:D_L:1;
N_t = length(t);    N_z = length(z);
pred = ones(N_z,1); C = ones(N_z,N_t);  

% =========================================================================
% Precomputation of known functions
% =========================================================================
W = w(z);
I_val = I(L);
e_g_val = eps_g(L);
e_diff_val = nu*eps_c(L) - eps_g(L);
% Dimensionless constants
xit = D_t*xi;   muz = mu*D_z;
% Columns represent the different time levels.
% The first column corresponds to the initial condition.

% =========================================================================
% Predictor-corrector time integration
% =========================================================================
for n = 1:N_t-1
    for j = 1:N_z
        % =================================================================
        % Predictor phase
        % =================================================================
        sum1 = j-1;
        sum2 = sum(C(1:j-1,n));
        alpha = exp(-muz * (e_g_val(1:end-1)*sum1 + ...
                e_diff_val(1:end-1)*sum2)) .* I_val(1:end-1);
        integrand_L = (2*alpha)./(1+alpha.^2);
        I_lambda = D_L*sum(integrand_L);
        pred(j) = C(j,n)*exp(-W(j)*xit*I_lambda);
        % =================================================================
        % Corrector phase
        % =================================================================
        sum2_trap = (pred(1) + 2*sum(pred(2:j-1)) + pred(j))*(j ~= 1);
        sum1_trap = (1 + 2*(j > 2)*(j-2) + 1)*(j ~= 1);
        alpha_trap = exp(-muz/2 * (e_g_val*sum1_trap + ...
                     e_diff_val*sum2_trap)) .* I_val;
        integrand_L_trap = (2*alpha_trap)./(1+alpha_trap.^2);
        I_lambda_trap = D_L/2 * (2*sum(integrand_L_trap(2:end-1)) + ...
                        integrand_L_trap(1) + integrand_L_trap(end));
        % =================================================================
        % Previous solution evaluated with the trapezoidal rule
        % =================================================================
        sum2_old_trap = (C(1,n) + 2*sum(C(2:j-1,n)) + C(j,n))*(j ~= 1);
        alpha = exp(-muz/2 * (e_g_val*sum1_trap + ...
                e_diff_val*sum2_old_trap)) .* I_val;
        integrand_L = (2*alpha)./(1+alpha.^2);
        I_lambda = D_L/2 * (2*sum(integrand_L(2:end-1)) + ...
                   integrand_L(1) + integrand_L(end));
        % =================================================================
        % Solution update
        % =================================================================
        C(j,n+1) = C(j,n)* exp(-W(j)*0.5*xit*(I_lambda_trap + I_lambda));
    end
end

end