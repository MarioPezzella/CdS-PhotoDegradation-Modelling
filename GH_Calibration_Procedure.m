% =========================================================================
% CdS MODEL CALIBRATION
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
% This script calibrates the parameters of the integrodifferential model
% for the photochemical degradation of cadmium yellow using preprocessed
% and smoothed experimental CdS concentration profiles.
%
% If calibration results are already available in Calibration_Results.mat,
% they are loaded directly. Otherwise, the script checks for the smoothed
% experimental data. If Smoothed_Data.mat is not available, the smoothing
% procedure is automatically performed by running GH_Data_Smoothing.m.
%
% The calibration combines a global simulated annealing optimization with
% a local fmincon refinement. The objective function is a weighted relative
% squared error that accounts for both spatial and temporal weights.
%
% The script also compares the calibrated model solution with both the
% smoothed experimental profiles and the original y-averaged profiles.
% Absolute calibration errors are plotted together with a 5% reference
% threshold.
% =========================================================================

%% Calibration
close all; clc

% =========================================================================
% Check for existing calibration results
% =========================================================================
resultsFile = fullfile(fileparts(mfilename('fullpath')), ...
                       'Calibration_Results.mat');
if isfile(resultsFile)
    fprintf('Loading previously computed calibration results.\n');
    load(resultsFile);
else
    fprintf('Calibration results not found. Starting model calibration.\n');
    % =====================================================================
    % Load or generate smoothed data
    % =====================================================================
    dataFile = fullfile(fileparts(mfilename('fullpath')), ...
                        'Smoothed_Data.mat');
    if isfile(dataFile)
        fprintf('Loading previously smoothed data.\n');
        load(dataFile);
    else
        fprintf('Smoothed data not found. Running GH_Data_Smoothing.m.\n');
        run(fullfile(fileparts(mfilename('fullpath')), ...
                     'GH_Data_Smoothing.m'));
        close all
        load(dataFile);
    end

    % =====================================================================
    % Experimental data and model parameters
    % =====================================================================
    fitted_value_C1 = min(1,fitted_value_C1);
    y_data = [fitted_value_C2, fitted_value_C1, ...
              fitted_value_D1, fitted_value_B2, fitted_value_B1];
    L = x_crop_all(end);    c_c_ref = 4.82/144.46;
    x0 = [3.5085, 5.4692, 1.3424];
    % =====================================================================
    % Model calibration
    % =====================================================================
    fprintf('Running model calibration.\n');
    [xi_opt,mu_opt,nu_opt,fval,exitflag,output] = Min_Sep(y_data,x0);
    fprintf('Model calibration completed.\n');
    % =====================================================================
    % Model simulation with calibrated parameters
    % =====================================================================
    fprintf('Computing model solution with calibrated parameters.\n');
    Sim = GH_Cal_CdS_PCTrap_Adim(xi_opt,mu_opt,nu_opt);
    % =====================================================================
    % Save calibration results
    % =====================================================================
    save(resultsFile,'-v7.3');
    fprintf('\n');
    fprintf('Calibration results saved in Calibration_Results.mat.\n');
end

%% Calibration error plots
% =========================================================================
% Compare the calibrated model solution with the smoothed and original
% experimental concentration profiles.
%
% The first row reports the absolute error with respect to the smoothed
% data used for calibration. The second row reports the absolute error with
% respect to the original y-averaged experimental data.
%
% A 5% reference threshold is included in both rows to facilitate the
% assessment of the calibration accuracy.
% =========================================================================
y_Sm_data = [fitted_value_C2,fitted_value_C1, ...
             fitted_value_D1,fitted_value_B2,fitted_value_B1];
y_original = [Adim_C2_mean,Adim_C1_mean, ...
              Adim_D1_mean,Adim_B2_mean,Adim_B1_mean];
figure('Position',[100 100 1400 600]);
t = tiledlayout(2,5,'TileSpacing','compact','Padding','compact');
% =========================================================================
% Plot settings
% =========================================================================
colors = {'y','g','r',[0.49,0.18,0.56],'b'};
labels = {'Sample t1','Sample t4','Sample t7','Sample t9','Sample t10'};
threshold = 5e-2;
% =========================================================================
% First row: errors with respect to smoothed data
% =========================================================================
for k = 1:5
    nexttile(k)
    err = abs(Sim(:,k)-y_Sm_data(:,k));
    semilogy(x_crop_all,err,'Color',colors{k},'LineWidth',2);
    hold on
    h = yline(threshold,'k:',LineWidth=2,Label='5% Threshold');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    xlabel('z (depth) [cm]')
    ylabel('$|c^*_i(z_j)-c_{l_j}^{k_i}|$','Interpreter','latex')
    title(labels{k})
    axis square tight
    ylim([1e-6 0.2])
end
% =========================================================================
% Second row: errors with respect to original averaged data
% =========================================================================
for k = 1:5
    nexttile(5+k)
    err_origin = abs(Sim(:,k)-y_original(:,k));
    semilogy(x_crop_all,err_origin,'Color',colors{k},'LineWidth',2);
    hold on
    h = yline(threshold,'k:',LineWidth=2,Label='5% Threshold');
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    xlabel('z (depth) [cm]')
    ylabel('$|\bar{c}(z_j,\tau_i)-c_{l_j}^{k_i}|$','Interpreter','latex')
    title(labels{k})
    axis square tight
    ylim([1e-5 0.2])
end
% =========================================================================
% Row titles
% =========================================================================
annotation('textbox',[0 0.95 1 0.05], 'String', ...
    'Calibration Absolute Errors with respect to Smoothed Data', ...
    'EdgeColor','none','HorizontalAlignment','center', ...
    'FontWeight','bold','FontSize',12);

annotation('textbox',[0 0.47 1 0.05], 'String', ...
    'Calibration Absolute Errors with respect to Cropped Averaged Data',...
    'EdgeColor','none','HorizontalAlignment','center', ...
    'FontWeight','bold','FontSize',12);

% =========================================================================
% Calibration function
% =========================================================================
function [xi_opt,mu_opt,nu_opt,fval,exitflag,output] = Min_Sep(y_data,x0)

    % =====================================================================
    % Spatial weights
    % =====================================================================
    m = size(y_data,1);
    n_times = size(y_data,2);
    alpha = [0.01,1,0.6,0.5,0.01];
    weights_space = zeros(m,n_times);
    for j = 1:n_times
        weights_space(:,j) = j*exp(-alpha(j)*(1:m)*0.25)';
    end
    % =====================================================================
    % Temporal weights
    % =====================================================================
    weights_time = [0.5,1,1.75,1,3];
    % =====================================================================
    % Objective function
    % =====================================================================
    objFun = @(x) weighted_objFun(x,y_data, weights_space,weights_time);
    % =====================================================================
    % Parameter bounds
    % =====================================================================
    x_lb = [0,0,0];
    L = 0.003775000000000;
    c_c_ref = 0.033365637546726;
    x_ub = [1e4, 2*L*c_c_ref*2.171393420017231e+04, 10];
    % =====================================================================
    % Simulated annealing
    % =====================================================================
    fprintf('  Global optimization: simulated annealing.\n');
    options_sa = optimoptions('simulannealbnd', ...
        'Display','iter', ...
        'MaxFunctionEvaluations',5e4, ...
        'MaxIterations',5e4, ...
        'FunctionTolerance',1e-13, ...
        'PlotFcn',{@saplotbestf,@saplottemperature, ...
                   @saplotf,@saplotstopping});
    [x_sa,~,~,~] = simulannealbnd(objFun,x0,x_lb,x_ub,options_sa);

    % =====================================================================
    % Local refinement with fmincon
    % =====================================================================
    fprintf('  Local refinement: fmincon.\n');
    options_local = optimoptions('fmincon', ...
        'Display','iter', ...
        'Algorithm','interior-point', ...
        'MaxFunctionEvaluations',1e4, ...
        'OptimalityTolerance',1e-13, ...
        'StepTolerance',1e-13, ...
        'UseParallel',true);
    problem = createOptimProblem('fmincon', ...
        'objective',objFun, ...
        'x0',x_sa, ...
        'lb',x_lb, ...
        'ub',x_ub, ...
        'options',options_local);
    [x_opt,fval,exitflag,output] = fmincon(problem);
    xi_opt = x_opt(1);
    mu_opt = x_opt(2);
    nu_opt = x_opt(3);
end

% =========================================================================
% Weighted objective function
% =========================================================================
function obj = weighted_objFun(x,y_data,weights_space,weights_time)
    model_output = GH_Cal_CdS_PCTrap_Adim(x(1),x(2),x(3));
    err_sq = (model_output-y_data).^2;
    dat_sq = y_data.^2;
    err_weighted = err_sq.*weights_space;
    dat_weighted = dat_sq.*weights_space;
    sum_err_per_time = sum(err_weighted,1);
    sum_data_per_time = sum(dat_weighted,1);
    num = sum(weights_time.*sum_err_per_time);
    den = sum(weights_time.*sum_data_per_time);
    obj = num/den;
end