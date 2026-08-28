% =========================================================================
% LONG-TERM CdS DEGRADATION SIMULATION
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
% This script computes the long-term evolution of the dimensionless CdS
% concentration using the calibrated parameters of the integrodifferential
% model.
%
% The dimensionless model is solved up to 30 years with increased spatial
% resolution. The dimensionless solution is then rescaled for each
% experimental sample to match the CdS concentration at the corresponding
% measurement time.
%
% The script also reconstructs dense two-dimensional concentration maps by
% interpolation of the experimental data along the y and z directions.
% These maps are used to visualize the predicted CdS concentration after
% 30 years of exposure.
%
% In addition, the script performs a quantitative analysis of CdSO4
% degradation for each experimental sample. The degradation front is
% identified from the mean relative CdSO4 concentration using a 5%
% threshold. The mean CdSO4 concentration within the degraded region and
% the total CdSO4 content over the entire sample are then computed at
% one and 30 years. The corresponding percentage increase in the total
% CdSO4 content is also reported.
%
% These quantities provide the degradation-front locations and CdSO4
% concentration values reported in Table 2 of the manuscript.
%
% Only the per-sample rescaling factors are cached, in
% Long_Simulations.mat. They are the sole quantity the 30-year integration
% is needed for and weigh about 1.2 MB, so the file can be distributed
% together with the code: a reviewer can reproduce every figure and every
% number of Table 2 without repeating the long simulation. If the file is
% absent, the simulation is performed and the factors are saved.
%
% The dense experimental maps are deliberately NOT cached. They are pure
% interpolation of the raw data, they amount to some 4.6 GB, and they are
% rebuilt at every run.
% =========================================================================

% =========================================================================
% Load calibration results
% =========================================================================
dataPath = fileparts(mfilename('fullpath'));
dataFile = fullfile(dataPath,'Calibration_Results.mat');
calibrationFile = fullfile(dataPath,'GH_Calibration_Procedure.m');
longSimulationFile = fullfile(dataPath,'Long_Simulations.mat');

if isfile(dataFile)
    fprintf('Loading calibration results.\n');
    load(dataFile);
else
    fprintf(['Calibration results not found. Running ', ...
        'GH_Calibration_Procedure.m.\n']);
    run(calibrationFile);
    close all
    load(dataFile);
end

% =========================================================================
% The MAT files are saved with the whole workspace, so they carry the path
% variables of the machine that produced them. Whenever a MAT file is
% shared rather than regenerated locally, those stale absolute paths would
% silently replace the ones computed above. They are therefore rebuilt here.
% =========================================================================
dataPath           = fileparts(mfilename('fullpath'));
dataFile           = fullfile(dataPath,'Calibration_Results.mat');
calibrationFile    = fullfile(dataPath,'GH_Calibration_Procedure.m');
longSimulationFile = fullfile(dataPath,'Long_Simulations.mat');


% =========================================================================
% Long-term dimensionless simulation
% =========================================================================
% Only the per-sample rescaling factors are cached. They are the sole
% quantity that the 30-year integration is needed for, they weigh about
% 1.2 MB in total, and everything else downstream (the dense experimental
% maps) is rebuilt from the raw data at every run in a matter of minutes.
%
% The variables are listed explicitly both on save and on load, so that the
% cache cannot carry anything else into the caller's workspace.
% =========================================================================
cachedVars = { ...
    'C_1Y_norm_C2','C_1Y_norm_C1','C_1Y_norm_D1', ...
    'C_1Y_norm_B2','C_1Y_norm_B1', ...
    'C_30Y_norm_C2','C_30Y_norm_C1','C_30Y_norm_D1', ...
    'C_30Y_norm_B2','C_30Y_norm_B1'};

if isfile(longSimulationFile)
    fprintf('Loading previously computed rescaling factors.\n');
    load(longSimulationFile,cachedVars{:});
else
    fprintf('Long-term simulation not found. Computing solution up to 30 years.\n');

    warning(['The long-term simulation is performed with high spatial ', ...
        'resolution. Consequently, the computation may require a ', ...
        'significant amount of time. Please wait until the simulation ', ...
        'is completed.']);

    % =====================================================================
    % Dimensionless 1D simulation with higher spatial resolution
    % =====================================================================
    [C_LONG,~] = GH_Cal_CdS_PCTrap_Adim_Denser(xi_opt,mu_opt,nu_opt);

    % =====================================================================
    % The dimensionless long-term solution is rescaled separately for each
    % sample. The normalization matches the model concentration at the
    % experimental measurement time of each sample.
    %
    % Columns of C_LONG are one day apart, so column 366 is the one-year
    % level and the last column the thirty-year one. Columns 2, 4, 8, 15
    % and 29 are the measurement times of t1, t4, t7, t9 and t10.
    %
    % This is the only place where C_LONG is needed: the rescaling is
    % performed here so that the 1.3 GB solution never has to be stored.
    % =====================================================================
    % After 1 Year
    C_1Y_norm_C2 = min(C_LONG(:,366)./C_LONG(1,2),1);
    C_1Y_norm_C1 = min(C_LONG(:,366)./C_LONG(1,4),1);
    C_1Y_norm_D1 = min(C_LONG(:,366)./C_LONG(1,8),1);
    C_1Y_norm_B2 = min(C_LONG(:,366)./C_LONG(1,15),1);
    C_1Y_norm_B1 = min(C_LONG(:,366)./C_LONG(1,29),1);
    % After 30 Years
    C_30Y_norm_C2 = min(C_LONG(:,end)./C_LONG(1,2),1);
    C_30Y_norm_C1 = min(C_LONG(:,end)./C_LONG(1,4),1);
    C_30Y_norm_D1 = min(C_LONG(:,end)./C_LONG(1,8),1);
    C_30Y_norm_B2 = min(C_LONG(:,end)./C_LONG(1,15),1);
    C_30Y_norm_B1 = min(C_LONG(:,end)./C_LONG(1,29),1);

    % =====================================================================
    % Save the rescaling factors
    % =====================================================================
    save(longSimulationFile,cachedVars{:},'-v7.3');

    fprintf('Rescaling factors saved in Long_Simulations.mat.\n');

    % The 30-year solution is no longer needed. Releasing it here keeps it
    % from coexisting with the dense maps allocated below.
    clear C_LONG
end

%==========================================================================
% The experimental concentration maps are interpolated along the y
% direction to obtain a denser representation of the measured profiles.
%==========================================================================
enlarge = 500;
y_C1_dense = y_C1(1):mean(diff(y_C1)/enlarge):y_C1(end);
y_C2_dense = y_C2(1):mean(diff(y_C2)/enlarge):y_C2(end);
y_D1_dense = y_D1(1):mean(diff(y_D1)/enlarge):y_D1(end);
y_B2_dense = y_B2(1):mean(diff(y_B2)/enlarge):y_B2(end);
y_B1_dense = y_B1(1):mean(diff(y_B1)/enlarge):y_B1(end);
z_all = numel(x_crop_all);

% Preallocate matrices
CdSC1_DS_dense = zeros(z_all,length(y_C1_dense));
CdSC2_DS_dense = zeros(z_all,length(y_C2_dense));
CdSD1_DS_dense = zeros(z_all,length(y_D1_dense));
CdSB1_DS_dense = zeros(z_all,length(y_B1_dense));
CdSB2_DS_dense = zeros(z_all,length(y_B2_dense));

for z = 1:z_all
    CdSC1_DS_dense(z,:) = min( interp1(y_C1,CdSC1_crop(z,:),y_C1_dense, ...
        'spline','extrap'),1);

    CdSC2_DS_dense(z,:) = min( interp1(y_C2,CdSC2_crop(z,:),y_C2_dense, ...
        'spline','extrap'),1);

    CdSD1_DS_dense(z,:) = min( interp1(y_D1,CdSD1_crop(z,:),y_D1_dense, ...
        'spline','extrap'),1);

    CdSB1_DS_dense(z,:) = min( interp1(y_B1,CdSB1_crop(z,:),y_B1_dense, ...
        'spline','extrap'),1);

    CdSB2_DS_dense(z,:) = min( interp1(y_B2,CdSB2_crop(z,:),y_B2_dense, ...
        'spline','extrap'),1);
end

%==========================================================================
% The experimental concentration maps are further interpolated along
% the z (depth) direction to obtain a denser spatial representation.
%==========================================================================
z_dense = (0:0.006622516556291/100:1).*x_crop_all(end);

% Preallocate matrices
CdSC1_DS_dense_MATX = zeros(length(z_dense),length(y_C1_dense));
CdSC2_DS_dense_MATX = zeros(length(z_dense),length(y_C2_dense));
CdSD1_DS_dense_MATX = zeros(length(z_dense),length(y_D1_dense));
CdSB1_DS_dense_MATX = zeros(length(z_dense),length(y_B1_dense));
CdSB2_DS_dense_MATX = zeros(length(z_dense),length(y_B2_dense));

for j = 1:length(y_C1_dense)
    CdSC1_DS_dense_MATX(:,j) = min( ...
        interp1(x_crop_all,CdSC1_DS_dense(:,j), ...
        z_dense,'spline','extrap'),1);
end

for j = 1:length(y_C2_dense)
    CdSC2_DS_dense_MATX(:,j) = min( ...
        interp1(x_crop_all,CdSC2_DS_dense(:,j), ...
        z_dense,'spline','extrap'),1);
end

for j = 1:length(y_B2_dense)
    CdSB2_DS_dense_MATX(:,j) = min( ...
        interp1(x_crop_all,CdSB2_DS_dense(:,j), ...
        z_dense,'spline','extrap'),1);
end

for j = 1:length(y_B1_dense)
    CdSB1_DS_dense_MATX(:,j) = min( ...
        interp1(x_crop_all,CdSB1_DS_dense(:,j), ...
        z_dense,'spline','extrap'),1);
end

for j = 1:length(y_D1_dense)
    CdSD1_DS_dense_MATX(:,j) = min( ...
        interp1(x_crop_all,CdSD1_DS_dense(:,j), ...
        z_dense,'spline','extrap'),1);
end

% =========================================================================
% PLOTS
% =========================================================================
% The first row reports the experimental CdS concentration maps at the
% measurements time. The second row reports the corresponding CdS
% concentration prediction after 30 years of exposure, obtained by
% multiplying the experimental spatial distribution by the rescaled
% dimensionless model solution.
% =========================================================================
figure()
t = tiledlayout(2,5,'TileSpacing','compact','Padding','compact');
% =========================================================================
% First row: experimental concentration maps
% =========================================================================
nexttile
imagesc(y_C2_dense*1e4,z_dense*1e4,CdSC2_DS_dense_MATX)
title('Sample t1')
colormap(jet(256))
ylabel('z (depth) [µm]'); xlabel('y [µm]')
clim([0.53 1])
axis square
ax = gca;
ax.FontSize = 11;

nexttile
imagesc(y_C1_dense*1e4,z_dense*1e4,CdSC1_DS_dense_MATX)
title('Sample t4')
colormap(jet(256))
ylabel('z (depth) [µm]'); xlabel('y [µm]')
clim([0.53 1])
axis square
ax = gca;
ax.FontSize = 11;

nexttile
imagesc(y_D1_dense*1e4,z_dense*1e4,CdSD1_DS_dense_MATX)
title('Sample t7')
colormap(jet(256))
ylabel('z (depth) [µm]'); xlabel('y [µm]')
clim([0.53 1])
axis square
ax = gca;
ax.FontSize = 11;

nexttile
imagesc(y_B2_dense*1e4,z_dense*1e4,CdSB2_DS_dense_MATX)
title('Sample t9')
colormap(jet(256))
ylabel('z (depth) [µm]'); xlabel('y [µm]')
clim([0.53 1])
axis square
ax = gca;
ax.FontSize = 11;

nexttile
imagesc(y_B1_dense*1e4,z_dense*1e4,CdSB1_DS_dense_MATX)
title('Sample t10')
colormap(jet(256))
ylabel('z (depth) [µm]'); xlabel('y [µm]')
clim([0.53 1])
axis square
ax = gca;
ax.FontSize = 11;

cb = colorbar('southoutside');
cb.Position = [0.06328125,0.533632286995516, ...
    0.88125,0.013841554559044];
cb.FontSize = 11;

% =========================================================================
% Second row: 30-year prediction
% =========================================================================
nexttile
imagesc(y_C2_dense*1e4,z_dense*1e4, ...
    CdSC2_DS_dense_MATX .* C_30Y_norm_C2(:,1))
title('Sample t1')
colormap(jet(256))
ylabel('z (depth) [µm]'); xlabel('y [µm]')
clim([0.53 1])
axis square
ax = gca;
ax.FontSize = 11;

nexttile
imagesc(y_C1_dense*1e4,z_dense*1e4, ...
    CdSC1_DS_dense_MATX .* C_30Y_norm_C1(:,1))
title('Sample t4')
colormap(jet(256))
ylabel('z (depth) [µm]'); xlabel('y [µm]')
clim([0.53 1])
axis square
ax = gca;
ax.FontSize = 11;

nexttile
imagesc(y_D1_dense*1e4,z_dense*1e4, ...
    CdSD1_DS_dense_MATX .* C_30Y_norm_D1(:,1))
title('Sample t7')
colormap(jet(256))
ylabel('z (depth) [µm]'); xlabel('y [µm]')
clim([0.53 1])
axis square
ax = gca;
ax.FontSize = 11;

nexttile
imagesc(y_B2_dense*1e4,z_dense*1e4, ...
    CdSB2_DS_dense_MATX .* C_30Y_norm_B2(:,1))
title('Sample t9')
colormap(jet(256))
ylabel('z (depth) [µm]'); xlabel('y [µm]')
clim([0.53 1])
axis square
ax = gca;
ax.FontSize = 11;

nexttile
imagesc(y_B1_dense*1e4,z_dense*1e4, ...
    CdSB1_DS_dense_MATX .* C_30Y_norm_B1(:,1))
title('Sample t10')
colormap(jet(256))
ylabel('z (depth) [µm]'); xlabel('y [µm]')
clim([0.53 1])
axis square
ax = gca;
ax.FontSize = 11;
% =========================================================================
% Global titles
% =========================================================================
annotation('textbox',[0,0.931494768310912,1,0.05], ...
    'String', ...
    'Simulated Relative CdS Concentrations - Samples Measurements (a)', ...
    'EdgeColor','none', ...
    'HorizontalAlignment','center', ...
    'FontWeight','bold', ...
    'FontSize',14);
annotation('textbox',[0,0.429641255605381,1,0.05], ...
    'String', ...
    'Relative CdS Concentrations Prediction After 30 Years of Exposure (b)', ...
    'EdgeColor','none', ...
    'HorizontalAlignment','center', ...
    'FontWeight','bold', ...
    'FontSize',14);


% =========================================================================
% Quantitative analysis of CdSO4 degradation
% =========================================================================
% For each sample, the relative CdSO4 concentration is computed as
% 1 - C_CdS. The degradation front is defined as the first depth at which
% the mean relative CdSO4 concentration falls below 5%.
%
% The mean CdSO4 concentration is computed over the region from the
% surface to the degradation front. The total CdSO4 content is computed
% over the entire two-dimensional sample map.
%
% The depth coordinate z_dense is expressed in cm and is converted to
% micrometres for the reported degradation-front positions.
% =========================================================================

% =========================================================================
% The rescaling factor is constant along y, so the two quantities needed
% below can be obtained from reductions of the measured maps rather than
% from the rescaled maps themselves:
%
%   mean_y( 1 - M.*v )  =  1 - mean_y(M).*v
%   sum_all( 1 - M.*v ) =  numel(M) - v' * sum_y(M)
%
% Both identities are exact. Working on the two reductions avoids building
% ten further arrays the size of the dense maps, which would otherwise take
% about 9 GB on top of the maps already in memory.
% =========================================================================

% Cell arrays of references: MATLAB copies on write, so these cost nothing
maps = {CdSC2_DS_dense_MATX, CdSC1_DS_dense_MATX, CdSD1_DS_dense_MATX, ...
        CdSB2_DS_dense_MATX, CdSB1_DS_dense_MATX};

norm_1Y = {C_1Y_norm_C2, C_1Y_norm_C1, C_1Y_norm_D1, ...
           C_1Y_norm_B2, C_1Y_norm_B1};

norm_30Y = {C_30Y_norm_C2, C_30Y_norm_C1, C_30Y_norm_D1, ...
            C_30Y_norm_B2, C_30Y_norm_B1};

sampleNames = {'t1','t4','t7','t9','t10'};

nSamples = numel(sampleNames);

front_1Y_um = zeros(nSamples,1);
front_30Y_um = zeros(nSamples,1);

mean_front_1Y = zeros(nSamples,1);
mean_front_30Y = zeros(nSamples,1);

total_CdSO4_1Y = zeros(nSamples,1);
total_CdSO4_30Y = zeros(nSamples,1);

ratio_CdSO4 = zeros(nSamples,1);
percentage_increase = zeros(nSamples,1);

for k = 1:nSamples

    % Reductions of the measured map along y, computed once per sample
    mean_CdS_y = mean(maps{k},2);
    sum_CdS_y  = sum(maps{k},2);
    nElements  = numel(maps{k});

    v_1Y  = norm_1Y{k}(:,1);
    v_30Y = norm_30Y{k}(:,1);

    % Mean relative CdSO4 concentration along the y direction
    mean_CdSO4_1Y = 1 - mean_CdS_y.*v_1Y;
    mean_CdSO4_30Y = 1 - mean_CdS_y.*v_30Y;

    % Degradation front defined by the 5% threshold
    idx_front_1Y = find(mean_CdSO4_1Y < 0.05,1,'first');
    idx_front_30Y = find(mean_CdSO4_30Y < 0.05,1,'first');

    % Front location converted from cm to micrometres
    if isempty(idx_front_1Y)
        front_1Y_um(k) = NaN;
        mean_front_1Y(k) = NaN;
    else
        front_1Y_um(k) = z_dense(idx_front_1Y)*1e4;
        mean_front_1Y(k) = mean(mean_CdSO4_1Y(1:idx_front_1Y));
    end

    if isempty(idx_front_30Y)
        front_30Y_um(k) = NaN;
        mean_front_30Y(k) = NaN;
    else
        front_30Y_um(k) = z_dense(idx_front_30Y)*1e4;
        mean_front_30Y(k) = mean(mean_CdSO4_30Y(1:idx_front_30Y));
    end

    % Total relative CdSO4 content over the entire sample
    total_CdSO4_1Y(k) = nElements - sum(sum_CdS_y.*v_1Y);
    total_CdSO4_30Y(k) = nElements - sum(sum_CdS_y.*v_30Y);

    % Ratio and percentage increase of the total CdSO4 content
    ratio_CdSO4(k) = total_CdSO4_30Y(k)/total_CdSO4_1Y(k);

    percentage_increase(k) = ...
        (total_CdSO4_30Y(k)-total_CdSO4_1Y(k)) / ...
        total_CdSO4_1Y(k) * 100;
end


% =========================================================================
% Display summary table
% =========================================================================
% The degradation-front locations are reported in micrometres.
% Concentrations are displayed in scientific notation with four decimal
% places.
% =========================================================================

fprintf('\n');
fprintf('==========================================================================\n');
fprintf('CdSO4 DEGRADATION SUMMARY\n');
fprintf('Degradation front locations are reported in micrometres (µm).\n');
fprintf('==========================================================================\n');

fprintf('%-8s %-18s %-18s %-22s %-22s %-18s\n', ...
    'Sample', ...
    'Front (1 year) [µm]', ...
    'Front (30 years) [µm]', ...
    'Mean CdSO4 (1 year)', ...
    'Mean CdSO4 (30 years)', ...
    'Increase [%]');

fprintf('--------------------------------------------------------------------------\n');

for k = 1:nSamples
    fprintf('%-8s %18.4f %18.4f %22.4e %22.4e %18.2f\n', ...
        sampleNames{k}, ...
        front_1Y_um(k), ...
        front_30Y_um(k), ...
        mean_front_1Y(k), ...
        mean_front_30Y(k), ...
        percentage_increase(k));
end

fprintf('==========================================================================\n');