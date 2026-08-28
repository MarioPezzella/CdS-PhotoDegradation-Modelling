% =========================================================================
% CdS DATA PROCESSING AND VISUALIZATION
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
% This script processes spatially resolved CdS concentration data obtained
% from synchrotron measurements of painted samples. The data are reshaped,
% cropped to the common positive concentration region, and saved as
% preprocessed data in PreProcessed_Data.mat. If the preprocessed data are
% already available, they are loaded directly to avoid repeating the data
% processing.
%
% The script also processes the experimental lamp emission spectra and
% UV-Vis reflectance data used in the mathematical model. The lamp spectra
% are truncated at the CdS band-gap wavelength, averaged, and
% adimensionalized to obtain the illumination profile. The reflectance data
% are used to compute the CdS molar absorptivity profile.
%
% The script also generates figures for each sample. Each figure contains
% the raw CdS relative concentration map and the corresponding y-averaged
% concentration profile with its standard deviation.
% =========================================================================

%% Data elaboration and preprocessing
clc; close all
% All paths are resolved relative to the location of this script, so that
% the pipeline behaves identically whether the script is executed directly
% or invoked through "run" by a downstream script.
scriptPath  = fileparts(mfilename('fullpath'));
rawDataPath = fullfile(scriptPath,'CdS_Data');
dataFile    = fullfile(scriptPath,'PreProcessed_Data.mat');
if isfile(dataFile)
    fprintf('Loading previously elaborated data.\n');
    load(dataFile);

else
    fprintf('Elaborated data not found. Processing raw data.\n');

    % =====================================================================
    % Sample t1 (C2)                                                      
    % =====================================================================
    MC2 = readmatrix(fullfile(rawDataPath,'Sample_t1.csv'));
    DatiC2 = MC2(:,3:4);
    CdSO4C2 = reshape(DatiC2(:,1),194,14);
    CdSC2 = reshape(DatiC2(:,2),194,14);
    % x is the depth, y the horizontal space
    dx = 0.25e-4;    % cm (0.25 micrometers)
    dy = 1e-4;       % cm (1 micrometer)
    [nRows,nCols] = size(CdSC2);
    x_C2 = dx*(0:nRows-1);
    y_C2 = dy*(0:nCols-1);
    % Crop the matrix to the deepest common positive region
    index_first_positive = arrayfun(@(j) ...
        find(CdSC2(:,j)>0,1,'first'),1:nCols);
    Pos_idx = max(index_first_positive);
    NRows_crop = nRows-Pos_idx+1;
    CdSC2_crop = ones(NRows_crop,nCols);
    x_crop_C2 = dx*(0:NRows_crop-1);
    for j = 1:nCols
        idx = index_first_positive(j);
        CdSC2_crop(:,j) = CdSC2(idx:NRows_crop+idx-1,j);
    end

    % =====================================================================
    % Sample t4 (C1)                                                      
    % =====================================================================
    MC = readmatrix(fullfile(rawDataPath,'Sample_t4.csv'));
    DatiC1 = MC(:,3:4);
    CdSO4C1 = reshape(DatiC1(:,1),284,17);
    CdSC1 = reshape(DatiC1(:,2),284,17);
    [nRows,nCols] = size(CdSC1);
    x_C1 = dx*(0:nRows-1);
    y_C1 = dy*(0:nCols-1);
    % Crop the matrix to the deepest common positive region
    index_first_positive = arrayfun(@(j) ...
        find(CdSC1(:,j)>0,1,'first'),1:nCols);
    Pos_idx = max(index_first_positive);
    NRows_crop = nRows-Pos_idx+1;
    CdSC1_crop = ones(NRows_crop,nCols);
    x_crop_C = dx*(0:NRows_crop-1);
    for j = 1:nCols
        idx = index_first_positive(j);
        CdSC1_crop(:,j) = CdSC1(idx:NRows_crop+idx-1,j);
    end

    % =====================================================================
    % Sample t7 (D1)                                                      
    % =====================================================================
    MD = readmatrix(fullfile(rawDataPath,'Sample_t7.csv'));
    DatiD1 = MD(:,4:5);
    CdSO4D1 = reshape(DatiD1(:,1),156,15);
    CdSD1 = reshape(DatiD1(:,2),156,15);
    [nRows,nCols] = size(CdSD1);
    x_D1 = dx*(0:nRows-1);
    y_D1 = dy*(0:nCols-1);
    % Crop the matrix to the deepest common positive region
    index_first_positive = arrayfun(@(j) ...
        find(CdSD1(:,j)>0,1,'first'),1:nCols);
    Pos_idx = max(index_first_positive);
    NRows_crop = nRows-Pos_idx+1;
    CdSD1_crop = ones(NRows_crop,nCols);
    x_crop_D = dx*(0:NRows_crop-1);
    for j = 1:nCols
        idx = index_first_positive(j);
        CdSD1_crop(:,j) = CdSD1(idx:NRows_crop+idx-1,j);
    end

    % =====================================================================
    % Sample t9 (B2)                                                      
    % =====================================================================
    MB2 = readmatrix(fullfile(rawDataPath,'Sample_t9.csv'));
    DatiB2 = MB2(:,3:4);
    CdSO4B2 = reshape(DatiB2(:,1),187,20);
    CdSB2 = reshape(DatiB2(:,2),187,20);
    [nRows,nCols] = size(CdSB2);
    x_B2 = dx*(0:nRows-1);
    y_B2 = dy*(0:nCols-1);
    % Crop the matrix to the deepest common positive region
    index_first_positive = arrayfun(@(j) ...
        find(CdSB2(:,j)>0,1,'first'),1:nCols);
    Pos_idx = max(index_first_positive);
    NRows_crop = nRows-Pos_idx+1;
    CdSB2_crop = ones(NRows_crop,nCols);
    x_crop_B2 = dx*(0:NRows_crop-1);
    for j = 1:nCols
        idx = index_first_positive(j);
        CdSB2_crop(:,j) = CdSB2(idx:NRows_crop+idx-1,j);
    end

    % =====================================================================
    % Sample t10 (B1)                                                     
    % =====================================================================
    MB = readmatrix(fullfile(rawDataPath,'Sample_t10.csv'));
    DatiB1 = MB(:,3:4);
    CdSO4B1 = reshape(DatiB1(:,1),234,15);
    CdSB1 = reshape(DatiB1(:,2),234,15);
    [nRows,nCols] = size(CdSB1);
    x_B1 = dx*(0:nRows-1);
    y_B1 = dy*(0:nCols-1);
    % Crop the matrix to the deepest common positive region
    index_first_positive = arrayfun(@(j) ...
        find(CdSB1(:,j)>0,1,'first'),1:nCols);
    Pos_idx = max(index_first_positive);
    NRows_crop = nRows-Pos_idx+1;
    CdSB1_crop = ones(NRows_crop,nCols);
    x_crop_B = dx*(0:NRows_crop-1);
    for j = 1:nCols
        idx = index_first_positive(j);
        CdSB1_crop(:,j) = CdSB1(idx:NRows_crop+idx-1,j);
    end
    
    % =====================================================================
    % Computation of y-Mean values and standard deviations                                     
    % =====================================================================
    % Sample t1 (C2)
    CdSC2_crop_mean = mean(CdSC2_crop,2);
    sigma_C2_std = std(CdSC2_crop,0,2);
    x_crop_C2_STD = [x_crop_C2 x_crop_C2(end:-1:1)];
    CdSC2_crop_CONF = [CdSC2_crop_mean+sigma_C2_std; ...
                      CdSC2_crop_mean(end:-1:1)-sigma_C2_std(end:-1:1)];
    % Sample t4 (C1)
    CdSC1_crop_mean = mean(CdSC1_crop,2);
    sigma_C1_std = std(CdSC1_crop,0,2);
    x_crop_C_STD = [x_crop_C x_crop_C(end:-1:1)];
    CdSC1_crop_CONF = [CdSC1_crop_mean+sigma_C1_std; ...
                       CdSC1_crop_mean(end:-1:1)-sigma_C1_std(end:-1:1)];
    % Sample t7 (D1)
    CdSD1_crop_mean = mean(CdSD1_crop,2);
    sigma_D1_std = std(CdSD1_crop,0,2);
    x_crop_D_STD = [x_crop_D x_crop_D(end:-1:1)];
    CdSD1_crop_CONF = [CdSD1_crop_mean+sigma_D1_std; ...
                       CdSD1_crop_mean(end:-1:1)-sigma_D1_std(end:-1:1)];
    % Sample t9 (B2)
    CdSB2_crop_mean = mean(CdSB2_crop,2);
    sigma_B2_std = std(CdSB2_crop,0,2);
    x_crop_B2_STD = [x_crop_B2 x_crop_B2(end:-1:1)];
    CdSB2_crop_CONF = [CdSB2_crop_mean+sigma_B2_std; ...
                       CdSB2_crop_mean(end:-1:1)-sigma_B2_std(end:-1:1)];
    % Sample t10 (B1)
    CdSB1_crop_mean = mean(CdSB1_crop,2);
    sigma_B1_std = std(CdSB1_crop,0,2);
    x_crop_B_STD = [x_crop_B x_crop_B(end:-1:1)];
    CdSB1_crop_CONF = [CdSB1_crop_mean+sigma_B1_std; ...
                       CdSB1_crop_mean(end:-1:1)-sigma_B1_std(end:-1:1)];

    % =====================================================================
    % Lamp emission spectra
    % =====================================================================
    % The three lamp spectra are truncated at the CdS band-gap wavelength
    % and then averaged and adimensionalized.

    lambda_g = 512;       % CdS band-gap wavelength [nm]

    % Sample t4
    C1_Spectrum = readmatrix(fullfile(rawDataPath,'Lamp_Spectrum_t4.txt'));
    lambda_C1 = C1_Spectrum(2:end-1,1);
    I_C1 = C1_Spectrum(2:end-1,2);
    cut_lambda_C1_index = find(lambda_C1 > lambda_g,1);
    lambda_C1_cut = lambda_C1(1:cut_lambda_C1_index);
    I_C1_cut = I_C1(1:cut_lambda_C1_index);

    % Sample t7
    D1_Spectrum = readmatrix(fullfile(rawDataPath,'Lamp_Spectrum_t7.txt'));
    lambda_D1 = D1_Spectrum(2:end-1,1);
    I_D1 = D1_Spectrum(2:end-1,2);
    cut_lambda_D1_index = find(lambda_D1 > lambda_g,1);
    lambda_D1_cut = lambda_D1(1:cut_lambda_D1_index);
    I_D1_cut = I_D1(1:cut_lambda_D1_index);

    % Sample t10
    B1_Spectrum = readmatrix(fullfile(rawDataPath,'Lamp_Spectrum_t10.txt'));
    lambda_B1 = B1_Spectrum(2:end-1,1);
    I_B1 = B1_Spectrum(2:end-1,2);
    cut_lambda_B1_index = find(lambda_B1 > lambda_g,1);
    lambda_B1_cut = lambda_B1(1:cut_lambda_B1_index);
    I_B1_cut = I_B1(1:cut_lambda_B1_index);

    % Mean lamp irradiance and Adimensionalization
    I_mean = (I_C1_cut + I_D1_cut + I_B1_cut)/3;
    lambda_max = lambda_C1_cut(end);
    lambda_adim = lambda_C1_cut/lambda_max;
    Iref = max(I_mean);
    I_mean_adim = subplus(I_mean/Iref);

    % =====================================================================
    % Reflectance data 
    % =====================================================================
    c_c_ref = 4.82/144.46; lambda_g = 512; lambda_0 = 380;
    R = readmatrix(fullfile(rawDataPath,'Reflectance_UV_Vis.txt'));
    WL = R(end:-1:1,1);
    Refl_un = abs(R(end:-1:1,2)/100);

    % =====================================================================
    % Water profile parameters
    % =====================================================================
    RH = 0.95;                   % Relative humidity [-]
    RHm = 0.45;                  % Lower relative humidity bound [-]
    Temp = 25;                   % Environment temperature [°C]
    npsi = 1;                    % Shape factor of water profile
    
    SatVap = 1e-6*(5.018 + 0.32321*Temp+8.1847e-3*Temp^2+3.1243e-4*Temp^3);
    w_ref = SatVap*RH/18.01528;  % Water reference value [mol/cm^3]
    wbm = SatVap*RHm/18.01528;  % Lower threshold of water [mol/cm^3]

    % =====================================================================
    % Save preprocessed data                                              
    % =====================================================================
    save(dataFile,'-v7.3');
end


%% Plots
% =========================================================================
% Common colour limits and colormap
cLim = [min(CdSB1_crop(:)) 1];

% Colours used for the mean profiles
profileColors = lines(5);

% =========================================================================
% RAW data before cropping
% =========================================================================
figure()
t = tiledlayout(1,5,'TileSpacing','compact','Padding','compact');

% Sample t1 (C2)
nexttile
imagesc(y_C2,x_C2,CdSC2)
title('Sample t1')
xlabel('y [cm]')
ylabel('z (depth) [cm]')
clim([min(CdSB1_crop(:)) 1])

% Sample t4 (C1)
nexttile
imagesc(y_C1,x_C1,CdSC1)
title('Sample t4')
xlabel('y [cm]')
ylabel('z (depth) [cm]')
clim([min(CdSB1_crop(:)) 1])

% Sample t7 (D1)
nexttile
imagesc(y_D1,x_D1,CdSD1)
title('Sample t7')
xlabel('y [cm]')
ylabel('z (depth) [cm]')
clim([min(CdSB1_crop(:)) 1])

% Sample t9 (B2)
nexttile
imagesc(y_B2,x_B2,CdSB2)
title('Sample t9')
xlabel('y [cm]')
ylabel('z (depth) [cm]')
clim([min(CdSB1_crop(:)) 1])

% Sample t10 (B1)
nexttile
imagesc(y_B1,x_B1,CdSB1)
title('Sample t10')
xlabel('y [cm]')
ylabel('z (depth) [cm]')
clim([min(CdSB1_crop(:)) 1])

colormap(jet)
sgtitle('CdS Relative Concentrations - RAW Data')
cb = colorbar; cb.Layout.Tile = 'south';

% =========================================================================
% Processed data and y-averaged profiles
% =========================================================================

% Sample t1 (C2)
figure()
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile
imagesc(y_C2,x_crop_C2,CdSC2_crop)
axis square tight
clim(cLim); colormap(jet);
title('Processed CdS Relative Concentration - Sample t1')
xlabel('y [cm]')
ylabel('z (depth) [cm]')
cb = colorbar;
cb.Label.String = 'Relative CdS concentration';
nexttile
hold on
hConf = fill(x_crop_C2_STD,CdSC2_crop_CONF, ...
    profileColors(1,:), ...
    'EdgeColor','none', ...
    'FaceAlpha',0.25);
hMean = plot(x_crop_C2,CdSC2_crop_mean, ...
    'Color',profileColors(1,:), ...
    'LineWidth',2);
axis square tight
xlabel('z (depth) [cm]')
ylabel('Relative CdS concentration')
title('y-Averaged CdS Relative Concentration - Sample t1')
legend([hConf hMean],'\mu \pm \sigma','Mean \mu', ...
    'Location','best')


% Sample t4 (C1)
figure()
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile
imagesc(y_C1,x_crop_C,CdSC1_crop)
axis square tight
clim(cLim); colormap(jet);
title('Processed CdS Relative Concentration - Sample t4')
xlabel('y [cm]')
ylabel('z (depth) [cm]')
cb = colorbar;
cb.Label.String = 'Relative CdS concentration';
nexttile
hold on
hConf = fill(x_crop_C_STD,CdSC1_crop_CONF, ...
    profileColors(2,:), ...
    'EdgeColor','none', ...
    'FaceAlpha',0.25);
hMean = plot(x_crop_C,CdSC1_crop_mean, ...
    'Color',profileColors(2,:), ...
    'LineWidth',2);
axis square tight
xlabel('z (depth) [cm]')
ylabel('Relative CdS concentration')
title('y-Averaged CdS Relative Concentration - Sample t4')
legend([hConf hMean],'\mu \pm \sigma','Mean \mu', ...
    'Location','best')

% Sample t7 (D1)
figure()
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile
imagesc(y_D1,x_crop_D,CdSD1_crop)
axis square tight
clim(cLim); colormap(jet);
title('Processed CdS Relative Concentration - Sample t7')
xlabel('y [cm]')
ylabel('z (depth) [cm]')
cb = colorbar;
cb.Label.String = 'Relative CdS concentration';
nexttile
hold on
hConf = fill(x_crop_D_STD,CdSD1_crop_CONF, ...
    profileColors(3,:), ...
    'EdgeColor','none', ...
    'FaceAlpha',0.25);
hMean = plot(x_crop_D,CdSD1_crop_mean, ...
    'Color',profileColors(3,:), ...
    'LineWidth',2);
axis square tight
xlabel('z (depth) [cm]')
ylabel('Relative CdS concentration')
title('y-Averaged CdS Relative Concentration - Sample t7')
legend([hConf hMean],'\mu \pm \sigma','Mean \mu', ...
    'Location','best')


% Sample t9 (B2)
figure()
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile
imagesc(y_B2,x_crop_B2,CdSB2_crop)
axis square tight
clim(cLim); colormap(jet);
title('Processed CdS Relative Concentration - Sample t9')
xlabel('y [cm]')
ylabel('z (depth) [cm]')
cb = colorbar;
cb.Label.String = 'Relative CdS concentration';
nexttile
hold on
hConf = fill(x_crop_B2_STD,CdSB2_crop_CONF, ...
    profileColors(4,:), ...
    'EdgeColor','none', ...
    'FaceAlpha',0.25);
hMean = plot(x_crop_B2,CdSB2_crop_mean, ...
    'Color',profileColors(4,:), ...
    'LineWidth',2);
axis square tight
xlabel('z (depth) [cm]')
ylabel('Relative CdS concentration')
title('y-Averaged CdS Relative Concentration - Sample t9')
legend([hConf hMean],'\mu \pm \sigma','Mean \mu', ...
    'Location','best')


% Sample t10 (B1)
figure()
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile
imagesc(y_B1,x_crop_B,CdSB1_crop)
axis square tight
clim(cLim); colormap(jet);
title('Processed CdS Relative Concentration - Sample t10')
xlabel('y [cm]')
ylabel('z (depth) [cm]')
cb = colorbar;
cb.Label.String = 'Relative CdS concentration';
nexttile
hold on
hConf = fill(x_crop_B_STD,CdSB1_crop_CONF, ...
    profileColors(5,:), ...
    'EdgeColor','none', ...
    'FaceAlpha',0.25);
hMean = plot(x_crop_B,CdSB1_crop_mean, ...
    'Color',profileColors(5,:), ...
    'LineWidth',2);
axis square tight
xlabel('z (depth) [cm]')
ylabel('Relative CdS concentration')
title('y-Averaged CdS Relative Concentration - Sample t10')
legend([hConf hMean],'\mu \pm \sigma','Mean \mu', ...
    'Location','best')