%% ukb_stats.m - Sex differences in subcortical volume and shape in mid-adulthood (UK Biobank)
%
% Lee HM, Liu S, Lalonde FM, et al. Sex Differences in Human Subcortical Anatomy
% are Robust and Replicable in Adulthood, Already Evident Prepubertally, and Aligned
% with Sex Chromosomal and Gonadal Influences. bioRxiv (2026).
% https://doi.org/10.64898/2026.10.05.756728
%
% Purpose:
%   Applies the HCP-YA models to the independent UK Biobank mid-adulthood
%   sample (ages 45-50 y), additionally controlling for scanning site.
%
% Inputs (participant-level data are not shared; see Data availability):
%   - <dataDir>/UKB.mat (struct UKB)
%   - <templateDir>/maget_template_surfaces_distance_area.mat
%
% Outputs:
%   - <dataDir>/UKB_results.mat (UKB struct with statistics)
%   - Figures in <figDir>
%
% Dependencies (not included in this repository):
%   - SurfStat (SurfStatLinMod, SurfStatT, term)
%   - fdr_bh (Benjamini-Hochberg FDR; MATLAB File Exchange)
%   - MATLAB Statistics and Machine Learning Toolbox
%   - daviolinplot (MATLAB File Exchange)
%   - Custom visualization helpers used only for figures (not included):
%     loadColorMaps, ThreshColorMap, HyoMagetView* (MAGeT-Brain surface viewers)
%
% Sections (%%) are designed to be run in order using "Run Section".

clear all; close all; clc;

%% Paths (edit to match your local setup)

toolboxDir         = '/path/to/matlab_toolboxes';          % SurfStat, fdr_bh, BrainSpace, plotting helpers
templateDir        = '/path/to/maget_templates';           % MAGeT-Brain template surfaces
codeDir            = '/path/to/this_repository';
dataDir            = '/path/to/processed_data';            % participant-level data and results (not shared)
corticalResultsDir = '/path/to/cortical_study_results';    % group-level cortical maps (Lee et al., 2026)
figDir             = '/path/to/figures/ukb';

addpath(genpath(toolboxDir));
addpath(genpath(codeDir));

%% Prep

% Load UKB data: preprocessed participant-level data
load(fullfile(dataDir, 'UKB.mat'));

% Load vertex-wise surface area maps for template surfaces
load(fullfile(templateDir, 'maget_template_surfaces_distance_area.mat'));

% Misc
loadColorMaps;
nPerm = 10000;
thresh = 0.05;

%% Define the terms and models

% Compute z-scores (struct)
STR_L = zscore(UKB.VOL.STR_L(UKB.QC.indSTR))';
STR_R = zscore(UKB.VOL.STR_R(UKB.QC.indSTR))';

GP_L = zscore(UKB.VOL.GP_L(UKB.QC.indGP))';
GP_R = zscore(UKB.VOL.GP_R(UKB.QC.indGP))';

THAL_L = zscore(UKB.VOL.THAL_L(UKB.QC.indTHAL))';
THAL_R = zscore(UKB.VOL.THAL_R(UKB.QC.indTHAL))';

HIPP_L = zscore(UKB.VOL.HIPP_L(UKB.QC.indHIPP))';
HIPP_R = zscore(UKB.VOL.HIPP_R(UKB.QC.indHIPP))';

AMYG_L = zscore(UKB.VOL.AMYG_L(UKB.QC.indAMYG))';
AMYG_R = zscore(UKB.VOL.AMYG_R(UKB.QC.indAMYG))';

% Compute z-scores (morph)
STR_MOR = zscore(UKB.MOR.STR(UKB.QC.indSTR,:));
GP_MOR = zscore(UKB.MOR.GP(UKB.QC.indGP,:));
THAL_MOR = zscore(UKB.MOR.THAL(UKB.QC.indTHAL,:));
HIPP_MOR = zscore(UKB.MOR.HIPP(UKB.QC.indHIPP,:));
AMYG_MOR = zscore(UKB.MOR.AMYG(UKB.QC.indAMYG,:));

% Prepare terms (struct)
SEX_STR = term(UKB.SEX(UKB.QC.indSTR));
AGE_STR = term(UKB.AGE(UKB.QC.indSTR) - mean(UKB.AGE(UKB.QC.indSTR)));
TBV_STR = term(UKB.TBV(UKB.QC.indSTR) - mean(UKB.TBV(UKB.QC.indSTR)));
SITE_STR = term(UKB.SITE(UKB.QC.indSTR,:));

SEX_GP = term(UKB.SEX(UKB.QC.indGP));
AGE_GP = term(UKB.AGE(UKB.QC.indGP) - mean(UKB.AGE(UKB.QC.indGP)));
TBV_GP = term(UKB.TBV(UKB.QC.indGP) - mean(UKB.TBV(UKB.QC.indGP)));
SITE_GP = term(UKB.SITE(UKB.QC.indGP,:));

SEX_THAL = term(UKB.SEX(UKB.QC.indTHAL));
AGE_THAL = term(UKB.AGE(UKB.QC.indTHAL) - mean(UKB.AGE(UKB.QC.indTHAL)));
TBV_THAL = term(UKB.TBV(UKB.QC.indTHAL) - mean(UKB.TBV(UKB.QC.indTHAL)));
SITE_THAL = term(UKB.SITE(UKB.QC.indTHAL,:));

SEX_HIPP = term(UKB.SEX(UKB.QC.indHIPP));
AGE_HIPP = term(UKB.AGE(UKB.QC.indHIPP) - mean(UKB.AGE(UKB.QC.indHIPP)));
TBV_HIPP = term(UKB.TBV(UKB.QC.indHIPP) - mean(UKB.TBV(UKB.QC.indHIPP)));
SITE_HIPP = term(UKB.SITE(UKB.QC.indHIPP,:));

SEX_AMYG = term(UKB.SEX(UKB.QC.indAMYG));
AGE_AMYG = term(UKB.AGE(UKB.QC.indAMYG) - mean(UKB.AGE(UKB.QC.indAMYG)));
TBV_AMYG = term(UKB.TBV(UKB.QC.indAMYG) - mean(UKB.TBV(UKB.QC.indAMYG)));
SITE_AMYG = term(UKB.SITE(UKB.QC.indAMYG,:));

% Prepare terms (morph)
TSA_STR = term(sum(UKB.MOR.STR(UKB.QC.indSTR,:),2) - mean(sum(UKB.MOR.STR(UKB.QC.indSTR,:),2)));
TSA_GP = term(sum(UKB.MOR.GP(UKB.QC.indGP,:),2) - mean(sum(UKB.MOR.GP(UKB.QC.indGP,:),2)));
TSA_THAL = term(sum(UKB.MOR.THAL(UKB.QC.indTHAL,:),2) - mean(sum(UKB.MOR.THAL(UKB.QC.indTHAL,:),2)));
TSA_HIPP = term(sum(UKB.MOR.HIPP(UKB.QC.indHIPP,:),2) - mean(sum(UKB.MOR.HIPP(UKB.QC.indHIPP,:),2)));
TSA_AMYG = term(sum(UKB.MOR.AMYG(UKB.QC.indAMYG,:),2) - mean(sum(UKB.MOR.AMYG(UKB.QC.indAMYG,:),2)));

% Define GLMs (struct)
UKB.STAT.VOL.TBV.M_STR = 1 + SEX_STR + AGE_STR + TBV_STR + SITE_STR;
UKB.STAT.VOL.TBV.M_GP = 1 + SEX_GP + AGE_GP + TBV_GP + SITE_GP;
UKB.STAT.VOL.TBV.M_THAL = 1 + SEX_THAL + AGE_THAL + TBV_THAL + SITE_THAL;
UKB.STAT.VOL.TBV.M_HIPP = 1 + SEX_HIPP + AGE_HIPP + TBV_HIPP + SITE_HIPP;
UKB.STAT.VOL.TBV.M_AMYG = 1 + SEX_AMYG + AGE_AMYG + TBV_AMYG + SITE_AMYG;

% Define GLMs (morph)
UKB.STAT.MOR.TSA.M_STR = 1 + SEX_STR + AGE_STR + TSA_STR + SITE_STR;
UKB.STAT.MOR.TSA.M_GP = 1 + SEX_GP + AGE_GP + TSA_GP + SITE_GP;
UKB.STAT.MOR.TSA.M_THAL = 1 + SEX_THAL + AGE_THAL + TSA_THAL + SITE_THAL;
UKB.STAT.MOR.TSA.M_HIPP = 1 + SEX_HIPP + AGE_HIPP + TSA_HIPP + SITE_HIPP;
UKB.STAT.MOR.TSA.M_AMYG = 1 + SEX_AMYG + AGE_AMYG + TSA_AMYG + SITE_AMYG;

%% Compute statistical parameters (TBV | TSA correction)

% Fit models (struct)
UKB.STAT.VOL.TBV.STR_L = SurfStatLinMod(STR_L, UKB.STAT.VOL.TBV.M_STR);
UKB.STAT.VOL.TBV.STR_L_SEX = SurfStatT(UKB.STAT.VOL.TBV.STR_L, SEX_STR);
STR_L_SEX_P = 2*(1 - tcdf(abs(UKB.STAT.VOL.TBV.STR_L_SEX.t), UKB.STAT.VOL.TBV.STR_L_SEX.df));

UKB.STAT.VOL.TBV.STR_R = SurfStatLinMod(STR_R, UKB.STAT.VOL.TBV.M_STR);
UKB.STAT.VOL.TBV.STR_R_SEX = SurfStatT(UKB.STAT.VOL.TBV.STR_R, SEX_STR);
STR_R_SEX_P = 2*(1 - tcdf(abs(UKB.STAT.VOL.TBV.STR_R_SEX.t), UKB.STAT.VOL.TBV.STR_R_SEX.df));

UKB.STAT.VOL.TBV.GP_L = SurfStatLinMod(GP_L, UKB.STAT.VOL.TBV.M_GP);
UKB.STAT.VOL.TBV.GP_L_SEX = SurfStatT(UKB.STAT.VOL.TBV.GP_L, SEX_GP);
GP_L_SEX_P = 2*(1 - tcdf(abs(UKB.STAT.VOL.TBV.GP_L_SEX.t), UKB.STAT.VOL.TBV.GP_L_SEX.df));

UKB.STAT.VOL.TBV.GP_R = SurfStatLinMod(GP_R, UKB.STAT.VOL.TBV.M_GP);
UKB.STAT.VOL.TBV.GP_R_SEX = SurfStatT(UKB.STAT.VOL.TBV.GP_R, SEX_GP);
GP_R_SEX_P = 2*(1 - tcdf(abs(UKB.STAT.VOL.TBV.GP_R_SEX.t), UKB.STAT.VOL.TBV.GP_R_SEX.df));

UKB.STAT.VOL.TBV.THAL_L = SurfStatLinMod(THAL_L, UKB.STAT.VOL.TBV.M_THAL);
UKB.STAT.VOL.TBV.THAL_L_SEX = SurfStatT(UKB.STAT.VOL.TBV.THAL_L, SEX_THAL);
THAL_L_SEX_P = 2*(1 - tcdf(abs(UKB.STAT.VOL.TBV.THAL_L_SEX.t), UKB.STAT.VOL.TBV.THAL_L_SEX.df));

UKB.STAT.VOL.TBV.THAL_R = SurfStatLinMod(THAL_R, UKB.STAT.VOL.TBV.M_THAL);
UKB.STAT.VOL.TBV.THAL_R_SEX = SurfStatT(UKB.STAT.VOL.TBV.THAL_R, SEX_THAL);
THAL_R_SEX_P = 2*(1 - tcdf(abs(UKB.STAT.VOL.TBV.THAL_R_SEX.t), UKB.STAT.VOL.TBV.THAL_R_SEX.df));

UKB.STAT.VOL.TBV.HIPP_L = SurfStatLinMod(HIPP_L, UKB.STAT.VOL.TBV.M_HIPP);
UKB.STAT.VOL.TBV.HIPP_L_SEX = SurfStatT(UKB.STAT.VOL.TBV.HIPP_L, SEX_HIPP);
HIPP_L_SEX_P = 2*(1 - tcdf(abs(UKB.STAT.VOL.TBV.HIPP_L_SEX.t), UKB.STAT.VOL.TBV.HIPP_L_SEX.df));

UKB.STAT.VOL.TBV.HIPP_R = SurfStatLinMod(HIPP_R, UKB.STAT.VOL.TBV.M_HIPP);
UKB.STAT.VOL.TBV.HIPP_R_SEX = SurfStatT(UKB.STAT.VOL.TBV.HIPP_R, SEX_HIPP);
HIPP_R_SEX_P = 2*(1 - tcdf(abs(UKB.STAT.VOL.TBV.HIPP_R_SEX.t), UKB.STAT.VOL.TBV.HIPP_R_SEX.df));

UKB.STAT.VOL.TBV.AMYG_L = SurfStatLinMod(AMYG_L, UKB.STAT.VOL.TBV.M_AMYG);
UKB.STAT.VOL.TBV.AMYG_L_SEX = SurfStatT(UKB.STAT.VOL.TBV.AMYG_L, SEX_AMYG);
AMYG_L_SEX_P = 2*(1 - tcdf(abs(UKB.STAT.VOL.TBV.AMYG_L_SEX.t), UKB.STAT.VOL.TBV.AMYG_L_SEX.df));

UKB.STAT.VOL.TBV.AMYG_R = SurfStatLinMod(AMYG_R, UKB.STAT.VOL.TBV.M_AMYG);
UKB.STAT.VOL.TBV.AMYG_R_SEX = SurfStatT(UKB.STAT.VOL.TBV.AMYG_R, SEX_AMYG);
AMYG_R_SEX_P = 2*(1 - tcdf(abs(UKB.STAT.VOL.TBV.AMYG_R_SEX.t), UKB.STAT.VOL.TBV.AMYG_R_SEX.df));

% Collate statistical values
UKB.STAT.VOL.TBV.SEX_B = [UKB.STAT.VOL.TBV.STR_L_SEX.coef(2) UKB.STAT.VOL.TBV.GP_L_SEX.coef(2) UKB.STAT.VOL.TBV.THAL_L_SEX.coef(2) UKB.STAT.VOL.TBV.HIPP_L_SEX.coef(2) UKB.STAT.VOL.TBV.AMYG_L_SEX.coef(2) ...
                          UKB.STAT.VOL.TBV.STR_R_SEX.coef(2) UKB.STAT.VOL.TBV.GP_R_SEX.coef(2) UKB.STAT.VOL.TBV.THAL_R_SEX.coef(2) UKB.STAT.VOL.TBV.HIPP_R_SEX.coef(2) UKB.STAT.VOL.TBV.AMYG_R_SEX.coef(2)];

UKB.STAT.VOL.TBV.SEX_P = [STR_L_SEX_P GP_L_SEX_P THAL_L_SEX_P HIPP_L_SEX_P AMYG_L_SEX_P ...
                          STR_R_SEX_P GP_R_SEX_P THAL_R_SEX_P HIPP_R_SEX_P AMYG_R_SEX_P];

UKB.STAT.VOL.TBV.SEX_Q = 10*UKB.STAT.VOL.TBV.SEX_P;

% Fit models (morph)
UKB.STAT.MOR.TSA.STR = SurfStatLinMod(STR_MOR, UKB.STAT.MOR.TSA.M_STR);
UKB.STAT.MOR.TSA.STR_SEX = SurfStatT(UKB.STAT.MOR.TSA.STR, SEX_STR);
UKB.STAT.MOR.TSA.STR_SEX_B = UKB.STAT.MOR.TSA.STR_SEX.coef(2,:);
UKB.STAT.MOR.TSA.STR_SEX_P = 2*(1 - tcdf(abs(UKB.STAT.MOR.TSA.STR_SEX.t), UKB.STAT.MOR.TSA.STR_SEX.df));
[~,~,~,UKB.STAT.MOR.TSA.STR_SEX_Q] = fdr_bh(UKB.STAT.MOR.TSA.STR_SEX_P);

UKB.STAT.MOR.TSA.GP = SurfStatLinMod(GP_MOR, UKB.STAT.MOR.TSA.M_GP);
UKB.STAT.MOR.TSA.GP_SEX = SurfStatT(UKB.STAT.MOR.TSA.GP, SEX_GP);
UKB.STAT.MOR.TSA.GP_SEX_B = UKB.STAT.MOR.TSA.GP_SEX.coef(2,:);
UKB.STAT.MOR.TSA.GP_SEX_P = 2*(1 - tcdf(abs(UKB.STAT.MOR.TSA.GP_SEX.t), UKB.STAT.MOR.TSA.GP_SEX.df));
[~,~,~,UKB.STAT.MOR.TSA.GP_SEX_Q] = fdr_bh(UKB.STAT.MOR.TSA.GP_SEX_P);

UKB.STAT.MOR.TSA.THAL = SurfStatLinMod(THAL_MOR, UKB.STAT.MOR.TSA.M_THAL);
UKB.STAT.MOR.TSA.THAL_SEX = SurfStatT(UKB.STAT.MOR.TSA.THAL, SEX_THAL);
UKB.STAT.MOR.TSA.THAL_SEX_B = UKB.STAT.MOR.TSA.THAL_SEX.coef(2,:);
UKB.STAT.MOR.TSA.THAL_SEX_P = 2*(1 - tcdf(abs(UKB.STAT.MOR.TSA.THAL_SEX.t), UKB.STAT.MOR.TSA.THAL_SEX.df));
[~,~,~,UKB.STAT.MOR.TSA.THAL_SEX_Q] = fdr_bh(UKB.STAT.MOR.TSA.THAL_SEX_P);

UKB.STAT.MOR.TSA.HIPP = SurfStatLinMod(HIPP_MOR, UKB.STAT.MOR.TSA.M_HIPP);
UKB.STAT.MOR.TSA.HIPP_SEX = SurfStatT(UKB.STAT.MOR.TSA.HIPP, SEX_HIPP);
UKB.STAT.MOR.TSA.HIPP_SEX_B = UKB.STAT.MOR.TSA.HIPP_SEX.coef(2,:);
UKB.STAT.MOR.TSA.HIPP_SEX_P = 2*(1 - tcdf(abs(UKB.STAT.MOR.TSA.HIPP_SEX.t), UKB.STAT.MOR.TSA.HIPP_SEX.df));
[~,~,~,UKB.STAT.MOR.TSA.HIPP_SEX_Q] = fdr_bh(UKB.STAT.MOR.TSA.HIPP_SEX_P);

UKB.STAT.MOR.TSA.AMYG = SurfStatLinMod(AMYG_MOR, UKB.STAT.MOR.TSA.M_AMYG);
UKB.STAT.MOR.TSA.AMYG_SEX = SurfStatT(UKB.STAT.MOR.TSA.AMYG, SEX_AMYG);
UKB.STAT.MOR.TSA.AMYG_SEX_B = UKB.STAT.MOR.TSA.AMYG_SEX.coef(2,:);
UKB.STAT.MOR.TSA.AMYG_SEX_P = 2*(1 - tcdf(abs(UKB.STAT.MOR.TSA.AMYG_SEX.t), UKB.STAT.MOR.TSA.AMYG_SEX.df));
[~,~,~,UKB.STAT.MOR.TSA.AMYG_SEX_Q] = fdr_bh(UKB.STAT.MOR.TSA.AMYG_SEX_P);

%% Save results

save(fullfile(dataDir, 'UKB_results.mat'), 'UKB', '-v7.3');

%% Visualize the results

% Effect size for regional shape
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(UKB.STAT.MOR.TSA.STR_SEX_B, [-0.5 0.5], cmap.PiYG); exportgraphics(h, fullfile(figDir, 'UKB_STR_TSA.png'), 'Resolution', 300); close(h);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(UKB.STAT.MOR.TSA.GP_SEX_B, [-0.5 0.5], cmap.PiYG); exportgraphics(h, fullfile(figDir, 'UKB_GP_TSA.png'), 'Resolution', 300); close(h);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(UKB.STAT.MOR.TSA.THAL_SEX_B, [-0.5 0.5], cmap.PiYG); exportgraphics(h, fullfile(figDir, 'UKB_THAL_TSA.png'), 'Resolution', 300); close(h);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(UKB.STAT.MOR.TSA.HIPP_SEX_B, [-0.5 0.5], cmap.PiYG); exportgraphics(h, fullfile(figDir, 'UKB_HIPP_TSA.png'), 'Resolution', 300); close(h);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(UKB.STAT.MOR.TSA.AMYG_SEX_B, [-0.5 0.5], cmap.PiYG); exportgraphics(h, fullfile(figDir, 'UKB_AMYG_TSA.png'), 'Resolution', 300); close(h);

% FDR thresholded with gray background
temp = UKB.STAT.MOR.TSA.STR_SEX_B .* (UKB.STAT.MOR.TSA.STR_SEX_Q < thresh); h = figure('WindowState', 'Maximized'); HyoMagetViewNL(temp, [-0.5 0.5], ThreshColorMap(temp, cmap.PiYG, [-0.5 0.5], 0.915)); exportgraphics(h, fullfile(figDir, 'UKB_STR_TSA_FDR.png'), 'Resolution', 300); close(h);
temp = UKB.STAT.MOR.TSA.GP_SEX_B .* (UKB.STAT.MOR.TSA.GP_SEX_Q < thresh); h = figure('WindowState', 'Maximized'); HyoMagetViewNL(temp, [-0.5 0.5], ThreshColorMap(temp, cmap.PiYG, [-0.5 0.5], 0.915)); exportgraphics(h, fullfile(figDir, 'UKB_GP_TSA_FDR.png'), 'Resolution', 300); close(h);
temp = UKB.STAT.MOR.TSA.THAL_SEX_B .* (UKB.STAT.MOR.TSA.THAL_SEX_Q < thresh); h = figure('WindowState', 'Maximized'); HyoMagetViewNL(temp, [-0.5 0.5], ThreshColorMap(temp, cmap.PiYG, [-0.5 0.5], 0.915)); exportgraphics(h, fullfile(figDir, 'UKB_THAL_TSA_FDR.png'), 'Resolution', 300); close(h);
temp = UKB.STAT.MOR.TSA.HIPP_SEX_B .* (UKB.STAT.MOR.TSA.HIPP_SEX_Q < thresh); h = figure('WindowState', 'Maximized'); HyoMagetViewNL(temp, [-0.5 0.5], ThreshColorMap(temp, cmap.PiYG, [-0.5 0.5], 0.915)); exportgraphics(h, fullfile(figDir, 'UKB_HIPP_TSA_FDR.png'), 'Resolution', 300); close(h);
temp = UKB.STAT.MOR.TSA.AMYG_SEX_B .* (UKB.STAT.MOR.TSA.AMYG_SEX_Q < thresh); h = figure('WindowState', 'Maximized'); HyoMagetViewNL(temp, [-0.5 0.5], ThreshColorMap(temp, cmap.PiYG, [-0.5 0.5], 0.915)); exportgraphics(h, fullfile(figDir, 'UKB_AMYG_TSA_FDR.png'), 'Resolution', 300); close(h);

% Display the areal proportion of significance
sum(S_STR_NA(UKB.STAT.MOR.TSA.STR_SEX_Q < thresh))
sum(S_GP_NA(UKB.STAT.MOR.TSA.GP_SEX_Q < thresh))
sum(S_THAL_NA(UKB.STAT.MOR.TSA.THAL_SEX_Q < thresh))
sum(S_HIPP_NA(UKB.STAT.MOR.TSA.HIPP_SEX_Q < thresh))
sum(S_AMYG_NA(UKB.STAT.MOR.TSA.AMYG_SEX_Q < thresh))

% Calculate the colorscale of thresholded maps
temp = UKB.STAT.MOR.TSA.STR_SEX_B .* (UKB.STAT.MOR.TSA.STR_SEX_Q < thresh);
[min(abs(temp(temp<0))) max(abs(temp(temp<0))) min(abs(temp(temp>0))) max(abs(temp(temp>0)))]

temp = UKB.STAT.MOR.TSA.GP_SEX_B .* (UKB.STAT.MOR.TSA.GP_SEX_Q < thresh);
[min(abs(temp(temp<0))) max(abs(temp(temp<0))) min(abs(temp(temp>0))) max(abs(temp(temp>0)))]

temp = UKB.STAT.MOR.TSA.THAL_SEX_B .* (UKB.STAT.MOR.TSA.THAL_SEX_Q < thresh);
[min(abs(temp(temp<0))) max(abs(temp(temp<0))) min(abs(temp(temp>0))) max(abs(temp(temp>0)))]

temp = UKB.STAT.MOR.TSA.HIPP_SEX_B .* (UKB.STAT.MOR.TSA.HIPP_SEX_Q < thresh);
[min(abs(temp(temp<0))) max(abs(temp(temp<0))) min(abs(temp(temp>0))) max(abs(temp(temp>0)))]

temp = UKB.STAT.MOR.TSA.AMYG_SEX_B .* (UKB.STAT.MOR.TSA.AMYG_SEX_Q < thresh);
[min(abs(temp(temp<0))) max(abs(temp(temp<0))) min(abs(temp(temp>0))) max(abs(temp(temp>0)))]

%% Rotating full view

% Collate into one vector
SUB_B = [UKB.STAT.MOR.TSA.STR_SEX_B UKB.STAT.MOR.TSA.GP_SEX_B UKB.STAT.MOR.TSA.THAL_SEX_B UKB.STAT.MOR.TSA.HIPP_SEX_B UKB.STAT.MOR.TSA.AMYG_SEX_B];
SUB_Q = [UKB.STAT.MOR.TSA.STR_SEX_Q UKB.STAT.MOR.TSA.GP_SEX_Q UKB.STAT.MOR.TSA.THAL_SEX_Q UKB.STAT.MOR.TSA.HIPP_SEX_Q UKB.STAT.MOR.TSA.AMYG_SEX_Q];

% Static
c = 1;
dA = 45;
for i =  0:dA:360-dA
    disp(['Rotation angle ' num2str(i) ' deg']);
    h = figure('WindowState', 'Maximized');
    HyoMagetViewSpin(SUB_B .* (SUB_Q < thresh), [-0.5 0.5], i, ThreshColorMap(SUB_B .* (SUB_Q < thresh), cmap.PiYG, [-0.5 0.5], 0.915), 'static');
    print(h, fullfile(figDir, 'static_ukb', ['UKB_TSA_ANG-' num2str(c) '.png']), '-dpng', '-r300');
    close(h);
    c = c + 1;
end

% Static (not thresholded)
c = 1;
dA = 45;
for i =  0:dA:360-dA
    disp(['Rotation angle ' num2str(i) ' deg']);
    h = figure('WindowState', 'Maximized');
    HyoMagetViewSpin(SUB_B, [-0.5 0.5], i, cmap.PiYG, 'static');
    print(h, fullfile(figDir, 'static_ukb', ['NOTHR_UKB_TSA_ANG-' num2str(c) '.png']), '-dpng', '-r300');
    close(h);
    c = c + 1;
end

%% Raincloud plot for bulk sex differences

% Extract absolute volume
STR = ((UKB.VOL.STR_L(UKB.QC.indSTR) + UKB.VOL.STR_R(UKB.QC.indSTR))');
GP = ((UKB.VOL.GP_L(UKB.QC.indGP) + UKB.VOL.GP_R(UKB.QC.indGP))');
THAL = ((UKB.VOL.THAL_L(UKB.QC.indTHAL) + UKB.VOL.THAL_R(UKB.QC.indTHAL))');
HIPP = ((UKB.VOL.HIPP_L(UKB.QC.indHIPP) + UKB.VOL.HIPP_R(UKB.QC.indHIPP))');
AMYG = ((UKB.VOL.AMYG_L(UKB.QC.indAMYG) + UKB.VOL.AMYG_R(UKB.QC.indAMYG))');

% Mean-center covariates
SEX_STR = UKB.SEX(UKB.QC.indSTR);
AGE_STR = UKB.AGE(UKB.QC.indSTR) - mean(UKB.AGE(UKB.QC.indSTR));
TBV_STR = UKB.TBV(UKB.QC.indSTR) - mean(UKB.TBV(UKB.QC.indSTR));
SITE_STR = UKB.SITE(UKB.QC.indSTR,2:end);

SEX_GP = UKB.SEX(UKB.QC.indGP);
AGE_GP = UKB.AGE(UKB.QC.indGP) - mean(UKB.AGE(UKB.QC.indGP));
TBV_GP = UKB.TBV(UKB.QC.indGP) - mean(UKB.TBV(UKB.QC.indGP));
SITE_GP = UKB.SITE(UKB.QC.indGP,2:end);

SEX_THAL = UKB.SEX(UKB.QC.indTHAL);
AGE_THAL = UKB.AGE(UKB.QC.indTHAL) - mean(UKB.AGE(UKB.QC.indTHAL));
TBV_THAL = UKB.TBV(UKB.QC.indTHAL) - mean(UKB.TBV(UKB.QC.indTHAL));
SITE_THAL = UKB.SITE(UKB.QC.indTHAL,2:end);

SEX_HIPP = UKB.SEX(UKB.QC.indHIPP);
AGE_HIPP = UKB.AGE(UKB.QC.indHIPP) - mean(UKB.AGE(UKB.QC.indHIPP));
TBV_HIPP = UKB.TBV(UKB.QC.indHIPP) - mean(UKB.TBV(UKB.QC.indHIPP));
SITE_HIPP = UKB.SITE(UKB.QC.indHIPP,2:end);

SEX_AMYG = UKB.SEX(UKB.QC.indAMYG);
AGE_AMYG = UKB.AGE(UKB.QC.indAMYG) - mean(UKB.AGE(UKB.QC.indAMYG));
TBV_AMYG = UKB.TBV(UKB.QC.indAMYG) - mean(UKB.TBV(UKB.QC.indAMYG));
SITE_AMYG = UKB.SITE(UKB.QC.indAMYG,2:end);

% Prepare terms (struct)
T_SEX_STR = term(SEX_STR);
T_AGE_STR = term(AGE_STR);
T_TBV_STR = term(TBV_STR);
T_SITE_STR = term(SITE_STR);
T_SEX_GP = term(SEX_GP);
T_AGE_GP = term(AGE_GP);
T_TBV_GP = term(TBV_GP);
T_SITE_GP = term(SITE_GP);
T_SEX_THAL = term(SEX_THAL);
T_AGE_THAL = term(AGE_THAL);
T_TBV_THAL = term(TBV_THAL);
T_SITE_THAL = term(SITE_THAL);
T_SEX_HIPP = term(SEX_HIPP);
T_AGE_HIPP = term(AGE_HIPP);
T_TBV_HIPP = term(TBV_HIPP);
T_SITE_HIPP = term(SITE_HIPP);
T_SEX_AMYG = term(SEX_AMYG);
T_AGE_AMYG = term(AGE_AMYG);
T_TBV_AMYG = term(TBV_AMYG);
T_SITE_AMYG = term(SITE_AMYG);

% Model
M_STR = 1 + T_SEX_STR + T_AGE_STR + T_TBV_STR + T_SITE_STR;
M_GP = 1 + T_SEX_GP + T_AGE_GP  + T_TBV_GP+ T_SITE_GP;
M_THAL = 1 + T_SEX_THAL + T_AGE_THAL + T_TBV_THAL+ T_SITE_THAL;
M_HIPP = 1 + T_SEX_HIPP + T_AGE_HIPP + T_TBV_HIPP+ T_SITE_HIPP;
M_AMYG = 1 + T_SEX_AMYG + T_AGE_AMYG + T_TBV_AMYG+ T_SITE_AMYG;

% Fit models
FM_STR = SurfStatLinMod(STR, M_STR);
FM_GP = SurfStatLinMod(GP, M_GP);
FM_THAL = SurfStatLinMod(THAL, M_THAL);
FM_HIPP = SurfStatLinMod(HIPP, M_HIPP);
FM_AMYG = SurfStatLinMod(AMYG, M_AMYG);

% Statistics
FM_STR_T = SurfStatT(FM_STR, SEX_STR);
FM_STR_P = 2*(1 - tcdf(abs(FM_STR_T.t), FM_STR_T.df));
FM_GP_T = SurfStatT(FM_GP, SEX_GP);
FM_GP_P = 2*(1 - tcdf(abs(FM_GP_T.t), FM_GP_T.df));
FM_THAL_T = SurfStatT(FM_THAL, SEX_THAL);
FM_THAL_P = 2*(1 - tcdf(abs(FM_THAL_T.t), FM_THAL_T.df));
FM_HIPP_T = SurfStatT(FM_HIPP, SEX_HIPP);
FM_HIPP_P = 2*(1 - tcdf(abs(FM_HIPP_T.t), FM_HIPP_T.df));
FM_AMYG_T = SurfStatT(FM_AMYG, SEX_AMYG);
FM_AMYG_P = 2*(1 - tcdf(abs(FM_AMYG_T.t), FM_AMYG_T.df));

% Multiple comparison
FM_STR_Q = FM_STR_P * 5;
FM_GP_Q = FM_GP_P * 5;
FM_THAL_Q = FM_THAL_P * 5;
FM_HIPP_Q = FM_HIPP_P * 5;
FM_AMYG_Q = FM_AMYG_P * 5;

[FM_STR_P FM_GP_P FM_THAL_P FM_HIPP_P FM_AMYG_P]
[FM_STR_Q FM_GP_Q FM_THAL_Q FM_HIPP_Q FM_AMYG_Q]

% Standardized effect sizes (z-scored inputs)
temp = SurfStatLinMod(zscore(STR), M_STR); temp.coef(2)
temp = SurfStatLinMod(zscore(GP), M_GP); temp.coef(2)
temp = SurfStatLinMod(zscore(THAL), M_THAL); temp.coef(2)
temp = SurfStatLinMod(zscore(HIPP), M_HIPP); temp.coef(2)
temp = SurfStatLinMod(zscore(AMYG), M_AMYG); temp.coef(2)

% Regress out age, TBV and site from the bulk volumes
for i = 1:size(STR,1)
    STR_SEX(i) = STR(i) - FM_STR.coef(3)*AGE_STR(i) - FM_STR.coef(4)*TBV_STR(i) - FM_STR.coef(5:6)'*SITE_STR(i,:)';
end

for i = 1:size(GP,1)
    GP_SEX(i) = GP(i) - FM_GP.coef(3)*AGE_GP(i) - FM_GP.coef(4)*TBV_GP(i) - FM_GP.coef(5:6)'*SITE_GP(i,:)';
end

for i = 1:size(THAL,1)
    THAL_SEX(i) = THAL(i) - FM_THAL.coef(3)*AGE_THAL(i) - FM_THAL.coef(4)*TBV_THAL(i) - FM_THAL.coef(5:6)'*SITE_THAL(i,:)';
end

for i = 1:size(HIPP,1)
    HIPP_SEX(i) = HIPP(i) - FM_HIPP.coef(3)*AGE_HIPP(i) - FM_HIPP.coef(4)*TBV_HIPP(i) - FM_HIPP.coef(5:6)'*SITE_HIPP(i,:)';
end

for i = 1:size(AMYG,1)
    AMYG_SEX(i) = AMYG(i) - FM_AMYG.coef(3)*AGE_AMYG(i) - FM_AMYG.coef(4)*TBV_AMYG(i) - FM_AMYG.coef(5:6)'*SITE_AMYG(i,:)';
end

% Separate the female and male values
STR_F = STR_SEX(SEX_STR==1);
GP_F = GP_SEX(SEX_GP==1);
THAL_F = THAL_SEX(SEX_THAL==1);
HIPP_F = HIPP_SEX(SEX_HIPP==1);
AMYG_F = AMYG_SEX(SEX_AMYG==1);

STR_M = STR_SEX(SEX_STR==0);
GP_M = GP_SEX(SEX_GP==0);
THAL_M = THAL_SEX(SEX_THAL==0);
HIPP_M = HIPP_SEX(SEX_HIPP==0);
AMYG_M = AMYG_SEX(SEX_AMYG==0);

% Display mean and SD
[mean(STR_F) std(STR_F) mean(STR_M) std(STR_M)]
[mean(GP_F) std(GP_F) mean(GP_M) std(GP_M)]
[mean(THAL_F) std(THAL_F) mean(THAL_M) std(THAL_M)]
[mean(HIPP_F) std(HIPP_F) mean(HIPP_M) std(HIPP_M)]
[mean(AMYG_F) std(AMYG_F) mean(AMYG_M) std(AMYG_M)]

% Add dummy elements to match the vector length
STR_F = [STR_F nan(1,1000-length(STR_F))];
GP_F = [GP_F nan(1,1000-length(GP_F))];
THAL_F = [THAL_F nan(1,1000-length(THAL_F))];
HIPP_F = [HIPP_F nan(1,1000-length(HIPP_F))];
AMYG_F = [AMYG_F nan(1,1000-length(AMYG_F))];

STR_M = [STR_M nan(1,1000-length(STR_M))];
GP_M = [GP_M nan(1,1000-length(GP_M))];
THAL_M = [THAL_M nan(1,1000-length(THAL_M))];
HIPP_M = [HIPP_M nan(1,1000-length(HIPP_M))];
AMYG_M = [AMYG_M nan(1,1000-length(AMYG_M))];

% STR
clear rainData;
rainData{1} = STR_F';
rainData{2} = STR_M';

group_inx = 1:2;
c =  [0.1529    0.3922    0.0980;...
      0.5569    0.0039    0.3216];
group_names = {'Females', 'Males'};
condition_names = {'Females', 'Males'};
h = figure; d = daviolinplot(rainData, 'groups', group_inx, 'colors', c, 'box', 3,...
    'boxcolor', 'w', 'scatter', 2, 'jitter', 1, 'scattercolor', 'same', 'smoothing', 300,...
    'scattersize', 10, 'scatteralpha', 0.4, 'outliers', 0, ...
    'xtlabels', condition_names, 'legend', group_names);
ylabel('Bulk volume after regressing out age and total brain volume');
set(d.sc,'MarkerEdgeColor','none');
set(gca,'FontSize',10);
xl = xlim; xlim([xl(1)-0.1, xl(2)+0.3]);
ylim([12000 28000]); yticks(12000:3200:28000);
ax = gca;
ax.YGrid = 'on'; ax.GridLineWidth = 2;
ax.TickDir = 'out';
h.Position = [0 0 1000 1000];
legend('off');
xlabel(''); ylabel('');
xticklabels([]); yticklabels([]);
exportgraphics(h, fullfile(figDir, 'raincloud_bulk_STR.png'), 'Resolution', 300);

% GP
rainData{1} = GP_F';
rainData{2} = GP_M';
group_inx = 1:2;
c =  [0.1529    0.3922    0.0980;...
      0.5569    0.0039    0.3216];
group_names = {'Females', 'Males'};
condition_names = {'Females', 'Males'};
h = figure; d = daviolinplot(rainData, 'groups', group_inx, 'colors', c, 'box', 3,...
    'boxcolor', 'w', 'scatter', 2, 'jitter', 1, 'scattercolor', 'same', 'smoothing', 50,...
    'scattersize', 10, 'scatteralpha', 0.4, 'outliers', 0, ...
    'xtlabels', condition_names, 'legend', group_names);
ylabel('Bulk volume after regressing out age and total brain volume');
set(d.sc,'MarkerEdgeColor','none');
set(gca,'FontSize',10);
xl = xlim; xlim([xl(1)-0.1, xl(2)+0.3]);
ylim([1800 4200]); yticks(1800:480:4200);
ax = gca;
ax.YGrid = 'on'; ax.GridLineWidth = 2;
ax.TickDir = 'out';
h.Position = [0 0 1000 1000];
legend('off');
xlabel(''); ylabel('');
xticklabels([]); yticklabels([]);
exportgraphics(h, fullfile(figDir, 'raincloud_bulk_GP.png'), 'Resolution', 300);

% THAL
rainData{1} = THAL_F';
rainData{2} = THAL_M';
group_inx = 1:2;
c =  [0.1529    0.3922    0.0980;...
      0.5569    0.0039    0.3216];
group_names = {'Females', 'Males'};
condition_names = {'Females', 'Males'};
h = figure; d = daviolinplot(rainData, 'groups', group_inx, 'colors', c, 'box', 3,...
    'boxcolor', 'w', 'scatter', 2, 'jitter', 1, 'scattercolor', 'same', 'smoothing', 175,...
    'scattersize', 10, 'scatteralpha', 0.4, 'outliers', 0, ...
    'xtlabels', condition_names, 'legend', group_names);
ylabel('Bulk volume after regressing out age and total brain volume');
set(d.sc,'MarkerEdgeColor','none');
set(gca,'FontSize',10);
xl = xlim; xlim([xl(1)-0.1, xl(2)+0.3]);
ylim([7750 16150]); yticks(7750:1680:16150);
ax = gca;
ax.YGrid = 'on'; ax.GridLineWidth = 2;
ax.TickDir = 'out';
h.Position = [0 0 1000 1000];
legend('off');
xlabel(''); ylabel('');
xticklabels([]); yticklabels([]);
exportgraphics(h, fullfile(figDir, 'raincloud_bulk_THAL.png'), 'Resolution', 300);

% HIPP
rainData{1} = HIPP_F';
rainData{2} = HIPP_M';
group_inx = 1:2;
c =  [0.1529    0.3922    0.0980;...
      0.5569    0.0039    0.3216];
group_names = {'Females', 'Males'};
condition_names = {'Females', 'Males'};
h = figure; d = daviolinplot(rainData, 'groups', group_inx, 'colors', c, 'box', 3,...
    'boxcolor', 'w', 'scatter', 2, 'jitter', 1, 'scattercolor', 'same', 'smoothing', 100,...
    'scattersize', 10, 'scatteralpha', 0.4, 'outliers', 0, ...
    'xtlabels', condition_names, 'legend', group_names);
ylabel('Bulk volume after regressing out age and total brain volume');
set(d.sc,'MarkerEdgeColor','none');
set(gca,'FontSize',10);
xl = xlim; xlim([xl(1)-0.1, xl(2)+0.3]);
ylim([1000 7000]); yticks(1000:1200:7000);
ax = gca;
ax.YGrid = 'on'; ax.GridLineWidth = 2;
ax.TickDir = 'out';
h.Position = [0 0 1000 1000];
legend('off');
xlabel(''); ylabel('');
xticklabels([]); yticklabels([]);
exportgraphics(h, fullfile(figDir, 'raincloud_bulk_HIPP.png'), 'Resolution', 300);

% AMYG
rainData{1} = AMYG_F';
rainData{2} = AMYG_M';
group_inx = 1:2;
c =  [0.1529    0.3922    0.0980;...
      0.5569    0.0039    0.3216];
group_names = {'Females', 'Males'};
condition_names = {'Females', 'Males'};
h = figure; d = daviolinplot(rainData, 'groups', group_inx, 'colors', c, 'box', 3,...
    'boxcolor', 'w', 'scatter', 2, 'jitter', 1, 'scattercolor', 'same', 'smoothing', 50,...
    'scattersize', 10, 'scatteralpha', 0.4, 'outliers', 0, ...
    'xtlabels', condition_names, 'legend', group_names);
ylabel('Bulk volume after regressing out age and total brain volume');
set(d.sc,'MarkerEdgeColor','none');
set(gca,'FontSize',10);
xl = xlim; xlim([xl(1)-0.1, xl(2)+0.3]);
ylim([1400 3600]); yticks(1400:440:3600);
ax = gca;
ax.YGrid = 'on'; ax.GridLineWidth = 2;
ax.TickDir = 'out';
h.Position = [0 0 1000 1000];
legend('off');
xlabel(''); ylabel('');
xticklabels([]); yticklabels([]);
exportgraphics(h, fullfile(figDir, 'raincloud_bulk_AMYG.png'), 'Resolution', 300);
