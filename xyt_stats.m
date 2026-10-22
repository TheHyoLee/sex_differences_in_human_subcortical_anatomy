%% xyt_stats.m - Effects of X- and Y-chromosome dosage and testicular hormones on subcortical shape
%
% Lee HM, Liu S, Lalonde FM, et al. Sex Differences in Human Subcortical Anatomy
% are Robust and Replicable in Adulthood, Already Evident Prepubertally, and Aligned
% with Sex Chromosomal and Gonadal Influences. bioRxiv (2026).
% https://doi.org/10.64898/2026.10.05.756728
%
% Purpose:
%   Estimates anatomical effects of X-chromosome dosage (XXY vs. XY),
%   Y-chromosome dosage (XYY vs. XY), and testicular hormone production
%   (IGD vs. XY controls), and tests their spatial congruence with HCP-YA
%   sex differences using Moran spectral randomization.
%
% Inputs (participant-level data are not shared; see Data availability):
%   - <dataDir>/XXY_*.mat, XYY_*.mat, IGD_*.mat (structs XXY, XYY, IGD)
%   - <dataDir>/HCP_results.mat (from b_hcp_stats.m)
%   - <templateDir>/maget_template_surfaces_distance_area.mat
%   - <corticalResultsDir>/sca_groupmaps_final.mat, igd_groupmaps_final.mat
%     (cortical maps from Lee et al., Nature Communications 2026)
%
% Outputs:
%   - <dataDir>/XXY_results_*.mat, XYY_results_*.mat, IGD_results_*.mat
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
figDir             = '/path/to/figures/xyt';

addpath(genpath(toolboxDir));
addpath(genpath(codeDir));

%% Prep

% Load data: preprocessed participant-level data
% load(fullfile(dataDir, 'XXY.mat'));
% load(fullfile(dataDir, 'XYY.mat'));
% load(fullfile(dataDir, 'IGD.mat'));

% ...or previously computed results
load(fullfile(dataDir, 'XXY_results.mat'));
load(fullfile(dataDir, 'XYY_results.mat'));
load(fullfile(dataDir, 'IGD_results.mat'));
load(fullfile(dataDir, 'HCP_results.mat')); % HCP-YA sex effects and spatial nulls (b_hcp_stats.m)

% Load surface area maps for the template surfaces
load(fullfile(templateDir, 'maget_template_surfaces_distance_area.mat'));

% Misc
loadColorMaps;
nPerm = 10000;
thresh = 0.05;

%% Define the terms and models for XXY

% Compute z-scores (struct)
XXY_STR_L = zscore(XXY.VOL.STR_L(XXY.QC.indSTR))';
XXY_STR_R = zscore(XXY.VOL.STR_R(XXY.QC.indSTR))';

XXY_GP_L = zscore(XXY.VOL.GP_L(XXY.QC.indGP))';
XXY_GP_R = zscore(XXY.VOL.GP_R(XXY.QC.indGP))';

XXY_THAL_L = zscore(XXY.VOL.THAL_L(XXY.QC.indTHAL))';
XXY_THAL_R = zscore(XXY.VOL.THAL_R(XXY.QC.indTHAL))';

XXY_HIPP_L = zscore(XXY.VOL.HIPP_L(XXY.QC.indHIPP))';
XXY_HIPP_R = zscore(XXY.VOL.HIPP_R(XXY.QC.indHIPP))';

XXY_AMYG_L = zscore(XXY.VOL.AMYG_L(XXY.QC.indAMYG))';
XXY_AMYG_R = zscore(XXY.VOL.AMYG_R(XXY.QC.indAMYG))';

% Compute z-scores (morph)
XXY_STR_MOR = zscore(XXY.MOR.STR(XXY.QC.indSTR,:));
XXY_GP_MOR = zscore(XXY.MOR.GP(XXY.QC.indGP,:));
XXY_THAL_MOR = zscore(XXY.MOR.THAL(XXY.QC.indTHAL,:));
XXY_HIPP_MOR = zscore(XXY.MOR.HIPP(XXY.QC.indHIPP,:));
XXY_AMYG_MOR = zscore(XXY.MOR.AMYG(XXY.QC.indAMYG,:));

% Prepare terms (struct)
XXY_SEX_STR = term(XXY.SEX(XXY.QC.indSTR));
XXY_AGE_STR = term(XXY.AGE(XXY.QC.indSTR) - mean(XXY.AGE(XXY.QC.indSTR)));
XXY_TBV_STR = term(XXY.TBV(XXY.QC.indSTR) - mean(XXY.TBV(XXY.QC.indSTR)));

XXY_SEX_GP = term(XXY.SEX(XXY.QC.indGP));
XXY_AGE_GP = term(XXY.AGE(XXY.QC.indGP) - mean(XXY.AGE(XXY.QC.indGP)));
XXY_TBV_GP = term(XXY.TBV(XXY.QC.indGP) - mean(XXY.TBV(XXY.QC.indGP)));

XXY_SEX_THAL = term(XXY.SEX(XXY.QC.indTHAL));
XXY_AGE_THAL = term(XXY.AGE(XXY.QC.indTHAL) - mean(XXY.AGE(XXY.QC.indTHAL)));
XXY_TBV_THAL = term(XXY.TBV(XXY.QC.indTHAL) - mean(XXY.TBV(XXY.QC.indTHAL)));

XXY_SEX_HIPP = term(XXY.SEX(XXY.QC.indHIPP));
XXY_AGE_HIPP = term(XXY.AGE(XXY.QC.indHIPP) - mean(XXY.AGE(XXY.QC.indHIPP)));
XXY_TBV_HIPP = term(XXY.TBV(XXY.QC.indHIPP) - mean(XXY.TBV(XXY.QC.indHIPP)));

XXY_SEX_AMYG = term(XXY.SEX(XXY.QC.indAMYG));
XXY_AGE_AMYG = term(XXY.AGE(XXY.QC.indAMYG) - mean(XXY.AGE(XXY.QC.indAMYG)));
XXY_TBV_AMYG = term(XXY.TBV(XXY.QC.indAMYG) - mean(XXY.TBV(XXY.QC.indAMYG)));

% Prepare terms (morph)
XXY_TSA_STR = term(sum(XXY.MOR.STR(XXY.QC.indSTR,:),2) - mean(sum(XXY.MOR.STR(XXY.QC.indSTR,:),2)));
XXY_TSA_GP = term(sum(XXY.MOR.GP(XXY.QC.indGP,:),2) - mean(sum(XXY.MOR.GP(XXY.QC.indGP,:),2)));
XXY_TSA_THAL = term(sum(XXY.MOR.THAL(XXY.QC.indTHAL,:),2) - mean(sum(XXY.MOR.THAL(XXY.QC.indTHAL,:),2)));
XXY_TSA_HIPP = term(sum(XXY.MOR.HIPP(XXY.QC.indHIPP,:),2) - mean(sum(XXY.MOR.HIPP(XXY.QC.indHIPP,:),2)));
XXY_TSA_AMYG = term(sum(XXY.MOR.AMYG(XXY.QC.indAMYG,:),2) - mean(sum(XXY.MOR.AMYG(XXY.QC.indAMYG,:),2)));

% Define GLMs (struct)
XXY.STAT.VOL.TBV.M_STR = 1 + XXY_SEX_STR + XXY_AGE_STR + XXY_TBV_STR;
XXY.STAT.VOL.TBV.M_GP = 1 + XXY_SEX_GP + XXY_AGE_GP + XXY_TBV_GP;
XXY.STAT.VOL.TBV.M_THAL = 1 + XXY_SEX_THAL + XXY_AGE_THAL + XXY_TBV_THAL;
XXY.STAT.VOL.TBV.M_HIPP = 1 + XXY_SEX_HIPP + XXY_AGE_HIPP + XXY_TBV_HIPP;
XXY.STAT.VOL.TBV.M_AMYG = 1 + XXY_SEX_AMYG + XXY_AGE_AMYG + XXY_TBV_AMYG;

% Define GLMs (morph)
XXY.STAT.MOR.TSA.M_STR = 1 + XXY_SEX_STR + XXY_AGE_STR + XXY_TSA_STR;
XXY.STAT.MOR.TSA.M_GP = 1 + XXY_SEX_GP + XXY_AGE_GP + XXY_TSA_GP;
XXY.STAT.MOR.TSA.M_THAL = 1 + XXY_SEX_THAL + XXY_AGE_THAL + XXY_TSA_THAL;
XXY.STAT.MOR.TSA.M_HIPP = 1 + XXY_SEX_HIPP + XXY_AGE_HIPP + XXY_TSA_HIPP;
XXY.STAT.MOR.TSA.M_AMYG = 1 + XXY_SEX_AMYG + XXY_AGE_AMYG + XXY_TSA_AMYG;

% Fit models (struct)
XXY.STAT.VOL.TBV.STR_L = SurfStatLinMod(XXY_STR_L, XXY.STAT.VOL.TBV.M_STR);
XXY.STAT.VOL.TBV.STR_L_SEX = SurfStatT(XXY.STAT.VOL.TBV.STR_L, XXY_SEX_STR);
XXY_STR_L_SEX_P = 2*(1 - tcdf(abs(XXY.STAT.VOL.TBV.STR_L_SEX.t), XXY.STAT.VOL.TBV.STR_L_SEX.df));

XXY.STAT.VOL.TBV.STR_R = SurfStatLinMod(XXY_STR_R, XXY.STAT.VOL.TBV.M_STR);
XXY.STAT.VOL.TBV.STR_R_SEX = SurfStatT(XXY.STAT.VOL.TBV.STR_R, XXY_SEX_STR);
XXY_STR_R_SEX_P = 2*(1 - tcdf(abs(XXY.STAT.VOL.TBV.STR_R_SEX.t), XXY.STAT.VOL.TBV.STR_R_SEX.df));

XXY.STAT.VOL.TBV.GP_L = SurfStatLinMod(XXY_GP_L, XXY.STAT.VOL.TBV.M_GP);
XXY.STAT.VOL.TBV.GP_L_SEX = SurfStatT(XXY.STAT.VOL.TBV.GP_L, XXY_SEX_GP);
XXY_GP_L_SEX_P = 2*(1 - tcdf(abs(XXY.STAT.VOL.TBV.GP_L_SEX.t), XXY.STAT.VOL.TBV.GP_L_SEX.df));

XXY.STAT.VOL.TBV.GP_R = SurfStatLinMod(XXY_GP_R, XXY.STAT.VOL.TBV.M_GP);
XXY.STAT.VOL.TBV.GP_R_SEX = SurfStatT(XXY.STAT.VOL.TBV.GP_R, XXY_SEX_GP);
XXY_GP_R_SEX_P = 2*(1 - tcdf(abs(XXY.STAT.VOL.TBV.GP_R_SEX.t), XXY.STAT.VOL.TBV.GP_R_SEX.df));

XXY.STAT.VOL.TBV.THAL_L = SurfStatLinMod(XXY_THAL_L, XXY.STAT.VOL.TBV.M_THAL);
XXY.STAT.VOL.TBV.THAL_L_SEX = SurfStatT(XXY.STAT.VOL.TBV.THAL_L, XXY_SEX_THAL);
XXY_THAL_L_SEX_P = 2*(1 - tcdf(abs(XXY.STAT.VOL.TBV.THAL_L_SEX.t), XXY.STAT.VOL.TBV.THAL_L_SEX.df));

XXY.STAT.VOL.TBV.THAL_R = SurfStatLinMod(XXY_THAL_R, XXY.STAT.VOL.TBV.M_THAL);
XXY.STAT.VOL.TBV.THAL_R_SEX = SurfStatT(XXY.STAT.VOL.TBV.THAL_R, XXY_SEX_THAL);
XXY_THAL_R_SEX_P = 2*(1 - tcdf(abs(XXY.STAT.VOL.TBV.THAL_R_SEX.t), XXY.STAT.VOL.TBV.THAL_R_SEX.df));

XXY.STAT.VOL.TBV.HIPP_L = SurfStatLinMod(XXY_HIPP_L, XXY.STAT.VOL.TBV.M_HIPP);
XXY.STAT.VOL.TBV.HIPP_L_SEX = SurfStatT(XXY.STAT.VOL.TBV.HIPP_L, XXY_SEX_HIPP);
XXY_HIPP_L_SEX_P = 2*(1 - tcdf(abs(XXY.STAT.VOL.TBV.HIPP_L_SEX.t), XXY.STAT.VOL.TBV.HIPP_L_SEX.df));

XXY.STAT.VOL.TBV.HIPP_R = SurfStatLinMod(XXY_HIPP_R, XXY.STAT.VOL.TBV.M_HIPP);
XXY.STAT.VOL.TBV.HIPP_R_SEX = SurfStatT(XXY.STAT.VOL.TBV.HIPP_R, XXY_SEX_HIPP);
XXY_HIPP_R_SEX_P = 2*(1 - tcdf(abs(XXY.STAT.VOL.TBV.HIPP_R_SEX.t), XXY.STAT.VOL.TBV.HIPP_R_SEX.df));

XXY.STAT.VOL.TBV.AMYG_L = SurfStatLinMod(XXY_AMYG_L, XXY.STAT.VOL.TBV.M_AMYG);
XXY.STAT.VOL.TBV.AMYG_L_SEX = SurfStatT(XXY.STAT.VOL.TBV.AMYG_L, XXY_SEX_AMYG);
XXY_AMYG_L_SEX_P = 2*(1 - tcdf(abs(XXY.STAT.VOL.TBV.AMYG_L_SEX.t), XXY.STAT.VOL.TBV.AMYG_L_SEX.df));

XXY.STAT.VOL.TBV.AMYG_R = SurfStatLinMod(XXY_AMYG_R, XXY.STAT.VOL.TBV.M_AMYG);
XXY.STAT.VOL.TBV.AMYG_R_SEX = SurfStatT(XXY.STAT.VOL.TBV.AMYG_R, XXY_SEX_AMYG);
XXY_AMYG_R_SEX_P = 2*(1 - tcdf(abs(XXY.STAT.VOL.TBV.AMYG_R_SEX.t), XXY.STAT.VOL.TBV.AMYG_R_SEX.df));

% Collate statistical values
XXY.STAT.VOL.TBV.SEX_B = [XXY.STAT.VOL.TBV.STR_L_SEX.coef(2) XXY.STAT.VOL.TBV.GP_L_SEX.coef(2) XXY.STAT.VOL.TBV.THAL_L_SEX.coef(2) XXY.STAT.VOL.TBV.HIPP_L_SEX.coef(2) XXY.STAT.VOL.TBV.AMYG_L_SEX.coef(2) ...
                          XXY.STAT.VOL.TBV.STR_R_SEX.coef(2) XXY.STAT.VOL.TBV.GP_R_SEX.coef(2) XXY.STAT.VOL.TBV.THAL_R_SEX.coef(2) XXY.STAT.VOL.TBV.HIPP_R_SEX.coef(2) XXY.STAT.VOL.TBV.AMYG_R_SEX.coef(2)];

XXY.STAT.VOL.TBV.SEX_P = [XXY_STR_L_SEX_P XXY_GP_L_SEX_P XXY_THAL_L_SEX_P XXY_HIPP_L_SEX_P XXY_AMYG_L_SEX_P ...
                          XXY_STR_R_SEX_P XXY_GP_R_SEX_P XXY_THAL_R_SEX_P XXY_HIPP_R_SEX_P XXY_AMYG_R_SEX_P];

XXY.STAT.VOL.TBV.SEX_Q = 10*XXY.STAT.VOL.TBV.SEX_P;

% Fit models (morph)
XXY.STAT.MOR.TSA.STR = SurfStatLinMod(XXY_STR_MOR, XXY.STAT.MOR.TSA.M_STR);
XXY.STAT.MOR.TSA.STR_SEX = SurfStatT(XXY.STAT.MOR.TSA.STR, XXY_SEX_STR);
XXY.STAT.MOR.TSA.STR_SEX_B = XXY.STAT.MOR.TSA.STR_SEX.coef(2,:);
XXY.STAT.MOR.TSA.STR_SEX_P = 2*(1 - tcdf(abs(XXY.STAT.MOR.TSA.STR_SEX.t), XXY.STAT.MOR.TSA.STR_SEX.df));
[~,~,~,XXY.STAT.MOR.TSA.STR_SEX_Q] = fdr_bh(XXY.STAT.MOR.TSA.STR_SEX_P);

XXY.STAT.MOR.TSA.GP = SurfStatLinMod(XXY_GP_MOR, XXY.STAT.MOR.TSA.M_GP);
XXY.STAT.MOR.TSA.GP_SEX = SurfStatT(XXY.STAT.MOR.TSA.GP, XXY_SEX_GP);
XXY.STAT.MOR.TSA.GP_SEX_B = XXY.STAT.MOR.TSA.GP_SEX.coef(2,:);
XXY.STAT.MOR.TSA.GP_SEX_P = 2*(1 - tcdf(abs(XXY.STAT.MOR.TSA.GP_SEX.t), XXY.STAT.MOR.TSA.GP_SEX.df));
[~,~,~,XXY.STAT.MOR.TSA.GP_SEX_Q] = fdr_bh(XXY.STAT.MOR.TSA.GP_SEX_P);

XXY.STAT.MOR.TSA.THAL = SurfStatLinMod(XXY_THAL_MOR, XXY.STAT.MOR.TSA.M_THAL);
XXY.STAT.MOR.TSA.THAL_SEX = SurfStatT(XXY.STAT.MOR.TSA.THAL, XXY_SEX_THAL);
XXY.STAT.MOR.TSA.THAL_SEX_B = XXY.STAT.MOR.TSA.THAL_SEX.coef(2,:);
XXY.STAT.MOR.TSA.THAL_SEX_P = 2*(1 - tcdf(abs(XXY.STAT.MOR.TSA.THAL_SEX.t), XXY.STAT.MOR.TSA.THAL_SEX.df));
[~,~,~,XXY.STAT.MOR.TSA.THAL_SEX_Q] = fdr_bh(XXY.STAT.MOR.TSA.THAL_SEX_P);

XXY.STAT.MOR.TSA.HIPP = SurfStatLinMod(XXY_HIPP_MOR, XXY.STAT.MOR.TSA.M_HIPP);
XXY.STAT.MOR.TSA.HIPP_SEX = SurfStatT(XXY.STAT.MOR.TSA.HIPP, XXY_SEX_HIPP);
XXY.STAT.MOR.TSA.HIPP_SEX_B = XXY.STAT.MOR.TSA.HIPP_SEX.coef(2,:);
XXY.STAT.MOR.TSA.HIPP_SEX_P = 2*(1 - tcdf(abs(XXY.STAT.MOR.TSA.HIPP_SEX.t), XXY.STAT.MOR.TSA.HIPP_SEX.df));
[~,~,~,XXY.STAT.MOR.TSA.HIPP_SEX_Q] = fdr_bh(XXY.STAT.MOR.TSA.HIPP_SEX_P);

XXY.STAT.MOR.TSA.AMYG = SurfStatLinMod(XXY_AMYG_MOR, XXY.STAT.MOR.TSA.M_AMYG);
XXY.STAT.MOR.TSA.AMYG_SEX = SurfStatT(XXY.STAT.MOR.TSA.AMYG, XXY_SEX_AMYG);
XXY.STAT.MOR.TSA.AMYG_SEX_B = XXY.STAT.MOR.TSA.AMYG_SEX.coef(2,:);
XXY.STAT.MOR.TSA.AMYG_SEX_P = 2*(1 - tcdf(abs(XXY.STAT.MOR.TSA.AMYG_SEX.t), XXY.STAT.MOR.TSA.AMYG_SEX.df));
[~,~,~,XXY.STAT.MOR.TSA.AMYG_SEX_Q] = fdr_bh(XXY.STAT.MOR.TSA.AMYG_SEX_P);

%% XYY

% Compute z-scores (struct)
XYY_STR_L = zscore(XYY.VOL.STR_L(XYY.QC.indSTR))';
XYY_STR_R = zscore(XYY.VOL.STR_R(XYY.QC.indSTR))';

XYY_GP_L = zscore(XYY.VOL.GP_L(XYY.QC.indGP))';
XYY_GP_R = zscore(XYY.VOL.GP_R(XYY.QC.indGP))';

XYY_THAL_L = zscore(XYY.VOL.THAL_L(XYY.QC.indTHAL))';
XYY_THAL_R = zscore(XYY.VOL.THAL_R(XYY.QC.indTHAL))';

XYY_HIPP_L = zscore(XYY.VOL.HIPP_L(XYY.QC.indHIPP))';
XYY_HIPP_R = zscore(XYY.VOL.HIPP_R(XYY.QC.indHIPP))';

XYY_AMYG_L = zscore(XYY.VOL.AMYG_L(XYY.QC.indAMYG))';
XYY_AMYG_R = zscore(XYY.VOL.AMYG_R(XYY.QC.indAMYG))';

% Compute z-scores (morph)
XYY_STR_MOR = zscore(XYY.MOR.STR(XYY.QC.indSTR,:));
XYY_GP_MOR = zscore(XYY.MOR.GP(XYY.QC.indGP,:));
XYY_THAL_MOR = zscore(XYY.MOR.THAL(XYY.QC.indTHAL,:));
XYY_HIPP_MOR = zscore(XYY.MOR.HIPP(XYY.QC.indHIPP,:));
XYY_AMYG_MOR = zscore(XYY.MOR.AMYG(XYY.QC.indAMYG,:));

% Prepare terms (struct)
XYY_SEX_STR = term(XYY.SEX(XYY.QC.indSTR));
XYY_AGE_STR = term(XYY.AGE(XYY.QC.indSTR) - mean(XYY.AGE(XYY.QC.indSTR)));
XYY_TBV_STR = term(XYY.TBV(XYY.QC.indSTR) - mean(XYY.TBV(XYY.QC.indSTR)));

XYY_SEX_GP = term(XYY.SEX(XYY.QC.indGP));
XYY_AGE_GP = term(XYY.AGE(XYY.QC.indGP) - mean(XYY.AGE(XYY.QC.indGP)));
XYY_TBV_GP = term(XYY.TBV(XYY.QC.indGP) - mean(XYY.TBV(XYY.QC.indGP)));

XYY_SEX_THAL = term(XYY.SEX(XYY.QC.indTHAL));
XYY_AGE_THAL = term(XYY.AGE(XYY.QC.indTHAL) - mean(XYY.AGE(XYY.QC.indTHAL)));
XYY_TBV_THAL = term(XYY.TBV(XYY.QC.indTHAL) - mean(XYY.TBV(XYY.QC.indTHAL)));

XYY_SEX_HIPP = term(XYY.SEX(XYY.QC.indHIPP));
XYY_AGE_HIPP = term(XYY.AGE(XYY.QC.indHIPP) - mean(XYY.AGE(XYY.QC.indHIPP)));
XYY_TBV_HIPP = term(XYY.TBV(XYY.QC.indHIPP) - mean(XYY.TBV(XYY.QC.indHIPP)));

XYY_SEX_AMYG = term(XYY.SEX(XYY.QC.indAMYG));
XYY_AGE_AMYG = term(XYY.AGE(XYY.QC.indAMYG) - mean(XYY.AGE(XYY.QC.indAMYG)));
XYY_TBV_AMYG = term(XYY.TBV(XYY.QC.indAMYG) - mean(XYY.TBV(XYY.QC.indAMYG)));

% Prepare terms (morph)
XYY_TSA_STR = term(sum(XYY.MOR.STR(XYY.QC.indSTR,:),2) - mean(sum(XYY.MOR.STR(XYY.QC.indSTR,:),2)));
XYY_TSA_GP = term(sum(XYY.MOR.GP(XYY.QC.indGP,:),2) - mean(sum(XYY.MOR.GP(XYY.QC.indGP,:),2)));
XYY_TSA_THAL = term(sum(XYY.MOR.THAL(XYY.QC.indTHAL,:),2) - mean(sum(XYY.MOR.THAL(XYY.QC.indTHAL,:),2)));
XYY_TSA_HIPP = term(sum(XYY.MOR.HIPP(XYY.QC.indHIPP,:),2) - mean(sum(XYY.MOR.HIPP(XYY.QC.indHIPP,:),2)));
XYY_TSA_AMYG = term(sum(XYY.MOR.AMYG(XYY.QC.indAMYG,:),2) - mean(sum(XYY.MOR.AMYG(XYY.QC.indAMYG,:),2)));

% Define GLMs (struct)
XYY.STAT.VOL.TBV.M_STR = 1 + XYY_SEX_STR + XYY_AGE_STR + XYY_TBV_STR;
XYY.STAT.VOL.TBV.M_GP = 1 + XYY_SEX_GP + XYY_AGE_GP + XYY_TBV_GP;
XYY.STAT.VOL.TBV.M_THAL = 1 + XYY_SEX_THAL + XYY_AGE_THAL + XYY_TBV_THAL;
XYY.STAT.VOL.TBV.M_HIPP = 1 + XYY_SEX_HIPP + XYY_AGE_HIPP + XYY_TBV_HIPP;
XYY.STAT.VOL.TBV.M_AMYG = 1 + XYY_SEX_AMYG + XYY_AGE_AMYG + XYY_TBV_AMYG;

% Define GLMs (morph)
XYY.STAT.MOR.TSA.M_STR = 1 + XYY_SEX_STR + XYY_AGE_STR + XYY_TSA_STR;
XYY.STAT.MOR.TSA.M_GP = 1 + XYY_SEX_GP + XYY_AGE_GP + XYY_TSA_GP;
XYY.STAT.MOR.TSA.M_THAL = 1 + XYY_SEX_THAL + XYY_AGE_THAL + XYY_TSA_THAL;
XYY.STAT.MOR.TSA.M_HIPP = 1 + XYY_SEX_HIPP + XYY_AGE_HIPP + XYY_TSA_HIPP;
XYY.STAT.MOR.TSA.M_AMYG = 1 + XYY_SEX_AMYG + XYY_AGE_AMYG + XYY_TSA_AMYG;

% Fit models (struct)
XYY.STAT.VOL.TBV.STR_L = SurfStatLinMod(XYY_STR_L, XYY.STAT.VOL.TBV.M_STR);
XYY.STAT.VOL.TBV.STR_L_SEX = SurfStatT(XYY.STAT.VOL.TBV.STR_L, XYY_SEX_STR);
XYY_STR_L_SEX_P = 2*(1 - tcdf(abs(XYY.STAT.VOL.TBV.STR_L_SEX.t), XYY.STAT.VOL.TBV.STR_L_SEX.df));

XYY.STAT.VOL.TBV.STR_R = SurfStatLinMod(XYY_STR_R, XYY.STAT.VOL.TBV.M_STR);
XYY.STAT.VOL.TBV.STR_R_SEX = SurfStatT(XYY.STAT.VOL.TBV.STR_R, XYY_SEX_STR);
XYY_STR_R_SEX_P = 2*(1 - tcdf(abs(XYY.STAT.VOL.TBV.STR_R_SEX.t), XYY.STAT.VOL.TBV.STR_R_SEX.df));

XYY.STAT.VOL.TBV.GP_L = SurfStatLinMod(XYY_GP_L, XYY.STAT.VOL.TBV.M_GP);
XYY.STAT.VOL.TBV.GP_L_SEX = SurfStatT(XYY.STAT.VOL.TBV.GP_L, XYY_SEX_GP);
XYY_GP_L_SEX_P = 2*(1 - tcdf(abs(XYY.STAT.VOL.TBV.GP_L_SEX.t), XYY.STAT.VOL.TBV.GP_L_SEX.df));

XYY.STAT.VOL.TBV.GP_R = SurfStatLinMod(XYY_GP_R, XYY.STAT.VOL.TBV.M_GP);
XYY.STAT.VOL.TBV.GP_R_SEX = SurfStatT(XYY.STAT.VOL.TBV.GP_R, XYY_SEX_GP);
XYY_GP_R_SEX_P = 2*(1 - tcdf(abs(XYY.STAT.VOL.TBV.GP_R_SEX.t), XYY.STAT.VOL.TBV.GP_R_SEX.df));

XYY.STAT.VOL.TBV.THAL_L = SurfStatLinMod(XYY_THAL_L, XYY.STAT.VOL.TBV.M_THAL);
XYY.STAT.VOL.TBV.THAL_L_SEX = SurfStatT(XYY.STAT.VOL.TBV.THAL_L, XYY_SEX_THAL);
XYY_THAL_L_SEX_P = 2*(1 - tcdf(abs(XYY.STAT.VOL.TBV.THAL_L_SEX.t), XYY.STAT.VOL.TBV.THAL_L_SEX.df));

XYY.STAT.VOL.TBV.THAL_R = SurfStatLinMod(XYY_THAL_R, XYY.STAT.VOL.TBV.M_THAL);
XYY.STAT.VOL.TBV.THAL_R_SEX = SurfStatT(XYY.STAT.VOL.TBV.THAL_R, XYY_SEX_THAL);
XYY_THAL_R_SEX_P = 2*(1 - tcdf(abs(XYY.STAT.VOL.TBV.THAL_R_SEX.t), XYY.STAT.VOL.TBV.THAL_R_SEX.df));

XYY.STAT.VOL.TBV.HIPP_L = SurfStatLinMod(XYY_HIPP_L, XYY.STAT.VOL.TBV.M_HIPP);
XYY.STAT.VOL.TBV.HIPP_L_SEX = SurfStatT(XYY.STAT.VOL.TBV.HIPP_L, XYY_SEX_HIPP);
XYY_HIPP_L_SEX_P = 2*(1 - tcdf(abs(XYY.STAT.VOL.TBV.HIPP_L_SEX.t), XYY.STAT.VOL.TBV.HIPP_L_SEX.df));

XYY.STAT.VOL.TBV.HIPP_R = SurfStatLinMod(XYY_HIPP_R, XYY.STAT.VOL.TBV.M_HIPP);
XYY.STAT.VOL.TBV.HIPP_R_SEX = SurfStatT(XYY.STAT.VOL.TBV.HIPP_R, XYY_SEX_HIPP);
XYY_HIPP_R_SEX_P = 2*(1 - tcdf(abs(XYY.STAT.VOL.TBV.HIPP_R_SEX.t), XYY.STAT.VOL.TBV.HIPP_R_SEX.df));

XYY.STAT.VOL.TBV.AMYG_L = SurfStatLinMod(XYY_AMYG_L, XYY.STAT.VOL.TBV.M_AMYG);
XYY.STAT.VOL.TBV.AMYG_L_SEX = SurfStatT(XYY.STAT.VOL.TBV.AMYG_L, XYY_SEX_AMYG);
XYY_AMYG_L_SEX_P = 2*(1 - tcdf(abs(XYY.STAT.VOL.TBV.AMYG_L_SEX.t), XYY.STAT.VOL.TBV.AMYG_L_SEX.df));

XYY.STAT.VOL.TBV.AMYG_R = SurfStatLinMod(XYY_AMYG_R, XYY.STAT.VOL.TBV.M_AMYG);
XYY.STAT.VOL.TBV.AMYG_R_SEX = SurfStatT(XYY.STAT.VOL.TBV.AMYG_R, XYY_SEX_AMYG);
XYY_AMYG_R_SEX_P = 2*(1 - tcdf(abs(XYY.STAT.VOL.TBV.AMYG_R_SEX.t), XYY.STAT.VOL.TBV.AMYG_R_SEX.df));

% Collate statistical values
XYY.STAT.VOL.TBV.SEX_B = [XYY.STAT.VOL.TBV.STR_L_SEX.coef(2) XYY.STAT.VOL.TBV.GP_L_SEX.coef(2) XYY.STAT.VOL.TBV.THAL_L_SEX.coef(2) XYY.STAT.VOL.TBV.HIPP_L_SEX.coef(2) XYY.STAT.VOL.TBV.AMYG_L_SEX.coef(2) ...
                          XYY.STAT.VOL.TBV.STR_R_SEX.coef(2) XYY.STAT.VOL.TBV.GP_R_SEX.coef(2) XYY.STAT.VOL.TBV.THAL_R_SEX.coef(2) XYY.STAT.VOL.TBV.HIPP_R_SEX.coef(2) XYY.STAT.VOL.TBV.AMYG_R_SEX.coef(2)];

XYY.STAT.VOL.TBV.SEX_P = [XYY_STR_L_SEX_P XYY_GP_L_SEX_P XYY_THAL_L_SEX_P XYY_HIPP_L_SEX_P XYY_AMYG_L_SEX_P ...
                          XYY_STR_R_SEX_P XYY_GP_R_SEX_P XYY_THAL_R_SEX_P XYY_HIPP_R_SEX_P XYY_AMYG_R_SEX_P];

XYY.STAT.VOL.TBV.SEX_Q = 10*XYY.STAT.VOL.TBV.SEX_P;

% Fit models (morph)
XYY.STAT.MOR.TSA.STR = SurfStatLinMod(XYY_STR_MOR, XYY.STAT.MOR.TSA.M_STR);
XYY.STAT.MOR.TSA.STR_SEX = SurfStatT(XYY.STAT.MOR.TSA.STR, XYY_SEX_STR);
XYY.STAT.MOR.TSA.STR_SEX_B = XYY.STAT.MOR.TSA.STR_SEX.coef(2,:);
XYY.STAT.MOR.TSA.STR_SEX_P = 2*(1 - tcdf(abs(XYY.STAT.MOR.TSA.STR_SEX.t), XYY.STAT.MOR.TSA.STR_SEX.df));
[~,~,~,XYY.STAT.MOR.TSA.STR_SEX_Q] = fdr_bh(XYY.STAT.MOR.TSA.STR_SEX_P);

XYY.STAT.MOR.TSA.GP = SurfStatLinMod(XYY_GP_MOR, XYY.STAT.MOR.TSA.M_GP);
XYY.STAT.MOR.TSA.GP_SEX = SurfStatT(XYY.STAT.MOR.TSA.GP, XYY_SEX_GP);
XYY.STAT.MOR.TSA.GP_SEX_B = XYY.STAT.MOR.TSA.GP_SEX.coef(2,:);
XYY.STAT.MOR.TSA.GP_SEX_P = 2*(1 - tcdf(abs(XYY.STAT.MOR.TSA.GP_SEX.t), XYY.STAT.MOR.TSA.GP_SEX.df));
[~,~,~,XYY.STAT.MOR.TSA.GP_SEX_Q] = fdr_bh(XYY.STAT.MOR.TSA.GP_SEX_P);

XYY.STAT.MOR.TSA.THAL = SurfStatLinMod(XYY_THAL_MOR, XYY.STAT.MOR.TSA.M_THAL);
XYY.STAT.MOR.TSA.THAL_SEX = SurfStatT(XYY.STAT.MOR.TSA.THAL, XYY_SEX_THAL);
XYY.STAT.MOR.TSA.THAL_SEX_B = XYY.STAT.MOR.TSA.THAL_SEX.coef(2,:);
XYY.STAT.MOR.TSA.THAL_SEX_P = 2*(1 - tcdf(abs(XYY.STAT.MOR.TSA.THAL_SEX.t), XYY.STAT.MOR.TSA.THAL_SEX.df));
[~,~,~,XYY.STAT.MOR.TSA.THAL_SEX_Q] = fdr_bh(XYY.STAT.MOR.TSA.THAL_SEX_P);

XYY.STAT.MOR.TSA.HIPP = SurfStatLinMod(XYY_HIPP_MOR, XYY.STAT.MOR.TSA.M_HIPP);
XYY.STAT.MOR.TSA.HIPP_SEX = SurfStatT(XYY.STAT.MOR.TSA.HIPP, XYY_SEX_HIPP);
XYY.STAT.MOR.TSA.HIPP_SEX_B = XYY.STAT.MOR.TSA.HIPP_SEX.coef(2,:);
XYY.STAT.MOR.TSA.HIPP_SEX_P = 2*(1 - tcdf(abs(XYY.STAT.MOR.TSA.HIPP_SEX.t), XYY.STAT.MOR.TSA.HIPP_SEX.df));
[~,~,~,XYY.STAT.MOR.TSA.HIPP_SEX_Q] = fdr_bh(XYY.STAT.MOR.TSA.HIPP_SEX_P);

XYY.STAT.MOR.TSA.AMYG = SurfStatLinMod(XYY_AMYG_MOR, XYY.STAT.MOR.TSA.M_AMYG);
XYY.STAT.MOR.TSA.AMYG_SEX = SurfStatT(XYY.STAT.MOR.TSA.AMYG, XYY_SEX_AMYG);
XYY.STAT.MOR.TSA.AMYG_SEX_B = XYY.STAT.MOR.TSA.AMYG_SEX.coef(2,:);
XYY.STAT.MOR.TSA.AMYG_SEX_P = 2*(1 - tcdf(abs(XYY.STAT.MOR.TSA.AMYG_SEX.t), XYY.STAT.MOR.TSA.AMYG_SEX.df));
[~,~,~,XYY.STAT.MOR.TSA.AMYG_SEX_Q] = fdr_bh(XYY.STAT.MOR.TSA.AMYG_SEX_P);

%% IGD

% Compute z-scores (struct)
IGD_STR_L = zscore(IGD.VOL.STR_L(IGD.QC.indSTR))';
IGD_STR_R = zscore(IGD.VOL.STR_R(IGD.QC.indSTR))';

IGD_GP_L = zscore(IGD.VOL.GP_L(IGD.QC.indGP))';
IGD_GP_R = zscore(IGD.VOL.GP_R(IGD.QC.indGP))';

IGD_THAL_L = zscore(IGD.VOL.THAL_L(IGD.QC.indTHAL))';
IGD_THAL_R = zscore(IGD.VOL.THAL_R(IGD.QC.indTHAL))';

IGD_HIPP_L = zscore(IGD.VOL.HIPP_L(IGD.QC.indHIPP))';
IGD_HIPP_R = zscore(IGD.VOL.HIPP_R(IGD.QC.indHIPP))';

IGD_AMYG_L = zscore(IGD.VOL.AMYG_L(IGD.QC.indAMYG))';
IGD_AMYG_R = zscore(IGD.VOL.AMYG_R(IGD.QC.indAMYG))';

% Compute z-scores (morph)
IGD_STR_MOR = zscore(IGD.MOR.STR(IGD.QC.indSTR,:));
IGD_GP_MOR = zscore(IGD.MOR.GP(IGD.QC.indGP,:));
IGD_THAL_MOR = zscore(IGD.MOR.THAL(IGD.QC.indTHAL,:));
IGD_HIPP_MOR = zscore(IGD.MOR.HIPP(IGD.QC.indHIPP,:));
IGD_AMYG_MOR = zscore(IGD.MOR.AMYG(IGD.QC.indAMYG,:));

% Prepare terms (struct)
IGD_GRP_STR = term(IGD.GRP(IGD.QC.indSTR));
IGD_AGE_STR = term(IGD.AGE(IGD.QC.indSTR) - mean(IGD.AGE(IGD.QC.indSTR)));
IGD_TBV_STR = term(IGD.TBV(IGD.QC.indSTR) - mean(IGD.TBV(IGD.QC.indSTR)));

IGD_GRP_GP = term(IGD.GRP(IGD.QC.indGP));
IGD_AGE_GP = term(IGD.AGE(IGD.QC.indGP) - mean(IGD.AGE(IGD.QC.indGP)));
IGD_TBV_GP = term(IGD.TBV(IGD.QC.indGP) - mean(IGD.TBV(IGD.QC.indGP)));

IGD_GRP_THAL = term(IGD.GRP(IGD.QC.indTHAL));
IGD_AGE_THAL = term(IGD.AGE(IGD.QC.indTHAL) - mean(IGD.AGE(IGD.QC.indTHAL)));
IGD_TBV_THAL = term(IGD.TBV(IGD.QC.indTHAL) - mean(IGD.TBV(IGD.QC.indTHAL)));

IGD_GRP_HIPP = term(IGD.GRP(IGD.QC.indHIPP));
IGD_AGE_HIPP = term(IGD.AGE(IGD.QC.indHIPP) - mean(IGD.AGE(IGD.QC.indHIPP)));
IGD_TBV_HIPP = term(IGD.TBV(IGD.QC.indHIPP) - mean(IGD.TBV(IGD.QC.indHIPP)));

IGD_GRP_AMYG = term(IGD.GRP(IGD.QC.indAMYG));
IGD_AGE_AMYG = term(IGD.AGE(IGD.QC.indAMYG) - mean(IGD.AGE(IGD.QC.indAMYG)));
IGD_TBV_AMYG = term(IGD.TBV(IGD.QC.indAMYG) - mean(IGD.TBV(IGD.QC.indAMYG)));

% Prepare terms (morph)
IGD_TSA_STR = term(sum(IGD.MOR.STR(IGD.QC.indSTR,:),2) - mean(sum(IGD.MOR.STR(IGD.QC.indSTR,:),2)));
IGD_TSA_GP = term(sum(IGD.MOR.GP(IGD.QC.indGP,:),2) - mean(sum(IGD.MOR.GP(IGD.QC.indGP,:),2)));
IGD_TSA_THAL = term(sum(IGD.MOR.THAL(IGD.QC.indTHAL,:),2) - mean(sum(IGD.MOR.THAL(IGD.QC.indTHAL,:),2)));
IGD_TSA_HIPP = term(sum(IGD.MOR.HIPP(IGD.QC.indHIPP,:),2) - mean(sum(IGD.MOR.HIPP(IGD.QC.indHIPP,:),2)));
IGD_TSA_AMYG = term(sum(IGD.MOR.AMYG(IGD.QC.indAMYG,:),2) - mean(sum(IGD.MOR.AMYG(IGD.QC.indAMYG,:),2)));

% Define GLMs (struct)
IGD.STAT.VOL.TBV.M_STR = 1 + IGD_GRP_STR + IGD_AGE_STR + IGD_TBV_STR;
IGD.STAT.VOL.TBV.M_GP = 1 + IGD_GRP_GP + IGD_AGE_GP + IGD_TBV_GP;
IGD.STAT.VOL.TBV.M_THAL = 1 + IGD_GRP_THAL + IGD_AGE_THAL + IGD_TBV_THAL;
IGD.STAT.VOL.TBV.M_HIPP = 1 + IGD_GRP_HIPP + IGD_AGE_HIPP + IGD_TBV_HIPP;
IGD.STAT.VOL.TBV.M_AMYG = 1 + IGD_GRP_AMYG + IGD_AGE_AMYG + IGD_TBV_AMYG;

% Define GLMs (morph)
IGD.STAT.MOR.TSA.M_STR = 1 + IGD_GRP_STR + IGD_AGE_STR + IGD_TSA_STR;
IGD.STAT.MOR.TSA.M_GP = 1 + IGD_GRP_GP + IGD_AGE_GP + IGD_TSA_GP;
IGD.STAT.MOR.TSA.M_THAL = 1 + IGD_GRP_THAL + IGD_AGE_THAL + IGD_TSA_THAL;
IGD.STAT.MOR.TSA.M_HIPP = 1 + IGD_GRP_HIPP + IGD_AGE_HIPP + IGD_TSA_HIPP;
IGD.STAT.MOR.TSA.M_AMYG = 1 + IGD_GRP_AMYG + IGD_AGE_AMYG + IGD_TSA_AMYG;

% Fit models (struct)
IGD.STAT.VOL.TBV.STR_L = SurfStatLinMod(IGD_STR_L, IGD.STAT.VOL.TBV.M_STR);
IGD.STAT.VOL.TBV.STR_L_GRP = SurfStatT(IGD.STAT.VOL.TBV.STR_L, IGD_GRP_STR);
IGD_STR_L_GRP_P = 2*(1 - tcdf(abs(IGD.STAT.VOL.TBV.STR_L_GRP.t), IGD.STAT.VOL.TBV.STR_L_GRP.df));

IGD.STAT.VOL.TBV.STR_R = SurfStatLinMod(IGD_STR_R, IGD.STAT.VOL.TBV.M_STR);
IGD.STAT.VOL.TBV.STR_R_GRP = SurfStatT(IGD.STAT.VOL.TBV.STR_R, IGD_GRP_STR);
IGD_STR_R_GRP_P = 2*(1 - tcdf(abs(IGD.STAT.VOL.TBV.STR_R_GRP.t), IGD.STAT.VOL.TBV.STR_R_GRP.df));

IGD.STAT.VOL.TBV.GP_L = SurfStatLinMod(IGD_GP_L, IGD.STAT.VOL.TBV.M_GP);
IGD.STAT.VOL.TBV.GP_L_GRP = SurfStatT(IGD.STAT.VOL.TBV.GP_L, IGD_GRP_GP);
IGD_GP_L_GRP_P = 2*(1 - tcdf(abs(IGD.STAT.VOL.TBV.GP_L_GRP.t), IGD.STAT.VOL.TBV.GP_L_GRP.df));

IGD.STAT.VOL.TBV.GP_R = SurfStatLinMod(IGD_GP_R, IGD.STAT.VOL.TBV.M_GP);
IGD.STAT.VOL.TBV.GP_R_GRP = SurfStatT(IGD.STAT.VOL.TBV.GP_R, IGD_GRP_GP);
IGD_GP_R_GRP_P = 2*(1 - tcdf(abs(IGD.STAT.VOL.TBV.GP_R_GRP.t), IGD.STAT.VOL.TBV.GP_R_GRP.df));

IGD.STAT.VOL.TBV.THAL_L = SurfStatLinMod(IGD_THAL_L, IGD.STAT.VOL.TBV.M_THAL);
IGD.STAT.VOL.TBV.THAL_L_GRP = SurfStatT(IGD.STAT.VOL.TBV.THAL_L, IGD_GRP_THAL);
IGD_THAL_L_GRP_P = 2*(1 - tcdf(abs(IGD.STAT.VOL.TBV.THAL_L_GRP.t), IGD.STAT.VOL.TBV.THAL_L_GRP.df));

IGD.STAT.VOL.TBV.THAL_R = SurfStatLinMod(IGD_THAL_R, IGD.STAT.VOL.TBV.M_THAL);
IGD.STAT.VOL.TBV.THAL_R_GRP = SurfStatT(IGD.STAT.VOL.TBV.THAL_R, IGD_GRP_THAL);
IGD_THAL_R_GRP_P = 2*(1 - tcdf(abs(IGD.STAT.VOL.TBV.THAL_R_GRP.t), IGD.STAT.VOL.TBV.THAL_R_GRP.df));

IGD.STAT.VOL.TBV.HIPP_L = SurfStatLinMod(IGD_HIPP_L, IGD.STAT.VOL.TBV.M_HIPP);
IGD.STAT.VOL.TBV.HIPP_L_GRP = SurfStatT(IGD.STAT.VOL.TBV.HIPP_L, IGD_GRP_HIPP);
IGD_HIPP_L_GRP_P = 2*(1 - tcdf(abs(IGD.STAT.VOL.TBV.HIPP_L_GRP.t), IGD.STAT.VOL.TBV.HIPP_L_GRP.df));

IGD.STAT.VOL.TBV.HIPP_R = SurfStatLinMod(IGD_HIPP_R, IGD.STAT.VOL.TBV.M_HIPP);
IGD.STAT.VOL.TBV.HIPP_R_GRP = SurfStatT(IGD.STAT.VOL.TBV.HIPP_R, IGD_GRP_HIPP);
IGD_HIPP_R_GRP_P = 2*(1 - tcdf(abs(IGD.STAT.VOL.TBV.HIPP_R_GRP.t), IGD.STAT.VOL.TBV.HIPP_R_GRP.df));

IGD.STAT.VOL.TBV.AMYG_L = SurfStatLinMod(IGD_AMYG_L, IGD.STAT.VOL.TBV.M_AMYG);
IGD.STAT.VOL.TBV.AMYG_L_GRP = SurfStatT(IGD.STAT.VOL.TBV.AMYG_L, IGD_GRP_AMYG);
IGD_AMYG_L_GRP_P = 2*(1 - tcdf(abs(IGD.STAT.VOL.TBV.AMYG_L_GRP.t), IGD.STAT.VOL.TBV.AMYG_L_GRP.df));

IGD.STAT.VOL.TBV.AMYG_R = SurfStatLinMod(IGD_AMYG_R, IGD.STAT.VOL.TBV.M_AMYG);
IGD.STAT.VOL.TBV.AMYG_R_GRP = SurfStatT(IGD.STAT.VOL.TBV.AMYG_R, IGD_GRP_AMYG);
IGD_AMYG_R_GRP_P = 2*(1 - tcdf(abs(IGD.STAT.VOL.TBV.AMYG_R_GRP.t), IGD.STAT.VOL.TBV.AMYG_R_GRP.df));

% Collate statistical values
IGD.STAT.VOL.TBV.GRP_B = [IGD.STAT.VOL.TBV.STR_L_GRP.coef(2) IGD.STAT.VOL.TBV.GP_L_GRP.coef(2) IGD.STAT.VOL.TBV.THAL_L_GRP.coef(2) IGD.STAT.VOL.TBV.HIPP_L_GRP.coef(2) IGD.STAT.VOL.TBV.AMYG_L_GRP.coef(2) ...
                          IGD.STAT.VOL.TBV.STR_R_GRP.coef(2) IGD.STAT.VOL.TBV.GP_R_GRP.coef(2) IGD.STAT.VOL.TBV.THAL_R_GRP.coef(2) IGD.STAT.VOL.TBV.HIPP_R_GRP.coef(2) IGD.STAT.VOL.TBV.AMYG_R_GRP.coef(2)];

IGD.STAT.VOL.TBV.GRP_P = [IGD_STR_L_GRP_P IGD_GP_L_GRP_P IGD_THAL_L_GRP_P IGD_HIPP_L_GRP_P IGD_AMYG_L_GRP_P ...
                          IGD_STR_R_GRP_P IGD_GP_R_GRP_P IGD_THAL_R_GRP_P IGD_HIPP_R_GRP_P IGD_AMYG_R_GRP_P];

IGD.STAT.VOL.TBV.GRP_Q = 10*IGD.STAT.VOL.TBV.GRP_P;

% Fit models (morph)
IGD.STAT.MOR.TSA.STR = SurfStatLinMod(IGD_STR_MOR, IGD.STAT.MOR.TSA.M_STR);
IGD.STAT.MOR.TSA.STR_GRP = SurfStatT(IGD.STAT.MOR.TSA.STR, IGD_GRP_STR);
IGD.STAT.MOR.TSA.STR_GRP_B = IGD.STAT.MOR.TSA.STR_GRP.coef(2,:);
IGD.STAT.MOR.TSA.STR_GRP_P = 2*(1 - tcdf(abs(IGD.STAT.MOR.TSA.STR_GRP.t), IGD.STAT.MOR.TSA.STR_GRP.df));
[~,~,~,IGD.STAT.MOR.TSA.STR_GRP_Q] = fdr_bh(IGD.STAT.MOR.TSA.STR_GRP_P);

IGD.STAT.MOR.TSA.GP = SurfStatLinMod(IGD_GP_MOR, IGD.STAT.MOR.TSA.M_GP);
IGD.STAT.MOR.TSA.GP_GRP = SurfStatT(IGD.STAT.MOR.TSA.GP, IGD_GRP_GP);
IGD.STAT.MOR.TSA.GP_GRP_B = IGD.STAT.MOR.TSA.GP_GRP.coef(2,:);
IGD.STAT.MOR.TSA.GP_GRP_P = 2*(1 - tcdf(abs(IGD.STAT.MOR.TSA.GP_GRP.t), IGD.STAT.MOR.TSA.GP_GRP.df));
[~,~,~,IGD.STAT.MOR.TSA.GP_GRP_Q] = fdr_bh(IGD.STAT.MOR.TSA.GP_GRP_P);

IGD.STAT.MOR.TSA.THAL = SurfStatLinMod(IGD_THAL_MOR, IGD.STAT.MOR.TSA.M_THAL);
IGD.STAT.MOR.TSA.THAL_GRP = SurfStatT(IGD.STAT.MOR.TSA.THAL, IGD_GRP_THAL);
IGD.STAT.MOR.TSA.THAL_GRP_B = IGD.STAT.MOR.TSA.THAL_GRP.coef(2,:);
IGD.STAT.MOR.TSA.THAL_GRP_P = 2*(1 - tcdf(abs(IGD.STAT.MOR.TSA.THAL_GRP.t), IGD.STAT.MOR.TSA.THAL_GRP.df));
[~,~,~,IGD.STAT.MOR.TSA.THAL_GRP_Q] = fdr_bh(IGD.STAT.MOR.TSA.THAL_GRP_P);

IGD.STAT.MOR.TSA.HIPP = SurfStatLinMod(IGD_HIPP_MOR, IGD.STAT.MOR.TSA.M_HIPP);
IGD.STAT.MOR.TSA.HIPP_GRP = SurfStatT(IGD.STAT.MOR.TSA.HIPP, IGD_GRP_HIPP);
IGD.STAT.MOR.TSA.HIPP_GRP_B = IGD.STAT.MOR.TSA.HIPP_GRP.coef(2,:);
IGD.STAT.MOR.TSA.HIPP_GRP_P = 2*(1 - tcdf(abs(IGD.STAT.MOR.TSA.HIPP_GRP.t), IGD.STAT.MOR.TSA.HIPP_GRP.df));
[~,~,~,IGD.STAT.MOR.TSA.HIPP_GRP_Q] = fdr_bh(IGD.STAT.MOR.TSA.HIPP_GRP_P);

IGD.STAT.MOR.TSA.AMYG = SurfStatLinMod(IGD_AMYG_MOR, IGD.STAT.MOR.TSA.M_AMYG);
IGD.STAT.MOR.TSA.AMYG_GRP = SurfStatT(IGD.STAT.MOR.TSA.AMYG, IGD_GRP_AMYG);
IGD.STAT.MOR.TSA.AMYG_GRP_B = IGD.STAT.MOR.TSA.AMYG_GRP.coef(2,:);
IGD.STAT.MOR.TSA.AMYG_GRP_P = 2*(1 - tcdf(abs(IGD.STAT.MOR.TSA.AMYG_GRP.t), IGD.STAT.MOR.TSA.AMYG_GRP.df));
[~,~,~,IGD.STAT.MOR.TSA.AMYG_GRP_Q] = fdr_bh(IGD.STAT.MOR.TSA.AMYG_GRP_P);

%% Save results

save(fullfile(dataDir, 'XXY_results.mat'), 'XXY', '-v7.3');
save(fullfile(dataDir, 'XYY_results.mat'), 'XYY', '-v7.3');
save(fullfile(dataDir, 'IGD_results.mat'), 'IGD', '-v7.3');

%% Append results from the cortex study

% Load cortical group maps for the clinical cohorts (Lee et al., 2026).
% Loaded into structs so that variables in these files do not overwrite the
% subcortical XXY/XYY/IGD structs in the workspace.
CORT_SCA = load(fullfile(corticalResultsDir, 'sca_groupmaps_final.mat'));
CORT_IGD = load(fullfile(corticalResultsDir, 'igd_groupmaps_final.mat'));

% Append cortical stats
XXY.STAT.COR.CV = CORT_SCA.XXY_CV_SEX;
XXY.STAT.COR.CV_B = CORT_SCA.XXY_CV_SEX.coef(2,:);
XXY.STAT.COR.CV_Q = CORT_SCA.XXY_CV_SEX_Q;
XXY.STAT.COR.SA = CORT_SCA.XXY_SA_SEX;
XXY.STAT.COR.SA_B = CORT_SCA.XXY_SA_SEX.coef(2,:);
XXY.STAT.COR.SA_Q = CORT_SCA.XXY_SA_SEX_Q;
XXY.STAT.COR.CT = CORT_SCA.XXY_CT_SEX;
XXY.STAT.COR.CT_B = CORT_SCA.XXY_CT_SEX.coef(2,:);
XXY.STAT.COR.CT_Q = CORT_SCA.XXY_CT_SEX_Q;

XYY.STAT.COR.CV = CORT_SCA.XYY_CV_SEX;
XYY.STAT.COR.CV_B = CORT_SCA.XYY_CV_SEX.coef(2,:);
XYY.STAT.COR.CV_Q = CORT_SCA.XYY_CV_SEX_Q;
XYY.STAT.COR.SA = CORT_SCA.XYY_SA_SEX;
XYY.STAT.COR.SA_B = CORT_SCA.XYY_SA_SEX.coef(2,:);
XYY.STAT.COR.SA_Q = CORT_SCA.XYY_SA_SEX_Q;
XYY.STAT.COR.CT = CORT_SCA.XYY_CT_SEX;
XYY.STAT.COR.CT_B = CORT_SCA.XYY_CT_SEX.coef(2,:);
XYY.STAT.COR.CT_Q = CORT_SCA.XYY_CT_SEX_Q;

IGD.STAT.COR.CV = CORT_IGD.IGD_CV_GRP;
IGD.STAT.COR.CV_B = CORT_IGD.IGD_CV_GRP.coef(2,:);
IGD.STAT.COR.CV_Q = CORT_IGD.IGD_CV_GRP_Q;
IGD.STAT.COR.SA = CORT_IGD.IGD_SA_GRP;
IGD.STAT.COR.SA_B = CORT_IGD.IGD_SA_GRP.coef(2,:);
IGD.STAT.COR.SA_Q = CORT_IGD.IGD_SA_GRP_Q;
IGD.STAT.COR.CT = CORT_IGD.IGD_CT_GRP;
IGD.STAT.COR.CT_B = CORT_IGD.IGD_CT_GRP.coef(2,:);
IGD.STAT.COR.CT_Q = CORT_IGD.IGD_CT_GRP_Q;

% Save appended results
save(fullfile(dataDir, 'XXY_results.mat'), 'XXY', '-v7.3');
save(fullfile(dataDir, 'XYY_results.mat'), 'XYY', '-v7.3');
save(fullfile(dataDir, 'IGD_results.mat'), 'IGD', '-v7.3');

%% Visualize the anatomical results

% Effect size maps
h = figure('WindowState', 'Maximized'); HyoMagetView(XXY.STAT.VOL.TBV.SEX_B, [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'xxy_vol_tbv_beta.png'), 'Resolution', 300);

% Thresholded
h = figure('WindowState', 'Maximized'); HyoMagetView(XXY.STAT.VOL.TBV.SEX_B .* (XXY.STAT.VOL.TBV.SEX_Q < thresh), [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'xxy_vol_tbv_beta_bonf_thr.png'), 'Resolution', 300);
h = figure('WindowState', 'Maximized'); HyoMagetView(XXY.STAT.VOL.TBV.SEX_B .* (XXY.STAT.VOL.TBV.SEX_P < thresh), [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'xxy_vol_tbv_beta_nom_thr.png'), 'Resolution', 300);

% Now the morph results
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(XXY.STAT.MOR.TSA.STR_SEX_B, [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'xxy_mor_str_tsa_beta.png'), 'Resolution', 300);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(XXY.STAT.MOR.TSA.GP_SEX_B, [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'xxy_mor_gp_tsa_beta.png'), 'Resolution', 300);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(XXY.STAT.MOR.TSA.THAL_SEX_B, [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'xxy_mor_thal_tsa_beta.png'), 'Resolution', 300);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(XXY.STAT.MOR.TSA.HIPP_SEX_B, [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'xxy_mor_hipp_tsa_beta.png'), 'Resolution', 300);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(XXY.STAT.MOR.TSA.AMYG_SEX_B, [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'xxy_mor_amyg_tsa_beta.png'), 'Resolution', 300);

h = figure('WindowState', 'Maximized'); HyoMagetViewNL(XXY.STAT.MOR.TSA.STR_SEX_B .* (XXY.STAT.MOR.TSA.STR_SEX_Q < thresh), [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'xxy_mor_str_tsa_beta_fdr_thr.png'), 'Resolution', 300);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(XXY.STAT.MOR.TSA.GP_SEX_B .* (XXY.STAT.MOR.TSA.GP_SEX_Q < thresh), [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'xxy_mor_gp_tsa_beta_fdr_thr.png'), 'Resolution', 300);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(XXY.STAT.MOR.TSA.THAL_SEX_B .* (XXY.STAT.MOR.TSA.THAL_SEX_Q < thresh), [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'xxy_mor_thal_tsa_beta_fdr_thr.png'), 'Resolution', 300);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(XXY.STAT.MOR.TSA.HIPP_SEX_B .* (XXY.STAT.MOR.TSA.HIPP_SEX_Q < thresh), [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'xxy_mor_hipp_tsa_beta_fdr_thr.png'), 'Resolution', 300);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(XXY.STAT.MOR.TSA.AMYG_SEX_B .* (XXY.STAT.MOR.TSA.AMYG_SEX_Q < thresh), [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'xxy_mor_amyg_tsa_beta_fdr_thr.png'), 'Resolution', 300);

close all;

% Effect size maps
h = figure('WindowState', 'Maximized'); HyoMagetView(XYY.STAT.VOL.TBV.SEX_B, [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'xyy_vol_tbv_beta.png'), 'Resolution', 300);

% Thresholded
h = figure('WindowState', 'Maximized'); HyoMagetView(XYY.STAT.VOL.TBV.SEX_B .* (XYY.STAT.VOL.TBV.SEX_Q < thresh), [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'xyy_vol_tbv_beta_bonf_thr.png'), 'Resolution', 300);
h = figure('WindowState', 'Maximized'); HyoMagetView(XYY.STAT.VOL.TBV.SEX_B .* (XYY.STAT.VOL.TBV.SEX_P < thresh), [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'xyy_vol_tbv_beta_nom_thr.png'), 'Resolution', 300);

% Now the morph results
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(XYY.STAT.MOR.TSA.STR_SEX_B, [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'xyy_mor_str_tsa_beta.png'), 'Resolution', 300);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(XYY.STAT.MOR.TSA.GP_SEX_B, [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'xyy_mor_gp_tsa_beta.png'), 'Resolution', 300);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(XYY.STAT.MOR.TSA.THAL_SEX_B, [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'xyy_mor_thal_tsa_beta.png'), 'Resolution', 300);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(XYY.STAT.MOR.TSA.HIPP_SEX_B, [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'xyy_mor_hipp_tsa_beta.png'), 'Resolution', 300);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(XYY.STAT.MOR.TSA.AMYG_SEX_B, [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'xyy_mor_amyg_tsa_beta.png'), 'Resolution', 300);

h = figure('WindowState', 'Maximized'); HyoMagetViewNL(XYY.STAT.MOR.TSA.STR_SEX_B .* (XYY.STAT.MOR.TSA.STR_SEX_Q < thresh), [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'xyy_mor_str_tsa_beta_fdr_thr.png'), 'Resolution', 300);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(XYY.STAT.MOR.TSA.GP_SEX_B .* (XYY.STAT.MOR.TSA.GP_SEX_Q < thresh), [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'xyy_mor_gp_tsa_beta_fdr_thr.png'), 'Resolution', 300);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(XYY.STAT.MOR.TSA.THAL_SEX_B .* (XYY.STAT.MOR.TSA.THAL_SEX_Q < thresh), [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'xyy_mor_thal_tsa_beta_fdr_thr.png'), 'Resolution', 300);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(XYY.STAT.MOR.TSA.HIPP_SEX_B .* (XYY.STAT.MOR.TSA.HIPP_SEX_Q < thresh), [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'xyy_mor_hipp_tsa_beta_fdr_thr.png'), 'Resolution', 300);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(XYY.STAT.MOR.TSA.AMYG_SEX_B .* (XYY.STAT.MOR.TSA.AMYG_SEX_Q < thresh), [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'xyy_mor_amyg_tsa_beta_fdr_thr.png'), 'Resolution', 300);

close all;

% Effect size maps
h = figure('WindowState', 'Maximized'); HyoMagetView(IGD.STAT.VOL.TBV.GRP_B, [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'igd_vol_tbv_beta.png'), 'Resolution', 300);

% Thresholded
h = figure('WindowState', 'Maximized'); HyoMagetView(IGD.STAT.VOL.TBV.GRP_B .* (IGD.STAT.VOL.TBV.GRP_Q < thresh), [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'igd_vol_tbv_beta_bonf_thr.png'), 'Resolution', 300);
h = figure('WindowState', 'Maximized'); HyoMagetView(IGD.STAT.VOL.TBV.GRP_B .* (IGD.STAT.VOL.TBV.GRP_P < thresh), [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'igd_vol_tbv_beta_nom_thr.png'), 'Resolution', 300);

% Now the morph results
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(IGD.STAT.MOR.TSA.STR_GRP_B, [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'igd_mor_str_tsa_beta.png'), 'Resolution', 300);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(IGD.STAT.MOR.TSA.GP_GRP_B, [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'igd_mor_gp_tsa_beta.png'), 'Resolution', 300);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(IGD.STAT.MOR.TSA.THAL_GRP_B, [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'igd_mor_thal_tsa_beta.png'), 'Resolution', 300);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(IGD.STAT.MOR.TSA.HIPP_GRP_B, [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'igd_mor_hipp_tsa_beta.png'), 'Resolution', 300);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(IGD.STAT.MOR.TSA.AMYG_GRP_B, [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'igd_mor_amyg_tsa_beta.png'), 'Resolution', 300);

h = figure('WindowState', 'Maximized'); HyoMagetViewNL(IGD.STAT.MOR.TSA.STR_GRP_B .* (IGD.STAT.MOR.TSA.STR_GRP_P < thresh), [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'igd_mor_str_tsa_beta_nominal_thr.png'), 'Resolution', 300);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(IGD.STAT.MOR.TSA.GP_GRP_B .* (IGD.STAT.MOR.TSA.GP_GRP_P < thresh), [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'igd_mor_gp_tsa_beta_nominal_thr.png'), 'Resolution', 300);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(IGD.STAT.MOR.TSA.THAL_GRP_B .* (IGD.STAT.MOR.TSA.THAL_GRP_P < thresh), [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'igd_mor_thal_tsa_beta_nominal_thr.png'), 'Resolution', 300);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(IGD.STAT.MOR.TSA.HIPP_GRP_B .* (IGD.STAT.MOR.TSA.HIPP_GRP_P < thresh), [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'igd_mor_hipp_tsa_beta_nominal_thr.png'), 'Resolution', 300);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(IGD.STAT.MOR.TSA.AMYG_GRP_B .* (IGD.STAT.MOR.TSA.AMYG_GRP_P < thresh), [-0.8 0.8], cmap.PuOr); exportgraphics(h, fullfile(figDir, 'igd_mor_amyg_tsa_beta_nominal_thr.png'), 'Resolution', 300);

close all;

% Correlations between effects
r_str_xy = corr(XXY.STAT.MOR.TSA.STR_SEX_B', XYY.STAT.MOR.TSA.STR_SEX_B')
r_str_xi = corr(XXY.STAT.MOR.TSA.STR_SEX_B', IGD.STAT.MOR.TSA.STR_GRP_B')
r_str_yi = corr(XYY.STAT.MOR.TSA.STR_SEX_B', IGD.STAT.MOR.TSA.STR_GRP_B')

r_gp_xy = corr(XXY.STAT.MOR.TSA.GP_SEX_B', XYY.STAT.MOR.TSA.GP_SEX_B')
r_gp_xi = corr(XXY.STAT.MOR.TSA.GP_SEX_B', IGD.STAT.MOR.TSA.GP_GRP_B')
r_gp_yi = corr(XYY.STAT.MOR.TSA.GP_SEX_B', IGD.STAT.MOR.TSA.GP_GRP_B')

r_thal_xy = corr(XXY.STAT.MOR.TSA.THAL_SEX_B', XYY.STAT.MOR.TSA.THAL_SEX_B')
r_thal_xi = corr(XXY.STAT.MOR.TSA.THAL_SEX_B', IGD.STAT.MOR.TSA.THAL_GRP_B')
r_thal_yi = corr(XYY.STAT.MOR.TSA.THAL_SEX_B', IGD.STAT.MOR.TSA.THAL_GRP_B')

r_hipp_xy = corr(XXY.STAT.MOR.TSA.HIPP_SEX_B', XYY.STAT.MOR.TSA.HIPP_SEX_B')
r_hipp_xi = corr(XXY.STAT.MOR.TSA.HIPP_SEX_B', IGD.STAT.MOR.TSA.HIPP_GRP_B')
r_hipp_yi = corr(XYY.STAT.MOR.TSA.HIPP_SEX_B', IGD.STAT.MOR.TSA.HIPP_GRP_B')

r_amyg_xy = corr(XXY.STAT.MOR.TSA.AMYG_SEX_B', XYY.STAT.MOR.TSA.AMYG_SEX_B')
r_amyg_xi = corr(XXY.STAT.MOR.TSA.AMYG_SEX_B', IGD.STAT.MOR.TSA.AMYG_GRP_B')
r_amyg_yi = corr(XYY.STAT.MOR.TSA.AMYG_SEX_B', IGD.STAT.MOR.TSA.AMYG_GRP_B')

% Calculate the colorscale of thresholded maps XXY
temp = XXY.STAT.MOR.TSA.STR_SEX_B .* (XXY.STAT.MOR.TSA.STR_SEX_Q < thresh);
[min(abs(temp(temp<0))) max(abs(temp(temp<0))) min(abs(temp(temp>0))) max(abs(temp(temp>0)))]

temp = XXY.STAT.MOR.TSA.GP_SEX_B .* (XXY.STAT.MOR.TSA.GP_SEX_Q < thresh);
[min(abs(temp(temp<0))) max(abs(temp(temp<0))) min(abs(temp(temp>0))) max(abs(temp(temp>0)))]

temp = XXY.STAT.MOR.TSA.THAL_SEX_B .* (XXY.STAT.MOR.TSA.THAL_SEX_Q < thresh);
[min(abs(temp(temp<0))) max(abs(temp(temp<0))) min(abs(temp(temp>0))) max(abs(temp(temp>0)))]

temp = XXY.STAT.MOR.TSA.HIPP_SEX_B .* (XXY.STAT.MOR.TSA.HIPP_SEX_Q < thresh);
[min(abs(temp(temp<0))) max(abs(temp(temp<0))) min(abs(temp(temp>0))) max(abs(temp(temp>0)))]

temp = XXY.STAT.MOR.TSA.AMYG_SEX_B .* (XXY.STAT.MOR.TSA.AMYG_SEX_Q < thresh);
[min(abs(temp(temp<0))) max(abs(temp(temp<0))) min(abs(temp(temp>0))) max(abs(temp(temp>0)))]

% Calculate the colorscale of thresholded maps XYY
temp = XYY.STAT.MOR.TSA.STR_SEX_B .* (XYY.STAT.MOR.TSA.STR_SEX_Q < thresh);
[min(abs(temp(temp<0))) max(abs(temp(temp<0))) min(abs(temp(temp>0))) max(abs(temp(temp>0)))]

temp = XYY.STAT.MOR.TSA.GP_SEX_B .* (XYY.STAT.MOR.TSA.GP_SEX_Q < thresh);
[min(abs(temp(temp<0))) max(abs(temp(temp<0))) min(abs(temp(temp>0))) max(abs(temp(temp>0)))]

temp = XYY.STAT.MOR.TSA.THAL_SEX_B .* (XYY.STAT.MOR.TSA.THAL_SEX_Q < thresh);
[min(abs(temp(temp<0))) max(abs(temp(temp<0))) min(abs(temp(temp>0))) max(abs(temp(temp>0)))]

temp = XYY.STAT.MOR.TSA.HIPP_SEX_B .* (XYY.STAT.MOR.TSA.HIPP_SEX_Q < thresh);
[min(abs(temp(temp<0))) max(abs(temp(temp<0))) min(abs(temp(temp>0))) max(abs(temp(temp>0)))]

temp = XYY.STAT.MOR.TSA.AMYG_SEX_B .* (XYY.STAT.MOR.TSA.AMYG_SEX_Q < thresh);
[min(abs(temp(temp<0))) max(abs(temp(temp<0))) min(abs(temp(temp>0))) max(abs(temp(temp>0)))]

% Calculate the colorscale of thresholded maps IGD
temp = IGD.STAT.MOR.TSA.STR_GRP_B .* (IGD.STAT.MOR.TSA.STR_GRP_P < thresh);
[min(abs(temp(temp<0))) max(abs(temp(temp<0))) min(abs(temp(temp>0))) max(abs(temp(temp>0)))]

temp = IGD.STAT.MOR.TSA.GP_GRP_B .* (IGD.STAT.MOR.TSA.GP_GRP_P < thresh);
[min(abs(temp(temp<0))) max(abs(temp(temp<0))) min(abs(temp(temp>0))) max(abs(temp(temp>0)))]

temp = IGD.STAT.MOR.TSA.THAL_GRP_B .* (IGD.STAT.MOR.TSA.THAL_GRP_P < thresh);
[min(abs(temp(temp<0))) max(abs(temp(temp<0))) min(abs(temp(temp>0))) max(abs(temp(temp>0)))]

temp = IGD.STAT.MOR.TSA.HIPP_GRP_B .* (IGD.STAT.MOR.TSA.HIPP_GRP_P < thresh);
[min(abs(temp(temp<0))) max(abs(temp(temp<0))) min(abs(temp(temp>0))) max(abs(temp(temp>0)))]

temp = IGD.STAT.MOR.TSA.AMYG_GRP_B .* (IGD.STAT.MOR.TSA.AMYG_GRP_P < thresh);
[min(abs(temp(temp<0))) max(abs(temp(temp<0))) min(abs(temp(temp>0))) max(abs(temp(temp>0)))]

%% Rotating full view

% XXY
SUB_B = [XXY.STAT.MOR.TSA.STR_SEX_B XXY.STAT.MOR.TSA.GP_SEX_B XXY.STAT.MOR.TSA.THAL_SEX_B XXY.STAT.MOR.TSA.HIPP_SEX_B XXY.STAT.MOR.TSA.AMYG_SEX_B];
SUB_Q = [XXY.STAT.MOR.TSA.STR_SEX_Q XXY.STAT.MOR.TSA.GP_SEX_Q XXY.STAT.MOR.TSA.THAL_SEX_Q XXY.STAT.MOR.TSA.HIPP_SEX_Q XXY.STAT.MOR.TSA.AMYG_SEX_Q];

% % XXY Video
% c = 1;
% dA = 0.2;
% for i = 0:dA:360-dA
%     disp(['Rotation angle ' num2str(i) ' deg']);
%     h = figure('WindowState', 'Maximized');
%     HyoMagetViewSpin(SUB_B .* (SUB_Q < thresh), [-0.8 0.8], i, cmap.PuOr);
%     print(h, fullfile(figDir, 'revolving_xxy', ['XXY_TSA_ANG-' num2str(c) '.png']), '-dpng', '-r300');
%     close(h);
%     c = c + 1;
% end

% XXY Static
c = 1;
dA = 45;
for i =  0:dA:360-dA
    disp(['Rotation angle ' num2str(i) ' deg']);
    h = figure('WindowState', 'Maximized');
    HyoMagetViewSpin(SUB_B .* (SUB_Q < thresh), [-0.8 0.8], i, ThreshColorMap(SUB_B .* (SUB_Q < thresh), cmap.PuOr, [-0.8 0.8], 0.915), 'static');
    print(h, fullfile(figDir, 'static_xxy', ['XXY_TSA_ANG-' num2str(c) '.png']), '-dpng', '-r300');
    close(h);
    c = c + 1;
end

% Static (not thresholded)
c = 1;
dA = 45;
for i =  0:dA:360-dA
    disp(['Rotation angle ' num2str(i) ' deg']);
    h = figure('WindowState', 'Maximized');
    HyoMagetViewSpin(SUB_B, [-0.8 0.8], i, cmap.PuOr, 'static');
    print(h, fullfile(figDir, 'static_xxy', ['NOTHR_XXY_TSA_ANG-' num2str(c) '.png']), '-dpng', '-r300');
    close(h);
    c = c + 1;
end

% XYY
SUB_B = [XYY.STAT.MOR.TSA.STR_SEX_B XYY.STAT.MOR.TSA.GP_SEX_B XYY.STAT.MOR.TSA.THAL_SEX_B XYY.STAT.MOR.TSA.HIPP_SEX_B XYY.STAT.MOR.TSA.AMYG_SEX_B];
SUB_Q = [XYY.STAT.MOR.TSA.STR_SEX_Q XYY.STAT.MOR.TSA.GP_SEX_Q XYY.STAT.MOR.TSA.THAL_SEX_Q XYY.STAT.MOR.TSA.HIPP_SEX_Q XYY.STAT.MOR.TSA.AMYG_SEX_Q];

% % XYY Video
% c = 1;
% dA = 0.2;
% for i = 0:dA:360-dA
%     disp(['Rotation angle ' num2str(i) ' deg']);
%     h = figure('WindowState', 'Maximized');
%     HyoMagetViewSpin(SUB_B .* (SUB_Q < thresh), [-0.8 0.8], i, cmap.PuOr);
%     print(h, fullfile(figDir, 'revolving_xyy', ['XYY_TSA_ANG-' num2str(c) '.png']), '-dpng', '-r300');
%     close(h);
%     c = c + 1;
% end

% XYY Static
c = 1;
dA = 45;
for i =  0:dA:360-dA
    disp(['Rotation angle ' num2str(i) ' deg']);
    h = figure('WindowState', 'Maximized');
    HyoMagetViewSpin(SUB_B .* (SUB_Q < thresh), [-0.8 0.8], i, ThreshColorMap(SUB_B .* (SUB_Q < thresh), cmap.PuOr, [-0.8 0.8], 0.915), 'static');
    print(h, fullfile(figDir, 'static_xyy', ['XYY_TSA_ANG-' num2str(c) '.png']), '-dpng', '-r300');
    close(h);
    c = c + 1;
end

% Static (not thresholded)
c = 1;
dA = 45;
for i =  0:dA:360-dA
    disp(['Rotation angle ' num2str(i) ' deg']);
    h = figure('WindowState', 'Maximized');
    HyoMagetViewSpin(SUB_B, [-0.8 0.8], i, cmap.PuOr, 'static');
    print(h, fullfile(figDir, 'static_xyy', ['NOTHR_XYY_TSA_ANG-' num2str(c) '.png']), '-dpng', '-r300');
    close(h);
    c = c + 1;
end

% IGD
SUB_B = [IGD.STAT.MOR.TSA.STR_GRP_B IGD.STAT.MOR.TSA.GP_GRP_B IGD.STAT.MOR.TSA.THAL_GRP_B IGD.STAT.MOR.TSA.HIPP_GRP_B IGD.STAT.MOR.TSA.AMYG_GRP_B];
SUB_Q = [IGD.STAT.MOR.TSA.STR_GRP_P IGD.STAT.MOR.TSA.GP_GRP_P IGD.STAT.MOR.TSA.THAL_GRP_P IGD.STAT.MOR.TSA.HIPP_GRP_P IGD.STAT.MOR.TSA.AMYG_GRP_P];

% % IGD Video
% c = 1;
% dA = 0.2;
% for i = 0:dA:360-dA
%     disp(['Rotation angle ' num2str(i) ' deg']);
%     h = figure('WindowState', 'Maximized');
%     HyoMagetViewSpin(SUB_B .* (SUB_Q < thresh), [-0.8 0.8], i, cmap.PuOr);
%     print(h, fullfile(figDir, 'revolving_igd', ['IGD_TSA_ANG-' num2str(c) '.png']), '-dpng', '-r300');
%     close(h);
%     c = c + 1;
% end

% IGD Static
c = 1;
dA = 45;
for i =  0:dA:360-dA
    disp(['Rotation angle ' num2str(i) ' deg']);
    h = figure('WindowState', 'Maximized');
    HyoMagetViewSpin(SUB_B .* (SUB_Q < thresh), [-0.8 0.8], i, ThreshColorMap(SUB_B .* (SUB_Q < thresh), cmap.PuOr, [-0.8 0.8], 0.915), 'static');
    print(h, fullfile(figDir, 'static_igd', ['IGD_TSA_ANG-' num2str(c) '.png']), '-dpng', '-r300');
    close(h);
    c = c + 1;
end

% Static (not thresholded)
c = 1;
dA = 45;
for i =  0:dA:360-dA
    disp(['Rotation angle ' num2str(i) ' deg']);
    h = figure('WindowState', 'Maximized');
    HyoMagetViewSpin(SUB_B, [-0.8 0.8], i, cmap.PuOr, 'static');
    print(h, fullfile(figDir, 'static_igd', ['NOTHR_IGD_TSA_ANG-' num2str(c) '.png']), '-dpng', '-r300');
    close(h);
    c = c + 1;
end

%% Regional congruence between HCP and XXY

% Combine significance and direction into a single map
XXY_STR_AD = zeros(1,26000);
HCP_STR_AD = zeros(1,26000);
STR_AD = zeros(1,26000);

XXY_STR_AD(XXY.STAT.MOR.TSA.STR_SEX_B > 0) =  1;
XXY_STR_AD(XXY.STAT.MOR.TSA.STR_SEX_B < 0) = -1;
XXY_STR_AD((XXY.STAT.MOR.TSA.STR_SEX_Q < thresh) & (XXY.STAT.MOR.TSA.STR_SEX_B > 0))  =  2;
XXY_STR_AD((XXY.STAT.MOR.TSA.STR_SEX_Q < thresh) & (XXY.STAT.MOR.TSA.STR_SEX_B < 0))  = -2;

HCP_STR_AD(HCP.STAT.MOR.TSA.STR_SEX_B > 0) =  1;
HCP_STR_AD(HCP.STAT.MOR.TSA.STR_SEX_B < 0) = -1;
HCP_STR_AD((HCP.STAT.MOR.TSA.STR_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.STR_SEX_B > 0))  =  2;
HCP_STR_AD((HCP.STAT.MOR.TSA.STR_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.STR_SEX_B < 0))  = -2;

% Generate agreement maps
STR_AD((HCP_STR_AD < 2) & (HCP_STR_AD > -2)) = 0;

STR_AD((HCP_STR_AD == 2) & (XXY_STR_AD == 2)) = -1.9;
STR_AD((HCP_STR_AD == 2) & (XXY_STR_AD == 1)) = -1.6;

STR_AD((HCP_STR_AD == -2) & (XXY_STR_AD == 2)) = -0.9;
STR_AD((HCP_STR_AD == -2) & (XXY_STR_AD == 1)) = -0.6;

STR_AD((HCP_STR_AD == 2) & (XXY_STR_AD == -2)) = 0.3;
STR_AD((HCP_STR_AD == 2) & (XXY_STR_AD == -1)) = 0.6;

STR_AD((HCP_STR_AD == -2) & (XXY_STR_AD == -2)) = 1.3;
STR_AD((HCP_STR_AD == -2) & (XXY_STR_AD == -1)) = 1.6;

% Visualize
h = figure('WindowState', 'Maximized'); HyoMagetViewFlat(STR_AD, [-2 2], cmap.quad_div);
exportgraphics(h, fullfile(figDir, 'xxy_conjunction_STR.png'), 'Resolution', 300);

% Bar graphs
STR_AD_Count = [sum(S_STR_NA(STR_AD == 1.3)) sum(S_STR_NA(STR_AD == 1.6)) 0 ...
               sum(S_STR_NA(STR_AD == -1.9)) sum(S_STR_NA(STR_AD == -1.6)) 0 ...
               sum(S_STR_NA(STR_AD == -0.9)) sum(S_STR_NA(STR_AD == -0.6)) 0 ...
               sum(S_STR_NA(STR_AD == 0.3)) sum(S_STR_NA(STR_AD == 0.6)) 0];

h = figure;
b = bar(1, STR_AD_Count, 'stacked', 'EdgeColor', 'none');
for i = 1:length(STR_AD_Count)
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_quad_div(i,:);
end
box off
axis off
exportgraphics(h, fullfile(figDir, 'xxy_conjunction_STR_bar.png'), 'Resolution', 300);

% Combine significance and direction into a single map
XXY_GP_AD = zeros(1,6000);
HCP_GP_AD = zeros(1,6000);
GP_AD = zeros(1,6000);

XXY_GP_AD(XXY.STAT.MOR.TSA.GP_SEX_B > 0) =  1;
XXY_GP_AD(XXY.STAT.MOR.TSA.GP_SEX_B < 0) = -1;
XXY_GP_AD((XXY.STAT.MOR.TSA.GP_SEX_Q < thresh) & (XXY.STAT.MOR.TSA.GP_SEX_B > 0))  =  2;
XXY_GP_AD((XXY.STAT.MOR.TSA.GP_SEX_Q < thresh) & (XXY.STAT.MOR.TSA.GP_SEX_B < 0))  = -2;

HCP_GP_AD(HCP.STAT.MOR.TSA.GP_SEX_B > 0) =  1;
HCP_GP_AD(HCP.STAT.MOR.TSA.GP_SEX_B < 0) = -1;
HCP_GP_AD((HCP.STAT.MOR.TSA.GP_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.GP_SEX_B > 0))  =  2;
HCP_GP_AD((HCP.STAT.MOR.TSA.GP_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.GP_SEX_B < 0))  = -2;

% Generate agreement maps
GP_AD((HCP_GP_AD < 2) & (HCP_GP_AD > -2)) = 0;

GP_AD((HCP_GP_AD == 2) & (XXY_GP_AD == 2)) = -1.9;
GP_AD((HCP_GP_AD == 2) & (XXY_GP_AD == 1)) = -1.6;

GP_AD((HCP_GP_AD == -2) & (XXY_GP_AD == 2)) = -0.9;
GP_AD((HCP_GP_AD == -2) & (XXY_GP_AD == 1)) = -0.6;

GP_AD((HCP_GP_AD == 2) & (XXY_GP_AD == -2)) = 0.3;
GP_AD((HCP_GP_AD == 2) & (XXY_GP_AD == -1)) = 0.6;

GP_AD((HCP_GP_AD == -2) & (XXY_GP_AD == -2)) = 1.3;
GP_AD((HCP_GP_AD == -2) & (XXY_GP_AD == -1)) = 1.6;

% Visualize
h = figure('WindowState', 'Maximized'); HyoMagetViewFlat(GP_AD, [-2 2], cmap.quad_div);
exportgraphics(h, fullfile(figDir, 'xxy_conjunction_GP.png'), 'Resolution', 300);

% Bar graphs
GP_AD_Count = [sum(S_GP_NA(GP_AD == 1.3)) sum(S_GP_NA(GP_AD == 1.6)) 0 ...
               sum(S_GP_NA(GP_AD == -1.9)) sum(S_GP_NA(GP_AD == -1.6)) 0 ...
               sum(S_GP_NA(GP_AD == -0.9)) sum(S_GP_NA(GP_AD == -0.6)) 0 ...
               sum(S_GP_NA(GP_AD == 0.3)) sum(S_GP_NA(GP_AD == 0.6)) 0];

h = figure;
b = bar(1, GP_AD_Count, 'stacked', 'EdgeColor', 'none');
for i = 1:length(GP_AD_Count)
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_quad_div(i,:);
end
box off
axis off
exportgraphics(h, fullfile(figDir, 'xxy_conjunction_GP_bar.png'), 'Resolution', 300);

% Combine significance and direction into a single map
XXY_THAL_AD = zeros(1,13000);
HCP_THAL_AD = zeros(1,13000);
THAL_AD = zeros(1,13000);

XXY_THAL_AD(XXY.STAT.MOR.TSA.THAL_SEX_B > 0) =  1;
XXY_THAL_AD(XXY.STAT.MOR.TSA.THAL_SEX_B < 0) = -1;
XXY_THAL_AD((XXY.STAT.MOR.TSA.THAL_SEX_Q < thresh) & (XXY.STAT.MOR.TSA.THAL_SEX_B > 0))  =  2;
XXY_THAL_AD((XXY.STAT.MOR.TSA.THAL_SEX_Q < thresh) & (XXY.STAT.MOR.TSA.THAL_SEX_B < 0))  = -2;

HCP_THAL_AD(HCP.STAT.MOR.TSA.THAL_SEX_B > 0) =  1;
HCP_THAL_AD(HCP.STAT.MOR.TSA.THAL_SEX_B < 0) = -1;
HCP_THAL_AD((HCP.STAT.MOR.TSA.THAL_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.THAL_SEX_B > 0))  =  2;
HCP_THAL_AD((HCP.STAT.MOR.TSA.THAL_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.THAL_SEX_B < 0))  = -2;

% Generate agreement maps
THAL_AD((HCP_THAL_AD < 2) & (HCP_THAL_AD > -2)) = 0;

THAL_AD((HCP_THAL_AD == 2) & (XXY_THAL_AD == 2)) = -1.9;
THAL_AD((HCP_THAL_AD == 2) & (XXY_THAL_AD == 1)) = -1.6;

THAL_AD((HCP_THAL_AD == -2) & (XXY_THAL_AD == 2)) = -0.9;
THAL_AD((HCP_THAL_AD == -2) & (XXY_THAL_AD == 1)) = -0.6;

THAL_AD((HCP_THAL_AD == 2) & (XXY_THAL_AD == -2)) = 0.3;
THAL_AD((HCP_THAL_AD == 2) & (XXY_THAL_AD == -1)) = 0.6;

THAL_AD((HCP_THAL_AD == -2) & (XXY_THAL_AD == -2)) = 1.3;
THAL_AD((HCP_THAL_AD == -2) & (XXY_THAL_AD == -1)) = 1.6;

% Visualize
h = figure('WindowState', 'Maximized'); HyoMagetViewFlat(THAL_AD, [-2 2], cmap.quad_div);
exportgraphics(h, fullfile(figDir, 'xxy_conjunction_THAL.png'), 'Resolution', 300);

% Bar graphs
THAL_AD_Count = [sum(S_THAL_NA(THAL_AD == 1.3)) sum(S_THAL_NA(THAL_AD == 1.6)) 0 ...
               sum(S_THAL_NA(THAL_AD == -1.9)) sum(S_THAL_NA(THAL_AD == -1.6)) 0 ...
               sum(S_THAL_NA(THAL_AD == -0.9)) sum(S_THAL_NA(THAL_AD == -0.6)) 0 ...
               sum(S_THAL_NA(THAL_AD == 0.3)) sum(S_THAL_NA(THAL_AD == 0.6)) 0];

h = figure;
b = bar(1, THAL_AD_Count, 'stacked', 'EdgeColor', 'none');
for i = 1:length(THAL_AD_Count)
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_quad_div(i,:);
end
box off
axis off
exportgraphics(h, fullfile(figDir, 'xxy_conjunction_THAL_bar.png'), 'Resolution', 300);

% Combine significance and direction into a single map
XXY_HIPP_AD = zeros(1,2367);
HCP_HIPP_AD = zeros(1,2367);
HIPP_AD = zeros(1,2367);

XXY_HIPP_AD(XXY.STAT.MOR.TSA.HIPP_SEX_B > 0) =  1;
XXY_HIPP_AD(XXY.STAT.MOR.TSA.HIPP_SEX_B < 0) = -1;
XXY_HIPP_AD((XXY.STAT.MOR.TSA.HIPP_SEX_Q < thresh) & (XXY.STAT.MOR.TSA.HIPP_SEX_B > 0))  =  2;
XXY_HIPP_AD((XXY.STAT.MOR.TSA.HIPP_SEX_Q < thresh) & (XXY.STAT.MOR.TSA.HIPP_SEX_B < 0))  = -2;

HCP_HIPP_AD(HCP.STAT.MOR.TSA.HIPP_SEX_B > 0) =  1;
HCP_HIPP_AD(HCP.STAT.MOR.TSA.HIPP_SEX_B < 0) = -1;
HCP_HIPP_AD((HCP.STAT.MOR.TSA.HIPP_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.HIPP_SEX_B > 0))  =  2;
HCP_HIPP_AD((HCP.STAT.MOR.TSA.HIPP_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.HIPP_SEX_B < 0))  = -2;

% Generate agreement maps
HIPP_AD((HCP_HIPP_AD < 2) & (HCP_HIPP_AD > -2)) = 0;

HIPP_AD((HCP_HIPP_AD == 2) & (XXY_HIPP_AD == 2)) = -1.9;
HIPP_AD((HCP_HIPP_AD == 2) & (XXY_HIPP_AD == 1)) = -1.6;

HIPP_AD((HCP_HIPP_AD == -2) & (XXY_HIPP_AD == 2)) = -0.9;
HIPP_AD((HCP_HIPP_AD == -2) & (XXY_HIPP_AD == 1)) = -0.6;

HIPP_AD((HCP_HIPP_AD == 2) & (XXY_HIPP_AD == -2)) = 0.3;
HIPP_AD((HCP_HIPP_AD == 2) & (XXY_HIPP_AD == -1)) = 0.6;

HIPP_AD((HCP_HIPP_AD == -2) & (XXY_HIPP_AD == -2)) = 1.3;
HIPP_AD((HCP_HIPP_AD == -2) & (XXY_HIPP_AD == -1)) = 1.6;

% Visualize
h = figure('WindowState', 'Maximized'); HyoMagetViewFlat(HIPP_AD, [-2 2], cmap.quad_div);
exportgraphics(h, fullfile(figDir, 'xxy_conjunction_HIPP.png'), 'Resolution', 300);

% Bar graphs
HIPP_AD_Count = [sum(S_HIPP_NA(HIPP_AD == 1.3)) sum(S_HIPP_NA(HIPP_AD == 1.6)) 0 ...
               sum(S_HIPP_NA(HIPP_AD == -1.9)) sum(S_HIPP_NA(HIPP_AD == -1.6)) 0 ...
               sum(S_HIPP_NA(HIPP_AD == -0.9)) sum(S_HIPP_NA(HIPP_AD == -0.6)) 0 ...
               sum(S_HIPP_NA(HIPP_AD == 0.3)) sum(S_HIPP_NA(HIPP_AD == 0.6)) 0];

h = figure;
b = bar(1, HIPP_AD_Count, 'stacked', 'EdgeColor', 'none');
for i = 1:length(HIPP_AD_Count)
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_quad_div(i,:);
end
box off
axis off
exportgraphics(h, fullfile(figDir, 'xxy_conjunction_HIPP_bar.png'), 'Resolution', 300);

% Combine significance and direction into a single map
XXY_AMYG_AD = zeros(1,2878);
HCP_AMYG_AD = zeros(1,2878);
AMYG_AD = zeros(1,2878);

XXY_AMYG_AD(XXY.STAT.MOR.TSA.AMYG_SEX_B > 0) =  1;
XXY_AMYG_AD(XXY.STAT.MOR.TSA.AMYG_SEX_B < 0) = -1;
XXY_AMYG_AD((XXY.STAT.MOR.TSA.AMYG_SEX_Q < thresh) & (XXY.STAT.MOR.TSA.AMYG_SEX_B > 0))  =  2;
XXY_AMYG_AD((XXY.STAT.MOR.TSA.AMYG_SEX_Q < thresh) & (XXY.STAT.MOR.TSA.AMYG_SEX_B < 0))  = -2;

HCP_AMYG_AD(HCP.STAT.MOR.TSA.AMYG_SEX_B > 0) =  1;
HCP_AMYG_AD(HCP.STAT.MOR.TSA.AMYG_SEX_B < 0) = -1;
HCP_AMYG_AD((HCP.STAT.MOR.TSA.AMYG_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.AMYG_SEX_B > 0))  =  2;
HCP_AMYG_AD((HCP.STAT.MOR.TSA.AMYG_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.AMYG_SEX_B < 0))  = -2;

% Generate agreement maps
AMYG_AD((HCP_AMYG_AD < 2) & (HCP_AMYG_AD > -2)) = 0;

AMYG_AD((HCP_AMYG_AD == 2) & (XXY_AMYG_AD == 2)) = -1.9;
AMYG_AD((HCP_AMYG_AD == 2) & (XXY_AMYG_AD == 1)) = -1.6;

AMYG_AD((HCP_AMYG_AD == -2) & (XXY_AMYG_AD == 2)) = -0.9;
AMYG_AD((HCP_AMYG_AD == -2) & (XXY_AMYG_AD == 1)) = -0.6;

AMYG_AD((HCP_AMYG_AD == 2) & (XXY_AMYG_AD == -2)) = 0.3;
AMYG_AD((HCP_AMYG_AD == 2) & (XXY_AMYG_AD == -1)) = 0.6;

AMYG_AD((HCP_AMYG_AD == -2) & (XXY_AMYG_AD == -2)) = 1.3;
AMYG_AD((HCP_AMYG_AD == -2) & (XXY_AMYG_AD == -1)) = 1.6;

% Visualize
h = figure('WindowState', 'Maximized'); HyoMagetViewFlat(AMYG_AD, [-2 2], cmap.quad_div);
exportgraphics(h, fullfile(figDir, 'xxy_conjunction_AMYG.png'), 'Resolution', 300);

% Bar graphs
AMYG_AD_Count = [sum(S_AMYG_NA(AMYG_AD == 1.3)) sum(S_AMYG_NA(AMYG_AD == 1.6)) 0 ...
               sum(S_AMYG_NA(AMYG_AD == -1.9)) sum(S_AMYG_NA(AMYG_AD == -1.6)) 0 ...
               sum(S_AMYG_NA(AMYG_AD == -0.9)) sum(S_AMYG_NA(AMYG_AD == -0.6)) 0 ...
               sum(S_AMYG_NA(AMYG_AD == 0.3)) sum(S_AMYG_NA(AMYG_AD == 0.6)) 0];

h = figure;
b = bar(1, AMYG_AD_Count, 'stacked', 'EdgeColor', 'none');
for i = 1:length(AMYG_AD_Count)
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_quad_div(i,:);
end
box off
axis off
exportgraphics(h, fullfile(figDir, 'xxy_conjunction_AMYG_bar.png'), 'Resolution', 300);

close all;

sum(STR_AD_Count(1:6)) ./ sum(STR_AD_Count)
sum(GP_AD_Count(1:6)) ./ sum(GP_AD_Count)
sum(THAL_AD_Count(1:6)) ./ sum(THAL_AD_Count)
sum(HIPP_AD_Count(1:6)) ./ sum(HIPP_AD_Count)
sum(AMYG_AD_Count(1:6)) ./ sum(AMYG_AD_Count)

%% Spatial permutation testing for HCP and XXY

% Categorize HCP sex difference maps
PERM_HCP_STR_AD = zeros(26000, nPerm);
PERM_HCP_STR_AD(HCP.STAT.MOR.MSR.STR_SEX_B > 0) =  1;
PERM_HCP_STR_AD(HCP.STAT.MOR.MSR.STR_SEX_B < 0) = -1;
PERM_HCP_STR_AD(HCP.STAT.MOR.MSR.STR_SEX_Q & (HCP.STAT.MOR.MSR.STR_SEX_B > 0))  =  2;
PERM_HCP_STR_AD(HCP.STAT.MOR.MSR.STR_SEX_Q & (HCP.STAT.MOR.MSR.STR_SEX_B < 0))  = -2;

PERM_HCP_GP_AD = zeros(6000, nPerm);
PERM_HCP_GP_AD(HCP.STAT.MOR.MSR.GP_SEX_B > 0) =  1;
PERM_HCP_GP_AD(HCP.STAT.MOR.MSR.GP_SEX_B < 0) = -1;
PERM_HCP_GP_AD(HCP.STAT.MOR.MSR.GP_SEX_Q & (HCP.STAT.MOR.MSR.GP_SEX_B > 0))  =  2;
PERM_HCP_GP_AD(HCP.STAT.MOR.MSR.GP_SEX_Q & (HCP.STAT.MOR.MSR.GP_SEX_B < 0))  = -2;

PERM_HCP_THAL_AD = zeros(13000, nPerm);
PERM_HCP_THAL_AD(HCP.STAT.MOR.MSR.THAL_SEX_B > 0) =  1;
PERM_HCP_THAL_AD(HCP.STAT.MOR.MSR.THAL_SEX_B < 0) = -1;
PERM_HCP_THAL_AD(HCP.STAT.MOR.MSR.THAL_SEX_Q & (HCP.STAT.MOR.MSR.THAL_SEX_B > 0))  =  2;
PERM_HCP_THAL_AD(HCP.STAT.MOR.MSR.THAL_SEX_Q & (HCP.STAT.MOR.MSR.THAL_SEX_B < 0))  = -2;

PERM_HCP_HIPP_AD = zeros(2367, nPerm);
PERM_HCP_HIPP_AD(HCP.STAT.MOR.MSR.HIPP_SEX_B > 0) =  1;
PERM_HCP_HIPP_AD(HCP.STAT.MOR.MSR.HIPP_SEX_B < 0) = -1;
PERM_HCP_HIPP_AD(HCP.STAT.MOR.MSR.HIPP_SEX_Q & (HCP.STAT.MOR.MSR.HIPP_SEX_B > 0))  =  2;
PERM_HCP_HIPP_AD(HCP.STAT.MOR.MSR.HIPP_SEX_Q & (HCP.STAT.MOR.MSR.HIPP_SEX_B < 0))  = -2;

PERM_HCP_AMYG_AD = zeros(2878, nPerm);
PERM_HCP_AMYG_AD(HCP.STAT.MOR.MSR.AMYG_SEX_B > 0) =  1;
PERM_HCP_AMYG_AD(HCP.STAT.MOR.MSR.AMYG_SEX_B < 0) = -1;
PERM_HCP_AMYG_AD(HCP.STAT.MOR.MSR.AMYG_SEX_Q & (HCP.STAT.MOR.MSR.AMYG_SEX_B > 0))  =  2;
PERM_HCP_AMYG_AD(HCP.STAT.MOR.MSR.AMYG_SEX_Q & (HCP.STAT.MOR.MSR.AMYG_SEX_B < 0))  =  -2;

% Compute agreement maps
PERM_STR_AD = zeros(26000, nPerm);
PERM_STR_AD((PERM_HCP_STR_AD == 2) & repmat(XXY_STR_AD == 2, nPerm, 1)') = -1.9;
PERM_STR_AD((PERM_HCP_STR_AD == 2) & repmat(XXY_STR_AD == 1, nPerm, 1)') = -1.6;
PERM_STR_AD((PERM_HCP_STR_AD == -2) & repmat(XXY_STR_AD == 2, nPerm, 1)') = -0.9;
PERM_STR_AD((PERM_HCP_STR_AD == -2) & repmat(XXY_STR_AD == 1, nPerm, 1)') = -0.6;
PERM_STR_AD((PERM_HCP_STR_AD == 2) & repmat(XXY_STR_AD == -2, nPerm, 1)') = 0.3;
PERM_STR_AD((PERM_HCP_STR_AD == 2) & repmat(XXY_STR_AD == -1, nPerm, 1)') = 0.6;
PERM_STR_AD((PERM_HCP_STR_AD == -2) & repmat(XXY_STR_AD == -2, nPerm, 1)') = 1.3;
PERM_STR_AD((PERM_HCP_STR_AD == -2) & repmat(XXY_STR_AD == -1, nPerm, 1)') = 1.6;

PERM_GP_AD = zeros(6000, nPerm);
PERM_GP_AD((PERM_HCP_GP_AD == 2) & repmat(XXY_GP_AD == 2, nPerm, 1)') = -1.9;
PERM_GP_AD((PERM_HCP_GP_AD == 2) & repmat(XXY_GP_AD == 1, nPerm, 1)') = -1.6;
PERM_GP_AD((PERM_HCP_GP_AD == -2) & repmat(XXY_GP_AD == 2, nPerm, 1)') = -0.9;
PERM_GP_AD((PERM_HCP_GP_AD == -2) & repmat(XXY_GP_AD == 1, nPerm, 1)') = -0.6;
PERM_GP_AD((PERM_HCP_GP_AD == 2) & repmat(XXY_GP_AD == -2, nPerm, 1)') = 0.3;
PERM_GP_AD((PERM_HCP_GP_AD == 2) & repmat(XXY_GP_AD == -1, nPerm, 1)') = 0.6;
PERM_GP_AD((PERM_HCP_GP_AD == -2) & repmat(XXY_GP_AD == -2, nPerm, 1)') = 1.3;
PERM_GP_AD((PERM_HCP_GP_AD == -2) & repmat(XXY_GP_AD == -1, nPerm, 1)') = 1.6;

PERM_THAL_AD = zeros(13000, nPerm);
PERM_THAL_AD((PERM_HCP_THAL_AD == 2) & repmat(XXY_THAL_AD == 2, nPerm, 1)') = -1.9;
PERM_THAL_AD((PERM_HCP_THAL_AD == 2) & repmat(XXY_THAL_AD == 1, nPerm, 1)') = -1.6;
PERM_THAL_AD((PERM_HCP_THAL_AD == -2) & repmat(XXY_THAL_AD == 2, nPerm, 1)') = -0.9;
PERM_THAL_AD((PERM_HCP_THAL_AD == -2) & repmat(XXY_THAL_AD == 1, nPerm, 1)') = -0.6;
PERM_THAL_AD((PERM_HCP_THAL_AD == 2) & repmat(XXY_THAL_AD == -2, nPerm, 1)') = 0.3;
PERM_THAL_AD((PERM_HCP_THAL_AD == 2) & repmat(XXY_THAL_AD == -1, nPerm, 1)') = 0.6;
PERM_THAL_AD((PERM_HCP_THAL_AD == -2) & repmat(XXY_THAL_AD == -2, nPerm, 1)') = 1.3;
PERM_THAL_AD((PERM_HCP_THAL_AD == -2) & repmat(XXY_THAL_AD == -1, nPerm, 1)') = 1.6;

PERM_HIPP_AD = zeros(2367, nPerm);
PERM_HIPP_AD((PERM_HCP_HIPP_AD == 2) & repmat(XXY_HIPP_AD == 2, nPerm, 1)') = -1.9;
PERM_HIPP_AD((PERM_HCP_HIPP_AD == 2) & repmat(XXY_HIPP_AD == 1, nPerm, 1)') = -1.6;
PERM_HIPP_AD((PERM_HCP_HIPP_AD == -2) & repmat(XXY_HIPP_AD == 2, nPerm, 1)') = -0.9;
PERM_HIPP_AD((PERM_HCP_HIPP_AD == -2) & repmat(XXY_HIPP_AD == 1, nPerm, 1)') = -0.6;
PERM_HIPP_AD((PERM_HCP_HIPP_AD == 2) & repmat(XXY_HIPP_AD == -2, nPerm, 1)') = 0.3;
PERM_HIPP_AD((PERM_HCP_HIPP_AD == 2) & repmat(XXY_HIPP_AD == -1, nPerm, 1)') = 0.6;
PERM_HIPP_AD((PERM_HCP_HIPP_AD == -2) & repmat(XXY_HIPP_AD == -2, nPerm, 1)') = 1.3;
PERM_HIPP_AD((PERM_HCP_HIPP_AD == -2) & repmat(XXY_HIPP_AD == -1, nPerm, 1)') = 1.6;

PERM_AMYG_AD = zeros(2878, nPerm);
PERM_AMYG_AD((PERM_HCP_AMYG_AD == 2) & repmat(XXY_AMYG_AD == 2, nPerm, 1)') = -1.9;
PERM_AMYG_AD((PERM_HCP_AMYG_AD == 2) & repmat(XXY_AMYG_AD == 1, nPerm, 1)') = -1.6;
PERM_AMYG_AD((PERM_HCP_AMYG_AD == -2) & repmat(XXY_AMYG_AD == 2, nPerm, 1)') = -0.9;
PERM_AMYG_AD((PERM_HCP_AMYG_AD == -2) & repmat(XXY_AMYG_AD == 1, nPerm, 1)') = -0.6;
PERM_AMYG_AD((PERM_HCP_AMYG_AD == 2) & repmat(XXY_AMYG_AD == -2, nPerm, 1)') = 0.3;
PERM_AMYG_AD((PERM_HCP_AMYG_AD == 2) & repmat(XXY_AMYG_AD == -1, nPerm, 1)') = 0.6;
PERM_AMYG_AD((PERM_HCP_AMYG_AD == -2) & repmat(XXY_AMYG_AD == -2, nPerm, 1)') = 1.3;
PERM_AMYG_AD((PERM_HCP_AMYG_AD == -2) & repmat(XXY_AMYG_AD == -1, nPerm, 1)') = 1.6;

% Calculate proportions
PERM_STR = sum(repmat(S_STR_NA, 1, nPerm) .* (abs(PERM_STR_AD) > 1));
PERM_GP = sum(repmat(S_GP_NA, 1, nPerm) .* (abs(PERM_GP_AD) > 1));
PERM_THAL = sum(repmat(S_THAL_NA, 1, nPerm) .* (abs(PERM_THAL_AD) > 1));
PERM_HIPP = sum(repmat(S_HIPP_NA, 1, nPerm) .* (abs(PERM_HIPP_AD) > 1));
PERM_AMYG = sum(repmat(S_AMYG_NA, 1, nPerm) .* (abs(PERM_AMYG_AD) > 1));

% Calculate p value
STR_PERM_P = sum(sum(S_STR_NA(abs(STR_AD) > 1)) < PERM_STR) ./ nPerm;
GP_PERM_P = sum(sum(S_GP_NA(abs(GP_AD) > 1)) < PERM_GP) ./ nPerm;
THAL_PERM_P = sum(sum(S_THAL_NA(abs(THAL_AD) > 1)) < PERM_THAL) ./ nPerm;
HIPP_PERM_P = sum(sum(S_HIPP_NA(abs(HIPP_AD) > 1)) < PERM_HIPP) ./ nPerm;
AMYG_PERM_P = sum(sum(S_AMYG_NA(abs(AMYG_AD) > 1)) < PERM_AMYG) ./ nPerm;

disp('MSR permutations');
[STR_PERM_P GP_PERM_P THAL_PERM_P HIPP_PERM_P AMYG_PERM_P]

%% Regional congruence between HCP and XYY

% Combine significance and direction into a single map
XYY_STR_AD = zeros(1,26000);
HCP_STR_AD = zeros(1,26000);
STR_AD = zeros(1,26000);

XYY_STR_AD(XYY.STAT.MOR.TSA.STR_SEX_B > 0) =  1;
XYY_STR_AD(XYY.STAT.MOR.TSA.STR_SEX_B < 0) = -1;
XYY_STR_AD((XYY.STAT.MOR.TSA.STR_SEX_Q < thresh) & (XYY.STAT.MOR.TSA.STR_SEX_B > 0))  =  2;
XYY_STR_AD((XYY.STAT.MOR.TSA.STR_SEX_Q < thresh) & (XYY.STAT.MOR.TSA.STR_SEX_B < 0))  = -2;

HCP_STR_AD(HCP.STAT.MOR.TSA.STR_SEX_B > 0) =  1;
HCP_STR_AD(HCP.STAT.MOR.TSA.STR_SEX_B < 0) = -1;
HCP_STR_AD((HCP.STAT.MOR.TSA.STR_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.STR_SEX_B > 0))  =  2;
HCP_STR_AD((HCP.STAT.MOR.TSA.STR_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.STR_SEX_B < 0))  = -2;

% Generate agreement maps
STR_AD((HCP_STR_AD < 2) & (HCP_STR_AD > -2)) = 0;

STR_AD((HCP_STR_AD == 2) & (XYY_STR_AD == -2)) = -1.9;
STR_AD((HCP_STR_AD == 2) & (XYY_STR_AD == -1)) = -1.6;

STR_AD((HCP_STR_AD == -2) & (XYY_STR_AD == -2)) = -0.9;
STR_AD((HCP_STR_AD == -2) & (XYY_STR_AD == -1)) = -0.6;

STR_AD((HCP_STR_AD == 2) & (XYY_STR_AD == 2)) = 0.3;
STR_AD((HCP_STR_AD == 2) & (XYY_STR_AD == 1)) = 0.6;

STR_AD((HCP_STR_AD == -2) & (XYY_STR_AD == 2)) = 1.3;
STR_AD((HCP_STR_AD == -2) & (XYY_STR_AD == 1)) = 1.6;

% Visualize
h = figure('WindowState', 'Maximized'); HyoMagetViewFlat(STR_AD, [-2 2], cmap.quad_div);
exportgraphics(h, fullfile(figDir, 'xyy_conjunction_STR.png'), 'Resolution', 300);

% Bar graphs
STR_AD_Count = [sum(S_STR_NA(STR_AD == 1.3)) sum(S_STR_NA(STR_AD == 1.6)) 0 ...
               sum(S_STR_NA(STR_AD == -1.9)) sum(S_STR_NA(STR_AD == -1.6)) 0 ...
               sum(S_STR_NA(STR_AD == -0.9)) sum(S_STR_NA(STR_AD == -0.6)) 0 ...
               sum(S_STR_NA(STR_AD == 0.3)) sum(S_STR_NA(STR_AD == 0.6)) 0];

h = figure;
b = bar(1, STR_AD_Count, 'stacked', 'EdgeColor', 'none');
for i = 1:length(STR_AD_Count)
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_quad_div(i,:);
end
box off
axis off
exportgraphics(h, fullfile(figDir, 'xyy_conjunction_STR_bar.png'), 'Resolution', 300);

% Combine significance and direction into a single map
XYY_GP_AD = zeros(1,6000);
HCP_GP_AD = zeros(1,6000);
GP_AD = zeros(1,6000);

XYY_GP_AD(XYY.STAT.MOR.TSA.GP_SEX_B > 0) =  1;
XYY_GP_AD(XYY.STAT.MOR.TSA.GP_SEX_B < 0) = -1;
XYY_GP_AD((XYY.STAT.MOR.TSA.GP_SEX_Q < thresh) & (XYY.STAT.MOR.TSA.GP_SEX_B > 0))  =  2;
XYY_GP_AD((XYY.STAT.MOR.TSA.GP_SEX_Q < thresh) & (XYY.STAT.MOR.TSA.GP_SEX_B < 0))  = -2;

HCP_GP_AD(HCP.STAT.MOR.TSA.GP_SEX_B > 0) =  1;
HCP_GP_AD(HCP.STAT.MOR.TSA.GP_SEX_B < 0) = -1;
HCP_GP_AD((HCP.STAT.MOR.TSA.GP_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.GP_SEX_B > 0))  =  2;
HCP_GP_AD((HCP.STAT.MOR.TSA.GP_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.GP_SEX_B < 0))  = -2;

% Generate agreement maps
GP_AD((HCP_GP_AD < 2) & (HCP_GP_AD > -2)) = 0;

GP_AD((HCP_GP_AD == 2) & (XYY_GP_AD == -2)) = -1.9;
GP_AD((HCP_GP_AD == 2) & (XYY_GP_AD == -1)) = -1.6;

GP_AD((HCP_GP_AD == -2) & (XYY_GP_AD == -2)) = -0.9;
GP_AD((HCP_GP_AD == -2) & (XYY_GP_AD == -1)) = -0.6;

GP_AD((HCP_GP_AD == 2) & (XYY_GP_AD == 2)) = 0.3;
GP_AD((HCP_GP_AD == 2) & (XYY_GP_AD == 1)) = 0.6;

GP_AD((HCP_GP_AD == -2) & (XYY_GP_AD == 2)) = 1.3;
GP_AD((HCP_GP_AD == -2) & (XYY_GP_AD == 1)) = 1.6;

% Visualize
h = figure('WindowState', 'Maximized'); HyoMagetViewFlat(GP_AD, [-2 2], cmap.quad_div);
exportgraphics(h, fullfile(figDir, 'xyy_conjunction_GP.png'), 'Resolution', 300);

% Bar graphs
GP_AD_Count = [sum(S_GP_NA(GP_AD == 1.3)) sum(S_GP_NA(GP_AD == 1.6)) 0 ...
               sum(S_GP_NA(GP_AD == -1.9)) sum(S_GP_NA(GP_AD == -1.6)) 0 ...
               sum(S_GP_NA(GP_AD == -0.9)) sum(S_GP_NA(GP_AD == -0.6)) 0 ...
               sum(S_GP_NA(GP_AD == 0.3)) sum(S_GP_NA(GP_AD == 0.6)) 0];

h = figure;
b = bar(1, GP_AD_Count, 'stacked', 'EdgeColor', 'none');
for i = 1:length(GP_AD_Count)
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_quad_div(i,:);
end
box off
axis off
exportgraphics(h, fullfile(figDir, 'xyy_conjunction_GP_bar.png'), 'Resolution', 300);

% Combine significance and direction into a single map
XYY_THAL_AD = zeros(1,13000);
HCP_THAL_AD = zeros(1,13000);
THAL_AD = zeros(1,13000);

XYY_THAL_AD(XYY.STAT.MOR.TSA.THAL_SEX_B > 0) =  1;
XYY_THAL_AD(XYY.STAT.MOR.TSA.THAL_SEX_B < 0) = -1;
XYY_THAL_AD((XYY.STAT.MOR.TSA.THAL_SEX_Q < thresh) & (XYY.STAT.MOR.TSA.THAL_SEX_B > 0))  =  2;
XYY_THAL_AD((XYY.STAT.MOR.TSA.THAL_SEX_Q < thresh) & (XYY.STAT.MOR.TSA.THAL_SEX_B < 0))  = -2;

HCP_THAL_AD(HCP.STAT.MOR.TSA.THAL_SEX_B > 0) =  1;
HCP_THAL_AD(HCP.STAT.MOR.TSA.THAL_SEX_B < 0) = -1;
HCP_THAL_AD((HCP.STAT.MOR.TSA.THAL_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.THAL_SEX_B > 0))  =  2;
HCP_THAL_AD((HCP.STAT.MOR.TSA.THAL_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.THAL_SEX_B < 0))  = -2;

% Generate agreement maps
THAL_AD((HCP_THAL_AD < 2) & (HCP_THAL_AD > -2)) = 0;

THAL_AD((HCP_THAL_AD == 2) & (XYY_THAL_AD == -2)) = -1.9;
THAL_AD((HCP_THAL_AD == 2) & (XYY_THAL_AD == -1)) = -1.6;

THAL_AD((HCP_THAL_AD == -2) & (XYY_THAL_AD == -2)) = -0.9;
THAL_AD((HCP_THAL_AD == -2) & (XYY_THAL_AD == -1)) = -0.6;

THAL_AD((HCP_THAL_AD == 2) & (XYY_THAL_AD == 2)) = 0.3;
THAL_AD((HCP_THAL_AD == 2) & (XYY_THAL_AD == 1)) = 0.6;

THAL_AD((HCP_THAL_AD == -2) & (XYY_THAL_AD == 2)) = 1.3;
THAL_AD((HCP_THAL_AD == -2) & (XYY_THAL_AD == 1)) = 1.6;

% Visualize
h = figure('WindowState', 'Maximized'); HyoMagetViewFlat(THAL_AD, [-2 2], cmap.quad_div);
exportgraphics(h, fullfile(figDir, 'xyy_conjunction_THAL.png'), 'Resolution', 300);

% Bar graphs
THAL_AD_Count = [sum(S_THAL_NA(THAL_AD == 1.3)) sum(S_THAL_NA(THAL_AD == 1.6)) 0 ...
               sum(S_THAL_NA(THAL_AD == -1.9)) sum(S_THAL_NA(THAL_AD == -1.6)) 0 ...
               sum(S_THAL_NA(THAL_AD == -0.9)) sum(S_THAL_NA(THAL_AD == -0.6)) 0 ...
               sum(S_THAL_NA(THAL_AD == 0.3)) sum(S_THAL_NA(THAL_AD == 0.6)) 0];

h = figure;
b = bar(1, THAL_AD_Count, 'stacked', 'EdgeColor', 'none');
for i = 1:length(THAL_AD_Count)
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_quad_div(i,:);
end
box off
axis off
exportgraphics(h, fullfile(figDir, 'xyy_conjunction_THAL_bar.png'), 'Resolution', 300);

% Combine significance and direction into a single map
XYY_HIPP_AD = zeros(1,2367);
HCP_HIPP_AD = zeros(1,2367);
HIPP_AD = zeros(1,2367);

XYY_HIPP_AD(XYY.STAT.MOR.TSA.HIPP_SEX_B > 0) =  1;
XYY_HIPP_AD(XYY.STAT.MOR.TSA.HIPP_SEX_B < 0) = -1;
XYY_HIPP_AD((XYY.STAT.MOR.TSA.HIPP_SEX_Q < thresh) & (XYY.STAT.MOR.TSA.HIPP_SEX_B > 0))  =  2;
XYY_HIPP_AD((XYY.STAT.MOR.TSA.HIPP_SEX_Q < thresh) & (XYY.STAT.MOR.TSA.HIPP_SEX_B < 0))  = -2;

HCP_HIPP_AD(HCP.STAT.MOR.TSA.HIPP_SEX_B > 0) =  1;
HCP_HIPP_AD(HCP.STAT.MOR.TSA.HIPP_SEX_B < 0) = -1;
HCP_HIPP_AD((HCP.STAT.MOR.TSA.HIPP_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.HIPP_SEX_B > 0))  =  2;
HCP_HIPP_AD((HCP.STAT.MOR.TSA.HIPP_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.HIPP_SEX_B < 0))  = -2;

% Generate agreement maps
HIPP_AD((HCP_HIPP_AD < 2) & (HCP_HIPP_AD > -2)) = 0;

HIPP_AD((HCP_HIPP_AD == 2) & (XYY_HIPP_AD == -2)) = -1.9;
HIPP_AD((HCP_HIPP_AD == 2) & (XYY_HIPP_AD == -1)) = -1.6;

HIPP_AD((HCP_HIPP_AD == -2) & (XYY_HIPP_AD == -2)) = -0.9;
HIPP_AD((HCP_HIPP_AD == -2) & (XYY_HIPP_AD == -1)) = -0.6;

HIPP_AD((HCP_HIPP_AD == 2) & (XYY_HIPP_AD == 2)) = 0.3;
HIPP_AD((HCP_HIPP_AD == 2) & (XYY_HIPP_AD == 1)) = 0.6;

HIPP_AD((HCP_HIPP_AD == -2) & (XYY_HIPP_AD == 2)) = 1.3;
HIPP_AD((HCP_HIPP_AD == -2) & (XYY_HIPP_AD == 1)) = 1.6;

% Visualize
h = figure('WindowState', 'Maximized'); HyoMagetViewFlat(HIPP_AD, [-2 2], cmap.quad_div);
exportgraphics(h, fullfile(figDir, 'xyy_conjunction_HIPP.png'), 'Resolution', 300);

% Bar graphs
HIPP_AD_Count = [sum(S_HIPP_NA(HIPP_AD == 1.3)) sum(S_HIPP_NA(HIPP_AD == 1.6)) 0 ...
               sum(S_HIPP_NA(HIPP_AD == -1.9)) sum(S_HIPP_NA(HIPP_AD == -1.6)) 0 ...
               sum(S_HIPP_NA(HIPP_AD == -0.9)) sum(S_HIPP_NA(HIPP_AD == -0.6)) 0 ...
               sum(S_HIPP_NA(HIPP_AD == 0.3)) sum(S_HIPP_NA(HIPP_AD == 0.6)) 0];

h = figure;
b = bar(1, HIPP_AD_Count, 'stacked', 'EdgeColor', 'none');
for i = 1:length(HIPP_AD_Count)
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_quad_div(i,:);
end
box off
axis off
exportgraphics(h, fullfile(figDir, 'xyy_conjunction_HIPP_bar.png'), 'Resolution', 300);

% Combine significance and direction into a single map
XYY_AMYG_AD = zeros(1,2878);
HCP_AMYG_AD = zeros(1,2878);
AMYG_AD = zeros(1,2878);

XYY_AMYG_AD(XYY.STAT.MOR.TSA.AMYG_SEX_B > 0) =  1;
XYY_AMYG_AD(XYY.STAT.MOR.TSA.AMYG_SEX_B < 0) = -1;
XYY_AMYG_AD((XYY.STAT.MOR.TSA.AMYG_SEX_Q < thresh) & (XYY.STAT.MOR.TSA.AMYG_SEX_B > 0))  =  2;
XYY_AMYG_AD((XYY.STAT.MOR.TSA.AMYG_SEX_Q < thresh) & (XYY.STAT.MOR.TSA.AMYG_SEX_B < 0))  = -2;

HCP_AMYG_AD(HCP.STAT.MOR.TSA.AMYG_SEX_B > 0) =  1;
HCP_AMYG_AD(HCP.STAT.MOR.TSA.AMYG_SEX_B < 0) = -1;
HCP_AMYG_AD((HCP.STAT.MOR.TSA.AMYG_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.AMYG_SEX_B > 0))  =  2;
HCP_AMYG_AD((HCP.STAT.MOR.TSA.AMYG_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.AMYG_SEX_B < 0))  = -2;

% Generate agreement maps
AMYG_AD((HCP_AMYG_AD < 2) & (HCP_AMYG_AD > -2)) = 0;

AMYG_AD((HCP_AMYG_AD == 2) & (XYY_AMYG_AD == -2)) = -1.9;
AMYG_AD((HCP_AMYG_AD == 2) & (XYY_AMYG_AD == -1)) = -1.6;

AMYG_AD((HCP_AMYG_AD == -2) & (XYY_AMYG_AD == -2)) = -0.9;
AMYG_AD((HCP_AMYG_AD == -2) & (XYY_AMYG_AD == -1)) = -0.6;

AMYG_AD((HCP_AMYG_AD == 2) & (XYY_AMYG_AD == 2)) = 0.3;
AMYG_AD((HCP_AMYG_AD == 2) & (XYY_AMYG_AD == 1)) = 0.6;

AMYG_AD((HCP_AMYG_AD == -2) & (XYY_AMYG_AD == 2)) = 1.3;
AMYG_AD((HCP_AMYG_AD == -2) & (XYY_AMYG_AD == 1)) = 1.6;

% Visualize
h = figure('WindowState', 'Maximized'); HyoMagetViewFlat(AMYG_AD, [-2 2], cmap.quad_div);
exportgraphics(h, fullfile(figDir, 'xyy_conjunction_AMYG.png'), 'Resolution', 300);

% Bar graphs
AMYG_AD_Count = [sum(S_AMYG_NA(AMYG_AD == 1.3)) sum(S_AMYG_NA(AMYG_AD == 1.6)) 0 ...
               sum(S_AMYG_NA(AMYG_AD == -1.9)) sum(S_AMYG_NA(AMYG_AD == -1.6)) 0 ...
               sum(S_AMYG_NA(AMYG_AD == -0.9)) sum(S_AMYG_NA(AMYG_AD == -0.6)) 0 ...
               sum(S_AMYG_NA(AMYG_AD == 0.3)) sum(S_AMYG_NA(AMYG_AD == 0.6)) 0];

h = figure;
b = bar(1, AMYG_AD_Count, 'stacked', 'EdgeColor', 'none');
for i = 1:length(AMYG_AD_Count)
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_quad_div(i,:);
end
box off
axis off
exportgraphics(h, fullfile(figDir, 'xyy_conjunction_AMYG_bar.png'), 'Resolution', 300);

close all;

sum(STR_AD_Count(1:6)) ./ sum(STR_AD_Count)
sum(GP_AD_Count(1:6)) ./ sum(GP_AD_Count)
sum(THAL_AD_Count(1:6)) ./ sum(THAL_AD_Count)
sum(HIPP_AD_Count(1:6)) ./ sum(HIPP_AD_Count)
sum(AMYG_AD_Count(1:6)) ./ sum(AMYG_AD_Count)

%% Spatial permutation testing for HCP and XYY

% Categorize HCP sex difference maps
PERM_HCP_STR_AD = zeros(26000, nPerm);
PERM_HCP_STR_AD(HCP.STAT.MOR.MSR.STR_SEX_B > 0) =  1;
PERM_HCP_STR_AD(HCP.STAT.MOR.MSR.STR_SEX_B < 0) = -1;
PERM_HCP_STR_AD(HCP.STAT.MOR.MSR.STR_SEX_Q & (HCP.STAT.MOR.MSR.STR_SEX_B > 0))  =  2;
PERM_HCP_STR_AD(HCP.STAT.MOR.MSR.STR_SEX_Q & (HCP.STAT.MOR.MSR.STR_SEX_B < 0))  = -2;

PERM_HCP_GP_AD = zeros(6000, nPerm);
PERM_HCP_GP_AD(HCP.STAT.MOR.MSR.GP_SEX_B > 0) =  1;
PERM_HCP_GP_AD(HCP.STAT.MOR.MSR.GP_SEX_B < 0) = -1;
PERM_HCP_GP_AD(HCP.STAT.MOR.MSR.GP_SEX_Q & (HCP.STAT.MOR.MSR.GP_SEX_B > 0))  =  2;
PERM_HCP_GP_AD(HCP.STAT.MOR.MSR.GP_SEX_Q & (HCP.STAT.MOR.MSR.GP_SEX_B < 0))  = -2;

PERM_HCP_THAL_AD = zeros(13000, nPerm);
PERM_HCP_THAL_AD(HCP.STAT.MOR.MSR.THAL_SEX_B > 0) =  1;
PERM_HCP_THAL_AD(HCP.STAT.MOR.MSR.THAL_SEX_B < 0) = -1;
PERM_HCP_THAL_AD(HCP.STAT.MOR.MSR.THAL_SEX_Q & (HCP.STAT.MOR.MSR.THAL_SEX_B > 0))  =  2;
PERM_HCP_THAL_AD(HCP.STAT.MOR.MSR.THAL_SEX_Q & (HCP.STAT.MOR.MSR.THAL_SEX_B < 0))  = -2;

PERM_HCP_HIPP_AD = zeros(2367, nPerm);
PERM_HCP_HIPP_AD(HCP.STAT.MOR.MSR.HIPP_SEX_B > 0) =  1;
PERM_HCP_HIPP_AD(HCP.STAT.MOR.MSR.HIPP_SEX_B < 0) = -1;
PERM_HCP_HIPP_AD(HCP.STAT.MOR.MSR.HIPP_SEX_Q & (HCP.STAT.MOR.MSR.HIPP_SEX_B > 0))  =  2;
PERM_HCP_HIPP_AD(HCP.STAT.MOR.MSR.HIPP_SEX_Q & (HCP.STAT.MOR.MSR.HIPP_SEX_B < 0))  = -2;

PERM_HCP_AMYG_AD = zeros(2878, nPerm);
PERM_HCP_AMYG_AD(HCP.STAT.MOR.MSR.AMYG_SEX_B > 0) =  1;
PERM_HCP_AMYG_AD(HCP.STAT.MOR.MSR.AMYG_SEX_B < 0) = -1;
PERM_HCP_AMYG_AD(HCP.STAT.MOR.MSR.AMYG_SEX_Q & (HCP.STAT.MOR.MSR.AMYG_SEX_B > 0))  =  2;
PERM_HCP_AMYG_AD(HCP.STAT.MOR.MSR.AMYG_SEX_Q & (HCP.STAT.MOR.MSR.AMYG_SEX_B < 0))  =  -2;

% Compute agreement maps
PERM_STR_AD = zeros(26000, nPerm);
PERM_STR_AD((PERM_HCP_STR_AD == 2) & repmat(XYY_STR_AD == -2, nPerm, 1)') = -1.9;
PERM_STR_AD((PERM_HCP_STR_AD == 2) & repmat(XYY_STR_AD == -1, nPerm, 1)') = -1.6;
PERM_STR_AD((PERM_HCP_STR_AD == -2) & repmat(XYY_STR_AD == -2, nPerm, 1)') = -0.9;
PERM_STR_AD((PERM_HCP_STR_AD == -2) & repmat(XYY_STR_AD == -1, nPerm, 1)') = -0.6;
PERM_STR_AD((PERM_HCP_STR_AD == 2) & repmat(XYY_STR_AD == 2, nPerm, 1)') = 0.3;
PERM_STR_AD((PERM_HCP_STR_AD == 2) & repmat(XYY_STR_AD == 1, nPerm, 1)') = 0.6;
PERM_STR_AD((PERM_HCP_STR_AD == -2) & repmat(XYY_STR_AD == 2, nPerm, 1)') = 1.3;
PERM_STR_AD((PERM_HCP_STR_AD == -2) & repmat(XYY_STR_AD == 1, nPerm, 1)') = 1.6;

PERM_GP_AD = zeros(6000, nPerm);
PERM_GP_AD((PERM_HCP_GP_AD == 2) & repmat(XYY_GP_AD == -2, nPerm, 1)') = -1.9;
PERM_GP_AD((PERM_HCP_GP_AD == 2) & repmat(XYY_GP_AD == -1, nPerm, 1)') = -1.6;
PERM_GP_AD((PERM_HCP_GP_AD == -2) & repmat(XYY_GP_AD == -2, nPerm, 1)') = -0.9;
PERM_GP_AD((PERM_HCP_GP_AD == -2) & repmat(XYY_GP_AD == -1, nPerm, 1)') = -0.6;
PERM_GP_AD((PERM_HCP_GP_AD == 2) & repmat(XYY_GP_AD == 2, nPerm, 1)') = 0.3;
PERM_GP_AD((PERM_HCP_GP_AD == 2) & repmat(XYY_GP_AD == 1, nPerm, 1)') = 0.6;
PERM_GP_AD((PERM_HCP_GP_AD == -2) & repmat(XYY_GP_AD == 2, nPerm, 1)') = 1.3;
PERM_GP_AD((PERM_HCP_GP_AD == -2) & repmat(XYY_GP_AD == 1, nPerm, 1)') = 1.6;

PERM_THAL_AD = zeros(13000, nPerm);
PERM_THAL_AD((PERM_HCP_THAL_AD == 2) & repmat(XYY_THAL_AD == -2, nPerm, 1)') = -1.9;
PERM_THAL_AD((PERM_HCP_THAL_AD == 2) & repmat(XYY_THAL_AD == -1, nPerm, 1)') = -1.6;
PERM_THAL_AD((PERM_HCP_THAL_AD == -2) & repmat(XYY_THAL_AD == -2, nPerm, 1)') = -0.9;
PERM_THAL_AD((PERM_HCP_THAL_AD == -2) & repmat(XYY_THAL_AD == -1, nPerm, 1)') = -0.6;
PERM_THAL_AD((PERM_HCP_THAL_AD == 2) & repmat(XYY_THAL_AD == 2, nPerm, 1)') = 0.3;
PERM_THAL_AD((PERM_HCP_THAL_AD == 2) & repmat(XYY_THAL_AD == 1, nPerm, 1)') = 0.6;
PERM_THAL_AD((PERM_HCP_THAL_AD == -2) & repmat(XYY_THAL_AD == 2, nPerm, 1)') = 1.3;
PERM_THAL_AD((PERM_HCP_THAL_AD == -2) & repmat(XYY_THAL_AD == 1, nPerm, 1)') = 1.6;

PERM_HIPP_AD = zeros(2367, nPerm);
PERM_HIPP_AD((PERM_HCP_HIPP_AD == 2) & repmat(XYY_HIPP_AD == -2, nPerm, 1)') = -1.9;
PERM_HIPP_AD((PERM_HCP_HIPP_AD == 2) & repmat(XYY_HIPP_AD == -1, nPerm, 1)') = -1.6;
PERM_HIPP_AD((PERM_HCP_HIPP_AD == -2) & repmat(XYY_HIPP_AD == -2, nPerm, 1)') = -0.9;
PERM_HIPP_AD((PERM_HCP_HIPP_AD == -2) & repmat(XYY_HIPP_AD == -1, nPerm, 1)') = -0.6;
PERM_HIPP_AD((PERM_HCP_HIPP_AD == 2) & repmat(XYY_HIPP_AD == 2, nPerm, 1)') = 0.3;
PERM_HIPP_AD((PERM_HCP_HIPP_AD == 2) & repmat(XYY_HIPP_AD == 1, nPerm, 1)') = 0.6;
PERM_HIPP_AD((PERM_HCP_HIPP_AD == -2) & repmat(XYY_HIPP_AD == 2, nPerm, 1)') = 1.3;
PERM_HIPP_AD((PERM_HCP_HIPP_AD == -2) & repmat(XYY_HIPP_AD == 1, nPerm, 1)') = 1.6;

PERM_AMYG_AD = zeros(2878, nPerm);
PERM_AMYG_AD((PERM_HCP_AMYG_AD == 2) & repmat(XYY_AMYG_AD == -2, nPerm, 1)') = -1.9;
PERM_AMYG_AD((PERM_HCP_AMYG_AD == 2) & repmat(XYY_AMYG_AD == -1, nPerm, 1)') = -1.6;
PERM_AMYG_AD((PERM_HCP_AMYG_AD == -2) & repmat(XYY_AMYG_AD == -2, nPerm, 1)') = -0.9;
PERM_AMYG_AD((PERM_HCP_AMYG_AD == -2) & repmat(XYY_AMYG_AD == -1, nPerm, 1)') = -0.6;
PERM_AMYG_AD((PERM_HCP_AMYG_AD == 2) & repmat(XYY_AMYG_AD == 2, nPerm, 1)') = 0.3;
PERM_AMYG_AD((PERM_HCP_AMYG_AD == 2) & repmat(XYY_AMYG_AD == 1, nPerm, 1)') = 0.6;
PERM_AMYG_AD((PERM_HCP_AMYG_AD == -2) & repmat(XYY_AMYG_AD == 2, nPerm, 1)') = 1.3;
PERM_AMYG_AD((PERM_HCP_AMYG_AD == -2) & repmat(XYY_AMYG_AD == 1, nPerm, 1)') = 1.6;

% Calculate proportions
PERM_STR = sum(repmat(S_STR_NA, 1, nPerm) .* (abs(PERM_STR_AD) > 1));
PERM_GP = sum(repmat(S_GP_NA, 1, nPerm) .* (abs(PERM_GP_AD) > 1));
PERM_THAL = sum(repmat(S_THAL_NA, 1, nPerm) .* (abs(PERM_THAL_AD) > 1));
PERM_HIPP = sum(repmat(S_HIPP_NA, 1, nPerm) .* (abs(PERM_HIPP_AD) > 1));
PERM_AMYG = sum(repmat(S_AMYG_NA, 1, nPerm) .* (abs(PERM_AMYG_AD) > 1));

% Calculate p value
STR_PERM_P = sum(sum(S_STR_NA(abs(STR_AD) > 1)) < PERM_STR) ./ nPerm;
GP_PERM_P = sum(sum(S_GP_NA(abs(GP_AD) > 1)) < PERM_GP) ./ nPerm;
THAL_PERM_P = sum(sum(S_THAL_NA(abs(THAL_AD) > 1)) < PERM_THAL) ./ nPerm;
HIPP_PERM_P = sum(sum(S_HIPP_NA(abs(HIPP_AD) > 1)) < PERM_HIPP) ./ nPerm;
AMYG_PERM_P = sum(sum(S_AMYG_NA(abs(AMYG_AD) > 1)) < PERM_AMYG) ./ nPerm;

disp('MSR permutations');
[STR_PERM_P GP_PERM_P THAL_PERM_P HIPP_PERM_P AMYG_PERM_P]

%% Regional congruence between HCP and IGD

% Combine significance and direction into a single map
IGD_STR_AD = zeros(1,26000);
HCP_STR_AD = zeros(1,26000);
STR_AD = zeros(1,26000);

IGD_STR_AD(IGD.STAT.MOR.TSA.STR_GRP_B > 0) =  1;
IGD_STR_AD(IGD.STAT.MOR.TSA.STR_GRP_B < 0) = -1;
IGD_STR_AD((IGD.STAT.MOR.TSA.STR_GRP_P < thresh) & (IGD.STAT.MOR.TSA.STR_GRP_B > 0))  =  2;
IGD_STR_AD((IGD.STAT.MOR.TSA.STR_GRP_P < thresh) & (IGD.STAT.MOR.TSA.STR_GRP_B < 0))  = -2;

HCP_STR_AD(HCP.STAT.MOR.TSA.STR_SEX_B > 0) =  1;
HCP_STR_AD(HCP.STAT.MOR.TSA.STR_SEX_B < 0) = -1;
HCP_STR_AD((HCP.STAT.MOR.TSA.STR_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.STR_SEX_B > 0))  =  2;
HCP_STR_AD((HCP.STAT.MOR.TSA.STR_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.STR_SEX_B < 0))  = -2;

% Generate agreement maps
STR_AD((HCP_STR_AD < 2) & (HCP_STR_AD > -2)) = 0;

STR_AD((HCP_STR_AD == 2) & (IGD_STR_AD == -2)) = -1.9;
STR_AD((HCP_STR_AD == 2) & (IGD_STR_AD == -1)) = -1.6;

STR_AD((HCP_STR_AD == -2) & (IGD_STR_AD == -2)) = -0.9;
STR_AD((HCP_STR_AD == -2) & (IGD_STR_AD == -1)) = -0.6;

STR_AD((HCP_STR_AD == 2) & (IGD_STR_AD == 2)) = 0.3;
STR_AD((HCP_STR_AD == 2) & (IGD_STR_AD == 1)) = 0.6;

STR_AD((HCP_STR_AD == -2) & (IGD_STR_AD == 2)) = 1.3;
STR_AD((HCP_STR_AD == -2) & (IGD_STR_AD == 1)) = 1.6;

% Visualize
h = figure('WindowState', 'Maximized'); HyoMagetViewFlat(STR_AD, [-2 2], cmap.quad_div);
exportgraphics(h, fullfile(figDir, 'igd_conjunction_STR.png'), 'Resolution', 300);

% Bar graphs
STR_AD_Count = [sum(S_STR_NA(STR_AD == 1.3)) sum(S_STR_NA(STR_AD == 1.6)) 0 ...
               sum(S_STR_NA(STR_AD == -1.9)) sum(S_STR_NA(STR_AD == -1.6)) 0 ...
               sum(S_STR_NA(STR_AD == -0.9)) sum(S_STR_NA(STR_AD == -0.6)) 0 ...
               sum(S_STR_NA(STR_AD == 0.3)) sum(S_STR_NA(STR_AD == 0.6)) 0];

h = figure;
b = bar(1, STR_AD_Count, 'stacked', 'EdgeColor', 'none');
for i = 1:length(STR_AD_Count)
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_quad_div(i,:);
end
box off
axis off
exportgraphics(h, fullfile(figDir, 'igd_conjunction_STR_bar.png'), 'Resolution', 300);

% Combine significance and direction into a single map
IGD_GP_AD = zeros(1,6000);
HCP_GP_AD = zeros(1,6000);
GP_AD = zeros(1,6000);

IGD_GP_AD(IGD.STAT.MOR.TSA.GP_GRP_B > 0) =  1;
IGD_GP_AD(IGD.STAT.MOR.TSA.GP_GRP_B < 0) = -1;
IGD_GP_AD((IGD.STAT.MOR.TSA.GP_GRP_P < thresh) & (IGD.STAT.MOR.TSA.GP_GRP_B > 0))  =  2;
IGD_GP_AD((IGD.STAT.MOR.TSA.GP_GRP_P < thresh) & (IGD.STAT.MOR.TSA.GP_GRP_B < 0))  = -2;

HCP_GP_AD(HCP.STAT.MOR.TSA.GP_SEX_B > 0) =  1;
HCP_GP_AD(HCP.STAT.MOR.TSA.GP_SEX_B < 0) = -1;
HCP_GP_AD((HCP.STAT.MOR.TSA.GP_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.GP_SEX_B > 0))  =  2;
HCP_GP_AD((HCP.STAT.MOR.TSA.GP_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.GP_SEX_B < 0))  = -2;

% Generate agreement maps
GP_AD((HCP_GP_AD < 2) & (HCP_GP_AD > -2)) = 0;

GP_AD((HCP_GP_AD == 2) & (IGD_GP_AD == -2)) = -1.9;
GP_AD((HCP_GP_AD == 2) & (IGD_GP_AD == -1)) = -1.6;

GP_AD((HCP_GP_AD == -2) & (IGD_GP_AD == -2)) = -0.9;
GP_AD((HCP_GP_AD == -2) & (IGD_GP_AD == -1)) = -0.6;

GP_AD((HCP_GP_AD == 2) & (IGD_GP_AD == 2)) = 0.3;
GP_AD((HCP_GP_AD == 2) & (IGD_GP_AD == 1)) = 0.6;

GP_AD((HCP_GP_AD == -2) & (IGD_GP_AD == 2)) = 1.3;
GP_AD((HCP_GP_AD == -2) & (IGD_GP_AD == 1)) = 1.6;

% Visualize
h = figure('WindowState', 'Maximized'); HyoMagetViewFlat(GP_AD, [-2 2], cmap.quad_div);
exportgraphics(h, fullfile(figDir, 'igd_conjunction_GP.png'), 'Resolution', 300);

% Bar graphs
GP_AD_Count = [sum(S_GP_NA(GP_AD == 1.3)) sum(S_GP_NA(GP_AD == 1.6)) 0 ...
               sum(S_GP_NA(GP_AD == -1.9)) sum(S_GP_NA(GP_AD == -1.6)) 0 ...
               sum(S_GP_NA(GP_AD == -0.9)) sum(S_GP_NA(GP_AD == -0.6)) 0 ...
               sum(S_GP_NA(GP_AD == 0.3)) sum(S_GP_NA(GP_AD == 0.6)) 0];

h = figure;
b = bar(1, GP_AD_Count, 'stacked', 'EdgeColor', 'none');
for i = 1:length(GP_AD_Count)
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_quad_div(i,:);
end
box off
axis off
exportgraphics(h, fullfile(figDir, 'igd_conjunction_GP_bar.png'), 'Resolution', 300);

% Combine significance and direction into a single map
IGD_THAL_AD = zeros(1,13000);
HCP_THAL_AD = zeros(1,13000);
THAL_AD = zeros(1,13000);

IGD_THAL_AD(IGD.STAT.MOR.TSA.THAL_GRP_B > 0) =  1;
IGD_THAL_AD(IGD.STAT.MOR.TSA.THAL_GRP_B < 0) = -1;
IGD_THAL_AD((IGD.STAT.MOR.TSA.THAL_GRP_P < thresh) & (IGD.STAT.MOR.TSA.THAL_GRP_B > 0))  =  2;
IGD_THAL_AD((IGD.STAT.MOR.TSA.THAL_GRP_P < thresh) & (IGD.STAT.MOR.TSA.THAL_GRP_B < 0))  = -2;

HCP_THAL_AD(HCP.STAT.MOR.TSA.THAL_SEX_B > 0) =  1;
HCP_THAL_AD(HCP.STAT.MOR.TSA.THAL_SEX_B < 0) = -1;
HCP_THAL_AD((HCP.STAT.MOR.TSA.THAL_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.THAL_SEX_B > 0))  =  2;
HCP_THAL_AD((HCP.STAT.MOR.TSA.THAL_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.THAL_SEX_B < 0))  = -2;

% Generate agreement maps
THAL_AD((HCP_THAL_AD < 2) & (HCP_THAL_AD > -2)) = 0;

THAL_AD((HCP_THAL_AD == 2) & (IGD_THAL_AD == -2)) = -1.9;
THAL_AD((HCP_THAL_AD == 2) & (IGD_THAL_AD == -1)) = -1.6;

THAL_AD((HCP_THAL_AD == -2) & (IGD_THAL_AD == -2)) = -0.9;
THAL_AD((HCP_THAL_AD == -2) & (IGD_THAL_AD == -1)) = -0.6;

THAL_AD((HCP_THAL_AD == 2) & (IGD_THAL_AD == 2)) = 0.3;
THAL_AD((HCP_THAL_AD == 2) & (IGD_THAL_AD == 1)) = 0.6;

THAL_AD((HCP_THAL_AD == -2) & (IGD_THAL_AD == 2)) = 1.3;
THAL_AD((HCP_THAL_AD == -2) & (IGD_THAL_AD == 1)) = 1.6;

% Visualize
h = figure('WindowState', 'Maximized'); HyoMagetViewFlat(THAL_AD, [-2 2], cmap.quad_div);
exportgraphics(h, fullfile(figDir, 'igd_conjunction_THAL.png'), 'Resolution', 300);

% Bar graphs
THAL_AD_Count = [sum(S_THAL_NA(THAL_AD == 1.3)) sum(S_THAL_NA(THAL_AD == 1.6)) 0 ...
               sum(S_THAL_NA(THAL_AD == -1.9)) sum(S_THAL_NA(THAL_AD == -1.6)) 0 ...
               sum(S_THAL_NA(THAL_AD == -0.9)) sum(S_THAL_NA(THAL_AD == -0.6)) 0 ...
               sum(S_THAL_NA(THAL_AD == 0.3)) sum(S_THAL_NA(THAL_AD == 0.6)) 0];

h = figure;
b = bar(1, THAL_AD_Count, 'stacked', 'EdgeColor', 'none');
for i = 1:length(THAL_AD_Count)
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_quad_div(i,:);
end
box off
axis off
exportgraphics(h, fullfile(figDir, 'igd_conjunction_THAL_bar.png'), 'Resolution', 300);

% Combine significance and direction into a single map
IGD_HIPP_AD = zeros(1,2367);
HCP_HIPP_AD = zeros(1,2367);
HIPP_AD = zeros(1,2367);

IGD_HIPP_AD(IGD.STAT.MOR.TSA.HIPP_GRP_B > 0) =  1;
IGD_HIPP_AD(IGD.STAT.MOR.TSA.HIPP_GRP_B < 0) = -1;
IGD_HIPP_AD((IGD.STAT.MOR.TSA.HIPP_GRP_P < thresh) & (IGD.STAT.MOR.TSA.HIPP_GRP_B > 0))  =  2;
IGD_HIPP_AD((IGD.STAT.MOR.TSA.HIPP_GRP_P < thresh) & (IGD.STAT.MOR.TSA.HIPP_GRP_B < 0))  = -2;

HCP_HIPP_AD(HCP.STAT.MOR.TSA.HIPP_SEX_B > 0) =  1;
HCP_HIPP_AD(HCP.STAT.MOR.TSA.HIPP_SEX_B < 0) = -1;
HCP_HIPP_AD((HCP.STAT.MOR.TSA.HIPP_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.HIPP_SEX_B > 0))  =  2;
HCP_HIPP_AD((HCP.STAT.MOR.TSA.HIPP_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.HIPP_SEX_B < 0))  = -2;

% Generate agreement maps
HIPP_AD((HCP_HIPP_AD < 2) & (HCP_HIPP_AD > -2)) = 0;

HIPP_AD((HCP_HIPP_AD == 2) & (IGD_HIPP_AD == -2)) = -1.9;
HIPP_AD((HCP_HIPP_AD == 2) & (IGD_HIPP_AD == -1)) = -1.6;

HIPP_AD((HCP_HIPP_AD == -2) & (IGD_HIPP_AD == -2)) = -0.9;
HIPP_AD((HCP_HIPP_AD == -2) & (IGD_HIPP_AD == -1)) = -0.6;

HIPP_AD((HCP_HIPP_AD == 2) & (IGD_HIPP_AD == 2)) = 0.3;
HIPP_AD((HCP_HIPP_AD == 2) & (IGD_HIPP_AD == 1)) = 0.6;

HIPP_AD((HCP_HIPP_AD == -2) & (IGD_HIPP_AD == 2)) = 1.3;
HIPP_AD((HCP_HIPP_AD == -2) & (IGD_HIPP_AD == 1)) = 1.6;

% Visualize
h = figure('WindowState', 'Maximized'); HyoMagetViewFlat(HIPP_AD, [-2 2], cmap.quad_div);
exportgraphics(h, fullfile(figDir, 'igd_conjunction_HIPP.png'), 'Resolution', 300);

% Bar graphs
HIPP_AD_Count = [sum(S_HIPP_NA(HIPP_AD == 1.3)) sum(S_HIPP_NA(HIPP_AD == 1.6)) 0 ...
               sum(S_HIPP_NA(HIPP_AD == -1.9)) sum(S_HIPP_NA(HIPP_AD == -1.6)) 0 ...
               sum(S_HIPP_NA(HIPP_AD == -0.9)) sum(S_HIPP_NA(HIPP_AD == -0.6)) 0 ...
               sum(S_HIPP_NA(HIPP_AD == 0.3)) sum(S_HIPP_NA(HIPP_AD == 0.6)) 0];

h = figure;
b = bar(1, HIPP_AD_Count, 'stacked', 'EdgeColor', 'none');
for i = 1:length(HIPP_AD_Count)
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_quad_div(i,:);
end
box off
axis off
exportgraphics(h, fullfile(figDir, 'igd_conjunction_HIPP_bar.png'), 'Resolution', 300);

% Combine significance and direction into a single map
IGD_AMYG_AD = zeros(1,2878);
HCP_AMYG_AD = zeros(1,2878);
AMYG_AD = zeros(1,2878);

IGD_AMYG_AD(IGD.STAT.MOR.TSA.AMYG_GRP_B > 0) =  1;
IGD_AMYG_AD(IGD.STAT.MOR.TSA.AMYG_GRP_B < 0) = -1;
IGD_AMYG_AD((IGD.STAT.MOR.TSA.AMYG_GRP_P < thresh) & (IGD.STAT.MOR.TSA.AMYG_GRP_B > 0))  =  2;
IGD_AMYG_AD((IGD.STAT.MOR.TSA.AMYG_GRP_P < thresh) & (IGD.STAT.MOR.TSA.AMYG_GRP_B < 0))  = -2;

HCP_AMYG_AD(HCP.STAT.MOR.TSA.AMYG_SEX_B > 0) =  1;
HCP_AMYG_AD(HCP.STAT.MOR.TSA.AMYG_SEX_B < 0) = -1;
HCP_AMYG_AD((HCP.STAT.MOR.TSA.AMYG_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.AMYG_SEX_B > 0))  =  2;
HCP_AMYG_AD((HCP.STAT.MOR.TSA.AMYG_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.AMYG_SEX_B < 0))  = -2;

% Generate agreement maps
AMYG_AD((HCP_AMYG_AD < 2) & (HCP_AMYG_AD > -2)) = 0;

AMYG_AD((HCP_AMYG_AD == 2) & (IGD_AMYG_AD == -2)) = -1.9;
AMYG_AD((HCP_AMYG_AD == 2) & (IGD_AMYG_AD == -1)) = -1.6;

AMYG_AD((HCP_AMYG_AD == -2) & (IGD_AMYG_AD == -2)) = -0.9;
AMYG_AD((HCP_AMYG_AD == -2) & (IGD_AMYG_AD == -1)) = -0.6;

AMYG_AD((HCP_AMYG_AD == 2) & (IGD_AMYG_AD == 2)) = 0.3;
AMYG_AD((HCP_AMYG_AD == 2) & (IGD_AMYG_AD == 1)) = 0.6;

AMYG_AD((HCP_AMYG_AD == -2) & (IGD_AMYG_AD == 2)) = 1.3;
AMYG_AD((HCP_AMYG_AD == -2) & (IGD_AMYG_AD == 1)) = 1.6;

% Visualize
h = figure('WindowState', 'Maximized'); HyoMagetViewFlat(AMYG_AD, [-2 2], cmap.quad_div);
exportgraphics(h, fullfile(figDir, 'igd_conjunction_AMYG.png'), 'Resolution', 300);

% Bar graphs
AMYG_AD_Count = [sum(S_AMYG_NA(AMYG_AD == 1.3)) sum(S_AMYG_NA(AMYG_AD == 1.6)) 0 ...
               sum(S_AMYG_NA(AMYG_AD == -1.9)) sum(S_AMYG_NA(AMYG_AD == -1.6)) 0 ...
               sum(S_AMYG_NA(AMYG_AD == -0.9)) sum(S_AMYG_NA(AMYG_AD == -0.6)) 0 ...
               sum(S_AMYG_NA(AMYG_AD == 0.3)) sum(S_AMYG_NA(AMYG_AD == 0.6)) 0];

h = figure;
b = bar(1, AMYG_AD_Count, 'stacked', 'EdgeColor', 'none');
for i = 1:length(AMYG_AD_Count)
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_quad_div(i,:);
end
box off
axis off
exportgraphics(h, fullfile(figDir, 'igd_conjunction_AMYG_bar.png'), 'Resolution', 300);

close all;

sum(STR_AD_Count(1:6)) ./ sum(STR_AD_Count)
sum(GP_AD_Count(1:6)) ./ sum(GP_AD_Count)
sum(THAL_AD_Count(1:6)) ./ sum(THAL_AD_Count)
sum(HIPP_AD_Count(1:6)) ./ sum(HIPP_AD_Count)
sum(AMYG_AD_Count(1:6)) ./ sum(AMYG_AD_Count)

%% Spatial permutation testing for HCP and IGD

% Categorize HCP sex difference maps
PERM_HCP_STR_AD = zeros(26000, nPerm);
PERM_HCP_STR_AD(HCP.STAT.MOR.MSR.STR_SEX_B > 0) =  1;
PERM_HCP_STR_AD(HCP.STAT.MOR.MSR.STR_SEX_B < 0) = -1;
PERM_HCP_STR_AD(HCP.STAT.MOR.MSR.STR_SEX_Q & (HCP.STAT.MOR.MSR.STR_SEX_B > 0))  =  2;
PERM_HCP_STR_AD(HCP.STAT.MOR.MSR.STR_SEX_Q & (HCP.STAT.MOR.MSR.STR_SEX_B < 0))  = -2;

PERM_HCP_GP_AD = zeros(6000, nPerm);
PERM_HCP_GP_AD(HCP.STAT.MOR.MSR.GP_SEX_B > 0) =  1;
PERM_HCP_GP_AD(HCP.STAT.MOR.MSR.GP_SEX_B < 0) = -1;
PERM_HCP_GP_AD(HCP.STAT.MOR.MSR.GP_SEX_Q & (HCP.STAT.MOR.MSR.GP_SEX_B > 0))  =  2;
PERM_HCP_GP_AD(HCP.STAT.MOR.MSR.GP_SEX_Q & (HCP.STAT.MOR.MSR.GP_SEX_B < 0))  = -2;

PERM_HCP_THAL_AD = zeros(13000, nPerm);
PERM_HCP_THAL_AD(HCP.STAT.MOR.MSR.THAL_SEX_B > 0) =  1;
PERM_HCP_THAL_AD(HCP.STAT.MOR.MSR.THAL_SEX_B < 0) = -1;
PERM_HCP_THAL_AD(HCP.STAT.MOR.MSR.THAL_SEX_Q & (HCP.STAT.MOR.MSR.THAL_SEX_B > 0))  =  2;
PERM_HCP_THAL_AD(HCP.STAT.MOR.MSR.THAL_SEX_Q & (HCP.STAT.MOR.MSR.THAL_SEX_B < 0))  = -2;

PERM_HCP_HIPP_AD = zeros(2367, nPerm);
PERM_HCP_HIPP_AD(HCP.STAT.MOR.MSR.HIPP_SEX_B > 0) =  1;
PERM_HCP_HIPP_AD(HCP.STAT.MOR.MSR.HIPP_SEX_B < 0) = -1;
PERM_HCP_HIPP_AD(HCP.STAT.MOR.MSR.HIPP_SEX_Q & (HCP.STAT.MOR.MSR.HIPP_SEX_B > 0))  =  2;
PERM_HCP_HIPP_AD(HCP.STAT.MOR.MSR.HIPP_SEX_Q & (HCP.STAT.MOR.MSR.HIPP_SEX_B < 0))  = -2;

PERM_HCP_AMYG_AD = zeros(2878, nPerm);
PERM_HCP_AMYG_AD(HCP.STAT.MOR.MSR.AMYG_SEX_B > 0) =  1;
PERM_HCP_AMYG_AD(HCP.STAT.MOR.MSR.AMYG_SEX_B < 0) = -1;
PERM_HCP_AMYG_AD(HCP.STAT.MOR.MSR.AMYG_SEX_Q & (HCP.STAT.MOR.MSR.AMYG_SEX_B > 0))  =  2;
PERM_HCP_AMYG_AD(HCP.STAT.MOR.MSR.AMYG_SEX_Q & (HCP.STAT.MOR.MSR.AMYG_SEX_B < 0))  =  -2;

% Compute agreement maps
PERM_STR_AD = zeros(26000, nPerm);
PERM_STR_AD((PERM_HCP_STR_AD == 2) & repmat(IGD_STR_AD == -2, nPerm, 1)') = -1.9;
PERM_STR_AD((PERM_HCP_STR_AD == 2) & repmat(IGD_STR_AD == -1, nPerm, 1)') = -1.6;
PERM_STR_AD((PERM_HCP_STR_AD == -2) & repmat(IGD_STR_AD == -2, nPerm, 1)') = -0.9;
PERM_STR_AD((PERM_HCP_STR_AD == -2) & repmat(IGD_STR_AD == -1, nPerm, 1)') = -0.6;
PERM_STR_AD((PERM_HCP_STR_AD == 2) & repmat(IGD_STR_AD == 2, nPerm, 1)') = 0.3;
PERM_STR_AD((PERM_HCP_STR_AD == 2) & repmat(IGD_STR_AD == 1, nPerm, 1)') = 0.6;
PERM_STR_AD((PERM_HCP_STR_AD == -2) & repmat(IGD_STR_AD == 2, nPerm, 1)') = 1.3;
PERM_STR_AD((PERM_HCP_STR_AD == -2) & repmat(IGD_STR_AD == 1, nPerm, 1)') = 1.6;

PERM_GP_AD = zeros(6000, nPerm);
PERM_GP_AD((PERM_HCP_GP_AD == 2) & repmat(IGD_GP_AD == -2, nPerm, 1)') = -1.9;
PERM_GP_AD((PERM_HCP_GP_AD == 2) & repmat(IGD_GP_AD == -1, nPerm, 1)') = -1.6;
PERM_GP_AD((PERM_HCP_GP_AD == -2) & repmat(IGD_GP_AD == -2, nPerm, 1)') = -0.9;
PERM_GP_AD((PERM_HCP_GP_AD == -2) & repmat(IGD_GP_AD == -1, nPerm, 1)') = -0.6;
PERM_GP_AD((PERM_HCP_GP_AD == 2) & repmat(IGD_GP_AD == 2, nPerm, 1)') = 0.3;
PERM_GP_AD((PERM_HCP_GP_AD == 2) & repmat(IGD_GP_AD == 1, nPerm, 1)') = 0.6;
PERM_GP_AD((PERM_HCP_GP_AD == -2) & repmat(IGD_GP_AD == 2, nPerm, 1)') = 1.3;
PERM_GP_AD((PERM_HCP_GP_AD == -2) & repmat(IGD_GP_AD == 1, nPerm, 1)') = 1.6;

PERM_THAL_AD = zeros(13000, nPerm);
PERM_THAL_AD((PERM_HCP_THAL_AD == 2) & repmat(IGD_THAL_AD == -2, nPerm, 1)') = -1.9;
PERM_THAL_AD((PERM_HCP_THAL_AD == 2) & repmat(IGD_THAL_AD == -1, nPerm, 1)') = -1.6;
PERM_THAL_AD((PERM_HCP_THAL_AD == -2) & repmat(IGD_THAL_AD == -2, nPerm, 1)') = -0.9;
PERM_THAL_AD((PERM_HCP_THAL_AD == -2) & repmat(IGD_THAL_AD == -1, nPerm, 1)') = -0.6;
PERM_THAL_AD((PERM_HCP_THAL_AD == 2) & repmat(IGD_THAL_AD == 2, nPerm, 1)') = 0.3;
PERM_THAL_AD((PERM_HCP_THAL_AD == 2) & repmat(IGD_THAL_AD == 1, nPerm, 1)') = 0.6;
PERM_THAL_AD((PERM_HCP_THAL_AD == -2) & repmat(IGD_THAL_AD == 2, nPerm, 1)') = 1.3;
PERM_THAL_AD((PERM_HCP_THAL_AD == -2) & repmat(IGD_THAL_AD == 1, nPerm, 1)') = 1.6;

PERM_HIPP_AD = zeros(2367, nPerm);
PERM_HIPP_AD((PERM_HCP_HIPP_AD == 2) & repmat(IGD_HIPP_AD == -2, nPerm, 1)') = -1.9;
PERM_HIPP_AD((PERM_HCP_HIPP_AD == 2) & repmat(IGD_HIPP_AD == -1, nPerm, 1)') = -1.6;
PERM_HIPP_AD((PERM_HCP_HIPP_AD == -2) & repmat(IGD_HIPP_AD == -2, nPerm, 1)') = -0.9;
PERM_HIPP_AD((PERM_HCP_HIPP_AD == -2) & repmat(IGD_HIPP_AD == -1, nPerm, 1)') = -0.6;
PERM_HIPP_AD((PERM_HCP_HIPP_AD == 2) & repmat(IGD_HIPP_AD == 2, nPerm, 1)') = 0.3;
PERM_HIPP_AD((PERM_HCP_HIPP_AD == 2) & repmat(IGD_HIPP_AD == 1, nPerm, 1)') = 0.6;
PERM_HIPP_AD((PERM_HCP_HIPP_AD == -2) & repmat(IGD_HIPP_AD == 2, nPerm, 1)') = 1.3;
PERM_HIPP_AD((PERM_HCP_HIPP_AD == -2) & repmat(IGD_HIPP_AD == 1, nPerm, 1)') = 1.6;

PERM_AMYG_AD = zeros(2878, nPerm);
PERM_AMYG_AD((PERM_HCP_AMYG_AD == 2) & repmat(IGD_AMYG_AD == -2, nPerm, 1)') = -1.9;
PERM_AMYG_AD((PERM_HCP_AMYG_AD == 2) & repmat(IGD_AMYG_AD == -1, nPerm, 1)') = -1.6;
PERM_AMYG_AD((PERM_HCP_AMYG_AD == -2) & repmat(IGD_AMYG_AD == -2, nPerm, 1)') = -0.9;
PERM_AMYG_AD((PERM_HCP_AMYG_AD == -2) & repmat(IGD_AMYG_AD == -1, nPerm, 1)') = -0.6;
PERM_AMYG_AD((PERM_HCP_AMYG_AD == 2) & repmat(IGD_AMYG_AD == 2, nPerm, 1)') = 0.3;
PERM_AMYG_AD((PERM_HCP_AMYG_AD == 2) & repmat(IGD_AMYG_AD == 1, nPerm, 1)') = 0.6;
PERM_AMYG_AD((PERM_HCP_AMYG_AD == -2) & repmat(IGD_AMYG_AD == 2, nPerm, 1)') = 1.3;
PERM_AMYG_AD((PERM_HCP_AMYG_AD == -2) & repmat(IGD_AMYG_AD == 1, nPerm, 1)') = 1.6;

% Calculate proportions
PERM_STR = sum(repmat(S_STR_NA, 1, nPerm) .* (abs(PERM_STR_AD) > 1));
PERM_GP = sum(repmat(S_GP_NA, 1, nPerm) .* (abs(PERM_GP_AD) > 1));
PERM_THAL = sum(repmat(S_THAL_NA, 1, nPerm) .* (abs(PERM_THAL_AD) > 1));
PERM_HIPP = sum(repmat(S_HIPP_NA, 1, nPerm) .* (abs(PERM_HIPP_AD) > 1));
PERM_AMYG = sum(repmat(S_AMYG_NA, 1, nPerm) .* (abs(PERM_AMYG_AD) > 1));

% Calculate p value
STR_PERM_P = sum(sum(S_STR_NA(abs(STR_AD) > 1)) < PERM_STR) ./ nPerm;
GP_PERM_P = sum(sum(S_GP_NA(abs(GP_AD) > 1)) < PERM_GP) ./ nPerm;
THAL_PERM_P = sum(sum(S_THAL_NA(abs(THAL_AD) > 1)) < PERM_THAL) ./ nPerm;
HIPP_PERM_P = sum(sum(S_HIPP_NA(abs(HIPP_AD) > 1)) < PERM_HIPP) ./ nPerm;
AMYG_PERM_P = sum(sum(S_AMYG_NA(abs(AMYG_AD) > 1)) < PERM_AMYG) ./ nPerm;

disp('MSR permutations');
[STR_PERM_P GP_PERM_P THAL_PERM_P HIPP_PERM_P AMYG_PERM_P]

%% Raincloud plot for bulk sex differences for XXY

% Extract absolute volume
STR = ((XXY.VOL.STR_L(XXY.QC.indSTR) + XXY.VOL.STR_R(XXY.QC.indSTR))');
GP = ((XXY.VOL.GP_L(XXY.QC.indGP) + XXY.VOL.GP_R(XXY.QC.indGP))');
THAL = ((XXY.VOL.THAL_L(XXY.QC.indTHAL) + XXY.VOL.THAL_R(XXY.QC.indTHAL))');
HIPP = ((XXY.VOL.HIPP_L(XXY.QC.indHIPP) + XXY.VOL.HIPP_R(XXY.QC.indHIPP))');
AMYG = ((XXY.VOL.AMYG_L(XXY.QC.indAMYG) + XXY.VOL.AMYG_R(XXY.QC.indAMYG))');

% Mean-center covariates
SEX_STR = XXY.SEX(XXY.QC.indSTR);
AGE_STR = XXY.AGE(XXY.QC.indSTR) - mean(XXY.AGE(XXY.QC.indSTR));
TBV_STR = XXY.TBV(XXY.QC.indSTR) - mean(XXY.TBV(XXY.QC.indSTR));
SEX_GP = XXY.SEX(XXY.QC.indGP);
AGE_GP = XXY.AGE(XXY.QC.indGP) - mean(XXY.AGE(XXY.QC.indGP));
TBV_GP = XXY.TBV(XXY.QC.indGP) - mean(XXY.TBV(XXY.QC.indGP));
SEX_THAL = XXY.SEX(XXY.QC.indTHAL);
AGE_THAL = XXY.AGE(XXY.QC.indTHAL) - mean(XXY.AGE(XXY.QC.indTHAL));
TBV_THAL = XXY.TBV(XXY.QC.indTHAL) - mean(XXY.TBV(XXY.QC.indTHAL));
SEX_HIPP = XXY.SEX(XXY.QC.indHIPP);
AGE_HIPP = XXY.AGE(XXY.QC.indHIPP) - mean(XXY.AGE(XXY.QC.indHIPP));
TBV_HIPP = XXY.TBV(XXY.QC.indHIPP) - mean(XXY.TBV(XXY.QC.indHIPP));
SEX_AMYG = XXY.SEX(XXY.QC.indAMYG);
AGE_AMYG = XXY.AGE(XXY.QC.indAMYG) - mean(XXY.AGE(XXY.QC.indAMYG));
TBV_AMYG = XXY.TBV(XXY.QC.indAMYG) - mean(XXY.TBV(XXY.QC.indAMYG));

% Prepare terms (struct)
T_SEX_STR = term(SEX_STR);
T_AGE_STR = term(AGE_STR);
T_TBV_STR = term(TBV_STR);
T_SEX_GP = term(SEX_GP);
T_AGE_GP = term(AGE_GP);
T_TBV_GP = term(TBV_GP);
T_SEX_THAL = term(SEX_THAL);
T_AGE_THAL = term(AGE_THAL);
T_TBV_THAL = term(TBV_THAL);
T_SEX_HIPP = term(SEX_HIPP);
T_AGE_HIPP = term(AGE_HIPP);
T_TBV_HIPP = term(TBV_HIPP);
T_SEX_AMYG = term(SEX_AMYG);
T_AGE_AMYG = term(AGE_AMYG);
T_TBV_AMYG = term(TBV_AMYG);

% Model
M_STR = 1 + T_SEX_STR + T_AGE_STR + T_TBV_STR;
M_GP = 1 + T_SEX_GP + T_AGE_GP  + T_TBV_GP;
M_THAL = 1 + T_SEX_THAL + T_AGE_THAL + T_TBV_THAL;
M_HIPP = 1 + T_SEX_HIPP + T_AGE_HIPP + T_TBV_HIPP;
M_AMYG = 1 + T_SEX_AMYG + T_AGE_AMYG + T_TBV_AMYG;

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

% Regress out the age and TBV from the bulk volumes
for i = 1:size(STR,1)
    STR_SEX(i) = STR(i) - FM_STR.coef(3)*AGE_STR(i) - FM_STR.coef(4)*TBV_STR(i);
end

for i = 1:size(GP,1)
    GP_SEX(i) = GP(i) - FM_GP.coef(3)*AGE_GP(i) - FM_GP.coef(4)*TBV_GP(i);
end

for i = 1:size(THAL,1)
    THAL_SEX(i) = THAL(i) - FM_THAL.coef(3)*AGE_THAL(i) - FM_THAL.coef(4)*TBV_THAL(i);
end

for i = 1:size(HIPP,1)
    HIPP_SEX(i) = HIPP(i) - FM_HIPP.coef(3)*AGE_HIPP(i) - FM_HIPP.coef(4)*TBV_HIPP(i);
end

for i = 1:size(AMYG,1)
    AMYG_SEX(i) = AMYG(i) - FM_AMYG.coef(3)*AGE_AMYG(i) - FM_AMYG.coef(4)*TBV_AMYG(i);
end

% Separate the groups
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
rainData{1} = STR_M';
rainData{2} = STR_F';

group_inx = 1:2;
c =  [0.1529    0.3922    0.0980;...
      0.5569    0.0039    0.3216];
group_names = {'XY', 'XXY'};
condition_names = {'XY', 'XXY'};
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
exportgraphics(h, fullfile(figDir, 'raincloud_XXY_bulk_STR.png'), 'Resolution', 300);

% GP
rainData{1} = GP_M';
rainData{2} = GP_F';
group_inx = 1:2;
c =  [0.1529    0.3922    0.0980;...
      0.5569    0.0039    0.3216];
group_names = {'XY', 'XXY'};
condition_names = {'XY', 'XXY'};
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
exportgraphics(h, fullfile(figDir, 'raincloud_XXY_bulk_GP.png'), 'Resolution', 300);

% THAL
rainData{1} = THAL_M';
rainData{2} = THAL_F';
group_inx = 1:2;
c =  [0.1529    0.3922    0.0980;...
      0.5569    0.0039    0.3216];
group_names = {'XY', 'XXY'};
condition_names = {'XY', 'XXY'};
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
exportgraphics(h, fullfile(figDir, 'raincloud_XXY_bulk_THAL.png'), 'Resolution', 300);

% HIPP
rainData{1} = HIPP_M';
rainData{2} = HIPP_F';
group_inx = 1:2;
c =  [0.1529    0.3922    0.0980;...
      0.5569    0.0039    0.3216];
group_names = {'XY', 'XXY'};
condition_names = {'XY', 'XXY'};
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
exportgraphics(h, fullfile(figDir, 'raincloud_XXY_bulk_HIPP.png'), 'Resolution', 300);

% AMYG
rainData{1} = AMYG_M';
rainData{2} = AMYG_F';
group_inx = 1:2;
c =  [0.1529    0.3922    0.0980;...
      0.5569    0.0039    0.3216];
group_names = {'XY', 'XXY'};
condition_names = {'XY', 'XXY'};
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
exportgraphics(h, fullfile(figDir, 'raincloud_XXY_bulk_AMYG.png'), 'Resolution', 300);

%% Raincloud plot for bulk sex differences for XYY

% Extract absolute volume
STR = ((XYY.VOL.STR_L(XYY.QC.indSTR) + XYY.VOL.STR_R(XYY.QC.indSTR))');
GP = ((XYY.VOL.GP_L(XYY.QC.indGP) + XYY.VOL.GP_R(XYY.QC.indGP))');
THAL = ((XYY.VOL.THAL_L(XYY.QC.indTHAL) + XYY.VOL.THAL_R(XYY.QC.indTHAL))');
HIPP = ((XYY.VOL.HIPP_L(XYY.QC.indHIPP) + XYY.VOL.HIPP_R(XYY.QC.indHIPP))');
AMYG = ((XYY.VOL.AMYG_L(XYY.QC.indAMYG) + XYY.VOL.AMYG_R(XYY.QC.indAMYG))');

% Mean-center covariates
SEX_STR = XYY.SEX(XYY.QC.indSTR);
AGE_STR = XYY.AGE(XYY.QC.indSTR) - mean(XYY.AGE(XYY.QC.indSTR));
TBV_STR = XYY.TBV(XYY.QC.indSTR) - mean(XYY.TBV(XYY.QC.indSTR));
SEX_GP = XYY.SEX(XYY.QC.indGP);
AGE_GP = XYY.AGE(XYY.QC.indGP) - mean(XYY.AGE(XYY.QC.indGP));
TBV_GP = XYY.TBV(XYY.QC.indGP) - mean(XYY.TBV(XYY.QC.indGP));
SEX_THAL = XYY.SEX(XYY.QC.indTHAL);
AGE_THAL = XYY.AGE(XYY.QC.indTHAL) - mean(XYY.AGE(XYY.QC.indTHAL));
TBV_THAL = XYY.TBV(XYY.QC.indTHAL) - mean(XYY.TBV(XYY.QC.indTHAL));
SEX_HIPP = XYY.SEX(XYY.QC.indHIPP);
AGE_HIPP = XYY.AGE(XYY.QC.indHIPP) - mean(XYY.AGE(XYY.QC.indHIPP));
TBV_HIPP = XYY.TBV(XYY.QC.indHIPP) - mean(XYY.TBV(XYY.QC.indHIPP));
SEX_AMYG = XYY.SEX(XYY.QC.indAMYG);
AGE_AMYG = XYY.AGE(XYY.QC.indAMYG) - mean(XYY.AGE(XYY.QC.indAMYG));
TBV_AMYG = XYY.TBV(XYY.QC.indAMYG) - mean(XYY.TBV(XYY.QC.indAMYG));

% Prepare terms (struct)
T_SEX_STR = term(SEX_STR);
T_AGE_STR = term(AGE_STR);
T_TBV_STR = term(TBV_STR);
T_SEX_GP = term(SEX_GP);
T_AGE_GP = term(AGE_GP);
T_TBV_GP = term(TBV_GP);
T_SEX_THAL = term(SEX_THAL);
T_AGE_THAL = term(AGE_THAL);
T_TBV_THAL = term(TBV_THAL);
T_SEX_HIPP = term(SEX_HIPP);
T_AGE_HIPP = term(AGE_HIPP);
T_TBV_HIPP = term(TBV_HIPP);
T_SEX_AMYG = term(SEX_AMYG);
T_AGE_AMYG = term(AGE_AMYG);
T_TBV_AMYG = term(TBV_AMYG);

% Model
M_STR = 1 + T_SEX_STR + T_AGE_STR + T_TBV_STR;
M_GP = 1 + T_SEX_GP + T_AGE_GP  + T_TBV_GP;
M_THAL = 1 + T_SEX_THAL + T_AGE_THAL + T_TBV_THAL;
M_HIPP = 1 + T_SEX_HIPP + T_AGE_HIPP + T_TBV_HIPP;
M_AMYG = 1 + T_SEX_AMYG + T_AGE_AMYG + T_TBV_AMYG;

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

% Regress out the age and TBV from the bulk volumes
for i = 1:size(STR,1)
    STR_SEX(i) = STR(i) - FM_STR.coef(3)*AGE_STR(i) - FM_STR.coef(4)*TBV_STR(i);
end

for i = 1:size(GP,1)
    GP_SEX(i) = GP(i) - FM_GP.coef(3)*AGE_GP(i) - FM_GP.coef(4)*TBV_GP(i);
end

for i = 1:size(THAL,1)
    THAL_SEX(i) = THAL(i) - FM_THAL.coef(3)*AGE_THAL(i) - FM_THAL.coef(4)*TBV_THAL(i);
end

for i = 1:size(HIPP,1)
    HIPP_SEX(i) = HIPP(i) - FM_HIPP.coef(3)*AGE_HIPP(i) - FM_HIPP.coef(4)*TBV_HIPP(i);
end

for i = 1:size(AMYG,1)
    AMYG_SEX(i) = AMYG(i) - FM_AMYG.coef(3)*AGE_AMYG(i) - FM_AMYG.coef(4)*TBV_AMYG(i);
end

% Separate the groups
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
rainData{1} = STR_M';
rainData{2} = STR_F';

group_inx = 1:2;
c =  [0.1529    0.3922    0.0980;...
      0.5569    0.0039    0.3216];
group_names = {'XY', 'XYY'};
condition_names = {'XY', 'XYY'};
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
exportgraphics(h, fullfile(figDir, 'raincloud_XYY_bulk_STR.png'), 'Resolution', 300);

% GP
rainData{1} = GP_M';
rainData{2} = GP_F';
group_inx = 1:2;
c =  [0.1529    0.3922    0.0980;...
      0.5569    0.0039    0.3216];
group_names = {'XY', 'XYY'};
condition_names = {'XY', 'XYY'};
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
exportgraphics(h, fullfile(figDir, 'raincloud_XYY_bulk_GP.png'), 'Resolution', 300);

% THAL
rainData{1} = THAL_M';
rainData{2} = THAL_F';
group_inx = 1:2;
c =  [0.1529    0.3922    0.0980;...
      0.5569    0.0039    0.3216];
group_names = {'XY', 'XYY'};
condition_names = {'XY', 'XYY'};
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
exportgraphics(h, fullfile(figDir, 'raincloud_XYY_bulk_THAL.png'), 'Resolution', 300);

% HIPP
rainData{1} = HIPP_M';
rainData{2} = HIPP_F';
group_inx = 1:2;
c =  [0.1529    0.3922    0.0980;...
      0.5569    0.0039    0.3216];
group_names = {'XY', 'XYY'};
condition_names = {'XY', 'XYY'};
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
exportgraphics(h, fullfile(figDir, 'raincloud_XYY_bulk_HIPP.png'), 'Resolution', 300);

% AMYG
rainData{1} = AMYG_M';
rainData{2} = AMYG_F';
group_inx = 1:2;
c =  [0.1529    0.3922    0.0980;...
      0.5569    0.0039    0.3216];
group_names = {'XY', 'XYY'};
condition_names = {'XY', 'XYY'};
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
exportgraphics(h, fullfile(figDir, 'raincloud_XYY_bulk_AMYG.png'), 'Resolution', 300);

%% Raincloud plot for bulk sex differences for IGD

% Extract absolute volume
STR = ((IGD.VOL.STR_L(IGD.QC.indSTR) + IGD.VOL.STR_R(IGD.QC.indSTR))');
GP = ((IGD.VOL.GP_L(IGD.QC.indGP) + IGD.VOL.GP_R(IGD.QC.indGP))');
THAL = ((IGD.VOL.THAL_L(IGD.QC.indTHAL) + IGD.VOL.THAL_R(IGD.QC.indTHAL))');
HIPP = ((IGD.VOL.HIPP_L(IGD.QC.indHIPP) + IGD.VOL.HIPP_R(IGD.QC.indHIPP))');
AMYG = ((IGD.VOL.AMYG_L(IGD.QC.indAMYG) + IGD.VOL.AMYG_R(IGD.QC.indAMYG))');

% Mean-center covariates
SEX_STR = IGD.GRP(IGD.QC.indSTR);
AGE_STR = IGD.AGE(IGD.QC.indSTR) - mean(IGD.AGE(IGD.QC.indSTR));
TBV_STR = IGD.TBV(IGD.QC.indSTR) - mean(IGD.TBV(IGD.QC.indSTR));
SEX_GP = IGD.GRP(IGD.QC.indGP);
AGE_GP = IGD.AGE(IGD.QC.indGP) - mean(IGD.AGE(IGD.QC.indGP));
TBV_GP = IGD.TBV(IGD.QC.indGP) - mean(IGD.TBV(IGD.QC.indGP));
SEX_THAL = IGD.GRP(IGD.QC.indTHAL);
AGE_THAL = IGD.AGE(IGD.QC.indTHAL) - mean(IGD.AGE(IGD.QC.indTHAL));
TBV_THAL = IGD.TBV(IGD.QC.indTHAL) - mean(IGD.TBV(IGD.QC.indTHAL));
SEX_HIPP = IGD.GRP(IGD.QC.indHIPP);
AGE_HIPP = IGD.AGE(IGD.QC.indHIPP) - mean(IGD.AGE(IGD.QC.indHIPP));
TBV_HIPP = IGD.TBV(IGD.QC.indHIPP) - mean(IGD.TBV(IGD.QC.indHIPP));
SEX_AMYG = IGD.GRP(IGD.QC.indAMYG);
AGE_AMYG = IGD.AGE(IGD.QC.indAMYG) - mean(IGD.AGE(IGD.QC.indAMYG));
TBV_AMYG = IGD.TBV(IGD.QC.indAMYG) - mean(IGD.TBV(IGD.QC.indAMYG));

% Prepare terms (struct)
T_SEX_STR = term(SEX_STR);
T_AGE_STR = term(AGE_STR);
T_TBV_STR = term(TBV_STR);
T_SEX_GP = term(SEX_GP);
T_AGE_GP = term(AGE_GP);
T_TBV_GP = term(TBV_GP);
T_SEX_THAL = term(SEX_THAL);
T_AGE_THAL = term(AGE_THAL);
T_TBV_THAL = term(TBV_THAL);
T_SEX_HIPP = term(SEX_HIPP);
T_AGE_HIPP = term(AGE_HIPP);
T_TBV_HIPP = term(TBV_HIPP);
T_SEX_AMYG = term(SEX_AMYG);
T_AGE_AMYG = term(AGE_AMYG);
T_TBV_AMYG = term(TBV_AMYG);

% Model
M_STR = 1 + T_SEX_STR + T_AGE_STR + T_TBV_STR;
M_GP = 1 + T_SEX_GP + T_AGE_GP  + T_TBV_GP;
M_THAL = 1 + T_SEX_THAL + T_AGE_THAL + T_TBV_THAL;
M_HIPP = 1 + T_SEX_HIPP + T_AGE_HIPP + T_TBV_HIPP;
M_AMYG = 1 + T_SEX_AMYG + T_AGE_AMYG + T_TBV_AMYG;

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

% Regress out the age and TBV from the bulk volumes
for i = 1:size(STR,1)
    STR_SEX(i) = STR(i) - FM_STR.coef(3)*AGE_STR(i) - FM_STR.coef(4)*TBV_STR(i);
end

for i = 1:size(GP,1)
    GP_SEX(i) = GP(i) - FM_GP.coef(3)*AGE_GP(i) - FM_GP.coef(4)*TBV_GP(i);
end

for i = 1:size(THAL,1)
    THAL_SEX(i) = THAL(i) - FM_THAL.coef(3)*AGE_THAL(i) - FM_THAL.coef(4)*TBV_THAL(i);
end

for i = 1:size(HIPP,1)
    HIPP_SEX(i) = HIPP(i) - FM_HIPP.coef(3)*AGE_HIPP(i) - FM_HIPP.coef(4)*TBV_HIPP(i);
end

for i = 1:size(AMYG,1)
    AMYG_SEX(i) = AMYG(i) - FM_AMYG.coef(3)*AGE_AMYG(i) - FM_AMYG.coef(4)*TBV_AMYG(i);
end

% Separate the groups
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
group_names = {'XY', 'IGD'};
condition_names = {'XY', 'IGD'};
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
exportgraphics(h, fullfile(figDir, 'raincloud_IGD_bulk_STR.png'), 'Resolution', 300);

% GP
rainData{1} = GP_F';
rainData{2} = GP_M';
group_inx = 1:2;
c =  [0.1529    0.3922    0.0980;...
      0.5569    0.0039    0.3216];
group_names = {'XY', 'IGD'};
condition_names = {'XY', 'IGD'};
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
exportgraphics(h, fullfile(figDir, 'raincloud_IGD_bulk_GP.png'), 'Resolution', 300);

% THAL
rainData{1} = THAL_F';
rainData{2} = THAL_M';
group_inx = 1:2;
c =  [0.1529    0.3922    0.0980;...
      0.5569    0.0039    0.3216];
group_names = {'XY', 'IGD'};
condition_names = {'XY', 'IGD'};
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
exportgraphics(h, fullfile(figDir, 'raincloud_IGD_bulk_THAL.png'), 'Resolution', 300);

% HIPP
rainData{1} = HIPP_F';
rainData{2} = HIPP_M';
group_inx = 1:2;
c =  [0.1529    0.3922    0.0980;...
      0.5569    0.0039    0.3216];
group_names = {'XY', 'IGD'};
condition_names = {'XY', 'IGD'};
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
exportgraphics(h, fullfile(figDir, 'raincloud_IGD_bulk_HIPP.png'), 'Resolution', 300);

% AMYG
rainData{1} = AMYG_F';
rainData{2} = AMYG_M';
group_inx = 1:2;
c =  [0.1529    0.3922    0.0980;...
      0.5569    0.0039    0.3216];
group_names = {'XY', 'IGD'};
condition_names = {'XY', 'IGD'};
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
exportgraphics(h, fullfile(figDir, 'raincloud_IGD_bulk_AMYG.png'), 'Resolution', 300);
