%% abcd_stats.m - Sex differences in pre/early puberty (ABCD-9Y) and comparison with adulthood
%
% Lee HM, Liu S, Lalonde FM, et al. Sex Differences in Human Subcortical Anatomy
% are Robust and Replicable in Adulthood, Already Evident Prepubertally, and Aligned
% with Sex Chromosomal and Gonadal Influences. bioRxiv (2026).
% https://doi.org/10.64898/2026.10.05.756728
%
% Purpose:
%   Estimates subcortical and cortical sex differences in the ABCD-9Y
%   sample, maps pubertal shifts (HCP-YA minus ABCD-9Y effect sizes), and
%   tests spatial similarity and conjunction between ages with Moran
%   spectral randomization (subcortex) and spin tests (cortex).
%
% Inputs (participant-level data are not shared; see Data availability):
%   - <dataDir>/ABCD.mat (struct ABCD)
%   - <dataDir>/HCP_results.mat (from b_hcp_stats.m)
%   - <templateDir>/maget_template_surfaces_distance_area.mat
%   - fsaverage_7.1.0.mat, glasser_fsaverage_7.1.0.mat, Vasa_G360.mat
%     (cortical surface, Glasser parcellation, and spin permutations; on the
%     MATLAB path via toolboxDir)
%
% Outputs:
%   - <dataDir>/ABCD_results.mat
%   - Figures in <figDir>
%
% Dependencies (not included in this repository):
%   - SurfStat (SurfStatLinMod, SurfStatT, term)
%   - fdr_bh (Benjamini-Hochberg FDR; MATLAB File Exchange)
%   - MATLAB Statistics and Machine Learning Toolbox
%   - daviolinplot, hexscatter (MATLAB File Exchange)
%   - Custom visualization helpers used only for figures (not included):
%     loadColorMaps, ThreshColorMap, HyoMagetView* (MAGeT-Brain surface viewers), HyoSurfacePlotFunction
%
% Sections (%%) are designed to be run in order using "Run Section".

clear all; close all; clc;

%% Paths (edit to match your local setup)

toolboxDir         = '/path/to/matlab_toolboxes';          % SurfStat, fdr_bh, BrainSpace, plotting helpers
templateDir        = '/path/to/maget_templates';           % MAGeT-Brain template surfaces
codeDir            = '/path/to/this_repository';
dataDir            = '/path/to/processed_data';            % participant-level data and results (not shared)
corticalResultsDir = '/path/to/cortical_study_results';    % group-level cortical maps (Lee et al., 2026)
figDir             = '/path/to/figures/abcd';

addpath(genpath(toolboxDir));
addpath(genpath(codeDir));

%% Prep

% Load ABCD data: preprocessed participant-level data
load(fullfile(dataDir, 'ABCD.mat'));

% HCP-YA results (from hcp_stats.m), required in either case
load(fullfile(dataDir, 'HCP_results.mat'));

% Load vertex-wise surface area maps for template surfaces
load(fullfile(templateDir, 'maget_template_surfaces_distance_area.mat'));

% Cortical surface template, Glasser parcellation and spin permutations (on the MATLAB path)
load('fsaverage_7.1.0.mat');
load('glasser_fsaverage_7.1.0.mat');
load('Vasa_G360.mat');

% Misc
loadColorMaps;
nPerm = 10000;
thresh = 0.05;

%% Compute sex effects for cortex

% Compute Z-score
CV = zscore(ABCD.COR.CV);
SA = zscore(ABCD.COR.SA);
CT = zscore(ABCD.COR.CT);

% Compute global phenotypes
CV_GS = sum(ABCD.COR.CV,2);
SA_GS = sum(ABCD.COR.SA,2);
CT_GS = sum(ABCD.COR.CT .* (ABCD.COR.SA ./ sum(ABCD.COR.SA,2)),2);

% Prepare terms
SEX = term(ABCD.SEX);
AGE = term(ABCD.AGE - mean(ABCD.AGE));
TBV = term(ABCD.TBV - mean(ABCD.TBV));
EUL = term(ABCD.EUL - mean(ABCD.EUL));
SITE = term(ABCD.SITE);
G_CV = term(CV_GS - mean(CV_GS));
G_SA = term(SA_GS - mean(SA_GS));
G_CT = term(CT_GS - mean(CT_GS));

% Define GLM
M_CV = 1 + SEX + AGE + EUL + G_CV + SITE;
M_SA = 1 + SEX + AGE + EUL + G_SA + SITE;
M_CT = 1 + SEX + AGE + EUL + G_CT + SITE;

% Compute stats
ABCD.STAT.COR.CV = SurfStatLinMod(CV, M_CV);
ABCD.STAT.COR.CV_SEX = SurfStatT(ABCD.STAT.COR.CV, SEX);
ABCD.STAT.COR.CV_B = ABCD.STAT.COR.CV_SEX.coef(2,:);
ABCD.STAT.COR.CV_P = 2*(1 - tcdf(abs(ABCD.STAT.COR.CV_SEX.t), ABCD.STAT.COR.CV_SEX.df));
[~,~,~,ABCD.STAT.COR.CV_Q] = fdr_bh(ABCD.STAT.COR.CV_P);

ABCD.STAT.COR.SA = SurfStatLinMod(SA, M_SA);
ABCD.STAT.COR.SA_SEX = SurfStatT(ABCD.STAT.COR.SA, SEX);
ABCD.STAT.COR.SA_B = ABCD.STAT.COR.SA_SEX.coef(2,:);
ABCD.STAT.COR.SA_P = 2*(1 - tcdf(abs(ABCD.STAT.COR.SA_SEX.t), ABCD.STAT.COR.SA_SEX.df));
[~,~,~,ABCD.STAT.COR.SA_Q] = fdr_bh(ABCD.STAT.COR.SA_P);

ABCD.STAT.COR.CT = SurfStatLinMod(CT, M_CT);
ABCD.STAT.COR.CT_SEX = SurfStatT(ABCD.STAT.COR.CT, SEX);
ABCD.STAT.COR.CT_B = ABCD.STAT.COR.CT_SEX.coef(2,:);
ABCD.STAT.COR.CT_P = 2*(1 - tcdf(abs(ABCD.STAT.COR.CT_SEX.t), ABCD.STAT.COR.CT_SEX.df));
[~,~,~,ABCD.STAT.COR.CT_Q] = fdr_bh(ABCD.STAT.COR.CT_P);

%% Compute pubertal shift maps for cortex

% Compute parameters for pubertal shift z-score mapping
ABCD_CV_XI = inv(ABCD.STAT.COR.CV.X' * ABCD.STAT.COR.CV.X);
ABCD_CV_SE = sqrt((ABCD.STAT.COR.CV.SSE / ABCD.STAT.COR.CV.df) * ABCD_CV_XI(2,2));
HCP_CV_XI = inv(HCP.STAT.COR.CV.X' * HCP.STAT.COR.CV.X);
HCP_CV_SE = sqrt((HCP.STAT.COR.CV.SSE / HCP.STAT.COR.CV.df) * HCP_CV_XI(2,2));
CV_DF = HCP.STAT.COR.CV.df + ABCD.STAT.COR.CV.df;

ABCD_SA_XI = inv(ABCD.STAT.COR.SA.X' * ABCD.STAT.COR.SA.X);
ABCD_SA_SE = sqrt((ABCD.STAT.COR.SA.SSE / ABCD.STAT.COR.SA.df) * ABCD_SA_XI(2,2));
HCP_SA_XI = inv(HCP.STAT.COR.SA.X' * HCP.STAT.COR.SA.X);
HCP_SA_SE = sqrt((HCP.STAT.COR.SA.SSE / HCP.STAT.COR.SA.df) * HCP_SA_XI(2,2));
SA_DF = HCP.STAT.COR.SA.df + ABCD.STAT.COR.SA.df;

ABCD_CT_XI = inv(ABCD.STAT.COR.CT.X' * ABCD.STAT.COR.CT.X);
ABCD_CT_SE = sqrt((ABCD.STAT.COR.CT.SSE / ABCD.STAT.COR.CT.df) * ABCD_CT_XI(2,2));
HCP_CT_XI = inv(HCP.STAT.COR.CT.X' * HCP.STAT.COR.CT.X);
HCP_CT_SE = sqrt((HCP.STAT.COR.CT.SSE / HCP.STAT.COR.CT.df) * HCP_CT_XI(2,2));
CT_DF = HCP.STAT.COR.CT.df + ABCD.STAT.COR.CT.df;

% Collect difference z score and significance maps
ABCD.STAT.COR.CV_PUB_Z = (HCP.STAT.COR.CV_B - ABCD.STAT.COR.CV_B) ./ sqrt(HCP_CV_SE.^2 + ABCD_CV_SE.^2);
ABCD.STAT.COR.CV_PUB_P = 2*(1 - tcdf(abs(ABCD.STAT.COR.CV_PUB_Z), CV_DF));
[~,~,~,ABCD.STAT.COR.CV_PUB_Q] = fdr_bh(ABCD.STAT.COR.CV_PUB_P);

ABCD.STAT.COR.SA_PUB_Z = (HCP.STAT.COR.SA_B - ABCD.STAT.COR.SA_B) ./ sqrt(HCP_SA_SE.^2 + ABCD_SA_SE.^2);
ABCD.STAT.COR.SA_PUB_P = 2*(1 - tcdf(abs(ABCD.STAT.COR.SA_PUB_Z), SA_DF));
[~,~,~,ABCD.STAT.COR.SA_PUB_Q] = fdr_bh(ABCD.STAT.COR.SA_PUB_P);

ABCD.STAT.COR.CT_PUB_Z = (HCP.STAT.COR.CT_B - ABCD.STAT.COR.CT_B) ./ sqrt(HCP_CT_SE.^2 + ABCD_CT_SE.^2);
ABCD.STAT.COR.CT_PUB_P = 2*(1 - tcdf(abs(ABCD.STAT.COR.CT_PUB_Z), CT_DF));
[~,~,~,ABCD.STAT.COR.CT_PUB_Q] = fdr_bh(ABCD.STAT.COR.CT_PUB_P);

%% Compute sex effects for subcortex

% Compute z-scores (struct)
STR_L = zscore(ABCD.SUB.VOL.STR_L(ABCD.QC.indSTR))';
STR_R = zscore(ABCD.SUB.VOL.STR_R(ABCD.QC.indSTR))';

GP_L = zscore(ABCD.SUB.VOL.GP_L(ABCD.QC.indGP))';
GP_R = zscore(ABCD.SUB.VOL.GP_R(ABCD.QC.indGP))';

THAL_L = zscore(ABCD.SUB.VOL.THAL_L(ABCD.QC.indTHAL))';
THAL_R = zscore(ABCD.SUB.VOL.THAL_R(ABCD.QC.indTHAL))';

HIPP_L = zscore(ABCD.SUB.VOL.HIPP_L(ABCD.QC.indHIPP))';
HIPP_R = zscore(ABCD.SUB.VOL.HIPP_R(ABCD.QC.indHIPP))';

AMYG_L = zscore(ABCD.SUB.VOL.AMYG_L(ABCD.QC.indAMYG))';
AMYG_R = zscore(ABCD.SUB.VOL.AMYG_R(ABCD.QC.indAMYG))';

% Compute z-scores (morph)
STR_MOR = zscore(ABCD.SUB.MOR.STR(ABCD.QC.indSTR,:));
GP_MOR = zscore(ABCD.SUB.MOR.GP(ABCD.QC.indGP,:));
THAL_MOR = zscore(ABCD.SUB.MOR.THAL(ABCD.QC.indTHAL,:));
HIPP_MOR = zscore(ABCD.SUB.MOR.HIPP(ABCD.QC.indHIPP,:));
AMYG_MOR = zscore(ABCD.SUB.MOR.AMYG(ABCD.QC.indAMYG,:));

% Prepare terms (struct)
SEX_STR = term(ABCD.SEX(ABCD.QC.indSTR));
AGE_STR = term(ABCD.AGE(ABCD.QC.indSTR) - mean(ABCD.AGE(ABCD.QC.indSTR)));
TBV_STR = term(ABCD.TBV(ABCD.QC.indSTR) - mean(ABCD.TBV(ABCD.QC.indSTR)));
SITE_STR = term(ABCD.SITE(ABCD.QC.indSTR));

SEX_GP = term(ABCD.SEX(ABCD.QC.indGP));
AGE_GP = term(ABCD.AGE(ABCD.QC.indGP) - mean(ABCD.AGE(ABCD.QC.indGP)));
TBV_GP = term(ABCD.TBV(ABCD.QC.indGP) - mean(ABCD.TBV(ABCD.QC.indGP)));
SITE_GP = term(ABCD.SITE(ABCD.QC.indGP));

SEX_THAL = term(ABCD.SEX(ABCD.QC.indTHAL));
AGE_THAL = term(ABCD.AGE(ABCD.QC.indTHAL) - mean(ABCD.AGE(ABCD.QC.indTHAL)));
TBV_THAL = term(ABCD.TBV(ABCD.QC.indTHAL) - mean(ABCD.TBV(ABCD.QC.indTHAL)));
SITE_THAL = term(ABCD.SITE(ABCD.QC.indTHAL));

SEX_HIPP = term(ABCD.SEX(ABCD.QC.indHIPP));
AGE_HIPP = term(ABCD.AGE(ABCD.QC.indHIPP) - mean(ABCD.AGE(ABCD.QC.indHIPP)));
TBV_HIPP = term(ABCD.TBV(ABCD.QC.indHIPP) - mean(ABCD.TBV(ABCD.QC.indHIPP)));
SITE_HIPP = term(ABCD.SITE(ABCD.QC.indHIPP));

SEX_AMYG = term(ABCD.SEX(ABCD.QC.indAMYG));
AGE_AMYG = term(ABCD.AGE(ABCD.QC.indAMYG) - mean(ABCD.AGE(ABCD.QC.indAMYG)));
TBV_AMYG = term(ABCD.TBV(ABCD.QC.indAMYG) - mean(ABCD.TBV(ABCD.QC.indAMYG)));
SITE_AMYG = term(ABCD.SITE(ABCD.QC.indAMYG));

% Prepare terms (morph)
TSA_STR = term(sum(ABCD.SUB.MOR.STR(ABCD.QC.indSTR,:),2) - mean(sum(ABCD.SUB.MOR.STR(ABCD.QC.indSTR,:),2)));
TSA_GP = term(sum(ABCD.SUB.MOR.GP(ABCD.QC.indGP,:),2) - mean(sum(ABCD.SUB.MOR.GP(ABCD.QC.indGP,:),2)));
TSA_THAL = term(sum(ABCD.SUB.MOR.THAL(ABCD.QC.indTHAL,:),2) - mean(sum(ABCD.SUB.MOR.THAL(ABCD.QC.indTHAL,:),2)));
TSA_HIPP = term(sum(ABCD.SUB.MOR.HIPP(ABCD.QC.indHIPP,:),2) - mean(sum(ABCD.SUB.MOR.HIPP(ABCD.QC.indHIPP,:),2)));
TSA_AMYG = term(sum(ABCD.SUB.MOR.AMYG(ABCD.QC.indAMYG,:),2) - mean(sum(ABCD.SUB.MOR.AMYG(ABCD.QC.indAMYG,:),2)));

% Define GLMs (struct)
ABCD.STAT.VOL.TBV.M_STR = 1 + SEX_STR + AGE_STR + TBV_STR + SITE_STR;
ABCD.STAT.VOL.TBV.M_GP = 1 + SEX_GP + AGE_GP + TBV_GP + SITE_GP;
ABCD.STAT.VOL.TBV.M_THAL = 1 + SEX_THAL + AGE_THAL + TBV_THAL + SITE_THAL;
ABCD.STAT.VOL.TBV.M_HIPP = 1 + SEX_HIPP + AGE_HIPP + TBV_HIPP + SITE_HIPP;
ABCD.STAT.VOL.TBV.M_AMYG = 1 + SEX_AMYG + AGE_AMYG + TBV_AMYG + SITE_AMYG;

% Define GLMs (morph)
ABCD.STAT.MOR.TSA.M_STR = 1 + SEX_STR + AGE_STR + TSA_STR + SITE_STR;
ABCD.STAT.MOR.TSA.M_GP = 1 + SEX_GP + AGE_GP + TSA_GP + SITE_GP;
ABCD.STAT.MOR.TSA.M_THAL = 1 + SEX_THAL + AGE_THAL + TSA_THAL + SITE_THAL;
ABCD.STAT.MOR.TSA.M_HIPP = 1 + SEX_HIPP + AGE_HIPP + TSA_HIPP + SITE_HIPP;
ABCD.STAT.MOR.TSA.M_AMYG = 1 + SEX_AMYG + AGE_AMYG + TSA_AMYG + SITE_AMYG;

% Fit models (struct)
ABCD.STAT.VOL.TBV.STR_L = SurfStatLinMod(STR_L, ABCD.STAT.VOL.TBV.M_STR);
ABCD.STAT.VOL.TBV.STR_L_SEX = SurfStatT(ABCD.STAT.VOL.TBV.STR_L, SEX_STR);
STR_L_SEX_P = 2*(1 - tcdf(abs(ABCD.STAT.VOL.TBV.STR_L_SEX.t), ABCD.STAT.VOL.TBV.STR_L_SEX.df));

ABCD.STAT.VOL.TBV.STR_R = SurfStatLinMod(STR_R, ABCD.STAT.VOL.TBV.M_STR);
ABCD.STAT.VOL.TBV.STR_R_SEX = SurfStatT(ABCD.STAT.VOL.TBV.STR_R, SEX_STR);
STR_R_SEX_P = 2*(1 - tcdf(abs(ABCD.STAT.VOL.TBV.STR_R_SEX.t), ABCD.STAT.VOL.TBV.STR_R_SEX.df));

ABCD.STAT.VOL.TBV.GP_L = SurfStatLinMod(GP_L, ABCD.STAT.VOL.TBV.M_GP);
ABCD.STAT.VOL.TBV.GP_L_SEX = SurfStatT(ABCD.STAT.VOL.TBV.GP_L, SEX_GP);
GP_L_SEX_P = 2*(1 - tcdf(abs(ABCD.STAT.VOL.TBV.GP_L_SEX.t), ABCD.STAT.VOL.TBV.GP_L_SEX.df));

ABCD.STAT.VOL.TBV.GP_R = SurfStatLinMod(GP_R, ABCD.STAT.VOL.TBV.M_GP);
ABCD.STAT.VOL.TBV.GP_R_SEX = SurfStatT(ABCD.STAT.VOL.TBV.GP_R, SEX_GP);
GP_R_SEX_P = 2*(1 - tcdf(abs(ABCD.STAT.VOL.TBV.GP_R_SEX.t), ABCD.STAT.VOL.TBV.GP_R_SEX.df));

ABCD.STAT.VOL.TBV.THAL_L = SurfStatLinMod(THAL_L, ABCD.STAT.VOL.TBV.M_THAL);
ABCD.STAT.VOL.TBV.THAL_L_SEX = SurfStatT(ABCD.STAT.VOL.TBV.THAL_L, SEX_THAL);
THAL_L_SEX_P = 2*(1 - tcdf(abs(ABCD.STAT.VOL.TBV.THAL_L_SEX.t), ABCD.STAT.VOL.TBV.THAL_L_SEX.df));

ABCD.STAT.VOL.TBV.THAL_R = SurfStatLinMod(THAL_R, ABCD.STAT.VOL.TBV.M_THAL);
ABCD.STAT.VOL.TBV.THAL_R_SEX = SurfStatT(ABCD.STAT.VOL.TBV.THAL_R, SEX_THAL);
THAL_R_SEX_P = 2*(1 - tcdf(abs(ABCD.STAT.VOL.TBV.THAL_R_SEX.t), ABCD.STAT.VOL.TBV.THAL_R_SEX.df));

ABCD.STAT.VOL.TBV.HIPP_L = SurfStatLinMod(HIPP_L, ABCD.STAT.VOL.TBV.M_HIPP);
ABCD.STAT.VOL.TBV.HIPP_L_SEX = SurfStatT(ABCD.STAT.VOL.TBV.HIPP_L, SEX_HIPP);
HIPP_L_SEX_P = 2*(1 - tcdf(abs(ABCD.STAT.VOL.TBV.HIPP_L_SEX.t), ABCD.STAT.VOL.TBV.HIPP_L_SEX.df));

ABCD.STAT.VOL.TBV.HIPP_R = SurfStatLinMod(HIPP_R, ABCD.STAT.VOL.TBV.M_HIPP);
ABCD.STAT.VOL.TBV.HIPP_R_SEX = SurfStatT(ABCD.STAT.VOL.TBV.HIPP_R, SEX_HIPP);
HIPP_R_SEX_P = 2*(1 - tcdf(abs(ABCD.STAT.VOL.TBV.HIPP_R_SEX.t), ABCD.STAT.VOL.TBV.HIPP_R_SEX.df));

ABCD.STAT.VOL.TBV.AMYG_L = SurfStatLinMod(AMYG_L, ABCD.STAT.VOL.TBV.M_AMYG);
ABCD.STAT.VOL.TBV.AMYG_L_SEX = SurfStatT(ABCD.STAT.VOL.TBV.AMYG_L, SEX_AMYG);
AMYG_L_SEX_P = 2*(1 - tcdf(abs(ABCD.STAT.VOL.TBV.AMYG_L_SEX.t), ABCD.STAT.VOL.TBV.AMYG_L_SEX.df));

ABCD.STAT.VOL.TBV.AMYG_R = SurfStatLinMod(AMYG_R, ABCD.STAT.VOL.TBV.M_AMYG);
ABCD.STAT.VOL.TBV.AMYG_R_SEX = SurfStatT(ABCD.STAT.VOL.TBV.AMYG_R, SEX_AMYG);
AMYG_R_SEX_P = 2*(1 - tcdf(abs(ABCD.STAT.VOL.TBV.AMYG_R_SEX.t), ABCD.STAT.VOL.TBV.AMYG_R_SEX.df));

% Collate statistical values
ABCD.STAT.VOL.TBV.SEX_B = [ABCD.STAT.VOL.TBV.STR_L_SEX.coef(2) ABCD.STAT.VOL.TBV.GP_L_SEX.coef(2) ABCD.STAT.VOL.TBV.THAL_L_SEX.coef(2) ABCD.STAT.VOL.TBV.HIPP_L_SEX.coef(2) ABCD.STAT.VOL.TBV.AMYG_L_SEX.coef(2) ...
                           ABCD.STAT.VOL.TBV.STR_R_SEX.coef(2) ABCD.STAT.VOL.TBV.GP_R_SEX.coef(2) ABCD.STAT.VOL.TBV.THAL_R_SEX.coef(2) ABCD.STAT.VOL.TBV.HIPP_R_SEX.coef(2) ABCD.STAT.VOL.TBV.AMYG_R_SEX.coef(2)];

ABCD.STAT.VOL.TBV.SEX_P = [STR_L_SEX_P GP_L_SEX_P THAL_L_SEX_P HIPP_L_SEX_P AMYG_L_SEX_P ...
                           STR_R_SEX_P GP_R_SEX_P THAL_R_SEX_P HIPP_R_SEX_P AMYG_R_SEX_P];

ABCD.STAT.VOL.TBV.SEX_Q = 10*ABCD.STAT.VOL.TBV.SEX_P;

% Fit models (morph)
ABCD.STAT.MOR.TSA.STR = SurfStatLinMod(STR_MOR, ABCD.STAT.MOR.TSA.M_STR);
ABCD.STAT.MOR.TSA.STR_SEX = SurfStatT(ABCD.STAT.MOR.TSA.STR, SEX_STR);
ABCD.STAT.MOR.TSA.STR_SEX_B = ABCD.STAT.MOR.TSA.STR_SEX.coef(2,:);
ABCD.STAT.MOR.TSA.STR_SEX_P = 2*(1 - tcdf(abs(ABCD.STAT.MOR.TSA.STR_SEX.t), ABCD.STAT.MOR.TSA.STR_SEX.df));
[~,~,~,ABCD.STAT.MOR.TSA.STR_SEX_Q] = fdr_bh(ABCD.STAT.MOR.TSA.STR_SEX_P);

ABCD.STAT.MOR.TSA.GP = SurfStatLinMod(GP_MOR, ABCD.STAT.MOR.TSA.M_GP);
ABCD.STAT.MOR.TSA.GP_SEX = SurfStatT(ABCD.STAT.MOR.TSA.GP, SEX_GP);
ABCD.STAT.MOR.TSA.GP_SEX_B = ABCD.STAT.MOR.TSA.GP_SEX.coef(2,:);
ABCD.STAT.MOR.TSA.GP_SEX_P = 2*(1 - tcdf(abs(ABCD.STAT.MOR.TSA.GP_SEX.t), ABCD.STAT.MOR.TSA.GP_SEX.df));
[~,~,~,ABCD.STAT.MOR.TSA.GP_SEX_Q] = fdr_bh(ABCD.STAT.MOR.TSA.GP_SEX_P);

ABCD.STAT.MOR.TSA.THAL = SurfStatLinMod(THAL_MOR, ABCD.STAT.MOR.TSA.M_THAL);
ABCD.STAT.MOR.TSA.THAL_SEX = SurfStatT(ABCD.STAT.MOR.TSA.THAL, SEX_THAL);
ABCD.STAT.MOR.TSA.THAL_SEX_B = ABCD.STAT.MOR.TSA.THAL_SEX.coef(2,:);
ABCD.STAT.MOR.TSA.THAL_SEX_P = 2*(1 - tcdf(abs(ABCD.STAT.MOR.TSA.THAL_SEX.t), ABCD.STAT.MOR.TSA.THAL_SEX.df));
[~,~,~,ABCD.STAT.MOR.TSA.THAL_SEX_Q] = fdr_bh(ABCD.STAT.MOR.TSA.THAL_SEX_P);

ABCD.STAT.MOR.TSA.HIPP = SurfStatLinMod(HIPP_MOR, ABCD.STAT.MOR.TSA.M_HIPP);
ABCD.STAT.MOR.TSA.HIPP_SEX = SurfStatT(ABCD.STAT.MOR.TSA.HIPP, SEX_HIPP);
ABCD.STAT.MOR.TSA.HIPP_SEX_B = ABCD.STAT.MOR.TSA.HIPP_SEX.coef(2,:);
ABCD.STAT.MOR.TSA.HIPP_SEX_P = 2*(1 - tcdf(abs(ABCD.STAT.MOR.TSA.HIPP_SEX.t), ABCD.STAT.MOR.TSA.HIPP_SEX.df));
[~,~,~,ABCD.STAT.MOR.TSA.HIPP_SEX_Q] = fdr_bh(ABCD.STAT.MOR.TSA.HIPP_SEX_P);

ABCD.STAT.MOR.TSA.AMYG = SurfStatLinMod(AMYG_MOR, ABCD.STAT.MOR.TSA.M_AMYG);
ABCD.STAT.MOR.TSA.AMYG_SEX = SurfStatT(ABCD.STAT.MOR.TSA.AMYG, SEX_AMYG);
ABCD.STAT.MOR.TSA.AMYG_SEX_B = ABCD.STAT.MOR.TSA.AMYG_SEX.coef(2,:);
ABCD.STAT.MOR.TSA.AMYG_SEX_P = 2*(1 - tcdf(abs(ABCD.STAT.MOR.TSA.AMYG_SEX.t), ABCD.STAT.MOR.TSA.AMYG_SEX.df));
[~,~,~,ABCD.STAT.MOR.TSA.AMYG_SEX_Q] = fdr_bh(ABCD.STAT.MOR.TSA.AMYG_SEX_P);

%% Pubertal shift mapping for subcortex

clear ABCD_VOL_XI ABCD_VOL_SE HCP_VOL_XI HCP_VOL_SE VOL_DF
ABCD_VOL_XI(:,:,1) = inv(ABCD.STAT.VOL.TBV.STR_L.X' * ABCD.STAT.VOL.TBV.STR_L.X);
ABCD_VOL_XI(:,:,2) = inv(ABCD.STAT.VOL.TBV.GP_L.X' * ABCD.STAT.VOL.TBV.GP_L.X);
ABCD_VOL_XI(:,:,3) = inv(ABCD.STAT.VOL.TBV.THAL_L.X' * ABCD.STAT.VOL.TBV.THAL_L.X);
ABCD_VOL_XI(:,:,4) = inv(ABCD.STAT.VOL.TBV.HIPP_L.X' * ABCD.STAT.VOL.TBV.HIPP_L.X);
ABCD_VOL_XI(:,:,5) = inv(ABCD.STAT.VOL.TBV.AMYG_L.X' * ABCD.STAT.VOL.TBV.AMYG_L.X);
ABCD_VOL_XI(:,:,6) = inv(ABCD.STAT.VOL.TBV.STR_R.X' * ABCD.STAT.VOL.TBV.STR_R.X);
ABCD_VOL_XI(:,:,7) = inv(ABCD.STAT.VOL.TBV.GP_R.X' * ABCD.STAT.VOL.TBV.GP_R.X);
ABCD_VOL_XI(:,:,8) = inv(ABCD.STAT.VOL.TBV.THAL_R.X' * ABCD.STAT.VOL.TBV.THAL_R.X);
ABCD_VOL_XI(:,:,9) = inv(ABCD.STAT.VOL.TBV.HIPP_R.X' * ABCD.STAT.VOL.TBV.HIPP_R.X);
ABCD_VOL_XI(:,:,10) = inv(ABCD.STAT.VOL.TBV.AMYG_R.X' * ABCD.STAT.VOL.TBV.AMYG_R.X);
ABCD_VOL_XI = squeeze(ABCD_VOL_XI(2,2,:));

ABCD_VOL_SE(1) = sqrt((ABCD.STAT.VOL.TBV.STR_L.SSE / ABCD.STAT.VOL.TBV.STR_L.df) * ABCD_VOL_XI(1));
ABCD_VOL_SE(2) = sqrt((ABCD.STAT.VOL.TBV.GP_L.SSE / ABCD.STAT.VOL.TBV.GP_L.df) * ABCD_VOL_XI(2));
ABCD_VOL_SE(3) = sqrt((ABCD.STAT.VOL.TBV.THAL_L.SSE / ABCD.STAT.VOL.TBV.THAL_L.df) * ABCD_VOL_XI(3));
ABCD_VOL_SE(4) = sqrt((ABCD.STAT.VOL.TBV.HIPP_L.SSE / ABCD.STAT.VOL.TBV.HIPP_L.df) * ABCD_VOL_XI(4));
ABCD_VOL_SE(5) = sqrt((ABCD.STAT.VOL.TBV.AMYG_L.SSE / ABCD.STAT.VOL.TBV.AMYG_L.df) * ABCD_VOL_XI(5));
ABCD_VOL_SE(6) = sqrt((ABCD.STAT.VOL.TBV.STR_R.SSE / ABCD.STAT.VOL.TBV.STR_R.df) * ABCD_VOL_XI(6));
ABCD_VOL_SE(7) = sqrt((ABCD.STAT.VOL.TBV.GP_R.SSE / ABCD.STAT.VOL.TBV.GP_R.df) * ABCD_VOL_XI(7));
ABCD_VOL_SE(8) = sqrt((ABCD.STAT.VOL.TBV.THAL_R.SSE / ABCD.STAT.VOL.TBV.THAL_R.df) * ABCD_VOL_XI(8));
ABCD_VOL_SE(9) = sqrt((ABCD.STAT.VOL.TBV.HIPP_R.SSE / ABCD.STAT.VOL.TBV.HIPP_R.df) * ABCD_VOL_XI(9));
ABCD_VOL_SE(10) = sqrt((ABCD.STAT.VOL.TBV.AMYG_R.SSE / ABCD.STAT.VOL.TBV.AMYG_R.df) * ABCD_VOL_XI(10));

HCP_VOL_XI(:,:,1) = inv(HCP.STAT.VOL.TBV.STR_L.X' * HCP.STAT.VOL.TBV.STR_L.X);
HCP_VOL_XI(:,:,2) = inv(HCP.STAT.VOL.TBV.GP_L.X' * HCP.STAT.VOL.TBV.GP_L.X);
HCP_VOL_XI(:,:,3) = inv(HCP.STAT.VOL.TBV.THAL_L.X' * HCP.STAT.VOL.TBV.THAL_L.X);
HCP_VOL_XI(:,:,4) = inv(HCP.STAT.VOL.TBV.HIPP_L.X' * HCP.STAT.VOL.TBV.HIPP_L.X);
HCP_VOL_XI(:,:,5) = inv(HCP.STAT.VOL.TBV.AMYG_L.X' * HCP.STAT.VOL.TBV.AMYG_L.X);
HCP_VOL_XI(:,:,6) = inv(HCP.STAT.VOL.TBV.STR_R.X' * HCP.STAT.VOL.TBV.STR_R.X);
HCP_VOL_XI(:,:,7) = inv(HCP.STAT.VOL.TBV.GP_R.X' * HCP.STAT.VOL.TBV.GP_R.X);
HCP_VOL_XI(:,:,8) = inv(HCP.STAT.VOL.TBV.THAL_R.X' * HCP.STAT.VOL.TBV.THAL_R.X);
HCP_VOL_XI(:,:,9) = inv(HCP.STAT.VOL.TBV.HIPP_R.X' * HCP.STAT.VOL.TBV.HIPP_R.X);
HCP_VOL_XI(:,:,10) = inv(HCP.STAT.VOL.TBV.AMYG_R.X' * HCP.STAT.VOL.TBV.AMYG_R.X);
HCP_VOL_XI = squeeze(HCP_VOL_XI(2,2,:));

HCP_VOL_SE(1) = sqrt((HCP.STAT.VOL.TBV.STR_L.SSE / HCP.STAT.VOL.TBV.STR_L.df) * HCP_VOL_XI(1));
HCP_VOL_SE(2) = sqrt((HCP.STAT.VOL.TBV.GP_L.SSE / HCP.STAT.VOL.TBV.GP_L.df) * HCP_VOL_XI(2));
HCP_VOL_SE(3) = sqrt((HCP.STAT.VOL.TBV.THAL_L.SSE / HCP.STAT.VOL.TBV.THAL_L.df) * HCP_VOL_XI(3));
HCP_VOL_SE(4) = sqrt((HCP.STAT.VOL.TBV.HIPP_L.SSE / HCP.STAT.VOL.TBV.HIPP_L.df) * HCP_VOL_XI(4));
HCP_VOL_SE(5) = sqrt((HCP.STAT.VOL.TBV.AMYG_L.SSE / HCP.STAT.VOL.TBV.AMYG_L.df) * HCP_VOL_XI(5));
HCP_VOL_SE(6) = sqrt((HCP.STAT.VOL.TBV.STR_R.SSE / HCP.STAT.VOL.TBV.STR_R.df) * HCP_VOL_XI(6));
HCP_VOL_SE(7) = sqrt((HCP.STAT.VOL.TBV.GP_R.SSE / HCP.STAT.VOL.TBV.GP_R.df) * HCP_VOL_XI(7));
HCP_VOL_SE(8) = sqrt((HCP.STAT.VOL.TBV.THAL_R.SSE / HCP.STAT.VOL.TBV.THAL_R.df) * HCP_VOL_XI(8));
HCP_VOL_SE(9) = sqrt((HCP.STAT.VOL.TBV.HIPP_R.SSE / HCP.STAT.VOL.TBV.HIPP_R.df) * HCP_VOL_XI(9));
HCP_VOL_SE(10) = sqrt((HCP.STAT.VOL.TBV.AMYG_R.SSE / HCP.STAT.VOL.TBV.AMYG_R.df) * HCP_VOL_XI(10));

VOL_DF(1) = HCP.STAT.VOL.TBV.STR_L.df + ABCD.STAT.VOL.TBV.STR_L.df;
VOL_DF(2) = HCP.STAT.VOL.TBV.GP_L.df + ABCD.STAT.VOL.TBV.GP_L.df;
VOL_DF(3) = HCP.STAT.VOL.TBV.THAL_L.df + ABCD.STAT.VOL.TBV.THAL_L.df;
VOL_DF(4) = HCP.STAT.VOL.TBV.HIPP_L.df + ABCD.STAT.VOL.TBV.HIPP_L.df;
VOL_DF(5) = HCP.STAT.VOL.TBV.AMYG_L.df + ABCD.STAT.VOL.TBV.AMYG_L.df;
VOL_DF(6) = HCP.STAT.VOL.TBV.STR_R.df + ABCD.STAT.VOL.TBV.STR_R.df;
VOL_DF(7) = HCP.STAT.VOL.TBV.GP_R.df + ABCD.STAT.VOL.TBV.GP_R.df;
VOL_DF(8) = HCP.STAT.VOL.TBV.THAL_R.df + ABCD.STAT.VOL.TBV.THAL_R.df;
VOL_DF(9) = HCP.STAT.VOL.TBV.HIPP_R.df + ABCD.STAT.VOL.TBV.HIPP_R.df;
VOL_DF(10) = HCP.STAT.VOL.TBV.AMYG_R.df + ABCD.STAT.VOL.TBV.AMYG_R.df;

ABCD.STAT.VOL.TBV.PUB_Z = (HCP.STAT.VOL.TBV.SEX_B - ABCD.STAT.VOL.TBV.SEX_B) ./ sqrt(HCP_VOL_SE.^2 + ABCD_VOL_SE.^2);
ABCD.STAT.VOL.TBV.PUB_P = 2*(1 - tcdf(abs(ABCD.STAT.VOL.TBV.PUB_Z), VOL_DF));
ABCD.STAT.VOL.TBV.PUB_Q = 10*(ABCD.STAT.VOL.TBV.PUB_P);

% Compute parameters for pubertal shift z-score mapping
ABCD_STR_XI = inv(ABCD.STAT.MOR.TSA.STR.X' * ABCD.STAT.MOR.TSA.STR.X);
ABCD_STR_SE = sqrt((ABCD.STAT.MOR.TSA.STR.SSE / ABCD.STAT.MOR.TSA.STR.df) * ABCD_STR_XI(2,2));
HCP_STR_XI = inv(HCP.STAT.MOR.TSA.STR.X' * HCP.STAT.MOR.TSA.STR.X);
HCP_STR_SE = sqrt((HCP.STAT.MOR.TSA.STR.SSE / HCP.STAT.MOR.TSA.STR.df) * HCP_STR_XI(2,2));
STR_DF = HCP.STAT.MOR.TSA.STR.df + ABCD.STAT.MOR.TSA.STR.df;

ABCD_GP_XI = inv(ABCD.STAT.MOR.TSA.GP.X' * ABCD.STAT.MOR.TSA.GP.X);
ABCD_GP_SE = sqrt((ABCD.STAT.MOR.TSA.GP.SSE / ABCD.STAT.MOR.TSA.GP.df) * ABCD_GP_XI(2,2));
HCP_GP_XI = inv(HCP.STAT.MOR.TSA.GP.X' * HCP.STAT.MOR.TSA.GP.X);
HCP_GP_SE = sqrt((HCP.STAT.MOR.TSA.GP.SSE / HCP.STAT.MOR.TSA.GP.df) * HCP_GP_XI(2,2));
GP_DF = HCP.STAT.MOR.TSA.GP.df + ABCD.STAT.MOR.TSA.GP.df;

ABCD_THAL_XI = inv(ABCD.STAT.MOR.TSA.THAL.X' * ABCD.STAT.MOR.TSA.THAL.X);
ABCD_THAL_SE = sqrt((ABCD.STAT.MOR.TSA.THAL.SSE / ABCD.STAT.MOR.TSA.THAL.df) * ABCD_THAL_XI(2,2));
HCP_THAL_XI = inv(HCP.STAT.MOR.TSA.THAL.X' * HCP.STAT.MOR.TSA.THAL.X);
HCP_THAL_SE = sqrt((HCP.STAT.MOR.TSA.THAL.SSE / HCP.STAT.MOR.TSA.THAL.df) * HCP_THAL_XI(2,2));
THAL_DF = HCP.STAT.MOR.TSA.THAL.df + ABCD.STAT.MOR.TSA.THAL.df;

ABCD_HIPP_XI = inv(ABCD.STAT.MOR.TSA.HIPP.X' * ABCD.STAT.MOR.TSA.HIPP.X);
ABCD_HIPP_SE = sqrt((ABCD.STAT.MOR.TSA.HIPP.SSE / ABCD.STAT.MOR.TSA.HIPP.df) * ABCD_HIPP_XI(2,2));
HCP_HIPP_XI = inv(HCP.STAT.MOR.TSA.HIPP.X' * HCP.STAT.MOR.TSA.HIPP.X);
HCP_HIPP_SE = sqrt((HCP.STAT.MOR.TSA.HIPP.SSE / HCP.STAT.MOR.TSA.HIPP.df) * HCP_HIPP_XI(2,2));
HIPP_DF = HCP.STAT.MOR.TSA.HIPP.df + ABCD.STAT.MOR.TSA.HIPP.df;

ABCD_AMYG_XI = inv(ABCD.STAT.MOR.TSA.AMYG.X' * ABCD.STAT.MOR.TSA.AMYG.X);
ABCD_AMYG_SE = sqrt((ABCD.STAT.MOR.TSA.AMYG.SSE / ABCD.STAT.MOR.TSA.AMYG.df) * ABCD_AMYG_XI(2,2));
HCP_AMYG_XI = inv(HCP.STAT.MOR.TSA.AMYG.X' * HCP.STAT.MOR.TSA.AMYG.X);
HCP_AMYG_SE = sqrt((HCP.STAT.MOR.TSA.AMYG.SSE / HCP.STAT.MOR.TSA.AMYG.df) * HCP_AMYG_XI(2,2));
AMYG_DF = HCP.STAT.MOR.TSA.AMYG.df + ABCD.STAT.MOR.TSA.AMYG.df;

% Collect difference z score and significance maps
ABCD.STAT.MOR.TSA.STR_PUB_Z = (HCP.STAT.MOR.TSA.STR_SEX_B - ABCD.STAT.MOR.TSA.STR_SEX_B) ./ sqrt(HCP_STR_SE.^2 + ABCD_STR_SE.^2);
ABCD.STAT.MOR.TSA.STR_PUB_P = 2*(1 - tcdf(abs(ABCD.STAT.MOR.TSA.STR_PUB_Z), STR_DF));
[~,~,~,ABCD.STAT.MOR.TSA.STR_PUB_Q] = fdr_bh(ABCD.STAT.MOR.TSA.STR_PUB_P);

ABCD.STAT.MOR.TSA.GP_PUB_Z = (HCP.STAT.MOR.TSA.GP_SEX_B - ABCD.STAT.MOR.TSA.GP_SEX_B) ./ sqrt(HCP_GP_SE.^2 + ABCD_GP_SE.^2);
ABCD.STAT.MOR.TSA.GP_PUB_P = 2*(1 - tcdf(abs(ABCD.STAT.MOR.TSA.GP_PUB_Z), GP_DF));
[~,~,~,ABCD.STAT.MOR.TSA.GP_PUB_Q] = fdr_bh(ABCD.STAT.MOR.TSA.GP_PUB_P);

ABCD.STAT.MOR.TSA.THAL_PUB_Z = (HCP.STAT.MOR.TSA.THAL_SEX_B - ABCD.STAT.MOR.TSA.THAL_SEX_B) ./ sqrt(HCP_THAL_SE.^2 + ABCD_THAL_SE.^2);
ABCD.STAT.MOR.TSA.THAL_PUB_P = 2*(1 - tcdf(abs(ABCD.STAT.MOR.TSA.THAL_PUB_Z), THAL_DF));
[~,~,~,ABCD.STAT.MOR.TSA.THAL_PUB_Q] = fdr_bh(ABCD.STAT.MOR.TSA.THAL_PUB_P);

ABCD.STAT.MOR.TSA.HIPP_PUB_Z = (HCP.STAT.MOR.TSA.HIPP_SEX_B - ABCD.STAT.MOR.TSA.HIPP_SEX_B) ./ sqrt(HCP_HIPP_SE.^2 + ABCD_HIPP_SE.^2);
ABCD.STAT.MOR.TSA.HIPP_PUB_P = 2*(1 - tcdf(abs(ABCD.STAT.MOR.TSA.HIPP_PUB_Z), HIPP_DF));
[~,~,~,ABCD.STAT.MOR.TSA.HIPP_PUB_Q] = fdr_bh(ABCD.STAT.MOR.TSA.HIPP_PUB_P);

ABCD.STAT.MOR.TSA.AMYG_PUB_Z = (HCP.STAT.MOR.TSA.AMYG_SEX_B - ABCD.STAT.MOR.TSA.AMYG_SEX_B) ./ sqrt(HCP_AMYG_SE.^2 + ABCD_AMYG_SE.^2);
ABCD.STAT.MOR.TSA.AMYG_PUB_P = 2*(1 - tcdf(abs(ABCD.STAT.MOR.TSA.AMYG_PUB_Z), AMYG_DF));
[~,~,~,ABCD.STAT.MOR.TSA.AMYG_PUB_Q] = fdr_bh(ABCD.STAT.MOR.TSA.AMYG_PUB_P);

[sum(S_STR_NA(ABCD.STAT.MOR.TSA.STR_PUB_Q < thresh)) ...
sum(S_GP_NA(ABCD.STAT.MOR.TSA.GP_PUB_Q < thresh)) ...
sum(S_THAL_NA(ABCD.STAT.MOR.TSA.THAL_PUB_Q < thresh)) ...
sum(S_HIPP_NA(ABCD.STAT.MOR.TSA.HIPP_PUB_Q < thresh)) ...
sum(S_AMYG_NA(ABCD.STAT.MOR.TSA.AMYG_PUB_Q < thresh))]

%% Save stats results

save(fullfile(dataDir, 'ABCD_results.mat'), 'ABCD', '-v7.3');

% Descriptive statistics: age by sex
[mean(ABCD.AGE(ABCD.SEX == 1)) std(ABCD.AGE(ABCD.SEX == 1))]
[mean(ABCD.AGE(ABCD.SEX == 0)) std(ABCD.AGE(ABCD.SEX == 0))]
[~,p] = ttest2(ABCD.AGE(ABCD.SEX == 1), ABCD.AGE(ABCD.SEX == 0))

%% Visualize cortical maps

% Standardized effect size maps (ABCD)
h = HyoSurfacePlotFunction(ABCD.STAT.COR.CV_B .* (ABCD.STAT.COR.CV_Q < thresh), IS, G360, cmap.PiYG, '', [-0.5 0.5], 'square');
exportgraphics(h, fullfile(figDir, 'abcd_cv_fdr.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(ABCD.STAT.COR.SA_B .* (ABCD.STAT.COR.SA_Q < thresh), IS, G360, cmap.PiYG, '', [-0.5 0.5], 'square');
exportgraphics(h, fullfile(figDir, 'abcd_sa_fdr.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(ABCD.STAT.COR.CT_B .* (ABCD.STAT.COR.CT_Q < thresh), IS, G360, cmap.PiYG, '', [-0.5 0.5], 'square');
exportgraphics(h, fullfile(figDir, 'abcd_ct_fdr.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(ABCD.STAT.COR.CV_B, IS, G360, cmap.PiYG, '', [-0.5 0.5], 'square');
exportgraphics(h, fullfile(figDir, 'abcd_cv.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(ABCD.STAT.COR.SA_B, IS, G360, cmap.PiYG, '', [-0.5 0.5], 'square');
exportgraphics(h, fullfile(figDir, 'abcd_sa.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(ABCD.STAT.COR.CT_B, IS, G360, cmap.PiYG, '', [-0.5 0.5], 'square');
exportgraphics(h, fullfile(figDir, 'abcd_ct.png'), 'Resolution', 300); close(h);

% Standardized effect size maps (HCP)
h = HyoSurfacePlotFunction(HCP.STAT.COR.CV_B .* (HCP.STAT.COR.CV_Q < thresh), IS, G360, cmap.PiYG, '', [-0.5 0.5], 'square');
exportgraphics(h, fullfile(figDir, 'hcp_cv_fdr.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(HCP.STAT.COR.SA_B .* (HCP.STAT.COR.SA_Q < thresh), IS, G360, cmap.PiYG, '', [-0.5 0.5], 'square');
exportgraphics(h, fullfile(figDir, 'hcp_sa_fdr.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(HCP.STAT.COR.CT_B .* (HCP.STAT.COR.CT_Q < thresh), IS, G360, cmap.PiYG, '', [-0.5 0.5], 'square');
exportgraphics(h, fullfile(figDir, 'hcp_ct_fdr.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(HCP.STAT.COR.CV_B, IS, G360, cmap.PiYG, '', [-0.5 0.5], 'square');
exportgraphics(h, fullfile(figDir, 'hcp_cv.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(HCP.STAT.COR.SA_B, IS, G360, cmap.PiYG, '', [-0.5 0.5], 'square');
exportgraphics(h, fullfile(figDir, 'hcp_sa.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction(HCP.STAT.COR.CT_B, IS, G360, cmap.PiYG, '', [-0.5 0.5], 'square');
exportgraphics(h, fullfile(figDir, 'hcp_ct.png'), 'Resolution', 300); close(h);

% Display proportion of cortical sheet with significant sex effect
[sum(ABCD.STAT.COR.CV_Q < thresh) sum(ABCD.STAT.COR.SA_Q < thresh) sum(ABCD.STAT.COR.CT_Q < thresh)] / 360

% Pubertal shift maps
h = HyoSurfacePlotFunction((HCP.STAT.COR.CV_B - ABCD.STAT.COR.CV_B) .* (ABCD.STAT.COR.CV_PUB_Q < thresh), IS, G360, cmap.PiYG, '', [-0.5 0.5], 'square');
exportgraphics(h, fullfile(figDir, 'pub_cv_fdr.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction((HCP.STAT.COR.SA_B - ABCD.STAT.COR.SA_B) .* (ABCD.STAT.COR.SA_PUB_Q < thresh), IS, G360, cmap.PiYG, '', [-0.5 0.5], 'square');
exportgraphics(h, fullfile(figDir, 'pub_sa_fdr.png'), 'Resolution', 300); close(h);

h = HyoSurfacePlotFunction((HCP.STAT.COR.CT_B - ABCD.STAT.COR.CT_B) .* (ABCD.STAT.COR.CT_PUB_Q < thresh), IS, G360, cmap.PiYG, '', [-0.5 0.5], 'square');
exportgraphics(h, fullfile(figDir, 'pub_ct_fdr.png'), 'Resolution', 300); close(h);

% Display proportion of cortical sheet with significant pubertal shift
[sum(ABCD.STAT.COR.CV_PUB_Q < thresh) sum(ABCD.STAT.COR.SA_PUB_Q < thresh) sum(ABCD.STAT.COR.CT_PUB_Q < thresh)] / 360

%% Visualize subcortical maps

% Compute FDR significance masks
indSTR_ABCD = ABCD.STAT.MOR.TSA.STR_SEX_Q < thresh;
indGP_ABCD = ABCD.STAT.MOR.TSA.GP_SEX_Q < thresh;
indTHAL_ABCD = ABCD.STAT.MOR.TSA.THAL_SEX_Q < thresh;
indHIPP_ABCD = ABCD.STAT.MOR.TSA.HIPP_SEX_Q < thresh;
indAMYG_ABCD = ABCD.STAT.MOR.TSA.AMYG_SEX_Q < thresh;

% ABCD results (not thresholded)
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(ABCD.STAT.VOL.TBV.SEX_B, [-0.5 0.5], cmap.PiYG); exportgraphics(h, fullfile(figDir, 'abcd_vol.png'), 'Resolution', 300); close(h);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(ABCD.STAT.MOR.TSA.STR_SEX_B, [-0.5 0.5], cmap.PiYG); exportgraphics(h, fullfile(figDir, 'abcd_str.png'), 'Resolution', 300); close(h);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(ABCD.STAT.MOR.TSA.GP_SEX_B, [-0.5 0.5], cmap.PiYG); exportgraphics(h, fullfile(figDir, 'abcd_gp.png'), 'Resolution', 300); close(h);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(ABCD.STAT.MOR.TSA.THAL_SEX_B, [-0.5 0.5], cmap.PiYG); exportgraphics(h, fullfile(figDir, 'abcd_thal.png'), 'Resolution', 300); close(h);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(ABCD.STAT.MOR.TSA.HIPP_SEX_B, [-0.5 0.5], cmap.PiYG); exportgraphics(h, fullfile(figDir, 'abcd_hipp.png'), 'Resolution', 300); close(h);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(ABCD.STAT.MOR.TSA.AMYG_SEX_B, [-0.5 0.5], cmap.PiYG); exportgraphics(h, fullfile(figDir, 'abcd_amyg.png'), 'Resolution', 300); close(h);

% ABCD FDR thresholded with gray background
temp = ABCD.STAT.MOR.TSA.STR_SEX_B .* (ABCD.STAT.MOR.TSA.STR_SEX_Q < thresh); h = figure('WindowState', 'Maximized'); HyoMagetViewNL(temp, [-0.5 0.5], ThreshColorMap(temp, cmap.PiYG, 0.915)); exportgraphics(h, fullfile(figDir, 'ABCD_STR_TSA_FDR.png'), 'Resolution', 300); close(h);
temp = ABCD.STAT.MOR.TSA.GP_SEX_B .* (ABCD.STAT.MOR.TSA.GP_SEX_Q < thresh); h = figure('WindowState', 'Maximized'); HyoMagetViewNL(temp, [-0.5 0.5], ThreshColorMap(temp, cmap.PiYG, 0.915)); exportgraphics(h, fullfile(figDir, 'ABCD_GP_TSA_FDR.png'), 'Resolution', 300); close(h);
temp = ABCD.STAT.MOR.TSA.THAL_SEX_B .* (ABCD.STAT.MOR.TSA.THAL_SEX_Q < thresh); h = figure('WindowState', 'Maximized'); HyoMagetViewNL(temp, [-0.5 0.5], ThreshColorMap(temp, cmap.PiYG, 0.915)); exportgraphics(h, fullfile(figDir, 'ABCD_THAL_TSA_FDR.png'), 'Resolution', 300); close(h);
temp = ABCD.STAT.MOR.TSA.HIPP_SEX_B .* (ABCD.STAT.MOR.TSA.HIPP_SEX_Q < thresh); h = figure('WindowState', 'Maximized'); HyoMagetViewNL(temp, [-0.5 0.5], ThreshColorMap(temp, cmap.PiYG, 0.915)); exportgraphics(h, fullfile(figDir, 'ABCD_HIPP_TSA_FDR.png'), 'Resolution', 300); close(h);
temp = ABCD.STAT.MOR.TSA.AMYG_SEX_B .* (ABCD.STAT.MOR.TSA.AMYG_SEX_Q < thresh); h = figure('WindowState', 'Maximized'); HyoMagetViewNL(temp, [-0.5 0.5], ThreshColorMap(temp, cmap.PiYG, 0.915)); exportgraphics(h, fullfile(figDir, 'ABCD_AMYG_TSA_FDR.png'), 'Resolution', 300); close(h);

% Display the areal proportion of significance
sum(S_STR_NA(ABCD.STAT.MOR.TSA.STR_SEX_Q < thresh))
sum(S_GP_NA(ABCD.STAT.MOR.TSA.GP_SEX_Q < thresh))
sum(S_THAL_NA(ABCD.STAT.MOR.TSA.THAL_SEX_Q < thresh))
sum(S_HIPP_NA(ABCD.STAT.MOR.TSA.HIPP_SEX_Q < thresh))
sum(S_AMYG_NA(ABCD.STAT.MOR.TSA.AMYG_SEX_Q < thresh))

% Calculate the colorscale of thresholded maps
temp = ABCD.STAT.MOR.TSA.STR_SEX_B .* (ABCD.STAT.MOR.TSA.STR_SEX_Q < thresh);
[min(abs(temp(temp<0))) max(abs(temp(temp<0))) min(abs(temp(temp>0))) max(abs(temp(temp>0)))]
temp = ABCD.STAT.MOR.TSA.GP_SEX_B .* (ABCD.STAT.MOR.TSA.GP_SEX_Q < thresh);
[min(abs(temp(temp<0))) max(abs(temp(temp<0))) min(abs(temp(temp>0))) max(abs(temp(temp>0)))]
temp = ABCD.STAT.MOR.TSA.THAL_SEX_B .* (ABCD.STAT.MOR.TSA.THAL_SEX_Q < thresh);
[min(abs(temp(temp<0))) max(abs(temp(temp<0))) min(abs(temp(temp>0))) max(abs(temp(temp>0)))]
temp = ABCD.STAT.MOR.TSA.HIPP_SEX_B .* (ABCD.STAT.MOR.TSA.HIPP_SEX_Q < thresh);
[min(abs(temp(temp<0))) max(abs(temp(temp<0))) min(abs(temp(temp>0))) max(abs(temp(temp>0)))]
temp = ABCD.STAT.MOR.TSA.AMYG_SEX_B .* (ABCD.STAT.MOR.TSA.AMYG_SEX_Q < thresh);
[min(abs(temp(temp<0))) max(abs(temp(temp<0))) min(abs(temp(temp>0))) max(abs(temp(temp>0)))]

% Visualize pubertal shift statistical maps (thresholded)
h = figure('WindowState', 'Maximized'); HyoMagetViewNL((HCP.STAT.VOL.TBV.SEX_B - ABCD.STAT.VOL.TBV.SEX_B) .* (ABCD.STAT.VOL.TBV.PUB_Q < thresh), [-0.5 0.5], cmap.PiYG); exportgraphics(h, fullfile(figDir, 'pub_vol_fdr.png'), 'Resolution', 300); close(h);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL((HCP.STAT.MOR.TSA.STR_SEX_B - ABCD.STAT.MOR.TSA.STR_SEX_B) .* (ABCD.STAT.MOR.TSA.STR_PUB_Q < thresh), [-0.5 0.5], cmap.PiYG); exportgraphics(h, fullfile(figDir, 'pub_str_fdr.png'), 'Resolution', 300); close(h);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL((HCP.STAT.MOR.TSA.GP_SEX_B - ABCD.STAT.MOR.TSA.GP_SEX_B) .* (ABCD.STAT.MOR.TSA.GP_PUB_Q < thresh), [-0.5 0.5], cmap.PiYG); exportgraphics(h, fullfile(figDir, 'pub_gp_fdr.png'), 'Resolution', 300); close(h);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL((HCP.STAT.MOR.TSA.THAL_SEX_B - ABCD.STAT.MOR.TSA.THAL_SEX_B) .* (ABCD.STAT.MOR.TSA.THAL_PUB_Q < thresh), [-0.5 0.5], cmap.PiYG); exportgraphics(h, fullfile(figDir, 'pub_thal_fdr.png'), 'Resolution', 300); close(h);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL((HCP.STAT.MOR.TSA.HIPP_SEX_B - ABCD.STAT.MOR.TSA.HIPP_SEX_B) .* (ABCD.STAT.MOR.TSA.HIPP_PUB_Q < thresh), [-0.5 0.5], cmap.PiYG); exportgraphics(h, fullfile(figDir, 'pub_hipp_fdr.png'), 'Resolution', 300); close(h);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL((HCP.STAT.MOR.TSA.AMYG_SEX_B - ABCD.STAT.MOR.TSA.AMYG_SEX_B) .* (ABCD.STAT.MOR.TSA.AMYG_PUB_Q < thresh), [-0.5 0.5], cmap.PiYG); exportgraphics(h, fullfile(figDir, 'pub_amyg_fdr.png'), 'Resolution', 300); close(h);

disp('ABCD Proportion of significantly sex-biased vertices (%)')
temp = [sum(S_STR_NA(indSTR_ABCD)) sum(S_GP_NA(indGP_ABCD)) sum(S_THAL_NA(indTHAL_ABCD)) sum(S_HIPP_NA(indHIPP_ABCD)) sum(S_AMYG_NA(indAMYG_ABCD))] * 100

disp('Proportion of significant differences in effect size (%)')
100*[sum(S_STR_NA(ABCD.STAT.MOR.TSA.STR_PUB_Q < thresh)) ...
sum(S_GP_NA(ABCD.STAT.MOR.TSA.GP_PUB_Q < thresh)) ...
sum(S_THAL_NA(ABCD.STAT.MOR.TSA.THAL_PUB_Q < thresh)) ...
sum(S_HIPP_NA(ABCD.STAT.MOR.TSA.HIPP_PUB_Q < thresh)) ...
sum(S_AMYG_NA(ABCD.STAT.MOR.TSA.AMYG_PUB_Q < thresh))]

%% Rotating full view

% Collate into one vector
SUB_B = [ABCD.STAT.MOR.TSA.STR_SEX_B ABCD.STAT.MOR.TSA.GP_SEX_B ABCD.STAT.MOR.TSA.THAL_SEX_B ABCD.STAT.MOR.TSA.HIPP_SEX_B ABCD.STAT.MOR.TSA.AMYG_SEX_B];
SUB_Q = [ABCD.STAT.MOR.TSA.STR_SEX_Q ABCD.STAT.MOR.TSA.GP_SEX_Q ABCD.STAT.MOR.TSA.THAL_SEX_Q ABCD.STAT.MOR.TSA.HIPP_SEX_Q ABCD.STAT.MOR.TSA.AMYG_SEX_Q];

% Static
c = 1;
dA = 45;
for i =  0:dA:360-dA
    disp(['Rotation angle ' num2str(i) ' deg']);
    h = figure('WindowState', 'Maximized');
    HyoMagetViewSpin(SUB_B .* (SUB_Q < thresh), [-0.5 0.5], i, ThreshColorMap(SUB_B .* (SUB_Q < thresh), cmap.PiYG, [-0.5 0.5], 0.915), 'static');
    print(h, fullfile(figDir, 'static_abcd', ['ABCD_TSA_ANG-' num2str(c) '.png']), '-dpng', '-r300');
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
    print(h, fullfile(figDir, 'static_abcd', ['NOTHR_ABCD_TSA_ANG-' num2str(c) '.png']), '-dpng', '-r300');
    close(h);
    c = c + 1;
end

%% Assess spatial correlations between HCP and ABCD (subcortex)

% Compute spatial correlations
R_STR = corr(HCP.STAT.MOR.TSA.STR_SEX_B', ABCD.STAT.MOR.TSA.STR_SEX_B');
R_GP = corr(HCP.STAT.MOR.TSA.GP_SEX_B', ABCD.STAT.MOR.TSA.GP_SEX_B');
R_THAL = corr(HCP.STAT.MOR.TSA.THAL_SEX_B', ABCD.STAT.MOR.TSA.THAL_SEX_B');
R_HIPP = corr(HCP.STAT.MOR.TSA.HIPP_SEX_B', ABCD.STAT.MOR.TSA.HIPP_SEX_B');
R_AMYG = corr(HCP.STAT.MOR.TSA.AMYG_SEX_B', ABCD.STAT.MOR.TSA.AMYG_SEX_B');

% Null distributions for spatial correlations
for i = 1:nPerm
    PERM_R_STR(i) = corr(HCP.STAT.MOR.MSR.STR_SEX_B(:,i), ABCD.STAT.MOR.TSA.STR_SEX_B');
    PERM_R_GP(i) = corr(HCP.STAT.MOR.MSR.GP_SEX_B(:,i), ABCD.STAT.MOR.TSA.GP_SEX_B');
    PERM_R_THAL(i) = corr(HCP.STAT.MOR.MSR.THAL_SEX_B(:,i), ABCD.STAT.MOR.TSA.THAL_SEX_B');
    PERM_R_HIPP(i) = corr(HCP.STAT.MOR.MSR.HIPP_SEX_B(:,i), ABCD.STAT.MOR.TSA.HIPP_SEX_B');
    PERM_R_AMYG(i) = corr(HCP.STAT.MOR.MSR.AMYG_SEX_B(:,i), ABCD.STAT.MOR.TSA.AMYG_SEX_B');
end

disp('Statistics for spatial correlations')
[R_STR R_GP R_THAL R_HIPP R_AMYG]
5*[sum(R_STR < PERM_R_STR) / nPerm sum(R_GP < PERM_R_GP) / nPerm sum(R_THAL < PERM_R_THAL) / nPerm sum(R_HIPP < PERM_R_HIPP) / nPerm sum(R_AMYG < PERM_R_AMYG) / nPerm]

% Test the significance of slope - "1"
H0 = [0 1];  c = 1;
mdl_STR = fitlm(ABCD.STAT.MOR.TSA.STR_SEX_B, HCP.STAT.MOR.TSA.STR_SEX_B);
mdl_GP = fitlm(ABCD.STAT.MOR.TSA.GP_SEX_B, HCP.STAT.MOR.TSA.GP_SEX_B);
mdl_THAL = fitlm(ABCD.STAT.MOR.TSA.THAL_SEX_B, HCP.STAT.MOR.TSA.THAL_SEX_B);
mdl_HIPP = fitlm(ABCD.STAT.MOR.TSA.HIPP_SEX_B, HCP.STAT.MOR.TSA.HIPP_SEX_B);
mdl_AMYG = fitlm(ABCD.STAT.MOR.TSA.AMYG_SEX_B, HCP.STAT.MOR.TSA.AMYG_SEX_B);

disp('Statistics for slope vs 1')
[mdl_STR.Coefficients.Estimate(2) mdl_GP.Coefficients.Estimate(2) mdl_THAL.Coefficients.Estimate(2) mdl_HIPP.Coefficients.Estimate(2) mdl_AMYG.Coefficients.Estimate(2)]
5*[coefTest(mdl_STR, H0, c) coefTest(mdl_GP, H0, c) coefTest(mdl_THAL, H0, c) coefTest(mdl_HIPP, H0, c) coefTest(mdl_AMYG, H0, c)]

% Generate scatter plots
h = figure('WindowState', 'maximized'); hold on;
hexscatter(ABCD.STAT.MOR.TSA.STR_SEX_B, HCP.STAT.MOR.TSA.STR_SEX_B, 'res', 50, 'xlim', [-1 1], 'ylim', [-1 1]);
scatter(ABCD.STAT.MOR.TSA.STR_SEX_B, HCP.STAT.MOR.TSA.STR_SEX_B, 5, [0 0.5 0.5], 'filled', 'MarkerEdgeColor', [0 0.5 0.5], 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0);
xlim([-1 1]); xticks(-1:0.25:1); ylim([-1 1]); yticks(-1:0.25:1); set(gca,'TickDir','out');
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--'; b.LineWidth = 2;
x = xline(0); x.LineStyle = '--'; x.LineWidth = 2; y = yline(0); y.LineStyle = '--'; y.LineWidth = 2;
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xlim([-1 1]); ylim([-1 1]);
axis square;
exportgraphics(h, fullfile(figDir, 'spatialcorr_abcd_str.png'), 'Resolution', 300);

h = figure('WindowState', 'maximized'); hold on;
hexscatter(ABCD.STAT.MOR.TSA.GP_SEX_B, HCP.STAT.MOR.TSA.GP_SEX_B, 'res', 50, 'xlim', [-1 1], 'ylim', [-1 1]);
scatter(ABCD.STAT.MOR.TSA.GP_SEX_B, HCP.STAT.MOR.TSA.GP_SEX_B, 5, [0 0.5 0.5], 'filled', 'MarkerEdgeColor', [0 0.5 0.5], 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0);
xlim([-1 1]); xticks(-1:0.25:1); ylim([-1 1]); yticks(-1:0.25:1); set(gca,'TickDir','out');
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--'; b.LineWidth = 2;
x = xline(0); x.LineStyle = '--'; x.LineWidth = 2; y = yline(0); y.LineStyle = '--'; y.LineWidth = 2;
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xlim([-1 1]); ylim([-1 1]);
axis square;
exportgraphics(h, fullfile(figDir, 'spatialcorr_abcd_gp.png'), 'Resolution', 300);

h = figure('WindowState', 'maximized'); hold on;
hexscatter(ABCD.STAT.MOR.TSA.THAL_SEX_B, HCP.STAT.MOR.TSA.THAL_SEX_B, 'res', 50, 'xlim', [-1 1], 'ylim', [-1 1]);
scatter(ABCD.STAT.MOR.TSA.THAL_SEX_B, HCP.STAT.MOR.TSA.THAL_SEX_B, 5, [0 0.5 0.5], 'filled', 'MarkerEdgeColor', [0 0.5 0.5], 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0);
xlim([-1 1]); xticks(-1:0.25:1); ylim([-1 1]); yticks(-1:0.25:1); set(gca,'TickDir','out');
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--'; b.LineWidth = 2;
x = xline(0); x.LineStyle = '--'; x.LineWidth = 2; y = yline(0); y.LineStyle = '--'; y.LineWidth = 2;
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xlim([-1 1]); ylim([-1 1]);
axis square;
exportgraphics(h, fullfile(figDir, 'spatialcorr_abcd_thal.png'), 'Resolution', 300);

h = figure('WindowState', 'maximized'); hold on;
hexscatter(ABCD.STAT.MOR.TSA.HIPP_SEX_B, HCP.STAT.MOR.TSA.HIPP_SEX_B, 'res', 50, 'xlim', [-1 1], 'ylim', [-1 1]);
scatter(ABCD.STAT.MOR.TSA.HIPP_SEX_B, HCP.STAT.MOR.TSA.HIPP_SEX_B, 5, [0 0.5 0.5], 'filled', 'MarkerEdgeColor', [0 0.5 0.5], 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0);
xlim([-1 1]); xticks(-1:0.25:1); ylim([-1 1]); yticks(-1:0.25:1); set(gca,'TickDir','out');
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--'; b.LineWidth = 2;
x = xline(0); x.LineStyle = '--'; x.LineWidth = 2; y = yline(0); y.LineStyle = '--'; y.LineWidth = 2;
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xlim([-1 1]); ylim([-1 1]);
axis square;
exportgraphics(h, fullfile(figDir, 'spatialcorr_abcd_hipp.png'), 'Resolution', 300);

h = figure('WindowState', 'maximized'); hold on;
hexscatter(ABCD.STAT.MOR.TSA.AMYG_SEX_B, HCP.STAT.MOR.TSA.AMYG_SEX_B, 'res', 50, 'xlim', [-1 1], 'ylim', [-1 1]);
scatter(ABCD.STAT.MOR.TSA.AMYG_SEX_B, HCP.STAT.MOR.TSA.AMYG_SEX_B, 5, [0 0.5 0.5], 'filled', 'MarkerEdgeColor', [0 0.5 0.5], 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0);
xlim([-1 1]); xticks(-1:0.25:1); ylim([-1 1]); yticks(-1:0.25:1); set(gca,'TickDir','out');
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--'; b.LineWidth = 2;
x = xline(0); x.LineStyle = '--'; x.LineWidth = 2; y = yline(0); y.LineStyle = '--'; y.LineWidth = 2;
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xlim([-1 1]); ylim([-1 1]);
axis square;
exportgraphics(h, fullfile(figDir, 'spatialcorr_abcd_amyg.png'), 'Resolution', 300);

close all;

%% Assess spatial correlations between HCP and ABCD (cortex)

% Compute spatial correlations
R_CV = corr(HCP.STAT.COR.CV_B', ABCD.STAT.COR.CV_B');
R_SA = corr(HCP.STAT.COR.SA_B', ABCD.STAT.COR.SA_B');
R_CT = corr(HCP.STAT.COR.CT_B', ABCD.STAT.COR.CT_B');

% Null distributions for spatial correlations
for i = 1:nPerm
    PERM_R_CV(i) = corr(HCP.STAT.COR.CV_B(Vasa(:,i))', ABCD.STAT.COR.CV_B');
    PERM_R_SA(i) = corr(HCP.STAT.COR.SA_B(Vasa(:,i))', ABCD.STAT.COR.SA_B');
    PERM_R_CT(i) = corr(HCP.STAT.COR.CT_B(Vasa(:,i))', ABCD.STAT.COR.CT_B');
end

disp('Statistics for spatial correlations')
[R_CV R_SA R_CT]
[sum(R_CV < PERM_R_CV) sum(R_SA < PERM_R_SA) sum(R_CT < PERM_R_CT)] ./ nPerm

% Test the significance of slope - "1"
H0 = [0 1];  c = 1;
mdl_CV = fitlm(ABCD.STAT.COR.CV_B, HCP.STAT.COR.CV_B);
mdl_SA = fitlm(ABCD.STAT.COR.SA_B, HCP.STAT.COR.SA_B);
mdl_CT = fitlm(ABCD.STAT.COR.CT_B, HCP.STAT.COR.CT_B);

disp('Statistics for slope vs 1')
[mdl_CV.Coefficients.Estimate(2) mdl_SA.Coefficients.Estimate(2) mdl_CT.Coefficients.Estimate(2)]
[coefTest(mdl_CV, H0, c) coefTest(mdl_SA, H0, c) coefTest(mdl_CT, H0, c)]

% Generate scatter plots
h = figure('WindowState', 'maximized'); hold on;
hexscatter(ABCD.STAT.COR.CV_B, HCP.STAT.COR.CV_B, 'res', 50, 'xlim', [-1 1], 'ylim', [-1 1]);
scatter(ABCD.STAT.COR.CV_B, HCP.STAT.COR.CV_B, 5, [0 0.5 0.5], 'filled', 'MarkerEdgeColor', [0 0.5 0.5], 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0);
xlim([-1 1]); xticks(-1:0.25:1); ylim([-1 1]); yticks(-1:0.25:1); set(gca,'TickDir','out');
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--'; b.LineWidth = 2;
x = xline(0); x.LineStyle = '--'; x.LineWidth = 2; y = yline(0); y.LineStyle = '--'; y.LineWidth = 2;
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xlim([-1 1]); ylim([-1 1]);
axis square;
exportgraphics(h, fullfile(figDir, 'spatialcorr_abcd_CV.png'), 'Resolution', 300);
close(h);

h = figure('WindowState', 'maximized'); hold on;
hexscatter(ABCD.STAT.COR.SA_B, HCP.STAT.COR.SA_B, 'res', 50, 'xlim', [-1 1], 'ylim', [-1 1]);
scatter(ABCD.STAT.COR.SA_B, HCP.STAT.COR.SA_B, 5, [0 0.5 0.5], 'filled', 'MarkerEdgeColor', [0 0.5 0.5], 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0);
xlim([-1 1]); xticks(-1:0.25:1); ylim([-1 1]); yticks(-1:0.25:1); set(gca,'TickDir','out');
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--'; b.LineWidth = 2;
x = xline(0); x.LineStyle = '--'; x.LineWidth = 2; y = yline(0); y.LineStyle = '--'; y.LineWidth = 2;
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xlim([-1 1]); ylim([-1 1]);
axis square;
exportgraphics(h, fullfile(figDir, 'spatialcorr_abcd_SA.png'), 'Resolution', 300);
close(h);

h = figure('WindowState', 'maximized'); hold on;
hexscatter(ABCD.STAT.COR.CT_B, HCP.STAT.COR.CT_B, 'res', 50, 'xlim', [-1 1], 'ylim', [-1 1]);
scatter(ABCD.STAT.COR.CT_B, HCP.STAT.COR.CT_B, 5, [0 0.5 0.5], 'filled', 'MarkerEdgeColor', [0 0.5 0.5], 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0);
xlim([-1 1]); xticks(-1:0.25:1); ylim([-1 1]); yticks(-1:0.25:1); set(gca,'TickDir','out');
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--'; b.LineWidth = 2;
x = xline(0); x.LineStyle = '--'; x.LineWidth = 2; y = yline(0); y.LineStyle = '--'; y.LineWidth = 2;
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xlim([-1 1]); ylim([-1 1]);
axis square;
exportgraphics(h, fullfile(figDir, 'spatialcorr_abcd_CT.png'), 'Resolution', 300);
close(h);

%% ABCD & HCP Conjunction

% STR
ABCD_STR_AD = zeros(1,26000);
HCP_STR_AD = zeros(1,26000);
STR_AD = zeros(1,26000);

ABCD_STR_AD(ABCD.STAT.MOR.TSA.STR_SEX_B > 0) =  1;
ABCD_STR_AD(ABCD.STAT.MOR.TSA.STR_SEX_B < 0) = -1;
ABCD_STR_AD((ABCD.STAT.MOR.TSA.STR_SEX_Q < thresh) & (ABCD.STAT.MOR.TSA.STR_SEX_B > 0))  =  2;
ABCD_STR_AD((ABCD.STAT.MOR.TSA.STR_SEX_Q < thresh) & (ABCD.STAT.MOR.TSA.STR_SEX_B < 0))  = -2;

HCP_STR_AD(HCP.STAT.MOR.TSA.STR_SEX_B > 0) =  1;
HCP_STR_AD(HCP.STAT.MOR.TSA.STR_SEX_B < 0) = -1;
HCP_STR_AD((HCP.STAT.MOR.TSA.STR_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.STR_SEX_B > 0))  =  2;
HCP_STR_AD((HCP.STAT.MOR.TSA.STR_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.STR_SEX_B < 0))  = -2;

% Established: FF, MM & stable (dark green or magenta)
STR_AD((ABCD_STR_AD == 2) & (HCP_STR_AD == 2)) = -1.9;
STR_AD((ABCD_STR_AD == -2) & (HCP_STR_AD == -2)) = 1.1;

% Emerging: fF, mM, mF, fM (medium green or magenta)
STR_AD((ABCD_STR_AD == 1) & (HCP_STR_AD == 2)) = -1.5;
STR_AD((ABCD_STR_AD == -1) & (HCP_STR_AD == -2)) = 1.5;
STR_AD((ABCD_STR_AD == -1) & (HCP_STR_AD == 2)) = -1.5;
STR_AD((ABCD_STR_AD == 1) & (HCP_STR_AD == -2)) = 1.5;

% Diminishing flip: Ff, Fm, Mm, Mf (light gray)
STR_AD((ABCD_STR_AD == 2) & (HCP_STR_AD == 1)) = 0.5;
STR_AD((ABCD_STR_AD == 2) & (HCP_STR_AD == -1)) = 0.5;
STR_AD((ABCD_STR_AD == -2) & (HCP_STR_AD == -1)) = 0.5;
STR_AD((ABCD_STR_AD == -2) & (HCP_STR_AD == 1)) = 0.5;

% Hard flip: MF, FM (dark gray)
STR_AD((ABCD_STR_AD == -2) & (HCP_STR_AD == 2)) = 0.1;
STR_AD((ABCD_STR_AD == 2) & (HCP_STR_AD == -2)) = 0.1;

% Insignificant sex-biases: ff, mf, mm, fm and insignificant delta (white)
STR_AD((ABCD_STR_AD == 1) & (HCP_STR_AD == 1)) = 0.9;
STR_AD((ABCD_STR_AD == -1) & (HCP_STR_AD == 1)) = 0.9;
STR_AD((ABCD_STR_AD == -1) & (HCP_STR_AD == -1)) = 0.9;
STR_AD((ABCD_STR_AD == 1) & (HCP_STR_AD == -1)) = 0.9;

h = figure('WindowState', 'Maximized'); HyoMagetViewFlat(STR_AD, [-2 2], cmap.tri_bi);
% exportgraphics(h, fullfile(figDir, 'abcd_hcp_overlap_STR.png'), 'Resolution', 300);

STR_AD_Count = [sum(S_STR_NA(STR_AD == 1.1)) sum(S_STR_NA(STR_AD == 1.5)) 0 ...
                sum(S_STR_NA(STR_AD == -1.9)) sum(S_STR_NA(STR_AD == -1.5)) 0 ...
                sum(S_STR_NA(STR_AD == 0.1)) sum(S_STR_NA(STR_AD == 0.5)) sum(S_STR_NA(STR_AD == 0.9))];

h = figure;
b = bar(1, STR_AD_Count, 'stacked', 'EdgeColor', 'none');
for i = 1:9
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_tri_bi(i,:);
end
box off
axis off
% exportgraphics(h, fullfile(figDir, 'abcd_hcp_overlap_bar_STR.png'), 'Resolution', 300);

% GP
ABCD_GP_AD = zeros(1,6000);
HCP_GP_AD = zeros(1,6000);
GP_AD = zeros(1,6000);

ABCD_GP_AD(ABCD.STAT.MOR.TSA.GP_SEX_B > 0) =  1;
ABCD_GP_AD(ABCD.STAT.MOR.TSA.GP_SEX_B < 0) = -1;
ABCD_GP_AD((ABCD.STAT.MOR.TSA.GP_SEX_Q < thresh) & (ABCD.STAT.MOR.TSA.GP_SEX_B > 0))  =  2;
ABCD_GP_AD((ABCD.STAT.MOR.TSA.GP_SEX_Q < thresh) & (ABCD.STAT.MOR.TSA.GP_SEX_B < 0))  = -2;

HCP_GP_AD(HCP.STAT.MOR.TSA.GP_SEX_B > 0) =  1;
HCP_GP_AD(HCP.STAT.MOR.TSA.GP_SEX_B < 0) = -1;
HCP_GP_AD((HCP.STAT.MOR.TSA.GP_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.GP_SEX_B > 0))  =  2;
HCP_GP_AD((HCP.STAT.MOR.TSA.GP_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.GP_SEX_B < 0))  = -2;

% Established: FF, MM & stable (dark green or magenta)
GP_AD((ABCD_GP_AD == 2) & (HCP_GP_AD == 2)) = -1.9;
GP_AD((ABCD_GP_AD == -2) & (HCP_GP_AD == -2)) = 1.1;

% Emerging: fF, mM, mF, fM (medium green or magenta)
GP_AD((ABCD_GP_AD == 1) & (HCP_GP_AD == 2)) = -1.5;
GP_AD((ABCD_GP_AD == -1) & (HCP_GP_AD == -2)) = 1.5;
GP_AD((ABCD_GP_AD == -1) & (HCP_GP_AD == 2)) = -1.5;
GP_AD((ABCD_GP_AD == 1) & (HCP_GP_AD == -2)) = 1.5;

% Diminishing flip: Ff, Fm, Mm, Mf (light gray)
GP_AD((ABCD_GP_AD == 2) & (HCP_GP_AD == 1)) = 0.5;
GP_AD((ABCD_GP_AD == 2) & (HCP_GP_AD == -1)) = 0.5;
GP_AD((ABCD_GP_AD == -2) & (HCP_GP_AD == -1)) = 0.5;
GP_AD((ABCD_GP_AD == -2) & (HCP_GP_AD == 1)) = 0.5;

% Hard flip: MF, FM (dark gray)
GP_AD((ABCD_GP_AD == -2) & (HCP_GP_AD == 2)) = 0.1;
GP_AD((ABCD_GP_AD == 2) & (HCP_GP_AD == -2)) = 0.1;

% Insignificant sex-biases: ff, mf, mm, fm and insignificant delta (white)
GP_AD((ABCD_GP_AD == 1) & (HCP_GP_AD == 1)) = 0.9;
GP_AD((ABCD_GP_AD == -1) & (HCP_GP_AD == 1)) = 0.9;
GP_AD((ABCD_GP_AD == -1) & (HCP_GP_AD == -1)) = 0.9;
GP_AD((ABCD_GP_AD == 1) & (HCP_GP_AD == -1)) = 0.9;

h = figure('WindowState', 'Maximized'); HyoMagetViewFlat(GP_AD, [-2 2], cmap.tri_bi);
% exportgraphics(h, fullfile(figDir, 'abcd_hcp_overlap_GP.png'), 'Resolution', 300);

GP_AD_Count = [sum(S_GP_NA(GP_AD == 1.1)) sum(S_GP_NA(GP_AD == 1.5)) 0 ...
                sum(S_GP_NA(GP_AD == -1.9)) sum(S_GP_NA(GP_AD == -1.5)) 0 ...
                sum(S_GP_NA(GP_AD == 0.1)) sum(S_GP_NA(GP_AD == 0.5)) sum(S_GP_NA(GP_AD == 0.9))];

h = figure;
b = bar(1, GP_AD_Count, 'stacked', 'EdgeColor', 'none');
for i = 1:9
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_tri_bi(i,:);
end
box off
axis off
% exportgraphics(h, fullfile(figDir, 'abcd_hcp_overlap_bar_GP.png'), 'Resolution', 300);

% THAL
ABCD_THAL_AD = zeros(1,13000);
HCP_THAL_AD = zeros(1,13000);
THAL_AD = zeros(1,13000);

ABCD_THAL_AD(ABCD.STAT.MOR.TSA.THAL_SEX_B > 0) =  1;
ABCD_THAL_AD(ABCD.STAT.MOR.TSA.THAL_SEX_B < 0) = -1;
ABCD_THAL_AD((ABCD.STAT.MOR.TSA.THAL_SEX_Q < thresh) & (ABCD.STAT.MOR.TSA.THAL_SEX_B > 0))  =  2;
ABCD_THAL_AD((ABCD.STAT.MOR.TSA.THAL_SEX_Q < thresh) & (ABCD.STAT.MOR.TSA.THAL_SEX_B < 0))  = -2;

HCP_THAL_AD(HCP.STAT.MOR.TSA.THAL_SEX_B > 0) =  1;
HCP_THAL_AD(HCP.STAT.MOR.TSA.THAL_SEX_B < 0) = -1;
HCP_THAL_AD((HCP.STAT.MOR.TSA.THAL_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.THAL_SEX_B > 0))  =  2;
HCP_THAL_AD((HCP.STAT.MOR.TSA.THAL_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.THAL_SEX_B < 0))  = -2;

% Established: FF, MM & stable (dark green or magenta)
THAL_AD((ABCD_THAL_AD == 2) & (HCP_THAL_AD == 2)) = -1.9;
THAL_AD((ABCD_THAL_AD == -2) & (HCP_THAL_AD == -2)) = 1.1;

% Emerging: fF, mM, mF, fM (medium green or magenta)
THAL_AD((ABCD_THAL_AD == 1) & (HCP_THAL_AD == 2)) = -1.5;
THAL_AD((ABCD_THAL_AD == -1) & (HCP_THAL_AD == -2)) = 1.5;
THAL_AD((ABCD_THAL_AD == -1) & (HCP_THAL_AD == 2)) = -1.5;
THAL_AD((ABCD_THAL_AD == 1) & (HCP_THAL_AD == -2)) = 1.5;

% Diminishing flip: Ff, Fm, Mm, Mf (light gray)
THAL_AD((ABCD_THAL_AD == 2) & (HCP_THAL_AD == 1)) = 0.5;
THAL_AD((ABCD_THAL_AD == 2) & (HCP_THAL_AD == -1)) = 0.5;
THAL_AD((ABCD_THAL_AD == -2) & (HCP_THAL_AD == -1)) = 0.5;
THAL_AD((ABCD_THAL_AD == -2) & (HCP_THAL_AD == 1)) = 0.5;

% Hard flip: MF, FM (dark gray)
THAL_AD((ABCD_THAL_AD == -2) & (HCP_THAL_AD == 2)) = 0.1;
THAL_AD((ABCD_THAL_AD == 2) & (HCP_THAL_AD == -2)) = 0.1;

% Insignificant sex-biases: ff, mf, mm, fm and insignificant delta (white)
THAL_AD((ABCD_THAL_AD == 1) & (HCP_THAL_AD == 1)) = 0.9;
THAL_AD((ABCD_THAL_AD == -1) & (HCP_THAL_AD == 1)) = 0.9;
THAL_AD((ABCD_THAL_AD == -1) & (HCP_THAL_AD == -1)) = 0.9;
THAL_AD((ABCD_THAL_AD == 1) & (HCP_THAL_AD == -1)) = 0.9;

h = figure('WindowState', 'Maximized'); HyoMagetViewFlat(THAL_AD, [-2 2], cmap.tri_bi);
% exportgraphics(h, fullfile(figDir, 'abcd_hcp_overlap_THAL.png'), 'Resolution', 300);

THAL_AD_Count = [sum(S_THAL_NA(THAL_AD == 1.1)) sum(S_THAL_NA(THAL_AD == 1.5)) 0 ...
                sum(S_THAL_NA(THAL_AD == -1.9)) sum(S_THAL_NA(THAL_AD == -1.5)) 0 ...
                sum(S_THAL_NA(THAL_AD == 0.1)) sum(S_THAL_NA(THAL_AD == 0.5)) sum(S_THAL_NA(THAL_AD == 0.9))];

h = figure;
b = bar(1, THAL_AD_Count, 'stacked', 'EdgeColor', 'none');
for i = 1:9
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_tri_bi(i,:);
end
box off
axis off
% exportgraphics(h, fullfile(figDir, 'abcd_hcp_overlap_bar_THAL.png'), 'Resolution', 300);

% HIPP
ABCD_HIPP_AD = zeros(1,2367);
HCP_HIPP_AD = zeros(1,2367);
HIPP_AD = zeros(1,2367);

ABCD_HIPP_AD(ABCD.STAT.MOR.TSA.HIPP_SEX_B > 0) =  1;
ABCD_HIPP_AD(ABCD.STAT.MOR.TSA.HIPP_SEX_B < 0) = -1;
ABCD_HIPP_AD((ABCD.STAT.MOR.TSA.HIPP_SEX_Q < thresh) & (ABCD.STAT.MOR.TSA.HIPP_SEX_B > 0))  =  2;
ABCD_HIPP_AD((ABCD.STAT.MOR.TSA.HIPP_SEX_Q < thresh) & (ABCD.STAT.MOR.TSA.HIPP_SEX_B < 0))  = -2;

HCP_HIPP_AD(HCP.STAT.MOR.TSA.HIPP_SEX_B > 0) =  1;
HCP_HIPP_AD(HCP.STAT.MOR.TSA.HIPP_SEX_B < 0) = -1;
HCP_HIPP_AD((HCP.STAT.MOR.TSA.HIPP_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.HIPP_SEX_B > 0))  =  2;
HCP_HIPP_AD((HCP.STAT.MOR.TSA.HIPP_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.HIPP_SEX_B < 0))  = -2;

% Established: FF, MM & stable (dark green or magenta)
HIPP_AD((ABCD_HIPP_AD == 2) & (HCP_HIPP_AD == 2)) = -1.9;
HIPP_AD((ABCD_HIPP_AD == -2) & (HCP_HIPP_AD == -2)) = 1.1;

% Emerging: fF, mM, mF, fM (medium green or magenta)
HIPP_AD((ABCD_HIPP_AD == 1) & (HCP_HIPP_AD == 2)) = -1.5;
HIPP_AD((ABCD_HIPP_AD == -1) & (HCP_HIPP_AD == -2)) = 1.5;
HIPP_AD((ABCD_HIPP_AD == -1) & (HCP_HIPP_AD == 2)) = -1.5;
HIPP_AD((ABCD_HIPP_AD == 1) & (HCP_HIPP_AD == -2)) = 1.5;

% Diminishing flip: Ff, Fm, Mm, Mf (light gray)
HIPP_AD((ABCD_HIPP_AD == 2) & (HCP_HIPP_AD == 1)) = 0.5;
HIPP_AD((ABCD_HIPP_AD == 2) & (HCP_HIPP_AD == -1)) = 0.5;
HIPP_AD((ABCD_HIPP_AD == -2) & (HCP_HIPP_AD == -1)) = 0.5;
HIPP_AD((ABCD_HIPP_AD == -2) & (HCP_HIPP_AD == 1)) = 0.5;

% Hard flip: MF, FM (dark gray)
HIPP_AD((ABCD_HIPP_AD == -2) & (HCP_HIPP_AD == 2)) = 0.1;
HIPP_AD((ABCD_HIPP_AD == 2) & (HCP_HIPP_AD == -2)) = 0.1;

% Insignificant sex-biases: ff, mf, mm, fm and insignificant delta (white)
HIPP_AD((ABCD_HIPP_AD == 1) & (HCP_HIPP_AD == 1)) = 0.9;
HIPP_AD((ABCD_HIPP_AD == -1) & (HCP_HIPP_AD == 1)) = 0.9;
HIPP_AD((ABCD_HIPP_AD == -1) & (HCP_HIPP_AD == -1)) = 0.9;
HIPP_AD((ABCD_HIPP_AD == 1) & (HCP_HIPP_AD == -1)) = 0.9;

h = figure('WindowState', 'Maximized'); HyoMagetViewFlat(HIPP_AD, [-2 2], cmap.tri_bi);
% exportgraphics(h, fullfile(figDir, 'abcd_hcp_overlap_HIPP.png'), 'Resolution', 300);

HIPP_AD_Count = [sum(S_HIPP_NA(HIPP_AD == 1.1)) sum(S_HIPP_NA(HIPP_AD == 1.5)) 0 ...
                sum(S_HIPP_NA(HIPP_AD == -1.9)) sum(S_HIPP_NA(HIPP_AD == -1.5)) 0 ...
                sum(S_HIPP_NA(HIPP_AD == 0.1)) sum(S_HIPP_NA(HIPP_AD == 0.5)) sum(S_HIPP_NA(HIPP_AD == 0.9))];

h = figure;
b = bar(1, HIPP_AD_Count, 'stacked', 'EdgeColor', 'none');
for i = 1:9
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_tri_bi(i,:);
end
box off
axis off
% exportgraphics(h, fullfile(figDir, 'abcd_hcp_overlap_bar_HIPP.png'), 'Resolution', 300);

% AMYG
ABCD_AMYG_AD = zeros(1,2878);
HCP_AMYG_AD = zeros(1,2878);
AMYG_AD = zeros(1,2878);

ABCD_AMYG_AD(ABCD.STAT.MOR.TSA.AMYG_SEX_B > 0) =  1;
ABCD_AMYG_AD(ABCD.STAT.MOR.TSA.AMYG_SEX_B < 0) = -1;
ABCD_AMYG_AD((ABCD.STAT.MOR.TSA.AMYG_SEX_Q < thresh) & (ABCD.STAT.MOR.TSA.AMYG_SEX_B > 0))  =  2;
ABCD_AMYG_AD((ABCD.STAT.MOR.TSA.AMYG_SEX_Q < thresh) & (ABCD.STAT.MOR.TSA.AMYG_SEX_B < 0))  = -2;

HCP_AMYG_AD(HCP.STAT.MOR.TSA.AMYG_SEX_B > 0) =  1;
HCP_AMYG_AD(HCP.STAT.MOR.TSA.AMYG_SEX_B < 0) = -1;
HCP_AMYG_AD((HCP.STAT.MOR.TSA.AMYG_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.AMYG_SEX_B > 0))  =  2;
HCP_AMYG_AD((HCP.STAT.MOR.TSA.AMYG_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.AMYG_SEX_B < 0))  = -2;

% Established: FF, MM & stable (dark green or magenta)
AMYG_AD((ABCD_AMYG_AD == 2) & (HCP_AMYG_AD == 2)) = -1.9;
AMYG_AD((ABCD_AMYG_AD == -2) & (HCP_AMYG_AD == -2)) = 1.1;

% Emerging: fF, mM, mF, fM (medium green or magenta)
AMYG_AD((ABCD_AMYG_AD == 1) & (HCP_AMYG_AD == 2)) = -1.5;
AMYG_AD((ABCD_AMYG_AD == -1) & (HCP_AMYG_AD == -2)) = 1.5;
AMYG_AD((ABCD_AMYG_AD == -1) & (HCP_AMYG_AD == 2)) = -1.5;
AMYG_AD((ABCD_AMYG_AD == 1) & (HCP_AMYG_AD == -2)) = 1.5;

% Diminishing flip: Ff, Fm, Mm, Mf (light gray)
AMYG_AD((ABCD_AMYG_AD == 2) & (HCP_AMYG_AD == 1)) = 0.5;
AMYG_AD((ABCD_AMYG_AD == 2) & (HCP_AMYG_AD == -1)) = 0.5;
AMYG_AD((ABCD_AMYG_AD == -2) & (HCP_AMYG_AD == -1)) = 0.5;
AMYG_AD((ABCD_AMYG_AD == -2) & (HCP_AMYG_AD == 1)) = 0.5;

% Hard flip: MF, FM (dark gray)
AMYG_AD((ABCD_AMYG_AD == -2) & (HCP_AMYG_AD == 2)) = 0.1;
AMYG_AD((ABCD_AMYG_AD == 2) & (HCP_AMYG_AD == -2)) = 0.1;

% Insignificant sex-biases: ff, mf, mm, fm and insignificant delta (white)
AMYG_AD((ABCD_AMYG_AD == 1) & (HCP_AMYG_AD == 1)) = 0.9;
AMYG_AD((ABCD_AMYG_AD == -1) & (HCP_AMYG_AD == 1)) = 0.9;
AMYG_AD((ABCD_AMYG_AD == -1) & (HCP_AMYG_AD == -1)) = 0.9;
AMYG_AD((ABCD_AMYG_AD == 1) & (HCP_AMYG_AD == -1)) = 0.9;

h = figure('WindowState', 'Maximized'); HyoMagetViewFlat(AMYG_AD, [-2 2], cmap.tri_bi);
% exportgraphics(h, fullfile(figDir, 'abcd_hcp_overlap_AMYG.png'), 'Resolution', 300);

AMYG_AD_Count = [sum(S_AMYG_NA(AMYG_AD == 1.1)) sum(S_AMYG_NA(AMYG_AD == 1.5)) 0 ...
                sum(S_AMYG_NA(AMYG_AD == -1.9)) sum(S_AMYG_NA(AMYG_AD == -1.5)) 0 ...
                sum(S_AMYG_NA(AMYG_AD == 0.1)) sum(S_AMYG_NA(AMYG_AD == 0.5)) sum(S_AMYG_NA(AMYG_AD == 0.9))];

h = figure;
b = bar(1, AMYG_AD_Count, 'stacked', 'EdgeColor', 'none');
for i = 1:9
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_tri_bi(i,:);
end
box off
axis off
% exportgraphics(h, fullfile(figDir, 'abcd_hcp_overlap_bar_AMYG.png'), 'Resolution', 300);

close all;

disp('Consistency')
[sum(STR_AD_Count([1:6 9])) / sum(STR_AD_Count)
sum(GP_AD_Count([1:6 9])) / sum(GP_AD_Count)
sum(THAL_AD_Count([1:6 9])) / sum(THAL_AD_Count)
sum(HIPP_AD_Count([1:6 9])) / sum(HIPP_AD_Count)
sum(AMYG_AD_Count([1:6 9])) / sum(AMYG_AD_Count)]

disp('Established')
[sum(STR_AD_Count([1 3])) / sum(STR_AD_Count)
sum(GP_AD_Count([1 3])) / sum(GP_AD_Count)
sum(THAL_AD_Count([1 3])) / sum(THAL_AD_Count)
sum(HIPP_AD_Count([1 3])) / sum(HIPP_AD_Count)
sum(AMYG_AD_Count([1 3])) / sum(AMYG_AD_Count)]

disp('Emerging')
[sum(STR_AD_Count([2 4])) / sum(STR_AD_Count)
sum(GP_AD_Count([2 4])) / sum(GP_AD_Count)
sum(THAL_AD_Count([2 4])) / sum(THAL_AD_Count)
sum(HIPP_AD_Count([2 4])) / sum(HIPP_AD_Count)
sum(AMYG_AD_Count([2 4])) / sum(AMYG_AD_Count)]

disp('Never sex-biased')
[sum(STR_AD_Count(9)) / sum(STR_AD_Count)
sum(GP_AD_Count(9)) / sum(GP_AD_Count)
sum(THAL_AD_Count(9)) / sum(THAL_AD_Count)
sum(HIPP_AD_Count(9)) / sum(HIPP_AD_Count)
sum(AMYG_AD_Count(9)) / sum(AMYG_AD_Count)]

%% Spin test for conjunction

% Discretize spatially permuted HCP sex difference maps
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
PERM_HCP_AMYG_AD(HCP.STAT.MOR.MSR.AMYG_SEX_Q & (HCP.STAT.MOR.MSR.AMYG_SEX_B < 0))  = -2;

% Compute agreement maps
PERM_STR_AD = zeros(26000, nPerm);
PERM_STR_AD(repmat(ABCD_STR_AD == 2, nPerm, 1)' & (PERM_HCP_STR_AD == 2)) = -1.9;
PERM_STR_AD(repmat(ABCD_STR_AD == -2, nPerm, 1)' & (PERM_HCP_STR_AD == -2)) = 1.1;
PERM_STR_AD(repmat(ABCD_STR_AD == 1, nPerm, 1)' & (PERM_HCP_STR_AD == 2)) = -1.5;
PERM_STR_AD(repmat(ABCD_STR_AD == -1, nPerm, 1)' & (PERM_HCP_STR_AD == 2)) = -1.5;
PERM_STR_AD(repmat(ABCD_STR_AD == 1, nPerm, 1)' & (PERM_HCP_STR_AD == -2)) = 1.5;
PERM_STR_AD(repmat(ABCD_STR_AD == -1, nPerm, 1)' & (PERM_HCP_STR_AD == -2)) = 1.5;
PERM_STR_AD(repmat(ABCD_STR_AD == -2, nPerm, 1)' & (PERM_HCP_STR_AD == 1)) = 0.5;
PERM_STR_AD(repmat(ABCD_STR_AD == 2, nPerm, 1)' & (PERM_HCP_STR_AD == 1)) = 0.5;
PERM_STR_AD(repmat(ABCD_STR_AD == -2, nPerm, 1)' & (PERM_HCP_STR_AD == -1)) = 0.5;
PERM_STR_AD(repmat(ABCD_STR_AD == 2, nPerm, 1)' & (PERM_HCP_STR_AD == -1)) = 0.5;
PERM_STR_AD(repmat(ABCD_STR_AD == -2, nPerm, 1)' & (PERM_HCP_STR_AD == 2)) = 0.1;
PERM_STR_AD(repmat(ABCD_STR_AD == 2, nPerm, 1)' & (PERM_HCP_STR_AD == -2)) = 0.1;
PERM_STR_AD(repmat(ABCD_STR_AD == 1, nPerm, 1)' & (PERM_HCP_STR_AD == 1)) = 0.9;
PERM_STR_AD(repmat(ABCD_STR_AD == -1, nPerm, 1)' & (PERM_HCP_STR_AD == 1)) = 0.9;
PERM_STR_AD(repmat(ABCD_STR_AD == -1, nPerm, 1)' & (PERM_HCP_STR_AD == -1)) = 0.9;
PERM_STR_AD(repmat(ABCD_STR_AD == 1, nPerm, 1)' & (PERM_HCP_STR_AD == -1)) = 0.9;

PERM_GP_AD = zeros(6000, nPerm);
PERM_GP_AD(repmat(ABCD_GP_AD == 2, nPerm, 1)' & (PERM_HCP_GP_AD == 2)) = -1.9;
PERM_GP_AD(repmat(ABCD_GP_AD == -2, nPerm, 1)' & (PERM_HCP_GP_AD == -2)) = 1.1;
PERM_GP_AD(repmat(ABCD_GP_AD == 1, nPerm, 1)' & (PERM_HCP_GP_AD == 2)) = -1.5;
PERM_GP_AD(repmat(ABCD_GP_AD == -1, nPerm, 1)' & (PERM_HCP_GP_AD == 2)) = -1.5;
PERM_GP_AD(repmat(ABCD_GP_AD == 1, nPerm, 1)' & (PERM_HCP_GP_AD == -2)) = 1.5;
PERM_GP_AD(repmat(ABCD_GP_AD == -1, nPerm, 1)' & (PERM_HCP_GP_AD == -2)) = 1.5;
PERM_GP_AD(repmat(ABCD_GP_AD == -2, nPerm, 1)' & (PERM_HCP_GP_AD == 1)) = 0.5;
PERM_GP_AD(repmat(ABCD_GP_AD == 2, nPerm, 1)' & (PERM_HCP_GP_AD == 1)) = 0.5;
PERM_GP_AD(repmat(ABCD_GP_AD == -2, nPerm, 1)' & (PERM_HCP_GP_AD == -1)) = 0.5;
PERM_GP_AD(repmat(ABCD_GP_AD == 2, nPerm, 1)' & (PERM_HCP_GP_AD == -1)) = 0.5;
PERM_GP_AD(repmat(ABCD_GP_AD == -2, nPerm, 1)' & (PERM_HCP_GP_AD == 2)) = 0.1;
PERM_GP_AD(repmat(ABCD_GP_AD == 2, nPerm, 1)' & (PERM_HCP_GP_AD == -2)) = 0.1;
PERM_GP_AD(repmat(ABCD_GP_AD == 1, nPerm, 1)' & (PERM_HCP_GP_AD == 1)) = 0.9;
PERM_GP_AD(repmat(ABCD_GP_AD == -1, nPerm, 1)' & (PERM_HCP_GP_AD == 1)) = 0.9;
PERM_GP_AD(repmat(ABCD_GP_AD == -1, nPerm, 1)' & (PERM_HCP_GP_AD == -1)) = 0.9;
PERM_GP_AD(repmat(ABCD_GP_AD == 1, nPerm, 1)' & (PERM_HCP_GP_AD == -1)) = 0.9;

PERM_THAL_AD = zeros(13000, nPerm);
PERM_THAL_AD(repmat(ABCD_THAL_AD == 2, nPerm, 1)' & (PERM_HCP_THAL_AD == 2)) = -1.9;
PERM_THAL_AD(repmat(ABCD_THAL_AD == -2, nPerm, 1)' & (PERM_HCP_THAL_AD == -2)) = 1.1;
PERM_THAL_AD(repmat(ABCD_THAL_AD == 1, nPerm, 1)' & (PERM_HCP_THAL_AD == 2)) = -1.5;
PERM_THAL_AD(repmat(ABCD_THAL_AD == -1, nPerm, 1)' & (PERM_HCP_THAL_AD == 2)) = -1.5;
PERM_THAL_AD(repmat(ABCD_THAL_AD == 1, nPerm, 1)' & (PERM_HCP_THAL_AD == -2)) = 1.5;
PERM_THAL_AD(repmat(ABCD_THAL_AD == -1, nPerm, 1)' & (PERM_HCP_THAL_AD == -2)) = 1.5;
PERM_THAL_AD(repmat(ABCD_THAL_AD == -2, nPerm, 1)' & (PERM_HCP_THAL_AD == 1)) = 0.5;
PERM_THAL_AD(repmat(ABCD_THAL_AD == 2, nPerm, 1)' & (PERM_HCP_THAL_AD == 1)) = 0.5;
PERM_THAL_AD(repmat(ABCD_THAL_AD == -2, nPerm, 1)' & (PERM_HCP_THAL_AD == -1)) = 0.5;
PERM_THAL_AD(repmat(ABCD_THAL_AD == 2, nPerm, 1)' & (PERM_HCP_THAL_AD == -1)) = 0.5;
PERM_THAL_AD(repmat(ABCD_THAL_AD == -2, nPerm, 1)' & (PERM_HCP_THAL_AD == 2)) = 0.1;
PERM_THAL_AD(repmat(ABCD_THAL_AD == 2, nPerm, 1)' & (PERM_HCP_THAL_AD == -2)) = 0.1;
PERM_THAL_AD(repmat(ABCD_THAL_AD == 1, nPerm, 1)' & (PERM_HCP_THAL_AD == 1)) = 0.9;
PERM_THAL_AD(repmat(ABCD_THAL_AD == -1, nPerm, 1)' & (PERM_HCP_THAL_AD == 1)) = 0.9;
PERM_THAL_AD(repmat(ABCD_THAL_AD == -1, nPerm, 1)' & (PERM_HCP_THAL_AD == -1)) = 0.9;
PERM_THAL_AD(repmat(ABCD_THAL_AD == 1, nPerm, 1)' & (PERM_HCP_THAL_AD == -1)) = 0.9;

PERM_HIPP_AD = zeros(2367, nPerm);
PERM_HIPP_AD(repmat(ABCD_HIPP_AD == 2, nPerm, 1)' & (PERM_HCP_HIPP_AD == 2)) = -1.9;
PERM_HIPP_AD(repmat(ABCD_HIPP_AD == -2, nPerm, 1)' & (PERM_HCP_HIPP_AD == -2)) = 1.1;
PERM_HIPP_AD(repmat(ABCD_HIPP_AD == 1, nPerm, 1)' & (PERM_HCP_HIPP_AD == 2)) = -1.5;
PERM_HIPP_AD(repmat(ABCD_HIPP_AD == -1, nPerm, 1)' & (PERM_HCP_HIPP_AD == 2)) = -1.5;
PERM_HIPP_AD(repmat(ABCD_HIPP_AD == 1, nPerm, 1)' & (PERM_HCP_HIPP_AD == -2)) = 1.5;
PERM_HIPP_AD(repmat(ABCD_HIPP_AD == -1, nPerm, 1)' & (PERM_HCP_HIPP_AD == -2)) = 1.5;
PERM_HIPP_AD(repmat(ABCD_HIPP_AD == -2, nPerm, 1)' & (PERM_HCP_HIPP_AD == 1)) = 0.5;
PERM_HIPP_AD(repmat(ABCD_HIPP_AD == 2, nPerm, 1)' & (PERM_HCP_HIPP_AD == 1)) = 0.5;
PERM_HIPP_AD(repmat(ABCD_HIPP_AD == -2, nPerm, 1)' & (PERM_HCP_HIPP_AD == -1)) = 0.5;
PERM_HIPP_AD(repmat(ABCD_HIPP_AD == 2, nPerm, 1)' & (PERM_HCP_HIPP_AD == -1)) = 0.5;
PERM_HIPP_AD(repmat(ABCD_HIPP_AD == -2, nPerm, 1)' & (PERM_HCP_HIPP_AD == 2)) = 0.1;
PERM_HIPP_AD(repmat(ABCD_HIPP_AD == 2, nPerm, 1)' & (PERM_HCP_HIPP_AD == -2)) = 0.1;
PERM_HIPP_AD(repmat(ABCD_HIPP_AD == 1, nPerm, 1)' & (PERM_HCP_HIPP_AD == 1)) = 0.9;
PERM_HIPP_AD(repmat(ABCD_HIPP_AD == -1, nPerm, 1)' & (PERM_HCP_HIPP_AD == 1)) = 0.9;
PERM_HIPP_AD(repmat(ABCD_HIPP_AD == -1, nPerm, 1)' & (PERM_HCP_HIPP_AD == -1)) = 0.9;
PERM_HIPP_AD(repmat(ABCD_HIPP_AD == 1, nPerm, 1)' & (PERM_HCP_HIPP_AD == -1)) = 0.9;

PERM_AMYG_AD = zeros(2878, nPerm);
PERM_AMYG_AD(repmat(ABCD_AMYG_AD == 2, nPerm, 1)' & (PERM_HCP_AMYG_AD == 2)) = -1.9;
PERM_AMYG_AD(repmat(ABCD_AMYG_AD == -2, nPerm, 1)' & (PERM_HCP_AMYG_AD == -2)) = 1.1;
PERM_AMYG_AD(repmat(ABCD_AMYG_AD == 1, nPerm, 1)' & (PERM_HCP_AMYG_AD == 2)) = -1.5;
PERM_AMYG_AD(repmat(ABCD_AMYG_AD == -1, nPerm, 1)' & (PERM_HCP_AMYG_AD == 2)) = -1.5;
PERM_AMYG_AD(repmat(ABCD_AMYG_AD == 1, nPerm, 1)' & (PERM_HCP_AMYG_AD == -2)) = 1.5;
PERM_AMYG_AD(repmat(ABCD_AMYG_AD == -1, nPerm, 1)' & (PERM_HCP_AMYG_AD == -2)) = 1.5;
PERM_AMYG_AD(repmat(ABCD_AMYG_AD == -2, nPerm, 1)' & (PERM_HCP_AMYG_AD == 1)) = 0.5;
PERM_AMYG_AD(repmat(ABCD_AMYG_AD == 2, nPerm, 1)' & (PERM_HCP_AMYG_AD == 1)) = 0.5;
PERM_AMYG_AD(repmat(ABCD_AMYG_AD == -2, nPerm, 1)' & (PERM_HCP_AMYG_AD == -1)) = 0.5;
PERM_AMYG_AD(repmat(ABCD_AMYG_AD == 2, nPerm, 1)' & (PERM_HCP_AMYG_AD == -1)) = 0.5;
PERM_AMYG_AD(repmat(ABCD_AMYG_AD == -2, nPerm, 1)' & (PERM_HCP_AMYG_AD == 2)) = 0.1;
PERM_AMYG_AD(repmat(ABCD_AMYG_AD == 2, nPerm, 1)' & (PERM_HCP_AMYG_AD == -2)) = 0.1;
PERM_AMYG_AD(repmat(ABCD_AMYG_AD == 1, nPerm, 1)' & (PERM_HCP_AMYG_AD == 1)) = 0.9;
PERM_AMYG_AD(repmat(ABCD_AMYG_AD == -1, nPerm, 1)' & (PERM_HCP_AMYG_AD == 1)) = 0.9;
PERM_AMYG_AD(repmat(ABCD_AMYG_AD == -1, nPerm, 1)' & (PERM_HCP_AMYG_AD == -1)) = 0.9;
PERM_AMYG_AD(repmat(ABCD_AMYG_AD == 1, nPerm, 1)' & (PERM_HCP_AMYG_AD == -1)) = 0.9;

% Calculate proportions
PERM_STR = sum(repmat(S_STR_NA, 1, nPerm) .* (abs(PERM_STR_AD) > 0.8));
PERM_GP = sum(repmat(S_GP_NA, 1, nPerm) .* (abs(PERM_GP_AD) > 0.8));
PERM_THAL = sum(repmat(S_THAL_NA, 1, nPerm) .* (abs(PERM_THAL_AD) > 0.8));
PERM_HIPP = sum(repmat(S_HIPP_NA, 1, nPerm) .* (abs(PERM_HIPP_AD) > 0.8));
PERM_AMYG = sum(repmat(S_AMYG_NA, 1, nPerm) .* (abs(PERM_AMYG_AD) > 0.8));

% Calculate p value
STR_PERM_P = sum(sum(S_STR_NA(abs(STR_AD) > 0.8)) < PERM_STR) ./ nPerm;
GP_PERM_P = sum(sum(S_GP_NA(abs(GP_AD) > 0.8)) < PERM_GP) ./ nPerm;
THAL_PERM_P = sum(sum(S_THAL_NA(abs(THAL_AD) > 0.8)) < PERM_THAL) ./ nPerm;
HIPP_PERM_P = sum(sum(S_HIPP_NA(abs(HIPP_AD) > 0.8)) < PERM_HIPP) ./ nPerm;
AMYG_PERM_P = sum(sum(S_AMYG_NA(abs(AMYG_AD) > 0.8)) < PERM_AMYG) ./ nPerm;

disp('MSR permutations');
[STR_PERM_P GP_PERM_P THAL_PERM_P HIPP_PERM_P AMYG_PERM_P]

% Raincloud
clear rainData;
sizeNorm =  [26000 6000 13000 2367 2878];
rainData{1} = [sum((PERM_STR_AD) == 1.1); sum((PERM_GP_AD) == 1.1); sum((PERM_THAL_AD) == 1.1); sum((PERM_HIPP_AD) == 1.1); sum((PERM_AMYG_AD) == 1.1)]' ./ sizeNorm;
rainData{2} = [sum((PERM_STR_AD) == 1.5); sum((PERM_GP_AD) == 1.5); sum((PERM_THAL_AD) == 1.5); sum((PERM_HIPP_AD) == 1.5); sum((PERM_AMYG_AD) == 1.5)]' ./ sizeNorm;
rainData{3} = [sum((PERM_STR_AD) == -1.9); sum((PERM_GP_AD) == -1.9); sum((PERM_THAL_AD) == -1.9); sum((PERM_HIPP_AD) == -1.9); sum((PERM_AMYG_AD) == -1.9)]' ./ sizeNorm;
rainData{4} = [sum((PERM_STR_AD) == -1.5); sum((PERM_GP_AD) == -1.5); sum((PERM_THAL_AD) == -1.5); sum((PERM_HIPP_AD) == -1.5); sum((PERM_AMYG_AD) == -1.5)]' ./ sizeNorm;
rainData{5} = [sum((PERM_STR_AD) == 0.1); sum((PERM_GP_AD) == 0.1); sum((PERM_THAL_AD) == 0.1); sum((PERM_HIPP_AD) == 0.1); sum((PERM_AMYG_AD) == 0.1)]' ./ sizeNorm;
rainData{6} = [sum((PERM_STR_AD) == 0.5); sum((PERM_GP_AD) == 0.5); sum((PERM_THAL_AD) == 0.5); sum((PERM_HIPP_AD) == 0.5); sum((PERM_AMYG_AD) == 0.5)]' ./ sizeNorm;
rainData{7} = [sum((PERM_STR_AD) == 0.9); sum((PERM_GP_AD) == 0.9); sum((PERM_THAL_AD) == 0.9); sum((PERM_HIPP_AD) == 0.9); sum((PERM_AMYG_AD) == 0.9)]' ./ sizeNorm;
group_inx = [ones(1,nPerm), 2.*ones(1,nPerm) 3.*ones(1,nPerm) 4.*ones(1,nPerm) 5.*ones(1,nPerm) 6.*ones(1,nPerm)];
c =  [0.5569    0.0039    0.3216
      0.8591    0.4207    0.6522
      0.1529    0.3922    0.0980
      0.4752    0.7183    0.2365
      0.2000    0.2000    0.2000
      0.5882    0.5882    0.5882
      1         1         1];
group_names = {'Established male-bias', 'Emerging male-bias', 'Established female-bias', 'Emerging female-bias', 'Hard flip', 'Diminishing', 'Never sex-biased'};
condition_names = {'Striatum', 'Globus pallidus', 'Thalamus', 'Hippocampus', 'Amygdala'};
h = figure; d = daviolinplot(rainData, 'groups', group_inx(1:90), 'colors', c, 'box', 3,...
    'boxcolor', 'w', 'scatter', 2, 'jitter', 1, 'scattercolor', 'same', 'smoothing', 1e-2,...
    'scattersize', 10, 'scatteralpha', 0.4, 'outliers', 0, ...
    'xtlabels', condition_names,...
    'legend', group_names);
ylabel('Null distribution of proportion overlap');
set(d.sc,'MarkerEdgeColor','none');
set(gca,'FontSize',10);
xl = xlim; xlim([xl(1)-0.1, xl(2)+0.3]);
ylim([0 1]); yticks(0:0.1:1);
ax = gca;
ax.YGrid = 'on'; ax.GridLineWidth = 2;
ax.TickDir = 'out';
h.Position = [0 0 2000 1000];
% legend('off');
saveas(h, fullfile(figDir, 'raincloud_null_overlap.tif'));

% Compute statistics
permData{1} = [sum((PERM_STR_AD) == 1.1); sum((PERM_STR_AD) == 1.5); sum((PERM_STR_AD) == -1.9); sum((PERM_STR_AD) == -1.5); sum((PERM_STR_AD) == 0.1); sum((PERM_STR_AD) == 0.5); sum((PERM_STR_AD) == 0.9)] ./ 26000;
permData{2} = [sum((PERM_GP_AD) == 1.1); sum((PERM_GP_AD) == 1.5); sum((PERM_GP_AD) == -1.9); sum((PERM_GP_AD) == -1.5); sum((PERM_GP_AD) == 0.1); sum((PERM_GP_AD) == 0.5); sum((PERM_GP_AD) == 0.9)] ./ 6000;
permData{3} = [sum((PERM_THAL_AD) == 1.1); sum((PERM_THAL_AD) == 1.5); sum((PERM_THAL_AD) == -1.9); sum((PERM_THAL_AD) == -1.5); sum((PERM_THAL_AD) == 0.1); sum((PERM_THAL_AD) == 0.5); sum((PERM_THAL_AD) == 0.9)] ./ 13000;
permData{4} = [sum((PERM_HIPP_AD) == 1.1); sum((PERM_HIPP_AD) == 1.5); sum((PERM_HIPP_AD) == -1.9); sum((PERM_HIPP_AD) == -1.5); sum((PERM_HIPP_AD) == 0.1); sum((PERM_HIPP_AD) == 0.5); sum((PERM_HIPP_AD) == 0.9)] ./ 2367;
permData{5} = [sum((PERM_AMYG_AD) == 1.1); sum((PERM_AMYG_AD) == 1.5); sum((PERM_AMYG_AD) == -1.9); sum((PERM_AMYG_AD) == -1.5); sum((PERM_AMYG_AD) == 0.1); sum((PERM_AMYG_AD) == 0.5); sum((PERM_AMYG_AD) == 0.9)] ./ 2878;

realData{1} = [sum((STR_AD) == 1.1); sum((STR_AD) == 1.5); sum((STR_AD) == -1.9); sum((STR_AD) == -1.5); sum((STR_AD) == 0.1); sum((STR_AD) == 0.5); sum((STR_AD) == 0.9)] ./ 26000;
realData{2} = [sum((GP_AD) == 1.1); sum((GP_AD) == 1.5); sum((GP_AD) == -1.9); sum((GP_AD) == -1.5); sum((GP_AD) == 0.1); sum((GP_AD) == 0.5); sum((GP_AD) == 0.9)] ./ 6000;
realData{3} = [sum((THAL_AD) == 1.1); sum((THAL_AD) == 1.5); sum((THAL_AD) == -1.9); sum((THAL_AD) == -1.5); sum((THAL_AD) == 0.1); sum((THAL_AD) == 0.5); sum((THAL_AD) == 0.9)] ./ 13000;
realData{4} = [sum((HIPP_AD) == 1.1); sum((HIPP_AD) == 1.5); sum((HIPP_AD) == -1.9); sum((HIPP_AD) == -1.5); sum((HIPP_AD) == 0.1); sum((HIPP_AD) == 0.5); sum((HIPP_AD) == 0.9)] ./ 2367;
realData{5} = [sum((AMYG_AD) == 1.1); sum((AMYG_AD) == 1.5); sum((AMYG_AD) == -1.9); sum((AMYG_AD) == -1.5); sum((AMYG_AD) == 0.1); sum((AMYG_AD) == 0.5); sum((AMYG_AD) == 0.9)] ./ 2878;

STR_PERM_P = sum(realData{1} < permData{1}, 2) ./ nPerm;
GP_PERM_P = sum(realData{2} < permData{2}, 2) ./ nPerm;
THAL_PERM_P = sum(realData{3} < permData{3}, 2) ./ nPerm;
HIPP_PERM_P = sum(realData{4} < permData{4}, 2) ./ nPerm;
AMYG_PERM_P = sum(realData{5} < permData{5}, 2) ./ nPerm;

disp('Proportion matrix');
% [mean(permData{1},2) mean(permData{2},2) mean(permData{3},2) mean(permData{4},2) mean(permData{5},2)]
cell2mat(realData)

disp('MSR permutations');
[STR_PERM_P GP_PERM_P THAL_PERM_P HIPP_PERM_P AMYG_PERM_P]

%% ABCD & HCP Conjunction for cortex

% CV
ABCD_CV_AD = zeros(1,360);
HCP_CV_AD = zeros(1,360);
CV_AD = zeros(1,360);

ABCD_CV_AD(ABCD.STAT.COR.CV_B > 0) =  1;
ABCD_CV_AD(ABCD.STAT.COR.CV_B < 0) = -1;
ABCD_CV_AD((ABCD.STAT.COR.CV_Q < thresh) & (ABCD.STAT.COR.CV_B > 0))  =  2;
ABCD_CV_AD((ABCD.STAT.COR.CV_Q < thresh) & (ABCD.STAT.COR.CV_B < 0))  = -2;

HCP_CV_AD(HCP.STAT.COR.CV_B > 0) =  1;
HCP_CV_AD(HCP.STAT.COR.CV_B < 0) = -1;
HCP_CV_AD((HCP.STAT.COR.CV_Q < thresh) & (HCP.STAT.COR.CV_B > 0))  =  2;
HCP_CV_AD((HCP.STAT.COR.CV_Q < thresh) & (HCP.STAT.COR.CV_B < 0))  = -2;

% Established: FF, MM & stable (dark green or magenta)
CV_AD((ABCD_CV_AD == 2) & (HCP_CV_AD == 2)) = -1.9;
CV_AD((ABCD_CV_AD == -2) & (HCP_CV_AD == -2)) = 1.1;

% Emerging: fF, mM, mF, fM (medium green or magenta)
CV_AD((ABCD_CV_AD == 1) & (HCP_CV_AD == 2)) = -1.5;
CV_AD((ABCD_CV_AD == -1) & (HCP_CV_AD == -2)) = 1.5;
CV_AD((ABCD_CV_AD == -1) & (HCP_CV_AD == 2)) = -1.5;
CV_AD((ABCD_CV_AD == 1) & (HCP_CV_AD == -2)) = 1.5;

% Diminishing flip: Ff, Fm, Mm, Mf (light gray)
CV_AD((ABCD_CV_AD == 2) & (HCP_CV_AD == 1)) = 0.5;
CV_AD((ABCD_CV_AD == 2) & (HCP_CV_AD == -1)) = 0.5;
CV_AD((ABCD_CV_AD == -2) & (HCP_CV_AD == -1)) = 0.5;
CV_AD((ABCD_CV_AD == -2) & (HCP_CV_AD == 1)) = 0.5;

% Hard flip: MF, FM (dark gray)
CV_AD((ABCD_CV_AD == -2) & (HCP_CV_AD == 2)) = 0.1;
CV_AD((ABCD_CV_AD == 2) & (HCP_CV_AD == -2)) = 0.1;

% Insignificant sex-biases: ff, mf, mm, fm and insignificant delta (white)
CV_AD((ABCD_CV_AD == 1) & (HCP_CV_AD == 1)) = 0.9;
CV_AD((ABCD_CV_AD == -1) & (HCP_CV_AD == 1)) = 0.9;
CV_AD((ABCD_CV_AD == -1) & (HCP_CV_AD == -1)) = 0.9;
CV_AD((ABCD_CV_AD == 1) & (HCP_CV_AD == -1)) = 0.9;

h = HyoSurfacePlotFunction(CV_AD, IS, G360, cmap.tri_bi, '', [-2 2], 'square');
% exportgraphics(h, fullfile(figDir, 'abcd_hcp_overlap_CV.png'), 'Resolution', 300);

CV_AD_Count = [sum(CV_AD == 1.1) sum(CV_AD == 1.5) 0 ...
                sum(CV_AD == -1.9) sum(CV_AD == -1.5) 0 ...
                sum(CV_AD == 0.1) sum(CV_AD == 0.5) sum(CV_AD == 0.9)];

h = figure;
b = bar(1, CV_AD_Count, 'stacked', 'EdgeColor', 'none');
for i = 1:9
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_tri_bi(i,:);
end
box off
axis off
% exportgraphics(h, fullfile(figDir, 'abcd_hcp_overlap_bar_CV.png'), 'Resolution', 300);

% SA
ABCD_SA_AD = zeros(1,360);
HCP_SA_AD = zeros(1,360);
SA_AD = zeros(1,360);

ABCD_SA_AD(ABCD.STAT.COR.SA_B > 0) =  1;
ABCD_SA_AD(ABCD.STAT.COR.SA_B < 0) = -1;
ABCD_SA_AD((ABCD.STAT.COR.SA_Q < thresh) & (ABCD.STAT.COR.SA_B > 0))  =  2;
ABCD_SA_AD((ABCD.STAT.COR.SA_Q < thresh) & (ABCD.STAT.COR.SA_B < 0))  = -2;

HCP_SA_AD(HCP.STAT.COR.SA_B > 0) =  1;
HCP_SA_AD(HCP.STAT.COR.SA_B < 0) = -1;
HCP_SA_AD((HCP.STAT.COR.SA_Q < thresh) & (HCP.STAT.COR.SA_B > 0))  =  2;
HCP_SA_AD((HCP.STAT.COR.SA_Q < thresh) & (HCP.STAT.COR.SA_B < 0))  = -2;

% Established: FF, MM & stable (dark green or magenta)
SA_AD((ABCD_SA_AD == 2) & (HCP_SA_AD == 2)) = -1.9;
SA_AD((ABCD_SA_AD == -2) & (HCP_SA_AD == -2)) = 1.1;

% Emerging: fF, mM, mF, fM (medium green or magenta)
SA_AD((ABCD_SA_AD == 1) & (HCP_SA_AD == 2)) = -1.5;
SA_AD((ABCD_SA_AD == -1) & (HCP_SA_AD == -2)) = 1.5;
SA_AD((ABCD_SA_AD == -1) & (HCP_SA_AD == 2)) = -1.5;
SA_AD((ABCD_SA_AD == 1) & (HCP_SA_AD == -2)) = 1.5;

% Diminishing flip: Ff, Fm, Mm, Mf (light gray)
SA_AD((ABCD_SA_AD == 2) & (HCP_SA_AD == 1)) = 0.5;
SA_AD((ABCD_SA_AD == 2) & (HCP_SA_AD == -1)) = 0.5;
SA_AD((ABCD_SA_AD == -2) & (HCP_SA_AD == -1)) = 0.5;
SA_AD((ABCD_SA_AD == -2) & (HCP_SA_AD == 1)) = 0.5;

% Hard flip: MF, FM (dark gray)
SA_AD((ABCD_SA_AD == -2) & (HCP_SA_AD == 2)) = 0.1;
SA_AD((ABCD_SA_AD == 2) & (HCP_SA_AD == -2)) = 0.1;

% Insignificant sex-biases: ff, mf, mm, fm and insignificant delta (white)
SA_AD((ABCD_SA_AD == 1) & (HCP_SA_AD == 1)) = 0.9;
SA_AD((ABCD_SA_AD == -1) & (HCP_SA_AD == 1)) = 0.9;
SA_AD((ABCD_SA_AD == -1) & (HCP_SA_AD == -1)) = 0.9;
SA_AD((ABCD_SA_AD == 1) & (HCP_SA_AD == -1)) = 0.9;

h = HyoSurfacePlotFunction(SA_AD, IS, G360, cmap.tri_bi, '', [-2 2], 'square');
% exportgraphics(h, fullfile(figDir, 'abcd_hcp_overlap_SA.png'), 'Resolution', 300);

SA_AD_Count = [sum(SA_AD == 1.1) sum(SA_AD == 1.5) 0 ...
                sum(SA_AD == -1.9) sum(SA_AD == -1.5) 0 ...
                sum(SA_AD == 0.1) sum(SA_AD == 0.5) sum(SA_AD == 0.9)];

h = figure;
b = bar(1, SA_AD_Count, 'stacked', 'EdgeColor', 'none');
for i = 1:9
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_tri_bi(i,:);
end
box off
axis off
% exportgraphics(h, fullfile(figDir, 'abcd_hcp_overlap_bar_SA.png'), 'Resolution', 300);

% CT
ABCD_CT_AD = zeros(1,360);
HCP_CT_AD = zeros(1,360);
CT_AD = zeros(1,360);

ABCD_CT_AD(ABCD.STAT.COR.CT_B > 0) =  1;
ABCD_CT_AD(ABCD.STAT.COR.CT_B < 0) = -1;
ABCD_CT_AD((ABCD.STAT.COR.CT_Q < thresh) & (ABCD.STAT.COR.CT_B > 0))  =  2;
ABCD_CT_AD((ABCD.STAT.COR.CT_Q < thresh) & (ABCD.STAT.COR.CT_B < 0))  = -2;

HCP_CT_AD(HCP.STAT.COR.CT_B > 0) =  1;
HCP_CT_AD(HCP.STAT.COR.CT_B < 0) = -1;
HCP_CT_AD((HCP.STAT.COR.CT_Q < thresh) & (HCP.STAT.COR.CT_B > 0))  =  2;
HCP_CT_AD((HCP.STAT.COR.CT_Q < thresh) & (HCP.STAT.COR.CT_B < 0))  = -2;

% Established: FF, MM & stable (dark green or magenta)
CT_AD((ABCD_CT_AD == 2) & (HCP_CT_AD == 2)) = -1.9;
CT_AD((ABCD_CT_AD == -2) & (HCP_CT_AD == -2)) = 1.1;

% Emerging: fF, mM, mF, fM (medium green or magenta)
CT_AD((ABCD_CT_AD == 1) & (HCP_CT_AD == 2)) = -1.5;
CT_AD((ABCD_CT_AD == -1) & (HCP_CT_AD == -2)) = 1.5;
CT_AD((ABCD_CT_AD == -1) & (HCP_CT_AD == 2)) = -1.5;
CT_AD((ABCD_CT_AD == 1) & (HCP_CT_AD == -2)) = 1.5;

% Diminishing flip: Ff, Fm, Mm, Mf (light gray)
CT_AD((ABCD_CT_AD == 2) & (HCP_CT_AD == 1)) = 0.5;
CT_AD((ABCD_CT_AD == 2) & (HCP_CT_AD == -1)) = 0.5;
CT_AD((ABCD_CT_AD == -2) & (HCP_CT_AD == -1)) = 0.5;
CT_AD((ABCD_CT_AD == -2) & (HCP_CT_AD == 1)) = 0.5;

% Hard flip: MF, FM (dark gray)
CT_AD((ABCD_CT_AD == -2) & (HCP_CT_AD == 2)) = 0.1;
CT_AD((ABCD_CT_AD == 2) & (HCP_CT_AD == -2)) = 0.1;

% Insignificant sex-biases: ff, mf, mm, fm and insignificant delta (white)
CT_AD((ABCD_CT_AD == 1) & (HCP_CT_AD == 1)) = 0.9;
CT_AD((ABCD_CT_AD == -1) & (HCP_CT_AD == 1)) = 0.9;
CT_AD((ABCD_CT_AD == -1) & (HCP_CT_AD == -1)) = 0.9;
CT_AD((ABCD_CT_AD == 1) & (HCP_CT_AD == -1)) = 0.9;

h = HyoSurfacePlotFunction(CT_AD, IS, G360, cmap.tri_bi, '', [-2 2], 'square');
% exportgraphics(h, fullfile(figDir, 'abcd_hcp_overlap_CT.png'), 'Resolution', 300);

CT_AD_Count = [sum(CT_AD == 1.1) sum(CT_AD == 1.5) 0 ...
                sum(CT_AD == -1.9) sum(CT_AD == -1.5) 0 ...
                sum(CT_AD == 0.1) sum(CT_AD == 0.5) sum(CT_AD == 0.9)];

h = figure;
b = bar(1, CT_AD_Count, 'stacked', 'EdgeColor', 'none');
for i = 1:9
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_tri_bi(i,:);
end
box off
axis off
% exportgraphics(h, fullfile(figDir, 'abcd_hcp_overlap_bar_CT.png'), 'Resolution', 300);

close all;

[sum(CV_AD_Count([1:6 9])) / sum(CV_AD_Count)]
[sum(SA_AD_Count([1:6 9])) / sum(SA_AD_Count)]
[sum(CT_AD_Count([1:6 9])) / sum(CT_AD_Count)]

%% Spin test (cortex)

% Spin permutations
for i = 1:nPerm
    HCP_CV_B_Spin(:,i) = HCP.STAT.COR.CV_B(Vasa(:,i));
    HCP_SA_B_Spin(:,i) = HCP.STAT.COR.SA_B(Vasa(:,i));
    HCP_CT_B_Spin(:,i) = HCP.STAT.COR.CT_B(Vasa(:,i));

    HCP_CV_Q_Spin(:,i) = HCP.STAT.COR.CV_Q(Vasa(:,i));
    HCP_SA_Q_Spin(:,i) = HCP.STAT.COR.SA_Q(Vasa(:,i));
    HCP_CT_Q_Spin(:,i) = HCP.STAT.COR.CT_Q(Vasa(:,i));
end

% Discretize spatially permuted HCP sex difference maps
PERM_HCP_CV_AD = zeros(360, nPerm);
PERM_HCP_CV_AD(HCP_CV_B_Spin > 0) =  1;
PERM_HCP_CV_AD(HCP_CV_B_Spin < 0) = -1;
PERM_HCP_CV_AD((HCP_CV_Q_Spin < thresh) & (HCP_CV_B_Spin > 0))  =  2;
PERM_HCP_CV_AD((HCP_CV_Q_Spin < thresh) & (HCP_CV_B_Spin < 0))  = -2;

PERM_HCP_SA_AD = zeros(360, nPerm);
PERM_HCP_SA_AD(HCP_CV_B_Spin > 0) =  1;
PERM_HCP_SA_AD(HCP_CV_B_Spin < 0) = -1;
PERM_HCP_SA_AD((HCP_CV_Q_Spin < thresh) & (HCP_CV_B_Spin > 0))  =  2;
PERM_HCP_SA_AD((HCP_CV_Q_Spin < thresh) & (HCP_CV_B_Spin < 0))  = -2;

PERM_HCP_CT_AD = zeros(360, nPerm);
PERM_HCP_CT_AD(HCP_CV_B_Spin > 0) =  1;
PERM_HCP_CT_AD(HCP_CV_B_Spin < 0) = -1;
PERM_HCP_CT_AD((HCP_CV_Q_Spin < thresh) & (HCP_CV_B_Spin > 0))  =  2;
PERM_HCP_CT_AD((HCP_CV_Q_Spin < thresh) & (HCP_CV_B_Spin < 0))  = -2;

% Compute agreement maps
PERM_CV_AD = zeros(360, nPerm);
PERM_CV_AD(repmat(ABCD_CV_AD == 2, nPerm, 1)' & (PERM_HCP_CV_AD == 2)) = -1.9;
PERM_CV_AD(repmat(ABCD_CV_AD == -2, nPerm, 1)' & (PERM_HCP_CV_AD == -2)) = 1.1;
PERM_CV_AD(repmat(ABCD_CV_AD == 1, nPerm, 1)' & (PERM_HCP_CV_AD == 2)) = -1.5;
PERM_CV_AD(repmat(ABCD_CV_AD == -1, nPerm, 1)' & (PERM_HCP_CV_AD == 2)) = -1.5;
PERM_CV_AD(repmat(ABCD_CV_AD == 1, nPerm, 1)' & (PERM_HCP_CV_AD == -2)) = 1.5;
PERM_CV_AD(repmat(ABCD_CV_AD == -1, nPerm, 1)' & (PERM_HCP_CV_AD == -2)) = 1.5;
PERM_CV_AD(repmat(ABCD_CV_AD == -2, nPerm, 1)' & (PERM_HCP_CV_AD == 1)) = 0.5;
PERM_CV_AD(repmat(ABCD_CV_AD == 2, nPerm, 1)' & (PERM_HCP_CV_AD == 1)) = 0.5;
PERM_CV_AD(repmat(ABCD_CV_AD == -2, nPerm, 1)' & (PERM_HCP_CV_AD == -1)) = 0.5;
PERM_CV_AD(repmat(ABCD_CV_AD == 2, nPerm, 1)' & (PERM_HCP_CV_AD == -1)) = 0.5;
PERM_CV_AD(repmat(ABCD_CV_AD == -2, nPerm, 1)' & (PERM_HCP_CV_AD == 2)) = 0.1;
PERM_CV_AD(repmat(ABCD_CV_AD == 2, nPerm, 1)' & (PERM_HCP_CV_AD == -2)) = 0.1;
PERM_CV_AD(repmat(ABCD_CV_AD == 1, nPerm, 1)' & (PERM_HCP_CV_AD == 1)) = 0.9;
PERM_CV_AD(repmat(ABCD_CV_AD == -1, nPerm, 1)' & (PERM_HCP_CV_AD == 1)) = 0.9;
PERM_CV_AD(repmat(ABCD_CV_AD == -1, nPerm, 1)' & (PERM_HCP_CV_AD == -1)) = 0.9;
PERM_CV_AD(repmat(ABCD_CV_AD == 1, nPerm, 1)' & (PERM_HCP_CV_AD == -1)) = 0.9;

PERM_SA_AD = zeros(360, nPerm);
PERM_SA_AD(repmat(ABCD_SA_AD == 2, nPerm, 1)' & (PERM_HCP_SA_AD == 2)) = -1.9;
PERM_SA_AD(repmat(ABCD_SA_AD == -2, nPerm, 1)' & (PERM_HCP_SA_AD == -2)) = 1.1;
PERM_SA_AD(repmat(ABCD_SA_AD == 1, nPerm, 1)' & (PERM_HCP_SA_AD == 2)) = -1.5;
PERM_SA_AD(repmat(ABCD_SA_AD == -1, nPerm, 1)' & (PERM_HCP_SA_AD == 2)) = -1.5;
PERM_SA_AD(repmat(ABCD_SA_AD == 1, nPerm, 1)' & (PERM_HCP_SA_AD == -2)) = 1.5;
PERM_SA_AD(repmat(ABCD_SA_AD == -1, nPerm, 1)' & (PERM_HCP_SA_AD == -2)) = 1.5;
PERM_SA_AD(repmat(ABCD_SA_AD == -2, nPerm, 1)' & (PERM_HCP_SA_AD == 1)) = 0.5;
PERM_SA_AD(repmat(ABCD_SA_AD == 2, nPerm, 1)' & (PERM_HCP_SA_AD == 1)) = 0.5;
PERM_SA_AD(repmat(ABCD_SA_AD == -2, nPerm, 1)' & (PERM_HCP_SA_AD == -1)) = 0.5;
PERM_SA_AD(repmat(ABCD_SA_AD == 2, nPerm, 1)' & (PERM_HCP_SA_AD == -1)) = 0.5;
PERM_SA_AD(repmat(ABCD_SA_AD == -2, nPerm, 1)' & (PERM_HCP_SA_AD == 2)) = 0.1;
PERM_SA_AD(repmat(ABCD_SA_AD == 2, nPerm, 1)' & (PERM_HCP_SA_AD == -2)) = 0.1;
PERM_SA_AD(repmat(ABCD_SA_AD == 1, nPerm, 1)' & (PERM_HCP_SA_AD == 1)) = 0.9;
PERM_SA_AD(repmat(ABCD_SA_AD == -1, nPerm, 1)' & (PERM_HCP_SA_AD == 1)) = 0.9;
PERM_SA_AD(repmat(ABCD_SA_AD == -1, nPerm, 1)' & (PERM_HCP_SA_AD == -1)) = 0.9;
PERM_SA_AD(repmat(ABCD_SA_AD == 1, nPerm, 1)' & (PERM_HCP_SA_AD == -1)) = 0.9;

PERM_CT_AD = zeros(360, nPerm);
PERM_CT_AD(repmat(ABCD_CT_AD == 2, nPerm, 1)' & (PERM_HCP_CT_AD == 2)) = -1.9;
PERM_CT_AD(repmat(ABCD_CT_AD == -2, nPerm, 1)' & (PERM_HCP_CT_AD == -2)) = 1.1;
PERM_CT_AD(repmat(ABCD_CT_AD == 1, nPerm, 1)' & (PERM_HCP_CT_AD == 2)) = -1.5;
PERM_CT_AD(repmat(ABCD_CT_AD == -1, nPerm, 1)' & (PERM_HCP_CT_AD == 2)) = -1.5;
PERM_CT_AD(repmat(ABCD_CT_AD == 1, nPerm, 1)' & (PERM_HCP_CT_AD == -2)) = 1.5;
PERM_CT_AD(repmat(ABCD_CT_AD == -1, nPerm, 1)' & (PERM_HCP_CT_AD == -2)) = 1.5;
PERM_CT_AD(repmat(ABCD_CT_AD == -2, nPerm, 1)' & (PERM_HCP_CT_AD == 1)) = 0.5;
PERM_CT_AD(repmat(ABCD_CT_AD == 2, nPerm, 1)' & (PERM_HCP_CT_AD == 1)) = 0.5;
PERM_CT_AD(repmat(ABCD_CT_AD == -2, nPerm, 1)' & (PERM_HCP_CT_AD == -1)) = 0.5;
PERM_CT_AD(repmat(ABCD_CT_AD == 2, nPerm, 1)' & (PERM_HCP_CT_AD == -1)) = 0.5;
PERM_CT_AD(repmat(ABCD_CT_AD == -2, nPerm, 1)' & (PERM_HCP_CT_AD == 2)) = 0.1;
PERM_CT_AD(repmat(ABCD_CT_AD == 2, nPerm, 1)' & (PERM_HCP_CT_AD == -2)) = 0.1;
PERM_CT_AD(repmat(ABCD_CT_AD == 1, nPerm, 1)' & (PERM_HCP_CT_AD == 1)) = 0.9;
PERM_CT_AD(repmat(ABCD_CT_AD == -1, nPerm, 1)' & (PERM_HCP_CT_AD == 1)) = 0.9;
PERM_CT_AD(repmat(ABCD_CT_AD == -1, nPerm, 1)' & (PERM_HCP_CT_AD == -1)) = 0.9;
PERM_CT_AD(repmat(ABCD_CT_AD == 1, nPerm, 1)' & (PERM_HCP_CT_AD == -1)) = 0.9;

% Calculate proportions
PERM_CV = sum(abs(PERM_CV_AD) > 0.8);
PERM_SA = sum(abs(PERM_SA_AD) > 0.8);
PERM_CT = sum(abs(PERM_CT_AD) > 0.8);

% Calculate p value
CV_PERM_P = sum(sum(abs(CV_AD) > 0.8) < PERM_CV) ./ nPerm;
SA_PERM_P = sum(sum(abs(SA_AD) > 0.8) < PERM_SA) ./ nPerm;
CT_PERM_P = sum(sum(abs(CT_AD) > 0.8) < PERM_CT) ./ nPerm;

disp('MSR permutations');
[CV_PERM_P SA_PERM_P CT_PERM_P]

%% Raincloud plot for bulk sex differences

% Extract absolute volume
STR = ((ABCD.SUB.VOL.STR_L(ABCD.QC.indSTR) + ABCD.SUB.VOL.STR_R(ABCD.QC.indSTR))');
GP = ((ABCD.SUB.VOL.GP_L(ABCD.QC.indGP) + ABCD.SUB.VOL.GP_R(ABCD.QC.indGP))');
THAL = ((ABCD.SUB.VOL.THAL_L(ABCD.QC.indTHAL) + ABCD.SUB.VOL.THAL_R(ABCD.QC.indTHAL))');
HIPP = ((ABCD.SUB.VOL.HIPP_L(ABCD.QC.indHIPP) + ABCD.SUB.VOL.HIPP_R(ABCD.QC.indHIPP))');
AMYG = ((ABCD.SUB.VOL.AMYG_L(ABCD.QC.indAMYG) + ABCD.SUB.VOL.AMYG_R(ABCD.QC.indAMYG))');

% Mean-center covariates
SEX_STR = ABCD.SEX(ABCD.QC.indSTR);
AGE_STR = ABCD.AGE(ABCD.QC.indSTR) - mean(ABCD.AGE(ABCD.QC.indSTR));
TBV_STR = ABCD.TBV(ABCD.QC.indSTR) - mean(ABCD.TBV(ABCD.QC.indSTR));
SITE_STR = dummyvar(categorical(ABCD.SITE(ABCD.QC.indSTR)));
SITE_STR = SITE_STR(:,2:end);

SEX_GP = ABCD.SEX(ABCD.QC.indGP);
AGE_GP = ABCD.AGE(ABCD.QC.indGP) - mean(ABCD.AGE(ABCD.QC.indGP));
TBV_GP = ABCD.TBV(ABCD.QC.indGP) - mean(ABCD.TBV(ABCD.QC.indGP));
SITE_GP = dummyvar(categorical(ABCD.SITE(ABCD.QC.indGP)));
SITE_GP = SITE_GP(:,2:end);

SEX_THAL = ABCD.SEX(ABCD.QC.indTHAL);
AGE_THAL = ABCD.AGE(ABCD.QC.indTHAL) - mean(ABCD.AGE(ABCD.QC.indTHAL));
TBV_THAL = ABCD.TBV(ABCD.QC.indTHAL) - mean(ABCD.TBV(ABCD.QC.indTHAL));
SITE_THAL = dummyvar(categorical(ABCD.SITE(ABCD.QC.indTHAL)));
SITE_THAL = SITE_THAL(:,2:end);

SEX_HIPP = ABCD.SEX(ABCD.QC.indHIPP);
AGE_HIPP = ABCD.AGE(ABCD.QC.indHIPP) - mean(ABCD.AGE(ABCD.QC.indHIPP));
TBV_HIPP = ABCD.TBV(ABCD.QC.indHIPP) - mean(ABCD.TBV(ABCD.QC.indHIPP));
SITE_HIPP = dummyvar(categorical(ABCD.SITE(ABCD.QC.indHIPP)));
SITE_HIPP = SITE_HIPP(:,2:end);

SEX_AMYG = ABCD.SEX(ABCD.QC.indAMYG);
AGE_AMYG = ABCD.AGE(ABCD.QC.indAMYG) - mean(ABCD.AGE(ABCD.QC.indAMYG));
TBV_AMYG = ABCD.TBV(ABCD.QC.indAMYG) - mean(ABCD.TBV(ABCD.QC.indAMYG));
SITE_AMYG = dummyvar(categorical(ABCD.SITE(ABCD.QC.indAMYG)));
SITE_AMYG = SITE_AMYG(:,2:end);

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

% Regress out the age and TBV from the bulk volumes
for i = 1:size(STR,1)
    STR_SEX(i) = STR(i) - FM_STR.coef(3)*AGE_STR(i) - FM_STR.coef(4)*TBV_STR(i) - FM_STR.coef(5:end)'*SITE_STR(i,:)';
end

for i = 1:size(GP,1)
    GP_SEX(i) = GP(i) - FM_GP.coef(3)*AGE_GP(i) - FM_GP.coef(4)*TBV_GP(i) - FM_GP.coef(5:end)'*SITE_GP(i,:)';
end

for i = 1:size(THAL,1)
    THAL_SEX(i) = THAL(i) - FM_THAL.coef(3)*AGE_THAL(i) - FM_THAL.coef(4)*TBV_THAL(i) - FM_THAL.coef(5:end)'*SITE_THAL(i,:)';
end

for i = 1:size(HIPP,1)
    HIPP_SEX(i) = HIPP(i) - FM_HIPP.coef(3)*AGE_HIPP(i) - FM_HIPP.coef(4)*TBV_HIPP(i) - FM_HIPP.coef(5:end)'*SITE_HIPP(i,:)';
end

for i = 1:size(AMYG,1)
    AMYG_SEX(i) = AMYG(i) - FM_AMYG.coef(3)*AGE_AMYG(i) - FM_AMYG.coef(4)*TBV_AMYG(i) - FM_AMYG.coef(5:end)'*SITE_AMYG(i,:)';
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

%% Exploratory: IGD effects vs. pubertal change in sex effects (HCP-YA minus ABCD-9Y)

load(fullfile(dataDir, 'IGD_results.mat')); % from e_xyt_stats.m

R_STR = corr(IGD.STAT.MOR.TSA.STR_GRP_B', (HCP.STAT.MOR.TSA.STR_SEX_B - ABCD.STAT.MOR.TSA.STR_SEX_B)')
R_GP = corr(IGD.STAT.MOR.TSA.GP_GRP_B', (HCP.STAT.MOR.TSA.GP_SEX_B - ABCD.STAT.MOR.TSA.GP_SEX_B)')
R_THAL = corr(IGD.STAT.MOR.TSA.THAL_GRP_B', (HCP.STAT.MOR.TSA.THAL_SEX_B - ABCD.STAT.MOR.TSA.THAL_SEX_B)')
R_HIPP = corr(IGD.STAT.MOR.TSA.HIPP_GRP_B', (HCP.STAT.MOR.TSA.HIPP_SEX_B - ABCD.STAT.MOR.TSA.HIPP_SEX_B)')
R_AMYG = corr(IGD.STAT.MOR.TSA.AMYG_GRP_B', (HCP.STAT.MOR.TSA.AMYG_SEX_B - ABCD.STAT.MOR.TSA.AMYG_SEX_B)')
