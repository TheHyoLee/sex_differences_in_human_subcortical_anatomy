%% hcp_stats.m - Sex differences in subcortical volume and shape in early adulthood (HCP-YA)
%
% Lee HM, Liu S, Lalonde FM, et al. Sex Differences in Human Subcortical Anatomy
% are Robust and Replicable in Adulthood, Already Evident Prepubertally, and Aligned
% with Sex Chromosomal and Gonadal Influences. bioRxiv (2026).
% https://doi.org/10.64898/2026.10.05.756728
%
% Purpose:
%   Fits vertex-wise and bulk-volume models of sex differences in the HCP-YA
%   sample (striatum, pallidum, thalamus, hippocampus, amygdala), compares
%   correction methods (raw, total size, allometric), generates spatial null
%   maps (Moran spectral randomization) used by downstream scripts, and
%   produces the corresponding figures.
%
% Inputs (participant-level data are not shared; see Data availability):
%   - <dataDir>/HCP.mat (struct HCP)
%   - <templateDir>/maget_template_surfaces_distance_area.mat
%   - <templateDir>/maget_template_surfaces_all.mat
%   - <corticalResultsDir>/hcp_groupmaps_final.mat (cortical sex-difference maps
%     from Lee et al., Nature Communications 2026)
%
% Outputs:
%   - <dataDir>/HCP_results.mat (HCP struct with statistics)
%   - Figures in <figDir>
%
% Dependencies (not included in this repository):
%   - SurfStat (SurfStatLinMod, SurfStatT, term)
%   - fdr_bh (Benjamini-Hochberg FDR; MATLAB File Exchange)
%   - MATLAB Statistics and Machine Learning Toolbox
%   - BrainSpace (surface_to_graph, compute_mem, moran_randomization)
%   - daviolinplot, hexscatter (MATLAB File Exchange)
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
figDir             = '/path/to/figures/hcp';

addpath(genpath(toolboxDir));
addpath(genpath(codeDir));

%% Prep

% Load data: preprocessed participant-level data
load(fullfile(dataDir, 'HCP.mat'));

% Load vertex-wise surface area maps for template surfaces
load(fullfile(templateDir, 'maget_template_surfaces_distance_area.mat'));

% Misc
loadColorMaps;
nPerm = 10000;
thresh = 0.05;

%% Define the terms and models

% Compute z-scores (struct)
STR_L = zscore(HCP.VOL.STR_L(HCP.QC.indSTR))';
STR_R = zscore(HCP.VOL.STR_R(HCP.QC.indSTR))';

GP_L = zscore(HCP.VOL.GP_L(HCP.QC.indGP))';
GP_R = zscore(HCP.VOL.GP_R(HCP.QC.indGP))';

THAL_L = zscore(HCP.VOL.THAL_L(HCP.QC.indTHAL))';
THAL_R = zscore(HCP.VOL.THAL_R(HCP.QC.indTHAL))';

HIPP_L = zscore(HCP.VOL.HIPP_L(HCP.QC.indHIPP))';
HIPP_R = zscore(HCP.VOL.HIPP_R(HCP.QC.indHIPP))';

AMYG_L = zscore(HCP.VOL.AMYG_L(HCP.QC.indAMYG))';
AMYG_R = zscore(HCP.VOL.AMYG_R(HCP.QC.indAMYG))';

% Compute z-scores (morph)
STR_MOR = zscore(HCP.MOR.STR(HCP.QC.indSTR,:));
GP_MOR = zscore(HCP.MOR.GP(HCP.QC.indGP,:));
THAL_MOR = zscore(HCP.MOR.THAL(HCP.QC.indTHAL,:));
HIPP_MOR = zscore(HCP.MOR.HIPP(HCP.QC.indHIPP,:));
AMYG_MOR = zscore(HCP.MOR.AMYG(HCP.QC.indAMYG,:));

% Compute log scores (struct)
STR_L_LOG = log10(HCP.VOL.STR_L(HCP.QC.indSTR))';
STR_R_LOG = log10(HCP.VOL.STR_R(HCP.QC.indSTR))';

GP_L_LOG = log10(HCP.VOL.GP_L(HCP.QC.indGP))';
GP_R_LOG = log10(HCP.VOL.GP_R(HCP.QC.indGP))';

THAL_L_LOG = log10(HCP.VOL.THAL_L(HCP.QC.indTHAL))';
THAL_R_LOG = log10(HCP.VOL.THAL_R(HCP.QC.indTHAL))';

HIPP_L_LOG = log10(HCP.VOL.HIPP_L(HCP.QC.indHIPP))';
HIPP_R_LOG = log10(HCP.VOL.HIPP_R(HCP.QC.indHIPP))';

AMYG_L_LOG = log10(HCP.VOL.AMYG_L(HCP.QC.indAMYG))';
AMYG_R_LOG = log10(HCP.VOL.AMYG_R(HCP.QC.indAMYG))';

% Compute log scores (morph)
STR_MOR_LOG = log10(HCP.MOR.STR(HCP.QC.indSTR,:));
GP_MOR_LOG = log10(HCP.MOR.GP(HCP.QC.indGP,:));
THAL_MOR_LOG = log10(HCP.MOR.THAL(HCP.QC.indTHAL,:));
HIPP_MOR_LOG = log10(HCP.MOR.HIPP(HCP.QC.indHIPP,:));
AMYG_MOR_LOG = log10(HCP.MOR.AMYG(HCP.QC.indAMYG,:));

% Prepare terms (struct)
SEX_STR = term(HCP.SEX(HCP.QC.indSTR));
AGE_STR = term(HCP.AGE(HCP.QC.indSTR) - mean(HCP.AGE(HCP.QC.indSTR)));
TBV_STR = term(HCP.TBV(HCP.QC.indSTR) - mean(HCP.TBV(HCP.QC.indSTR)));
TBV_STR_LOG = term(log10(HCP.TBV(HCP.QC.indSTR)) - mean(log10(HCP.TBV(HCP.QC.indSTR))));

SEX_GP = term(HCP.SEX(HCP.QC.indGP));
AGE_GP = term(HCP.AGE(HCP.QC.indGP) - mean(HCP.AGE(HCP.QC.indGP)));
TBV_GP = term(HCP.TBV(HCP.QC.indGP) - mean(HCP.TBV(HCP.QC.indGP)));
TBV_GP_LOG = term(log10(HCP.TBV(HCP.QC.indGP)) - mean(log10(HCP.TBV(HCP.QC.indGP))));

SEX_THAL = term(HCP.SEX(HCP.QC.indTHAL));
AGE_THAL = term(HCP.AGE(HCP.QC.indTHAL) - mean(HCP.AGE(HCP.QC.indTHAL)));
TBV_THAL = term(HCP.TBV(HCP.QC.indTHAL) - mean(HCP.TBV(HCP.QC.indTHAL)));
TBV_THAL_LOG = term(log10(HCP.TBV(HCP.QC.indTHAL)) - mean(log10(HCP.TBV(HCP.QC.indTHAL))));

SEX_HIPP = term(HCP.SEX(HCP.QC.indHIPP));
AGE_HIPP = term(HCP.AGE(HCP.QC.indHIPP) - mean(HCP.AGE(HCP.QC.indHIPP)));
TBV_HIPP = term(HCP.TBV(HCP.QC.indHIPP) - mean(HCP.TBV(HCP.QC.indHIPP)));
TBV_HIPP_LOG = term(log10(HCP.TBV(HCP.QC.indHIPP)) - mean(log10(HCP.TBV(HCP.QC.indHIPP))));

SEX_AMYG = term(HCP.SEX(HCP.QC.indAMYG));
AGE_AMYG = term(HCP.AGE(HCP.QC.indAMYG) - mean(HCP.AGE(HCP.QC.indAMYG)));
TBV_AMYG = term(HCP.TBV(HCP.QC.indAMYG) - mean(HCP.TBV(HCP.QC.indAMYG)));
TBV_AMYG_LOG = term(log10(HCP.TBV(HCP.QC.indAMYG)) - mean(log10(HCP.TBV(HCP.QC.indAMYG))));

% Prepare terms (morph)
TSA_STR = term(sum(HCP.MOR.STR(HCP.QC.indSTR,:),2) - mean(sum(HCP.MOR.STR(HCP.QC.indSTR,:),2)));
TSA_STR_LOG = term(log10(sum(HCP.MOR.STR(HCP.QC.indSTR,:),2)) - mean(log10(sum(HCP.MOR.STR(HCP.QC.indSTR,:),2))));

TSA_GP = term(sum(HCP.MOR.GP(HCP.QC.indGP,:),2) - mean(sum(HCP.MOR.GP(HCP.QC.indGP,:),2)));
TSA_GP_LOG = term(log10(sum(HCP.MOR.GP(HCP.QC.indGP,:),2)) - mean(log10(sum(HCP.MOR.GP(HCP.QC.indGP,:),2))));

TSA_THAL = term(sum(HCP.MOR.THAL(HCP.QC.indTHAL,:),2) - mean(sum(HCP.MOR.THAL(HCP.QC.indTHAL,:),2)));
TSA_THAL_LOG = term(log10(sum(HCP.MOR.THAL(HCP.QC.indTHAL,:),2)) - mean(log10(sum(HCP.MOR.THAL(HCP.QC.indTHAL,:),2))));

TSA_HIPP = term(sum(HCP.MOR.HIPP(HCP.QC.indHIPP,:),2) - mean(sum(HCP.MOR.HIPP(HCP.QC.indHIPP,:),2)));
TSA_HIPP_LOG = term(log10(sum(HCP.MOR.HIPP(HCP.QC.indHIPP,:),2)) - mean(log10(sum(HCP.MOR.HIPP(HCP.QC.indHIPP,:),2))));

TSA_AMYG = term(sum(HCP.MOR.AMYG(HCP.QC.indAMYG,:),2) - mean(sum(HCP.MOR.AMYG(HCP.QC.indAMYG,:),2)));
TSA_AMYG_LOG = term(log10(sum(HCP.MOR.AMYG(HCP.QC.indAMYG,:),2)) - mean(log10(sum(HCP.MOR.AMYG(HCP.QC.indAMYG,:),2))));

% Define GLMs (struct)
HCP.STAT.VOL.RAW.M_STR = 1 + SEX_STR + AGE_STR;
HCP.STAT.VOL.RAW.M_GP = 1 + SEX_GP + AGE_GP;
HCP.STAT.VOL.RAW.M_THAL = 1 + SEX_THAL + AGE_THAL;
HCP.STAT.VOL.RAW.M_HIPP = 1 + SEX_HIPP + AGE_HIPP;
HCP.STAT.VOL.RAW.M_AMYG = 1 + SEX_AMYG + AGE_AMYG;

HCP.STAT.VOL.TBV.M_STR = 1 + SEX_STR + AGE_STR + TBV_STR;
HCP.STAT.VOL.TBV.M_GP = 1 + SEX_GP + AGE_GP + TBV_GP;
HCP.STAT.VOL.TBV.M_THAL = 1 + SEX_THAL + AGE_THAL + TBV_THAL;
HCP.STAT.VOL.TBV.M_HIPP = 1 + SEX_HIPP + AGE_HIPP + TBV_HIPP;
HCP.STAT.VOL.TBV.M_AMYG = 1 + SEX_AMYG + AGE_AMYG + TBV_AMYG;

HCP.STAT.VOL.LOG.M_STR = 1 + SEX_STR + AGE_STR + TBV_STR_LOG;
HCP.STAT.VOL.LOG.M_GP = 1 + SEX_GP + AGE_GP + TBV_GP_LOG;
HCP.STAT.VOL.LOG.M_THAL = 1 + SEX_THAL + AGE_THAL + TBV_THAL_LOG;
HCP.STAT.VOL.LOG.M_HIPP = 1 + SEX_HIPP + AGE_HIPP + TBV_HIPP_LOG;
HCP.STAT.VOL.LOG.M_AMYG = 1 + SEX_AMYG + AGE_AMYG + TBV_AMYG_LOG;

% Define GLMs (morph)
HCP.STAT.MOR.RAW.M_STR = 1 + SEX_STR + AGE_STR;
HCP.STAT.MOR.RAW.M_GP = 1 + SEX_GP + AGE_GP;
HCP.STAT.MOR.RAW.M_THAL = 1 + SEX_THAL + AGE_THAL;
HCP.STAT.MOR.RAW.M_HIPP = 1 + SEX_HIPP + AGE_HIPP;
HCP.STAT.MOR.RAW.M_AMYG = 1 + SEX_AMYG + AGE_AMYG;

HCP.STAT.MOR.TSA.M_STR = 1 + SEX_STR + AGE_STR + TSA_STR;
HCP.STAT.MOR.TSA.M_GP = 1 + SEX_GP + AGE_GP + TSA_GP;
HCP.STAT.MOR.TSA.M_THAL = 1 + SEX_THAL + AGE_THAL + TSA_THAL;
HCP.STAT.MOR.TSA.M_HIPP = 1 + SEX_HIPP + AGE_HIPP + TSA_HIPP;
HCP.STAT.MOR.TSA.M_AMYG = 1 + SEX_AMYG + AGE_AMYG + TSA_AMYG;

HCP.STAT.MOR.LOG.M_STR = 1 + SEX_STR + AGE_STR + TSA_STR_LOG;
HCP.STAT.MOR.LOG.M_GP = 1 + SEX_GP + AGE_GP + TSA_GP_LOG;
HCP.STAT.MOR.LOG.M_THAL = 1 + SEX_THAL + AGE_THAL + TSA_THAL_LOG;
HCP.STAT.MOR.LOG.M_HIPP = 1 + SEX_HIPP + AGE_HIPP + TSA_HIPP_LOG;
HCP.STAT.MOR.LOG.M_AMYG = 1 + SEX_AMYG + AGE_AMYG + TSA_AMYG_LOG;

%% Compute statistical parameters (TBV | TSA correction)

% Fit models (struct)
HCP.STAT.VOL.TBV.STR_L = SurfStatLinMod(STR_L, HCP.STAT.VOL.TBV.M_STR);
HCP.STAT.VOL.TBV.STR_L_SEX = SurfStatT(HCP.STAT.VOL.TBV.STR_L, SEX_STR);
STR_L_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.VOL.TBV.STR_L_SEX.t), HCP.STAT.VOL.TBV.STR_L_SEX.df));

HCP.STAT.VOL.TBV.STR_R = SurfStatLinMod(STR_R, HCP.STAT.VOL.TBV.M_STR);
HCP.STAT.VOL.TBV.STR_R_SEX = SurfStatT(HCP.STAT.VOL.TBV.STR_R, SEX_STR);
STR_R_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.VOL.TBV.STR_R_SEX.t), HCP.STAT.VOL.TBV.STR_R_SEX.df));

HCP.STAT.VOL.TBV.GP_L = SurfStatLinMod(GP_L, HCP.STAT.VOL.TBV.M_GP);
HCP.STAT.VOL.TBV.GP_L_SEX = SurfStatT(HCP.STAT.VOL.TBV.GP_L, SEX_GP);
GP_L_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.VOL.TBV.GP_L_SEX.t), HCP.STAT.VOL.TBV.GP_L_SEX.df));

HCP.STAT.VOL.TBV.GP_R = SurfStatLinMod(GP_R, HCP.STAT.VOL.TBV.M_GP);
HCP.STAT.VOL.TBV.GP_R_SEX = SurfStatT(HCP.STAT.VOL.TBV.GP_R, SEX_GP);
GP_R_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.VOL.TBV.GP_R_SEX.t), HCP.STAT.VOL.TBV.GP_R_SEX.df));

HCP.STAT.VOL.TBV.THAL_L = SurfStatLinMod(THAL_L, HCP.STAT.VOL.TBV.M_THAL);
HCP.STAT.VOL.TBV.THAL_L_SEX = SurfStatT(HCP.STAT.VOL.TBV.THAL_L, SEX_THAL);
THAL_L_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.VOL.TBV.THAL_L_SEX.t), HCP.STAT.VOL.TBV.THAL_L_SEX.df));

HCP.STAT.VOL.TBV.THAL_R = SurfStatLinMod(THAL_R, HCP.STAT.VOL.TBV.M_THAL);
HCP.STAT.VOL.TBV.THAL_R_SEX = SurfStatT(HCP.STAT.VOL.TBV.THAL_R, SEX_THAL);
THAL_R_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.VOL.TBV.THAL_R_SEX.t), HCP.STAT.VOL.TBV.THAL_R_SEX.df));

HCP.STAT.VOL.TBV.HIPP_L = SurfStatLinMod(HIPP_L, HCP.STAT.VOL.TBV.M_HIPP);
HCP.STAT.VOL.TBV.HIPP_L_SEX = SurfStatT(HCP.STAT.VOL.TBV.HIPP_L, SEX_HIPP);
HIPP_L_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.VOL.TBV.HIPP_L_SEX.t), HCP.STAT.VOL.TBV.HIPP_L_SEX.df));

HCP.STAT.VOL.TBV.HIPP_R = SurfStatLinMod(HIPP_R, HCP.STAT.VOL.TBV.M_HIPP);
HCP.STAT.VOL.TBV.HIPP_R_SEX = SurfStatT(HCP.STAT.VOL.TBV.HIPP_R, SEX_HIPP);
HIPP_R_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.VOL.TBV.HIPP_R_SEX.t), HCP.STAT.VOL.TBV.HIPP_R_SEX.df));

HCP.STAT.VOL.TBV.AMYG_L = SurfStatLinMod(AMYG_L, HCP.STAT.VOL.TBV.M_AMYG);
HCP.STAT.VOL.TBV.AMYG_L_SEX = SurfStatT(HCP.STAT.VOL.TBV.AMYG_L, SEX_AMYG);
AMYG_L_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.VOL.TBV.AMYG_L_SEX.t), HCP.STAT.VOL.TBV.AMYG_L_SEX.df));

HCP.STAT.VOL.TBV.AMYG_R = SurfStatLinMod(AMYG_R, HCP.STAT.VOL.TBV.M_AMYG);
HCP.STAT.VOL.TBV.AMYG_R_SEX = SurfStatT(HCP.STAT.VOL.TBV.AMYG_R, SEX_AMYG);
AMYG_R_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.VOL.TBV.AMYG_R_SEX.t), HCP.STAT.VOL.TBV.AMYG_R_SEX.df));

% Collate statistical values
HCP.STAT.VOL.TBV.SEX_B = [HCP.STAT.VOL.TBV.STR_L_SEX.coef(2) HCP.STAT.VOL.TBV.GP_L_SEX.coef(2) HCP.STAT.VOL.TBV.THAL_L_SEX.coef(2) HCP.STAT.VOL.TBV.HIPP_L_SEX.coef(2) HCP.STAT.VOL.TBV.AMYG_L_SEX.coef(2) ...
                          HCP.STAT.VOL.TBV.STR_R_SEX.coef(2) HCP.STAT.VOL.TBV.GP_R_SEX.coef(2) HCP.STAT.VOL.TBV.THAL_R_SEX.coef(2) HCP.STAT.VOL.TBV.HIPP_R_SEX.coef(2) HCP.STAT.VOL.TBV.AMYG_R_SEX.coef(2)];

HCP.STAT.VOL.TBV.SEX_P = [STR_L_SEX_P GP_L_SEX_P THAL_L_SEX_P HIPP_L_SEX_P AMYG_L_SEX_P ...
                          STR_R_SEX_P GP_R_SEX_P THAL_R_SEX_P HIPP_R_SEX_P AMYG_R_SEX_P];

HCP.STAT.VOL.TBV.SEX_Q = 10*HCP.STAT.VOL.TBV.SEX_P;

% Fit models (morph)
HCP.STAT.MOR.TSA.STR = SurfStatLinMod(STR_MOR, HCP.STAT.MOR.TSA.M_STR);
HCP.STAT.MOR.TSA.STR_SEX = SurfStatT(HCP.STAT.MOR.TSA.STR, SEX_STR);
HCP.STAT.MOR.TSA.STR_SEX_B = HCP.STAT.MOR.TSA.STR_SEX.coef(2,:);
HCP.STAT.MOR.TSA.STR_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.MOR.TSA.STR_SEX.t), HCP.STAT.MOR.TSA.STR_SEX.df));
[~,~,~,HCP.STAT.MOR.TSA.STR_SEX_Q] = fdr_bh(HCP.STAT.MOR.TSA.STR_SEX_P);

HCP.STAT.MOR.TSA.GP = SurfStatLinMod(GP_MOR, HCP.STAT.MOR.TSA.M_GP);
HCP.STAT.MOR.TSA.GP_SEX = SurfStatT(HCP.STAT.MOR.TSA.GP, SEX_GP);
HCP.STAT.MOR.TSA.GP_SEX_B = HCP.STAT.MOR.TSA.GP_SEX.coef(2,:);
HCP.STAT.MOR.TSA.GP_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.MOR.TSA.GP_SEX.t), HCP.STAT.MOR.TSA.GP_SEX.df));
[~,~,~,HCP.STAT.MOR.TSA.GP_SEX_Q] = fdr_bh(HCP.STAT.MOR.TSA.GP_SEX_P);

HCP.STAT.MOR.TSA.THAL = SurfStatLinMod(THAL_MOR, HCP.STAT.MOR.TSA.M_THAL);
HCP.STAT.MOR.TSA.THAL_SEX = SurfStatT(HCP.STAT.MOR.TSA.THAL, SEX_THAL);
HCP.STAT.MOR.TSA.THAL_SEX_B = HCP.STAT.MOR.TSA.THAL_SEX.coef(2,:);
HCP.STAT.MOR.TSA.THAL_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.MOR.TSA.THAL_SEX.t), HCP.STAT.MOR.TSA.THAL_SEX.df));
[~,~,~,HCP.STAT.MOR.TSA.THAL_SEX_Q] = fdr_bh(HCP.STAT.MOR.TSA.THAL_SEX_P);

HCP.STAT.MOR.TSA.HIPP = SurfStatLinMod(HIPP_MOR, HCP.STAT.MOR.TSA.M_HIPP);
HCP.STAT.MOR.TSA.HIPP_SEX = SurfStatT(HCP.STAT.MOR.TSA.HIPP, SEX_HIPP);
HCP.STAT.MOR.TSA.HIPP_SEX_B = HCP.STAT.MOR.TSA.HIPP_SEX.coef(2,:);
HCP.STAT.MOR.TSA.HIPP_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.MOR.TSA.HIPP_SEX.t), HCP.STAT.MOR.TSA.HIPP_SEX.df));
[~,~,~,HCP.STAT.MOR.TSA.HIPP_SEX_Q] = fdr_bh(HCP.STAT.MOR.TSA.HIPP_SEX_P);

HCP.STAT.MOR.TSA.AMYG = SurfStatLinMod(AMYG_MOR, HCP.STAT.MOR.TSA.M_AMYG);
HCP.STAT.MOR.TSA.AMYG_SEX = SurfStatT(HCP.STAT.MOR.TSA.AMYG, SEX_AMYG);
HCP.STAT.MOR.TSA.AMYG_SEX_B = HCP.STAT.MOR.TSA.AMYG_SEX.coef(2,:);
HCP.STAT.MOR.TSA.AMYG_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.MOR.TSA.AMYG_SEX.t), HCP.STAT.MOR.TSA.AMYG_SEX.df));
[~,~,~,HCP.STAT.MOR.TSA.AMYG_SEX_Q] = fdr_bh(HCP.STAT.MOR.TSA.AMYG_SEX_P);

%% Compute statistical parameters (raw)

% Fit models (struct)
HCP.STAT.VOL.RAW.STR_L = SurfStatLinMod(STR_L, HCP.STAT.VOL.RAW.M_STR);
HCP.STAT.VOL.RAW.STR_L_SEX = SurfStatT(HCP.STAT.VOL.RAW.STR_L, SEX_STR);
STR_L_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.VOL.RAW.STR_L_SEX.t), HCP.STAT.VOL.RAW.STR_L_SEX.df));

HCP.STAT.VOL.RAW.STR_R = SurfStatLinMod(STR_R, HCP.STAT.VOL.RAW.M_STR);
HCP.STAT.VOL.RAW.STR_R_SEX = SurfStatT(HCP.STAT.VOL.RAW.STR_R, SEX_STR);
STR_R_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.VOL.RAW.STR_R_SEX.t), HCP.STAT.VOL.RAW.STR_R_SEX.df));

HCP.STAT.VOL.RAW.GP_L = SurfStatLinMod(GP_L, HCP.STAT.VOL.RAW.M_GP);
HCP.STAT.VOL.RAW.GP_L_SEX = SurfStatT(HCP.STAT.VOL.RAW.GP_L, SEX_GP);
GP_L_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.VOL.RAW.GP_L_SEX.t), HCP.STAT.VOL.RAW.GP_L_SEX.df));

HCP.STAT.VOL.RAW.GP_R = SurfStatLinMod(GP_R, HCP.STAT.VOL.RAW.M_GP);
HCP.STAT.VOL.RAW.GP_R_SEX = SurfStatT(HCP.STAT.VOL.RAW.GP_R, SEX_GP);
GP_R_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.VOL.RAW.GP_R_SEX.t), HCP.STAT.VOL.RAW.GP_R_SEX.df));

HCP.STAT.VOL.RAW.THAL_L = SurfStatLinMod(THAL_L, HCP.STAT.VOL.RAW.M_THAL);
HCP.STAT.VOL.RAW.THAL_L_SEX = SurfStatT(HCP.STAT.VOL.RAW.THAL_L, SEX_THAL);
THAL_L_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.VOL.RAW.THAL_L_SEX.t), HCP.STAT.VOL.RAW.THAL_L_SEX.df));

HCP.STAT.VOL.RAW.THAL_R = SurfStatLinMod(THAL_R, HCP.STAT.VOL.RAW.M_THAL);
HCP.STAT.VOL.RAW.THAL_R_SEX = SurfStatT(HCP.STAT.VOL.RAW.THAL_R, SEX_THAL);
THAL_R_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.VOL.RAW.THAL_R_SEX.t), HCP.STAT.VOL.RAW.THAL_R_SEX.df));

HCP.STAT.VOL.RAW.HIPP_L = SurfStatLinMod(HIPP_L, HCP.STAT.VOL.RAW.M_HIPP);
HCP.STAT.VOL.RAW.HIPP_L_SEX = SurfStatT(HCP.STAT.VOL.RAW.HIPP_L, SEX_HIPP);
HIPP_L_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.VOL.RAW.HIPP_L_SEX.t), HCP.STAT.VOL.RAW.HIPP_L_SEX.df));

HCP.STAT.VOL.RAW.HIPP_R = SurfStatLinMod(HIPP_R, HCP.STAT.VOL.RAW.M_HIPP);
HCP.STAT.VOL.RAW.HIPP_R_SEX = SurfStatT(HCP.STAT.VOL.RAW.HIPP_R, SEX_HIPP);
HIPP_R_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.VOL.RAW.HIPP_R_SEX.t), HCP.STAT.VOL.RAW.HIPP_R_SEX.df));

HCP.STAT.VOL.RAW.AMYG_L = SurfStatLinMod(AMYG_L, HCP.STAT.VOL.RAW.M_AMYG);
HCP.STAT.VOL.RAW.AMYG_L_SEX = SurfStatT(HCP.STAT.VOL.RAW.AMYG_L, SEX_AMYG);
AMYG_L_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.VOL.RAW.AMYG_L_SEX.t), HCP.STAT.VOL.RAW.AMYG_L_SEX.df));

HCP.STAT.VOL.RAW.AMYG_R = SurfStatLinMod(AMYG_R, HCP.STAT.VOL.RAW.M_AMYG);
HCP.STAT.VOL.RAW.AMYG_R_SEX = SurfStatT(HCP.STAT.VOL.RAW.AMYG_R, SEX_AMYG);
AMYG_R_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.VOL.RAW.AMYG_R_SEX.t), HCP.STAT.VOL.RAW.AMYG_R_SEX.df));

% Collate statistical values
HCP.STAT.VOL.RAW.SEX_B = [HCP.STAT.VOL.RAW.STR_L_SEX.coef(2) HCP.STAT.VOL.RAW.GP_L_SEX.coef(2) HCP.STAT.VOL.RAW.THAL_L_SEX.coef(2) HCP.STAT.VOL.RAW.HIPP_L_SEX.coef(2) HCP.STAT.VOL.RAW.AMYG_L_SEX.coef(2) ...
                          HCP.STAT.VOL.RAW.STR_R_SEX.coef(2) HCP.STAT.VOL.RAW.GP_R_SEX.coef(2) HCP.STAT.VOL.RAW.THAL_R_SEX.coef(2) HCP.STAT.VOL.RAW.HIPP_R_SEX.coef(2) HCP.STAT.VOL.RAW.AMYG_R_SEX.coef(2)];

HCP.STAT.VOL.RAW.SEX_P = [STR_L_SEX_P GP_L_SEX_P THAL_L_SEX_P HIPP_L_SEX_P AMYG_L_SEX_P ...
                          STR_R_SEX_P GP_R_SEX_P THAL_R_SEX_P HIPP_R_SEX_P AMYG_R_SEX_P];

HCP.STAT.VOL.RAW.SEX_Q = 10*HCP.STAT.VOL.RAW.SEX_P;

% Fit models (morph)
HCP.STAT.MOR.RAW.STR = SurfStatLinMod(STR_MOR, HCP.STAT.MOR.RAW.M_STR);
HCP.STAT.MOR.RAW.STR_SEX = SurfStatT(HCP.STAT.MOR.RAW.STR, SEX_STR);
HCP.STAT.MOR.RAW.STR_SEX_B = HCP.STAT.MOR.RAW.STR_SEX.coef(2,:);
HCP.STAT.MOR.RAW.STR_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.MOR.RAW.STR_SEX.t), HCP.STAT.MOR.RAW.STR_SEX.df));
[~,~,~,HCP.STAT.MOR.RAW.STR_SEX_Q] = fdr_bh(HCP.STAT.MOR.RAW.STR_SEX_P);

HCP.STAT.MOR.RAW.GP = SurfStatLinMod(GP_MOR, HCP.STAT.MOR.RAW.M_GP);
HCP.STAT.MOR.RAW.GP_SEX = SurfStatT(HCP.STAT.MOR.RAW.GP, SEX_GP);
HCP.STAT.MOR.RAW.GP_SEX_B = HCP.STAT.MOR.RAW.GP_SEX.coef(2,:);
HCP.STAT.MOR.RAW.GP_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.MOR.RAW.GP_SEX.t), HCP.STAT.MOR.RAW.GP_SEX.df));
[~,~,~,HCP.STAT.MOR.RAW.GP_SEX_Q] = fdr_bh(HCP.STAT.MOR.RAW.GP_SEX_P);

HCP.STAT.MOR.RAW.THAL = SurfStatLinMod(THAL_MOR, HCP.STAT.MOR.RAW.M_THAL);
HCP.STAT.MOR.RAW.THAL_SEX = SurfStatT(HCP.STAT.MOR.RAW.THAL, SEX_THAL);
HCP.STAT.MOR.RAW.THAL_SEX_B = HCP.STAT.MOR.RAW.THAL_SEX.coef(2,:);
HCP.STAT.MOR.RAW.THAL_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.MOR.RAW.THAL_SEX.t), HCP.STAT.MOR.RAW.THAL_SEX.df));
[~,~,~,HCP.STAT.MOR.RAW.THAL_SEX_Q] = fdr_bh(HCP.STAT.MOR.RAW.THAL_SEX_P);

HCP.STAT.MOR.RAW.HIPP = SurfStatLinMod(HIPP_MOR, HCP.STAT.MOR.RAW.M_HIPP);
HCP.STAT.MOR.RAW.HIPP_SEX = SurfStatT(HCP.STAT.MOR.RAW.HIPP, SEX_HIPP);
HCP.STAT.MOR.RAW.HIPP_SEX_B = HCP.STAT.MOR.RAW.HIPP_SEX.coef(2,:);
HCP.STAT.MOR.RAW.HIPP_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.MOR.RAW.HIPP_SEX.t), HCP.STAT.MOR.RAW.HIPP_SEX.df));
[~,~,~,HCP.STAT.MOR.RAW.HIPP_SEX_Q] = fdr_bh(HCP.STAT.MOR.RAW.HIPP_SEX_P);

HCP.STAT.MOR.RAW.AMYG = SurfStatLinMod(AMYG_MOR, HCP.STAT.MOR.RAW.M_AMYG);
HCP.STAT.MOR.RAW.AMYG_SEX = SurfStatT(HCP.STAT.MOR.RAW.AMYG, SEX_AMYG);
HCP.STAT.MOR.RAW.AMYG_SEX_B = HCP.STAT.MOR.RAW.AMYG_SEX.coef(2,:);
HCP.STAT.MOR.RAW.AMYG_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.MOR.RAW.AMYG_SEX.t), HCP.STAT.MOR.RAW.AMYG_SEX.df));
[~,~,~,HCP.STAT.MOR.RAW.AMYG_SEX_Q] = fdr_bh(HCP.STAT.MOR.RAW.AMYG_SEX_P);

%% Compute statistical parameters (allometric correction)

% Fit models (struct)
HCP.STAT.VOL.LOG.STR_L = SurfStatLinMod(STR_L_LOG, HCP.STAT.VOL.LOG.M_STR);
HCP.STAT.VOL.LOG.STR_L_SEX = SurfStatT(HCP.STAT.VOL.LOG.STR_L, SEX_STR);
STR_L_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.VOL.LOG.STR_L_SEX.t), HCP.STAT.VOL.LOG.STR_L_SEX.df));

HCP.STAT.VOL.LOG.STR_R = SurfStatLinMod(STR_R_LOG, HCP.STAT.VOL.LOG.M_STR);
HCP.STAT.VOL.LOG.STR_R_SEX = SurfStatT(HCP.STAT.VOL.LOG.STR_R, SEX_STR);
STR_R_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.VOL.LOG.STR_R_SEX.t), HCP.STAT.VOL.LOG.STR_R_SEX.df));

HCP.STAT.VOL.LOG.GP_L = SurfStatLinMod(GP_L_LOG, HCP.STAT.VOL.LOG.M_GP);
HCP.STAT.VOL.LOG.GP_L_SEX = SurfStatT(HCP.STAT.VOL.LOG.GP_L, SEX_GP);
GP_L_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.VOL.LOG.GP_L_SEX.t), HCP.STAT.VOL.LOG.GP_L_SEX.df));

HCP.STAT.VOL.LOG.GP_R = SurfStatLinMod(GP_R_LOG, HCP.STAT.VOL.LOG.M_GP);
HCP.STAT.VOL.LOG.GP_R_SEX = SurfStatT(HCP.STAT.VOL.LOG.GP_R, SEX_GP);
GP_R_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.VOL.LOG.GP_R_SEX.t), HCP.STAT.VOL.LOG.GP_R_SEX.df));

HCP.STAT.VOL.LOG.THAL_L = SurfStatLinMod(THAL_L_LOG, HCP.STAT.VOL.LOG.M_THAL);
HCP.STAT.VOL.LOG.THAL_L_SEX = SurfStatT(HCP.STAT.VOL.LOG.THAL_L, SEX_THAL);
THAL_L_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.VOL.LOG.THAL_L_SEX.t), HCP.STAT.VOL.LOG.THAL_L_SEX.df));

HCP.STAT.VOL.LOG.THAL_R = SurfStatLinMod(THAL_R_LOG, HCP.STAT.VOL.LOG.M_THAL);
HCP.STAT.VOL.LOG.THAL_R_SEX = SurfStatT(HCP.STAT.VOL.LOG.THAL_R, SEX_THAL);
THAL_R_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.VOL.LOG.THAL_R_SEX.t), HCP.STAT.VOL.LOG.THAL_R_SEX.df));

HCP.STAT.VOL.LOG.HIPP_L = SurfStatLinMod(HIPP_L_LOG, HCP.STAT.VOL.LOG.M_HIPP);
HCP.STAT.VOL.LOG.HIPP_L_SEX = SurfStatT(HCP.STAT.VOL.LOG.HIPP_L, SEX_HIPP);
HIPP_L_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.VOL.LOG.HIPP_L_SEX.t), HCP.STAT.VOL.LOG.HIPP_L_SEX.df));

HCP.STAT.VOL.LOG.HIPP_R = SurfStatLinMod(HIPP_R_LOG, HCP.STAT.VOL.LOG.M_HIPP);
HCP.STAT.VOL.LOG.HIPP_R_SEX = SurfStatT(HCP.STAT.VOL.LOG.HIPP_R, SEX_HIPP);
HIPP_R_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.VOL.LOG.HIPP_R_SEX.t), HCP.STAT.VOL.LOG.HIPP_R_SEX.df));

HCP.STAT.VOL.LOG.AMYG_L = SurfStatLinMod(AMYG_L_LOG, HCP.STAT.VOL.LOG.M_AMYG);
HCP.STAT.VOL.LOG.AMYG_L_SEX = SurfStatT(HCP.STAT.VOL.LOG.AMYG_L, SEX_AMYG);
AMYG_L_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.VOL.LOG.AMYG_L_SEX.t), HCP.STAT.VOL.LOG.AMYG_L_SEX.df));

HCP.STAT.VOL.LOG.AMYG_R = SurfStatLinMod(AMYG_R_LOG, HCP.STAT.VOL.LOG.M_AMYG);
HCP.STAT.VOL.LOG.AMYG_R_SEX = SurfStatT(HCP.STAT.VOL.LOG.AMYG_R, SEX_AMYG);
AMYG_R_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.VOL.LOG.AMYG_R_SEX.t), HCP.STAT.VOL.LOG.AMYG_R_SEX.df));

% Collate statistical values
HCP.STAT.VOL.LOG.SEX_B = [HCP.STAT.VOL.LOG.STR_L_SEX.coef(2) HCP.STAT.VOL.LOG.GP_L_SEX.coef(2) HCP.STAT.VOL.LOG.THAL_L_SEX.coef(2) HCP.STAT.VOL.LOG.HIPP_L_SEX.coef(2) HCP.STAT.VOL.LOG.AMYG_L_SEX.coef(2) ...
                          HCP.STAT.VOL.LOG.STR_R_SEX.coef(2) HCP.STAT.VOL.LOG.GP_R_SEX.coef(2) HCP.STAT.VOL.LOG.THAL_R_SEX.coef(2) HCP.STAT.VOL.LOG.HIPP_R_SEX.coef(2) HCP.STAT.VOL.LOG.AMYG_R_SEX.coef(2)];

HCP.STAT.VOL.LOG.SEX_P = [STR_L_SEX_P GP_L_SEX_P THAL_L_SEX_P HIPP_L_SEX_P AMYG_L_SEX_P ...
                          STR_R_SEX_P GP_R_SEX_P THAL_R_SEX_P HIPP_R_SEX_P AMYG_R_SEX_P];

HCP.STAT.VOL.LOG.SEX_Q = 10*HCP.STAT.VOL.LOG.SEX_P;

% Derive standardized effect size
HCP.STAT.VOL.LOG.SEX_SD_N = [mean(HCP.SEX(HCP.QC.indSTR)) mean(HCP.SEX(HCP.QC.indGP)) mean(HCP.SEX(HCP.QC.indTHAL)) mean(HCP.SEX(HCP.QC.indHIPP)) mean(HCP.SEX(HCP.QC.indAMYG)) ...
                             mean(HCP.SEX(HCP.QC.indSTR)) mean(HCP.SEX(HCP.QC.indGP)) mean(HCP.SEX(HCP.QC.indTHAL)) mean(HCP.SEX(HCP.QC.indHIPP)) mean(HCP.SEX(HCP.QC.indAMYG))];

HCP.STAT.VOL.LOG.SEX_SD_D = [std(log10(HCP.VOL.STR_L(HCP.QC.indSTR))) std(log10(HCP.VOL.GP_L(HCP.QC.indGP))) std(log10(HCP.VOL.THAL_L(HCP.QC.indTHAL))) std(log10(HCP.VOL.HIPP_L(HCP.QC.indHIPP))) std(log10(HCP.VOL.AMYG_L(HCP.QC.indAMYG))) ...
                             std(log10(HCP.VOL.STR_R(HCP.QC.indSTR))) std(log10(HCP.VOL.GP_R(HCP.QC.indGP))) std(log10(HCP.VOL.THAL_R(HCP.QC.indTHAL))) std(log10(HCP.VOL.HIPP_R(HCP.QC.indHIPP))) std(log10(HCP.VOL.AMYG_R(HCP.QC.indAMYG)))];

HCP.STAT.VOL.LOG.SEX_B_SD = HCP.STAT.VOL.LOG.SEX_B .* HCP.STAT.VOL.LOG.SEX_SD_N ./ HCP.STAT.VOL.LOG.SEX_SD_D;

% Fit models (morph)
HCP.STAT.MOR.LOG.STR = SurfStatLinMod(STR_MOR_LOG, HCP.STAT.MOR.LOG.M_STR);
HCP.STAT.MOR.LOG.STR_SEX = SurfStatT(HCP.STAT.MOR.LOG.STR, SEX_STR);
HCP.STAT.MOR.LOG.STR_SEX_B = HCP.STAT.MOR.LOG.STR_SEX.coef(2,:);
HCP.STAT.MOR.LOG.STR_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.MOR.LOG.STR_SEX.t), HCP.STAT.MOR.LOG.STR_SEX.df));
[~,~,~,HCP.STAT.MOR.LOG.STR_SEX_Q] = fdr_bh(HCP.STAT.MOR.LOG.STR_SEX_P);

HCP.STAT.MOR.LOG.GP = SurfStatLinMod(GP_MOR_LOG, HCP.STAT.MOR.LOG.M_GP);
HCP.STAT.MOR.LOG.GP_SEX = SurfStatT(HCP.STAT.MOR.LOG.GP, SEX_GP);
HCP.STAT.MOR.LOG.GP_SEX_B = HCP.STAT.MOR.LOG.GP_SEX.coef(2,:);
HCP.STAT.MOR.LOG.GP_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.MOR.LOG.GP_SEX.t), HCP.STAT.MOR.LOG.GP_SEX.df));
[~,~,~,HCP.STAT.MOR.LOG.GP_SEX_Q] = fdr_bh(HCP.STAT.MOR.LOG.GP_SEX_P);

HCP.STAT.MOR.LOG.THAL = SurfStatLinMod(THAL_MOR_LOG, HCP.STAT.MOR.LOG.M_THAL);
HCP.STAT.MOR.LOG.THAL_SEX = SurfStatT(HCP.STAT.MOR.LOG.THAL, SEX_THAL);
HCP.STAT.MOR.LOG.THAL_SEX_B = HCP.STAT.MOR.LOG.THAL_SEX.coef(2,:);
HCP.STAT.MOR.LOG.THAL_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.MOR.LOG.THAL_SEX.t), HCP.STAT.MOR.LOG.THAL_SEX.df));
[~,~,~,HCP.STAT.MOR.LOG.THAL_SEX_Q] = fdr_bh(HCP.STAT.MOR.LOG.THAL_SEX_P);

HCP.STAT.MOR.LOG.HIPP = SurfStatLinMod(HIPP_MOR_LOG, HCP.STAT.MOR.LOG.M_HIPP);
HCP.STAT.MOR.LOG.HIPP_SEX = SurfStatT(HCP.STAT.MOR.LOG.HIPP, SEX_HIPP);
HCP.STAT.MOR.LOG.HIPP_SEX_B = HCP.STAT.MOR.LOG.HIPP_SEX.coef(2,:);
HCP.STAT.MOR.LOG.HIPP_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.MOR.LOG.HIPP_SEX.t), HCP.STAT.MOR.LOG.HIPP_SEX.df));
[~,~,~,HCP.STAT.MOR.LOG.HIPP_SEX_Q] = fdr_bh(HCP.STAT.MOR.LOG.HIPP_SEX_P);

HCP.STAT.MOR.LOG.AMYG = SurfStatLinMod(AMYG_MOR_LOG, HCP.STAT.MOR.LOG.M_AMYG);
HCP.STAT.MOR.LOG.AMYG_SEX = SurfStatT(HCP.STAT.MOR.LOG.AMYG, SEX_AMYG);
HCP.STAT.MOR.LOG.AMYG_SEX_B = HCP.STAT.MOR.LOG.AMYG_SEX.coef(2,:);
HCP.STAT.MOR.LOG.AMYG_SEX_P = 2*(1 - tcdf(abs(HCP.STAT.MOR.LOG.AMYG_SEX.t), HCP.STAT.MOR.LOG.AMYG_SEX.df));
[~,~,~,HCP.STAT.MOR.LOG.AMYG_SEX_Q] = fdr_bh(HCP.STAT.MOR.LOG.AMYG_SEX_P);

% Derive standardized effect size
HCP.STAT.MOR.LOG.STR_SEX_B_SD = HCP.STAT.MOR.LOG.STR_SEX_B .* mean(HCP.SEX(HCP.QC.indSTR)) ./ std(STR_MOR_LOG);
HCP.STAT.MOR.LOG.GP_SEX_B_SD = HCP.STAT.MOR.LOG.GP_SEX_B .* mean(HCP.SEX(HCP.QC.indGP)) ./ std(GP_MOR_LOG);
HCP.STAT.MOR.LOG.THAL_SEX_B_SD = HCP.STAT.MOR.LOG.THAL_SEX_B .* mean(HCP.SEX(HCP.QC.indTHAL)) ./ std(THAL_MOR_LOG);
HCP.STAT.MOR.LOG.HIPP_SEX_B_SD = HCP.STAT.MOR.LOG.HIPP_SEX_B .* mean(HCP.SEX(HCP.QC.indHIPP)) ./ std(HIPP_MOR_LOG);
HCP.STAT.MOR.LOG.AMYG_SEX_B_SD = HCP.STAT.MOR.LOG.AMYG_SEX_B .* mean(HCP.SEX(HCP.QC.indAMYG)) ./ std(AMYG_MOR_LOG);

%% General null distribution

% Load subcortical templates
load(fullfile(templateDir, 'maget_template_surfaces_all.mat'));

% Compute geodesic distance
S_STR_L_GD = distances(surface_to_graph(S_STR_L, 'geodesic'));
S_STR_R_GD = distances(surface_to_graph(S_STR_R, 'geodesic'));

S_GP_L_GD = distances(surface_to_graph(S_GP_L, 'geodesic'));
S_GP_R_GD = distances(surface_to_graph(S_GP_R, 'geodesic'));

S_THAL_L_GD = distances(surface_to_graph(S_THAL_L, 'geodesic'));
S_THAL_R_GD = distances(surface_to_graph(S_THAL_R, 'geodesic'));

S_HIPP_L_GD = distances(surface_to_graph(S_HIPP_L, 'geodesic'));
S_HIPP_R_GD = distances(surface_to_graph(S_HIPP_R, 'geodesic'));

S_AMYG_L_GD = distances(surface_to_graph(S_AMYG_L, 'geodesic'));
S_AMYG_R_GD = distances(surface_to_graph(S_AMYG_R, 'geodesic'));

% Compute the mean eigenvectors
S_STR_L_MEM = compute_mem(S_STR_L_GD);
S_STR_R_MEM = compute_mem(S_STR_R_GD);

S_GP_L_MEM = compute_mem(S_GP_L_GD);
S_GP_R_MEM = compute_mem(S_GP_R_GD);

S_THAL_L_MEM = compute_mem(S_THAL_L_GD);
S_THAL_R_MEM = compute_mem(S_THAL_R_GD);

S_HIPP_L_MEM = compute_mem(S_HIPP_L_GD);
S_HIPP_R_MEM = compute_mem(S_HIPP_R_GD);

S_AMYG_L_MEM = compute_mem(S_AMYG_L_GD);
S_AMYG_R_MEM = compute_mem(S_AMYG_R_GD);

% Compute Moran spectral randomization
n_moran = 10000;
HCP.STAT.MOR.MSR.STR_SEX_B = [squeeze(moran_randomization(HCP.STAT.MOR.TSA.STR_SEX_B(1:end/2)', S_STR_L_MEM, n_moran, 'procedure', 'singleton', 'joint', false, 'random_state', 1));
                              squeeze(moran_randomization(HCP.STAT.MOR.TSA.STR_SEX_B(end/2+1:end)', S_STR_R_MEM, n_moran, 'procedure', 'singleton', 'joint', false, 'random_state', 1))];
HCP.STAT.MOR.MSR.STR_SEX_Q = [squeeze(moran_randomization(HCP.STAT.MOR.TSA.STR_SEX_Q(1:end/2)' < thresh, S_STR_L_MEM, n_moran, 'procedure', 'singleton', 'joint', false, 'random_state', 1));
                                squeeze(moran_randomization(HCP.STAT.MOR.TSA.STR_SEX_Q(end/2+1:end)' < thresh, S_STR_R_MEM, n_moran, 'procedure', 'singleton', 'joint', false, 'random_state', 1))] >= 0.5;

HCP.STAT.MOR.MSR.GP_SEX_B = [squeeze(moran_randomization(HCP.STAT.MOR.TSA.GP_SEX_B(1:end/2)', S_GP_L_MEM, n_moran, 'procedure', 'singleton', 'joint', false, 'random_state', 1));
                              squeeze(moran_randomization(HCP.STAT.MOR.TSA.GP_SEX_B(end/2+1:end)', S_GP_R_MEM, n_moran, 'procedure', 'singleton', 'joint', false, 'random_state', 1))];
HCP.STAT.MOR.MSR.GP_SEX_Q = [squeeze(moran_randomization(HCP.STAT.MOR.TSA.GP_SEX_Q(1:end/2)' < thresh, S_GP_L_MEM, n_moran, 'procedure', 'singleton', 'joint', false, 'random_state', 1));
                                squeeze(moran_randomization(HCP.STAT.MOR.TSA.GP_SEX_Q(end/2+1:end)' < thresh, S_GP_R_MEM, n_moran, 'procedure', 'singleton', 'joint', false, 'random_state', 1))] >= 0.5;

HCP.STAT.MOR.MSR.THAL_SEX_B = [squeeze(moran_randomization(HCP.STAT.MOR.TSA.THAL_SEX_B(1:end/2)', S_THAL_L_MEM, n_moran, 'procedure', 'singleton', 'joint', false, 'random_state', 1));
                              squeeze(moran_randomization(HCP.STAT.MOR.TSA.THAL_SEX_B(end/2+1:end)', S_THAL_R_MEM, n_moran, 'procedure', 'singleton', 'joint', false, 'random_state', 1))];
HCP.STAT.MOR.MSR.THAL_SEX_Q = [squeeze(moran_randomization(HCP.STAT.MOR.TSA.THAL_SEX_Q(1:end/2)' < thresh, S_THAL_L_MEM, n_moran, 'procedure', 'singleton', 'joint', false, 'random_state', 1));
                                squeeze(moran_randomization(HCP.STAT.MOR.TSA.THAL_SEX_Q(end/2+1:end)' < thresh, S_THAL_R_MEM, n_moran, 'procedure', 'singleton', 'joint', false, 'random_state', 1))] >= 0.5;

HCP.STAT.MOR.MSR.HIPP_SEX_B = [squeeze(moran_randomization(HCP.STAT.MOR.TSA.HIPP_SEX_B(1:1152)', S_HIPP_L_MEM, n_moran, 'procedure', 'singleton', 'joint', false, 'random_state', 1));
                              squeeze(moran_randomization(HCP.STAT.MOR.TSA.HIPP_SEX_B(1153:end)', S_HIPP_R_MEM, n_moran, 'procedure', 'singleton', 'joint', false, 'random_state', 1))];
HCP.STAT.MOR.MSR.HIPP_SEX_Q = [squeeze(moran_randomization(HCP.STAT.MOR.TSA.HIPP_SEX_Q(1:1152)' < thresh, S_HIPP_L_MEM, n_moran, 'procedure', 'singleton', 'joint', false, 'random_state', 1));
                                squeeze(moran_randomization(HCP.STAT.MOR.TSA.HIPP_SEX_Q(1153:end)' < thresh, S_HIPP_R_MEM, n_moran, 'procedure', 'singleton', 'joint', false, 'random_state', 1))] >= 0.5;

HCP.STAT.MOR.MSR.AMYG_SEX_B = [squeeze(moran_randomization(HCP.STAT.MOR.TSA.AMYG_SEX_B(1:1473)', S_AMYG_L_MEM, n_moran, 'procedure', 'singleton', 'joint', false, 'random_state', 1));
                              squeeze(moran_randomization(HCP.STAT.MOR.TSA.AMYG_SEX_B(1474:end)', S_AMYG_R_MEM, n_moran, 'procedure', 'singleton', 'joint', false, 'random_state', 1))];
HCP.STAT.MOR.MSR.AMYG_SEX_Q = [squeeze(moran_randomization(HCP.STAT.MOR.TSA.AMYG_SEX_Q(1:1473)' < thresh, S_AMYG_L_MEM, n_moran, 'procedure', 'singleton', 'joint', false, 'random_state', 1));
                                squeeze(moran_randomization(HCP.STAT.MOR.TSA.AMYG_SEX_Q(1474:end)' < thresh, S_AMYG_R_MEM, n_moran, 'procedure', 'singleton', 'joint', false, 'random_state', 1))] >= 0.5;

%% Append results from the study on cortical sex differences

% Load HCP results from cortical sex differences
load(fullfile(corticalResultsDir, 'hcp_groupmaps_final.mat'));
HCP.STAT.COR.CV = HCP_CV_SEX;
HCP.STAT.COR.CV_B = HCP_CV_SEX.coef(2,:);
HCP.STAT.COR.CV_Q = HCP_CV_SEX_Q;

HCP.STAT.COR.SA = HCP_SA_SEX;
HCP.STAT.COR.SA_B = HCP_SA_SEX.coef(2,:);
HCP.STAT.COR.SA_Q = HCP_SA_SEX_Q;

HCP.STAT.COR.CT = HCP_CT_SEX;
HCP.STAT.COR.CT_B = HCP_CT_SEX.coef(2,:);
HCP.STAT.COR.CT_Q = HCP_CT_SEX_Q;

HCP.QC.indCOR = ones(1,1085);

%% Save results

save(fullfile(dataDir, 'HCP_results.mat'), 'HCP', '-v7.3');

%% Visualize effect sizes

close all;

h = figure('WindowState', 'Maximized'); HyoMagetView([zscore(S_STR_SA); zscore(S_GP_SA); zscore(S_THAL_SA); zscore(S_HIPP_SA); zscore(S_AMYG_SA)], [-2.5 4], cmap.Purples); exportgraphics(h, fullfile(figDir, 'Surface_Area.png'), 'Resolution', 300); close(h);

% FDR thresholded with gray background
temp = HCP.STAT.MOR.TSA.STR_SEX_B .* (HCP.STAT.MOR.TSA.STR_SEX_Q < thresh); h = figure('WindowState', 'Maximized'); HyoMagetViewNL(temp, [-0.5 0.5], ThreshColorMap(temp, cmap.PiYG, 0.915)); exportgraphics(h, fullfile(figDir, 'HCP_STR_TSA_FDR.png'), 'Resolution', 300); close(h);
temp = HCP.STAT.MOR.TSA.GP_SEX_B .* (HCP.STAT.MOR.TSA.GP_SEX_Q < thresh); h = figure('WindowState', 'Maximized'); HyoMagetViewNL(temp, [-0.5 0.5], ThreshColorMap(temp, cmap.PiYG, 0.915)); exportgraphics(h, fullfile(figDir, 'HCP_GP_TSA_FDR.png'), 'Resolution', 300); close(h);
temp = HCP.STAT.MOR.TSA.THAL_SEX_B .* (HCP.STAT.MOR.TSA.THAL_SEX_Q < thresh); h = figure('WindowState', 'Maximized'); HyoMagetViewNL(temp, [-0.5 0.5], ThreshColorMap(temp, cmap.PiYG, 0.915)); exportgraphics(h, fullfile(figDir, 'HCP_THAL_TSA_FDR.png'), 'Resolution', 300); close(h);
temp = HCP.STAT.MOR.TSA.HIPP_SEX_B .* (HCP.STAT.MOR.TSA.HIPP_SEX_Q < thresh); h = figure('WindowState', 'Maximized'); HyoMagetViewNL(temp, [-0.5 0.5], ThreshColorMap(temp, cmap.PiYG, 0.915)); exportgraphics(h, fullfile(figDir, 'HCP_HIPP_TSA_FDR.png'), 'Resolution', 300); close(h);
temp = HCP.STAT.MOR.TSA.AMYG_SEX_B .* (HCP.STAT.MOR.TSA.AMYG_SEX_Q < thresh); h = figure('WindowState', 'Maximized'); HyoMagetViewNL(temp, [-0.5 0.5], ThreshColorMap(temp, cmap.PiYG, 0.915)); exportgraphics(h, fullfile(figDir, 'HCP_AMYG_TSA_FDR.png'), 'Resolution', 300); close(h);

% FDR thresholded (flat view)
h = figure('WindowState', 'Maximized'); HyoMagetViewNLFlat(HCP.STAT.MOR.TSA.STR_SEX_B .* (HCP.STAT.MOR.TSA.STR_SEX_Q < thresh), [-0.5 0.5], cmap.PiYG_thr); exportgraphics(h, fullfile(figDir, 'HCP_STR_TSA_FDR_flat.png'), 'Resolution', 300); close(h);
h = figure('WindowState', 'Maximized'); HyoMagetViewNLFlat(HCP.STAT.MOR.TSA.GP_SEX_B .* (HCP.STAT.MOR.TSA.GP_SEX_Q < thresh), [-0.5 0.5], cmap.PiYG_thr); exportgraphics(h, fullfile(figDir, 'HCP_GP_TSA_FDR_flat.png'), 'Resolution', 300); close(h);
h = figure('WindowState', 'Maximized'); HyoMagetViewNLFlat(HCP.STAT.MOR.TSA.THAL_SEX_B .* (HCP.STAT.MOR.TSA.THAL_SEX_Q < thresh), [-0.5 0.5], cmap.PiYG_thr); exportgraphics(h, fullfile(figDir, 'HCP_THAL_TSA_FDR_flat.png'), 'Resolution', 300); close(h);
h = figure('WindowState', 'Maximized'); HyoMagetViewNLFlat(HCP.STAT.MOR.TSA.HIPP_SEX_B .* (HCP.STAT.MOR.TSA.HIPP_SEX_Q < thresh), [-0.5 0.5], cmap.PiYG_thr); exportgraphics(h, fullfile(figDir, 'HCP_HIPP_TSA_FDR_flat.png'), 'Resolution', 300); close(h);
h = figure('WindowState', 'Maximized'); HyoMagetViewNLFlat(HCP.STAT.MOR.TSA.AMYG_SEX_B .* (HCP.STAT.MOR.TSA.AMYG_SEX_Q < thresh), [-0.5 0.5], cmap.PiYG_thr); exportgraphics(h, fullfile(figDir, 'HCP_AMYG_TSA_FDR_flat.png'), 'Resolution', 300); close(h);

% Not thresholded
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(HCP.STAT.MOR.TSA.STR_SEX_B, [-0.5 0.5], cmap.PiYG); exportgraphics(h, fullfile(figDir, 'HCP_STR_TSA.png'), 'Resolution', 300); close(h);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(HCP.STAT.MOR.TSA.GP_SEX_B, [-0.5 0.5], cmap.PiYG); exportgraphics(h, fullfile(figDir, 'HCP_GP_TSA.png'), 'Resolution', 300); close(h);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(HCP.STAT.MOR.TSA.THAL_SEX_B, [-0.5 0.5], cmap.PiYG); exportgraphics(h, fullfile(figDir, 'HCP_THAL_TSA.png'), 'Resolution', 300); close(h);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(HCP.STAT.MOR.TSA.HIPP_SEX_B, [-0.5 0.5], cmap.PiYG); exportgraphics(h, fullfile(figDir, 'HCP_HIPP_TSA.png'), 'Resolution', 300); close(h);
h = figure('WindowState', 'Maximized'); HyoMagetViewNL(HCP.STAT.MOR.TSA.AMYG_SEX_B, [-0.5 0.5], cmap.PiYG); exportgraphics(h, fullfile(figDir, 'HCP_AMYG_TSA.png'), 'Resolution', 300); close(h);

% Display the areal proportion of significance
sum(S_STR_NA(HCP.STAT.MOR.TSA.STR_SEX_Q < thresh))
sum(S_GP_NA(HCP.STAT.MOR.TSA.GP_SEX_Q < thresh))
sum(S_THAL_NA(HCP.STAT.MOR.TSA.THAL_SEX_Q < thresh))
sum(S_HIPP_NA(HCP.STAT.MOR.TSA.HIPP_SEX_Q < thresh))
sum(S_AMYG_NA(HCP.STAT.MOR.TSA.AMYG_SEX_Q < thresh))

% Calculate the colorscale of thresholded maps
temp = HCP.STAT.MOR.TSA.STR_SEX_B .* (HCP.STAT.MOR.TSA.STR_SEX_Q < thresh);
[min(abs(temp(temp<0))) max(abs(temp(temp<0))) min(abs(temp(temp>0))) max(abs(temp(temp>0)))]

temp = HCP.STAT.MOR.TSA.GP_SEX_B .* (HCP.STAT.MOR.TSA.GP_SEX_Q < thresh);
[min(abs(temp(temp<0))) max(abs(temp(temp<0))) min(abs(temp(temp>0))) max(abs(temp(temp>0)))]

temp = HCP.STAT.MOR.TSA.THAL_SEX_B .* (HCP.STAT.MOR.TSA.THAL_SEX_Q < thresh);
[min(abs(temp(temp<0))) max(abs(temp(temp<0))) min(abs(temp(temp>0))) max(abs(temp(temp>0)))]

temp = HCP.STAT.MOR.TSA.HIPP_SEX_B .* (HCP.STAT.MOR.TSA.HIPP_SEX_Q < thresh);
[min(abs(temp(temp<0))) max(abs(temp(temp<0))) min(abs(temp(temp>0))) max(abs(temp(temp>0)))]

temp = HCP.STAT.MOR.TSA.AMYG_SEX_B .* (HCP.STAT.MOR.TSA.AMYG_SEX_Q < thresh);
[min(abs(temp(temp<0))) max(abs(temp(temp<0))) min(abs(temp(temp>0))) max(abs(temp(temp>0)))]

%% Rotating full view

% Collate into one vector
SUB_B = [HCP.STAT.MOR.TSA.STR_SEX_B HCP.STAT.MOR.TSA.GP_SEX_B HCP.STAT.MOR.TSA.THAL_SEX_B HCP.STAT.MOR.TSA.HIPP_SEX_B HCP.STAT.MOR.TSA.AMYG_SEX_B];
SUB_Q = [HCP.STAT.MOR.TSA.STR_SEX_Q HCP.STAT.MOR.TSA.GP_SEX_Q HCP.STAT.MOR.TSA.THAL_SEX_Q HCP.STAT.MOR.TSA.HIPP_SEX_Q HCP.STAT.MOR.TSA.AMYG_SEX_Q];

% Full video
% c = 1;
% dA = 0.2;
% for i = 0:dA:360-dA
%     disp(['Rotation angle ' num2str(i) ' deg']);
%     h = figure('WindowState', 'Maximized');
%     HyoMagetViewSpin(SUB_B .* (SUB_Q < thresh), [-0.5 0.5], i, cmap.PiYG);
%     print(h, fullfile(figDir, 'revolving', ['HCP_TSA_ANG-' num2str(c) '.png']), '-dpng', '-r300');
%     close(h);
%     c = c + 1;
% end

% Static
c = 1;
dA = 45;
for i =  0:dA:360-dA
    disp(['Rotation angle ' num2str(i) ' deg']);
    h = figure('WindowState', 'Maximized');
    HyoMagetViewSpin(SUB_B .* (SUB_Q < thresh), [-0.5 0.5], i, ThreshColorMap(SUB_B .* (SUB_Q < thresh), cmap.PiYG, [-0.5 0.5], 0.915), 'static');
    print(h, fullfile(figDir, 'static_hcp', ['HCP_TSA_ANG-' num2str(c) '.png']), '-dpng', '-r300');
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
    print(h, fullfile(figDir, 'static_hcp', ['NOTHR_HCP_TSA_ANG-' num2str(c) '.png']), '-dpng', '-r300');
    close(h);
    c = c + 1;
end

%% Raincloud plot for bulk sex differences

% Extract absolute volume
STR = ((HCP.VOL.STR_L(HCP.QC.indSTR) + HCP.VOL.STR_R(HCP.QC.indSTR))');
GP = ((HCP.VOL.GP_L(HCP.QC.indGP) + HCP.VOL.GP_R(HCP.QC.indGP))');
THAL = ((HCP.VOL.THAL_L(HCP.QC.indTHAL) + HCP.VOL.THAL_R(HCP.QC.indTHAL))');
HIPP = ((HCP.VOL.HIPP_L(HCP.QC.indHIPP) + HCP.VOL.HIPP_R(HCP.QC.indHIPP))');
AMYG = ((HCP.VOL.AMYG_L(HCP.QC.indAMYG) + HCP.VOL.AMYG_R(HCP.QC.indAMYG))');

% Mean-center covariates
SEX_STR = HCP.SEX(HCP.QC.indSTR);
AGE_STR = HCP.AGE(HCP.QC.indSTR) - mean(HCP.AGE(HCP.QC.indSTR));
TBV_STR = HCP.TBV(HCP.QC.indSTR) - mean(HCP.TBV(HCP.QC.indSTR));
SEX_GP = HCP.SEX(HCP.QC.indGP);
AGE_GP = HCP.AGE(HCP.QC.indGP) - mean(HCP.AGE(HCP.QC.indGP));
TBV_GP = HCP.TBV(HCP.QC.indGP) - mean(HCP.TBV(HCP.QC.indGP));
SEX_THAL = HCP.SEX(HCP.QC.indTHAL);
AGE_THAL = HCP.AGE(HCP.QC.indTHAL) - mean(HCP.AGE(HCP.QC.indTHAL));
TBV_THAL = HCP.TBV(HCP.QC.indTHAL) - mean(HCP.TBV(HCP.QC.indTHAL));
SEX_HIPP = HCP.SEX(HCP.QC.indHIPP);
AGE_HIPP = HCP.AGE(HCP.QC.indHIPP) - mean(HCP.AGE(HCP.QC.indHIPP));
TBV_HIPP = HCP.TBV(HCP.QC.indHIPP) - mean(HCP.TBV(HCP.QC.indHIPP));
SEX_AMYG = HCP.SEX(HCP.QC.indAMYG);
AGE_AMYG = HCP.AGE(HCP.QC.indAMYG) - mean(HCP.AGE(HCP.QC.indAMYG));
TBV_AMYG = HCP.TBV(HCP.QC.indAMYG) - mean(HCP.TBV(HCP.QC.indAMYG));

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

% % Display mean and SD
% [mean(STR_F) std(STR_F) mean(STR_M) std(STR_M)]
% [mean(GP_F) std(GP_F) mean(GP_M) std(GP_M)]
% [mean(THAL_F) std(THAL_F) mean(THAL_M) std(THAL_M)]
% [mean(HIPP_F) std(HIPP_F) mean(HIPP_M) std(HIPP_M)]
% [mean(AMYG_F) std(AMYG_F) mean(AMYG_M) std(AMYG_M)]

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

%% Scatter between correction methods (raw, tbv and allo)

% Calculate correlations
R_STR_TSA_RAW = corr(HCP.STAT.MOR.TSA.STR_SEX_B', HCP.STAT.MOR.RAW.STR_SEX_B');
R_STR_TSA_LOG = corr(HCP.STAT.MOR.TSA.STR_SEX_B', HCP.STAT.MOR.LOG.STR_SEX_B');
R_STR_LOG_RAW = corr(HCP.STAT.MOR.LOG.STR_SEX_B', HCP.STAT.MOR.RAW.STR_SEX_B');
[R_STR_TSA_RAW R_STR_TSA_LOG R_STR_LOG_RAW]

R_GP_TSA_RAW = corr(HCP.STAT.MOR.TSA.GP_SEX_B', HCP.STAT.MOR.RAW.GP_SEX_B');
R_GP_TSA_LOG = corr(HCP.STAT.MOR.TSA.GP_SEX_B', HCP.STAT.MOR.LOG.GP_SEX_B');
R_GP_LOG_RAW = corr(HCP.STAT.MOR.LOG.GP_SEX_B', HCP.STAT.MOR.RAW.GP_SEX_B');
[R_GP_TSA_RAW R_GP_TSA_LOG R_GP_LOG_RAW]

R_THAL_TSA_RAW = corr(HCP.STAT.MOR.TSA.THAL_SEX_B', HCP.STAT.MOR.RAW.THAL_SEX_B');
R_THAL_TSA_LOG = corr(HCP.STAT.MOR.TSA.THAL_SEX_B', HCP.STAT.MOR.LOG.THAL_SEX_B');
R_THAL_LOG_RAW = corr(HCP.STAT.MOR.LOG.THAL_SEX_B', HCP.STAT.MOR.RAW.THAL_SEX_B');
[R_THAL_TSA_RAW R_THAL_TSA_LOG R_THAL_LOG_RAW]

R_HIPP_TSA_RAW = corr(HCP.STAT.MOR.TSA.HIPP_SEX_B', HCP.STAT.MOR.RAW.HIPP_SEX_B');
R_HIPP_TSA_LOG = corr(HCP.STAT.MOR.TSA.HIPP_SEX_B', HCP.STAT.MOR.LOG.HIPP_SEX_B');
R_HIPP_LOG_RAW = corr(HCP.STAT.MOR.LOG.HIPP_SEX_B', HCP.STAT.MOR.RAW.HIPP_SEX_B');
[R_HIPP_TSA_RAW R_HIPP_TSA_LOG R_HIPP_LOG_RAW]

R_AMYG_TSA_RAW = corr(HCP.STAT.MOR.TSA.AMYG_SEX_B', HCP.STAT.MOR.RAW.AMYG_SEX_B');
R_AMYG_TSA_LOG = corr(HCP.STAT.MOR.TSA.AMYG_SEX_B', HCP.STAT.MOR.LOG.AMYG_SEX_B');
R_AMYG_LOG_RAW = corr(HCP.STAT.MOR.LOG.AMYG_SEX_B', HCP.STAT.MOR.RAW.AMYG_SEX_B');
[R_AMYG_TSA_RAW R_AMYG_TSA_LOG R_AMYG_LOG_RAW]

% Test the significance of slope - "1"
H0 = [0 1];  c = 1;
mdl_STR_TSA_RAW = fitlm(HCP.STAT.MOR.TSA.STR_SEX_B, HCP.STAT.MOR.RAW.STR_SEX_B);
mdl_STR_TSA_LOG = fitlm(HCP.STAT.MOR.TSA.STR_SEX_B, HCP.STAT.MOR.LOG.STR_SEX_B);
mdl_STR_LOG_RAW = fitlm(HCP.STAT.MOR.LOG.STR_SEX_B, HCP.STAT.MOR.RAW.STR_SEX_B);

disp('Statistics for slope vs 1')
[mdl_STR_TSA_RAW.Coefficients.Estimate(2) mdl_STR_TSA_LOG.Coefficients.Estimate(2) mdl_STR_LOG_RAW.Coefficients.Estimate(2)]
5*[coefTest(mdl_STR_TSA_RAW, H0, c) coefTest(mdl_STR_TSA_LOG, H0, c) coefTest(mdl_STR_LOG_RAW, H0, c)]

% Plots for STR
h = figure('WindowState', 'maximized'); hold on;
hexscatter(HCP.STAT.MOR.RAW.STR_SEX_B, HCP.STAT.MOR.TSA.STR_SEX_B, 'res', 50, 'xlim', [-1.5 1.5], 'ylim', [-1.5 1.5]);
scatter(HCP.STAT.MOR.RAW.STR_SEX_B, HCP.STAT.MOR.TSA.STR_SEX_B, 5, [0 0.5 0.5], 'filled', 'MarkerEdgeColor', [0 0.5 0.5], 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0);
xlim([-1.5 1.5]); xticks(-1.5:0.5:1.5); ylim([-1.5 1.5]); yticks(-1.5:0.5:1.5); set(gca,'TickDir','out');
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--';
x = xline(0); x.LineStyle = '--'; y = yline(0); y.LineStyle = '--';
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xlim([-1.5 1.5]); ylim([-1.5 1.5]);
axis square;
exportgraphics(h, fullfile(figDir, 'spatialcorr_str_tsa_raw.png'), 'Resolution', 300);
close(h);

h = figure('WindowState', 'maximized'); hold on;
hexscatter(HCP.STAT.MOR.LOG.STR_SEX_B, HCP.STAT.MOR.TSA.STR_SEX_B, 'res', 50, 'xlim', [-0.05 0.05], 'ylim', [-1 1]);
scatter(HCP.STAT.MOR.LOG.STR_SEX_B, HCP.STAT.MOR.TSA.STR_SEX_B, 5, [0 0.5 0.5], 'filled', 'MarkerEdgeColor', [0 0.5 0.5], 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0);
xlim([-0.05 0.05]); xticks(-0.05:0.025:0.05); ylim([-1 1]); yticks(-1:0.5:1); set(gca,'TickDir','out');
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--';
x = xline(0); x.LineStyle = '--'; y = yline(0); y.LineStyle = '--';
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xlim([-0.05 0.05]); ylim([-1 1]);
axis square;
exportgraphics(h, fullfile(figDir, 'spatialcorr_str_tsa_log.png'), 'Resolution', 300);
close(h);

h = figure('WindowState', 'maximized'); hold on;
hexscatter(HCP.STAT.MOR.RAW.STR_SEX_B, HCP.STAT.MOR.LOG.STR_SEX_B, 'res', 50, 'xlim', [-1.5 1.5], 'ylim', [-0.05 0.05]);
scatter(HCP.STAT.MOR.RAW.STR_SEX_B, HCP.STAT.MOR.LOG.STR_SEX_B, 5, [0 0.5 0.5], 'filled', 'MarkerEdgeColor', [0 0.5 0.5], 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0);
xlim([-1.5 1.5]); xticks(-1.5:0.5:1.5); ylim([-0.05 0.05]); yticks(-0.05:0.025:0.05); set(gca,'TickDir','out');
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--';
x = xline(0); x.LineStyle = '--'; y = yline(0); y.LineStyle = '--';
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xlim([-1.5 1.5]); ylim([-0.05 0.05]);
axis square;
exportgraphics(h, fullfile(figDir, 'spatialcorr_str_log_raw.png'), 'Resolution', 300);
close(h);

h = figure('WindowState', 'maximized'); hold on;
hexscatter(HCP.STAT.MOR.LOG.STR_SEX_B_SD, HCP.STAT.MOR.TSA.STR_SEX_B, 'res', 50, 'xlim', [-1 1], 'ylim', [-1 1]);
scatter(HCP.STAT.MOR.LOG.STR_SEX_B_SD, HCP.STAT.MOR.TSA.STR_SEX_B, 5, [0 0.5 0.5], 'filled', 'MarkerEdgeColor', [0 0.5 0.5], 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0);
xlim([-1 1]); xticks(-1:0.5:1); ylim([-1 1]); yticks(-1:0.5:1); set(gca,'TickDir','out');
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--';
x = xline(0); x.LineStyle = '--'; y = yline(0); y.LineStyle = '--';
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xlim([-1 1]); ylim([-1 1]);
axis square;
exportgraphics(h, fullfile(figDir, 'spatialcorr_str_tsa_log-sd.png'), 'Resolution', 300);
close(h);

% Plots for GP
h = figure('WindowState', 'maximized'); hold on;
hexscatter(HCP.STAT.MOR.RAW.GP_SEX_B, HCP.STAT.MOR.TSA.GP_SEX_B, 'res', 50, 'xlim', [-1.5 1.5], 'ylim', [-1.5 1.5]);
scatter(HCP.STAT.MOR.RAW.GP_SEX_B, HCP.STAT.MOR.TSA.GP_SEX_B, 5, [0 0.5 0.5], 'filled', 'MarkerEdgeColor', [0 0.5 0.5], 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0);
xlim([-1.5 1.5]); xticks(-1.5:0.5:1.5); ylim([-1.5 1.5]); yticks(-1.5:0.5:1.5); set(gca,'TickDir','out');
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--';
x = xline(0); x.LineStyle = '--'; y = yline(0); y.LineStyle = '--';
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xlim([-1.5 1.5]); ylim([-1.5 1.5]);
axis square;
exportgraphics(h, fullfile(figDir, 'spatialcorr_gp_tsa_raw.png'), 'Resolution', 300);
close(h);

h = figure('WindowState', 'maximized'); hold on;
hexscatter(HCP.STAT.MOR.LOG.GP_SEX_B, HCP.STAT.MOR.TSA.GP_SEX_B, 'res', 50, 'xlim', [-0.05 0.05], 'ylim', [-1 1]);
scatter(HCP.STAT.MOR.LOG.GP_SEX_B, HCP.STAT.MOR.TSA.GP_SEX_B, 5, [0 0.5 0.5], 'filled', 'MarkerEdgeColor', [0 0.5 0.5], 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0);
xlim([-0.05 0.05]); xticks(-0.05:0.025:0.05); ylim([-1 1]); yticks(-1:0.5:1); set(gca,'TickDir','out');
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--';
x = xline(0); x.LineStyle = '--'; y = yline(0); y.LineStyle = '--';
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xlim([-0.05 0.05]); ylim([-1 1]);
axis square;
exportgraphics(h, fullfile(figDir, 'spatialcorr_gp_tsa_log.png'), 'Resolution', 300);
close(h);

h = figure('WindowState', 'maximized'); hold on;
hexscatter(HCP.STAT.MOR.RAW.GP_SEX_B, HCP.STAT.MOR.LOG.GP_SEX_B, 'res', 50, 'xlim', [-1.5 1.5], 'ylim', [-0.05 0.05]);
scatter(HCP.STAT.MOR.RAW.GP_SEX_B, HCP.STAT.MOR.LOG.GP_SEX_B, 5, [0 0.5 0.5], 'filled', 'MarkerEdgeColor', [0 0.5 0.5], 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0);
xlim([-1.5 1.5]); xticks(-1.5:0.5:1.5); ylim([-0.05 0.05]); yticks(-0.05:0.025:0.05); set(gca,'TickDir','out');
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--';
x = xline(0); x.LineStyle = '--'; y = yline(0); y.LineStyle = '--';
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xlim([-1.5 1.5]); ylim([-0.05 0.05]);
axis square;
exportgraphics(h, fullfile(figDir, 'spatialcorr_gp_log_raw.png'), 'Resolution', 300);
close(h);

h = figure('WindowState', 'maximized'); hold on;
hexscatter(HCP.STAT.MOR.LOG.GP_SEX_B_SD, HCP.STAT.MOR.TSA.GP_SEX_B, 'res', 50, 'xlim', [-1 1], 'ylim', [-1 1]);
scatter(HCP.STAT.MOR.LOG.GP_SEX_B_SD, HCP.STAT.MOR.TSA.GP_SEX_B, 5, [0 0.5 0.5], 'filled', 'MarkerEdgeColor', [0 0.5 0.5], 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0);
xlim([-1 1]); xticks(-1:0.5:1); ylim([-1 1]); yticks(-1:0.5:1); set(gca,'TickDir','out');
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--';
x = xline(0); x.LineStyle = '--'; y = yline(0); y.LineStyle = '--';
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xlim([-1 1]); ylim([-1 1]);
axis square;
exportgraphics(h, fullfile(figDir, 'spatialcorr_gp_tsa_log-sd.png'), 'Resolution', 300);
close(h);

% Plots for THAL
h = figure('WindowState', 'maximized'); hold on;
hexscatter(HCP.STAT.MOR.RAW.THAL_SEX_B, HCP.STAT.MOR.TSA.THAL_SEX_B, 'res', 50, 'xlim', [-1.5 1.5], 'ylim', [-1.5 1.5]);
scatter(HCP.STAT.MOR.RAW.THAL_SEX_B, HCP.STAT.MOR.TSA.THAL_SEX_B, 5, [0 0.5 0.5], 'filled', 'MarkerEdgeColor', [0 0.5 0.5], 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0);
xlim([-1.5 1.5]); xticks(-1.5:0.5:1.5); ylim([-1.5 1.5]); yticks(-1.5:0.5:1.5); set(gca,'TickDir','out');
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--';
x = xline(0); x.LineStyle = '--'; y = yline(0); y.LineStyle = '--';
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xlim([-1.5 1.5]); ylim([-1.5 1.5]);
axis square;
exportgraphics(h, fullfile(figDir, 'spatialcorr_thal_tsa_raw.png'), 'Resolution', 300);
close(h);

h = figure('WindowState', 'maximized'); hold on;
hexscatter(HCP.STAT.MOR.LOG.THAL_SEX_B, HCP.STAT.MOR.TSA.THAL_SEX_B, 'res', 50, 'xlim', [-0.05 0.05], 'ylim', [-1 1]);
scatter(HCP.STAT.MOR.LOG.THAL_SEX_B, HCP.STAT.MOR.TSA.THAL_SEX_B, 5, [0 0.5 0.5], 'filled', 'MarkerEdgeColor', [0 0.5 0.5], 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0);
xlim([-0.075 0.075]); xticks(-0.075:0.025:0.075); ylim([-1.5 1.5]); yticks(-1.5:0.5:1.5); set(gca,'TickDir','out');
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--';
x = xline(0); x.LineStyle = '--'; y = yline(0); y.LineStyle = '--';
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xlim([-0.075 0.075]); ylim([-1.5 1.5]);
axis square;
exportgraphics(h, fullfile(figDir, 'spatialcorr_thal_tsa_log.png'), 'Resolution', 300);
close(h);

h = figure('WindowState', 'maximized'); hold on;
hexscatter(HCP.STAT.MOR.RAW.THAL_SEX_B, HCP.STAT.MOR.LOG.THAL_SEX_B, 'res', 50, 'xlim', [-1.5 1.5], 'ylim', [-0.05 0.05]);
scatter(HCP.STAT.MOR.RAW.THAL_SEX_B, HCP.STAT.MOR.LOG.THAL_SEX_B, 5, [0 0.5 0.5], 'filled', 'MarkerEdgeColor', [0 0.5 0.5], 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0);
xlim([-1.5 1.5]); xticks(-1.5:0.5:1.5); ylim([-0.075 0.075]); yticks(-0.075:0.025:0.075); set(gca,'TickDir','out');
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--';
x = xline(0); x.LineStyle = '--'; y = yline(0); y.LineStyle = '--';
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xlim([-1.5 1.5]); ylim([-0.075 0.075]);
axis square;
exportgraphics(h, fullfile(figDir, 'spatialcorr_thal_log_raw.png'), 'Resolution', 300);
close(h);

h = figure('WindowState', 'maximized'); hold on;
hexscatter(HCP.STAT.MOR.LOG.THAL_SEX_B_SD, HCP.STAT.MOR.TSA.THAL_SEX_B, 'res', 50, 'xlim', [-1 1], 'ylim', [-1 1]);
scatter(HCP.STAT.MOR.LOG.THAL_SEX_B_SD, HCP.STAT.MOR.TSA.THAL_SEX_B, 5, [0 0.5 0.5], 'filled', 'MarkerEdgeColor', [0 0.5 0.5], 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0);
xlim([-1 1]); xticks(-1:0.5:1); ylim([-1 1]); yticks(-1:0.5:1); set(gca,'TickDir','out');
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--';
x = xline(0); x.LineStyle = '--'; y = yline(0); y.LineStyle = '--';
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xlim([-1 1]); ylim([-1 1]);
axis square;
exportgraphics(h, fullfile(figDir, 'spatialcorr_thal_tsa_log-sd.png'), 'Resolution', 300);
close(h);

% Plots for HIPP
h = figure('WindowState', 'maximized'); hold on;
hexscatter(HCP.STAT.MOR.RAW.HIPP_SEX_B, HCP.STAT.MOR.TSA.HIPP_SEX_B, 'res', 50, 'xlim', [-1.5 1.5], 'ylim', [-1.5 1.5]);
scatter(HCP.STAT.MOR.RAW.HIPP_SEX_B, HCP.STAT.MOR.TSA.HIPP_SEX_B, 5, [0 0.5 0.5], 'filled', 'MarkerEdgeColor', [0 0.5 0.5], 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0);
xlim([-1.5 1.5]); xticks(-1.5:0.5:1.5); ylim([-1.5 1.5]); yticks(-1.5:0.5:1.5); set(gca,'TickDir','out');
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--';
x = xline(0); x.LineStyle = '--'; y = yline(0); y.LineStyle = '--';
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xlim([-1.5 1.5]); ylim([-1.5 1.5]);
axis square;
exportgraphics(h, fullfile(figDir, 'spatialcorr_hipp_tsa_raw.png'), 'Resolution', 300);
close(h);

h = figure('WindowState', 'maximized'); hold on;
hexscatter(HCP.STAT.MOR.LOG.HIPP_SEX_B, HCP.STAT.MOR.TSA.HIPP_SEX_B, 'res', 50, 'xlim', [-0.05 0.05], 'ylim', [-1 1]);
scatter(HCP.STAT.MOR.LOG.HIPP_SEX_B, HCP.STAT.MOR.TSA.HIPP_SEX_B, 5, [0 0.5 0.5], 'filled', 'MarkerEdgeColor', [0 0.5 0.5], 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0);
xlim([-0.05 0.05]); xticks(-0.05:0.025:0.05); ylim([-1 1]); yticks(-1:0.5:1); set(gca,'TickDir','out');
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--';
x = xline(0); x.LineStyle = '--'; y = yline(0); y.LineStyle = '--';
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xlim([-0.05 0.05]); ylim([-1 1]);
axis square;
exportgraphics(h, fullfile(figDir, 'spatialcorr_hipp_tsa_log.png'), 'Resolution', 300);
close(h);

h = figure('WindowState', 'maximized'); hold on;
hexscatter(HCP.STAT.MOR.RAW.HIPP_SEX_B, HCP.STAT.MOR.LOG.HIPP_SEX_B, 'res', 50, 'xlim', [-1.5 1.5], 'ylim', [-0.05 0.05]);
scatter(HCP.STAT.MOR.RAW.HIPP_SEX_B, HCP.STAT.MOR.LOG.HIPP_SEX_B, 5, [0 0.5 0.5], 'filled', 'MarkerEdgeColor', [0 0.5 0.5], 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0);
xlim([-1.5 1.5]); xticks(-1.5:0.5:1.5); ylim([-0.05 0.05]); yticks(-0.05:0.025:0.05); set(gca,'TickDir','out');
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--';
x = xline(0); x.LineStyle = '--'; y = yline(0); y.LineStyle = '--';
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xlim([-1.5 1.5]); ylim([-0.05 0.05]);
axis square;
exportgraphics(h, fullfile(figDir, 'spatialcorr_hipp_log_raw.png'), 'Resolution', 300);
close(h);

h = figure('WindowState', 'maximized'); hold on;
hexscatter(HCP.STAT.MOR.LOG.HIPP_SEX_B_SD, HCP.STAT.MOR.TSA.HIPP_SEX_B, 'res', 50, 'xlim', [-1 1], 'ylim', [-1 1]);
scatter(HCP.STAT.MOR.LOG.HIPP_SEX_B_SD, HCP.STAT.MOR.TSA.HIPP_SEX_B, 5, [0 0.5 0.5], 'filled', 'MarkerEdgeColor', [0 0.5 0.5], 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0);
xlim([-1 1]); xticks(-1:0.5:1); ylim([-1 1]); yticks(-1:0.5:1); set(gca,'TickDir','out');
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--';
x = xline(0); x.LineStyle = '--'; y = yline(0); y.LineStyle = '--';
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xlim([-1 1]); ylim([-1 1]);
axis square;
exportgraphics(h, fullfile(figDir, 'spatialcorr_hipp_tsa_log-sd.png'), 'Resolution', 300);
close(h);

% Plots for AMYG
h = figure('WindowState', 'maximized'); hold on;
hexscatter(HCP.STAT.MOR.RAW.AMYG_SEX_B, HCP.STAT.MOR.TSA.AMYG_SEX_B, 'res', 50, 'xlim', [-1.5 1.5], 'ylim', [-1.5 1.5]);
scatter(HCP.STAT.MOR.RAW.AMYG_SEX_B, HCP.STAT.MOR.TSA.AMYG_SEX_B, 5, [0 0.5 0.5], 'filled', 'MarkerEdgeColor', [0 0.5 0.5], 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0);
xlim([-1.5 1.5]); xticks(-1.5:0.5:1.5); ylim([-1.5 1.5]); yticks(-1.5:0.5:1.5); set(gca,'TickDir','out');
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--';
x = xline(0); x.LineStyle = '--'; y = yline(0); y.LineStyle = '--';
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xlim([-1.5 1.5]); ylim([-1.5 1.5]);
axis square;
exportgraphics(h, fullfile(figDir, 'spatialcorr_amyg_tsa_raw.png'), 'Resolution', 300);
close(h);

h = figure('WindowState', 'maximized'); hold on;
hexscatter(HCP.STAT.MOR.LOG.AMYG_SEX_B, HCP.STAT.MOR.TSA.AMYG_SEX_B, 'res', 50, 'xlim', [-0.05 0.05], 'ylim', [-1 1]);
scatter(HCP.STAT.MOR.LOG.AMYG_SEX_B, HCP.STAT.MOR.TSA.AMYG_SEX_B, 5, [0 0.5 0.5], 'filled', 'MarkerEdgeColor', [0 0.5 0.5], 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0);
xlim([-0.05 0.05]); xticks(-0.05:0.025:0.05); ylim([-1 1]); yticks(-1:0.5:1); set(gca,'TickDir','out');
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--';
x = xline(0); x.LineStyle = '--'; y = yline(0); y.LineStyle = '--';
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xlim([-0.05 0.05]); ylim([-1 1]);
axis square;
exportgraphics(h, fullfile(figDir, 'spatialcorr_amyg_tsa_log.png'), 'Resolution', 300);
close(h);

h = figure('WindowState', 'maximized'); hold on;
hexscatter(HCP.STAT.MOR.RAW.AMYG_SEX_B, HCP.STAT.MOR.LOG.AMYG_SEX_B, 'res', 50, 'xlim', [-1.5 1.5], 'ylim', [-0.05 0.05]);
scatter(HCP.STAT.MOR.RAW.AMYG_SEX_B, HCP.STAT.MOR.LOG.AMYG_SEX_B, 5, [0 0.5 0.5], 'filled', 'MarkerEdgeColor', [0 0.5 0.5], 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0);
xlim([-1.5 1.5]); xticks(-1.5:0.5:1.5); ylim([-0.05 0.05]); yticks(-0.05:0.025:0.05); set(gca,'TickDir','out');
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--';
x = xline(0); x.LineStyle = '--'; y = yline(0); y.LineStyle = '--';
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xlim([-1.5 1.5]); ylim([-0.05 0.05]);
axis square;
exportgraphics(h, fullfile(figDir, 'spatialcorr_amyg_log_raw.png'), 'Resolution', 300);
close(h);

h = figure('WindowState', 'maximized'); hold on;
hexscatter(HCP.STAT.MOR.LOG.AMYG_SEX_B_SD, HCP.STAT.MOR.TSA.AMYG_SEX_B, 'res', 50, 'xlim', [-1 1], 'ylim', [-1 1]);
scatter(HCP.STAT.MOR.LOG.AMYG_SEX_B_SD, HCP.STAT.MOR.TSA.AMYG_SEX_B, 5, [0 0.5 0.5], 'filled', 'MarkerEdgeColor', [0 0.5 0.5], 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0);
xlim([-1 1]); xticks(-1:0.5:1); ylim([-1 1]); yticks(-1:0.5:1); set(gca,'TickDir','out');
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--';
x = xline(0); x.LineStyle = '--'; y = yline(0); y.LineStyle = '--';
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xlim([-1 1]); ylim([-1 1]);
axis square;
exportgraphics(h, fullfile(figDir, 'spatialcorr_amyg_tsa_log-sd.png'), 'Resolution', 300);
close(h);
