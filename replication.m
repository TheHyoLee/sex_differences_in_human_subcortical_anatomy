%% replication.m - Within-sample reproducibility and between-sample replicability
%
% Lee HM, Liu S, Lalonde FM, et al. Sex Differences in Human Subcortical Anatomy
% are Robust and Replicable in Adulthood, Already Evident Prepubertally, and Aligned
% with Sex Chromosomal and Gonadal Influences. bioRxiv (2026).
% https://doi.org/10.64898/2026.10.05.756728
%
% Purpose:
%   (1) Within-sample reproducibility of HCP-YA sex-difference maps using
%       10,000 sex-balanced split-half subsamples (n = 50 to 450 per half).
%   (2) Between-sample replicability (HCP-YA vs. UKB): spatial correlations,
%       conjunction (agreement) maps, and Moran spectral randomization tests.
%   Requires the Parallel Computing Toolbox for parfor (or replace with for).
%
% Inputs (participant-level data are not shared; see Data availability):
%   - <dataDir>/HCP_results.mat (from hcp_stats.m)
%   - <dataDir>/UKB_results.mat (from ukb_stats.m)
%   - <templateDir>/maget_template_surfaces_distance_area.mat
%
% Outputs:
%   - <dataDir>/HCP_replication.mat
%   - Figures in <figDir>
%
% Dependencies (not included in this repository):
%   - SurfStat (SurfStatLinMod, SurfStatT, term)
%   - fdr_bh (Benjamini-Hochberg FDR; MATLAB File Exchange)
%   - MATLAB Statistics and Machine Learning Toolbox
%   - Parallel Computing Toolbox (parfor)
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
figDir             = '/path/to/figures/replication';

addpath(genpath(toolboxDir));
addpath(genpath(codeDir));

%% Prep

% Load results from hcp_stats.m and ukb_stats.m
load(fullfile(dataDir, 'HCP_results.mat'));
load(fullfile(dataDir, 'UKB_results.mat'));

% Load surface area maps for the template surfaces
load(fullfile(templateDir, 'maget_template_surfaces_distance_area.mat'));

% Misc
loadColorMaps;
nPerm = 10000;
thresh = 0.05;

%% Within-sample reproducibility: sex-balanced split-half subsampling

% (struct)
STR_L = HCP.VOL.STR_L(HCP.QC.indSTR)';
STR_R = HCP.VOL.STR_R(HCP.QC.indSTR)';
GP_L = HCP.VOL.GP_L(HCP.QC.indGP)';
GP_R = HCP.VOL.GP_R(HCP.QC.indGP)';
THAL_L = HCP.VOL.THAL_L(HCP.QC.indTHAL)';
THAL_R = HCP.VOL.THAL_R(HCP.QC.indTHAL)';
HIPP_L = HCP.VOL.HIPP_L(HCP.QC.indHIPP)';
HIPP_R = HCP.VOL.HIPP_R(HCP.QC.indHIPP)';
AMYG_L = HCP.VOL.AMYG_L(HCP.QC.indAMYG)';
AMYG_R = HCP.VOL.AMYG_R(HCP.QC.indAMYG)';

% (morph)
STR_MOR = HCP.MOR.STR(HCP.QC.indSTR,:);
GP_MOR = HCP.MOR.GP(HCP.QC.indGP,:);
THAL_MOR = HCP.MOR.THAL(HCP.QC.indTHAL,:);
HIPP_MOR = HCP.MOR.HIPP(HCP.QC.indHIPP,:);
AMYG_MOR = HCP.MOR.AMYG(HCP.QC.indAMYG,:);

% (TBV)
STR_TBV = HCP.TBV(HCP.QC.indSTR);
GP_TBV = HCP.TBV(HCP.QC.indGP);
THAL_TBV = HCP.TBV(HCP.QC.indTHAL);
HIPP_TBV = HCP.TBV(HCP.QC.indHIPP);
AMYG_TBV = HCP.TBV(HCP.QC.indAMYG);

% (TSA)
STR_TSA = sum(HCP.MOR.STR(HCP.QC.indSTR,:),2);
GP_TSA = sum(HCP.MOR.GP(HCP.QC.indGP,:),2);
THAL_TSA = sum(HCP.MOR.THAL(HCP.QC.indTHAL,:),2);
HIPP_TSA = sum(HCP.MOR.HIPP(HCP.QC.indHIPP,:),2);
AMYG_TSA = sum(HCP.MOR.AMYG(HCP.QC.indAMYG,:),2);

% (Age)
STR_AGE = HCP.AGE(HCP.QC.indSTR);
GP_AGE = HCP.AGE(HCP.QC.indGP);
THAL_AGE = HCP.AGE(HCP.QC.indTHAL);
HIPP_AGE = HCP.AGE(HCP.QC.indHIPP);
AMYG_AGE = HCP.AGE(HCP.QC.indAMYG);

% (Sex)
STR_SEX = HCP.SEX(HCP.QC.indSTR);
GP_SEX = HCP.SEX(HCP.QC.indGP);
THAL_SEX = HCP.SEX(HCP.QC.indTHAL);
HIPP_SEX = HCP.SEX(HCP.QC.indHIPP);
AMYG_SEX = HCP.SEX(HCP.QC.indAMYG);

% (ID)
STR_ID = HCP.ID(HCP.QC.indSTR);
GP_ID = HCP.ID(HCP.QC.indGP);
THAL_ID = HCP.ID(HCP.QC.indTHAL);
HIPP_ID = HCP.ID(HCP.QC.indHIPP);
AMYG_ID = HCP.ID(HCP.QC.indAMYG);

% Subsampling parameters
nBS = 10000;
SizeBS = 50:50:450;

% Split-half subsampling (without replacement)
for iS = 1:length(SizeBS)

    disp(['[Computing stats] Current sample size: ' num2str(SizeBS(iS))]);

    parfor iBS = 1:nBS % Iteration

        % PART 1
        % Sex-stratified sampling without replacement, split into two halves
        [~, indBS_F_STR] = datasample(STR_ID(STR_SEX==1), SizeBS(iS), 'Replace', false);
        [~, indBS_M_STR] = datasample(STR_ID(STR_SEX==0), SizeBS(iS), 'Replace', false);
        tempF_STR = find(STR_SEX==1);
        tempM_STR = find(STR_SEX==0);
        indBS1_STR = [tempF_STR(indBS_F_STR(1:end/2)); tempM_STR(indBS_M_STR(1:end/2))];
        indBS2_STR = [tempF_STR(indBS_F_STR(end/2+1:end)); tempM_STR(indBS_M_STR(end/2+1:end))];

        [~, indBS_F_GP] = datasample(GP_ID(GP_SEX==1), SizeBS(iS), 'Replace', false);
        [~, indBS_M_GP] = datasample(GP_ID(GP_SEX==0), SizeBS(iS), 'Replace', false);
        tempF_GP = find(GP_SEX==1);
        tempM_GP = find(GP_SEX==0);
        indBS1_GP = [tempF_GP(indBS_F_GP(1:end/2)); tempM_GP(indBS_M_GP(1:end/2))];
        indBS2_GP = [tempF_GP(indBS_F_GP(end/2+1:end)); tempM_GP(indBS_M_GP(end/2+1:end))];

        [~, indBS_F_THAL] = datasample(THAL_ID(THAL_SEX==1), SizeBS(iS), 'Replace', false);
        [~, indBS_M_THAL] = datasample(THAL_ID(THAL_SEX==0), SizeBS(iS), 'Replace', false);
        tempF_THAL = find(THAL_SEX==1);
        tempM_THAL = find(THAL_SEX==0);
        indBS1_THAL = [tempF_THAL(indBS_F_THAL(1:end/2)); tempM_THAL(indBS_M_THAL(1:end/2))];
        indBS2_THAL = [tempF_THAL(indBS_F_THAL(end/2+1:end)); tempM_THAL(indBS_M_THAL(end/2+1:end))];

        [~, indBS_F_HIPP] = datasample(HIPP_ID(HIPP_SEX==1), SizeBS(iS), 'Replace', false);
        [~, indBS_M_HIPP] = datasample(HIPP_ID(HIPP_SEX==0), SizeBS(iS), 'Replace', false);
        tempF_HIPP = find(HIPP_SEX==1);
        tempM_HIPP = find(HIPP_SEX==0);
        indBS1_HIPP = [tempF_HIPP(indBS_F_HIPP(1:end/2)); tempM_HIPP(indBS_M_HIPP(1:end/2))];
        indBS2_HIPP = [tempF_HIPP(indBS_F_HIPP(end/2+1:end)); tempM_HIPP(indBS_M_HIPP(end/2+1:end))];

        [~, indBS_F_AMYG] = datasample(AMYG_ID(AMYG_SEX==1), SizeBS(iS), 'Replace', false);
        [~, indBS_M_AMYG] = datasample(AMYG_ID(AMYG_SEX==0), SizeBS(iS), 'Replace', false);
        tempF_AMYG = find(AMYG_SEX==1);
        tempM_AMYG = find(AMYG_SEX==0);
        indBS1_AMYG = [tempF_AMYG(indBS_F_AMYG(1:end/2)); tempM_AMYG(indBS_M_AMYG(1:end/2))];
        indBS2_AMYG = [tempF_AMYG(indBS_F_AMYG(end/2+1:end)); tempM_AMYG(indBS_M_AMYG(end/2+1:end))];

        % PART 2
        % Prepare the terms for first half sample
        iSTR_VOL_L_1 = zscore(STR_L(indBS1_STR));
        iSTR_VOL_R_1 = zscore(STR_R(indBS1_STR));
        iSTR_MOR_1 = zscore(STR_MOR(indBS1_STR,:));
        iSTR_AGE_1 = term(STR_AGE(indBS1_STR) - mean(STR_AGE(indBS1_STR)));
        iSTR_SEX_1 = term(STR_SEX(indBS1_STR));
        iSTR_TBV_1 = term(STR_TBV(indBS1_STR) - mean(STR_TBV(indBS1_STR)));
        iSTR_TSA_1 = term(STR_TSA(indBS1_STR) - mean(STR_TSA(indBS1_STR)));

        iGP_VOL_L_1 = zscore(GP_L(indBS1_GP));
        iGP_VOL_R_1 = zscore(GP_R(indBS1_GP));
        iGP_MOR_1 = zscore(GP_MOR(indBS1_GP,:));
        iGP_AGE_1 = term(GP_AGE(indBS1_GP) - mean(GP_AGE(indBS1_GP)));
        iGP_SEX_1 = term(GP_SEX(indBS1_GP));
        iGP_TBV_1 = term(GP_TBV(indBS1_GP) - mean(GP_TBV(indBS1_GP)));
        iGP_TSA_1 = term(GP_TSA(indBS1_GP) - mean(GP_TSA(indBS1_GP)));

        iTHAL_VOL_L_1 = zscore(THAL_L(indBS1_THAL));
        iTHAL_VOL_R_1 = zscore(THAL_R(indBS1_THAL));
        iTHAL_MOR_1 = zscore(THAL_MOR(indBS1_THAL,:));
        iTHAL_AGE_1 = term(THAL_AGE(indBS1_THAL) - mean(THAL_AGE(indBS1_THAL)));
        iTHAL_SEX_1 = term(THAL_SEX(indBS1_THAL));
        iTHAL_TBV_1 = term(THAL_TBV(indBS1_THAL) - mean(THAL_TBV(indBS1_THAL)));
        iTHAL_TSA_1 = term(THAL_TSA(indBS1_THAL) - mean(THAL_TSA(indBS1_THAL)));

        iHIPP_VOL_L_1 = zscore(HIPP_L(indBS1_HIPP));
        iHIPP_VOL_R_1 = zscore(HIPP_R(indBS1_HIPP));
        iHIPP_MOR_1 = zscore(HIPP_MOR(indBS1_HIPP,:));
        iHIPP_AGE_1 = term(HIPP_AGE(indBS1_HIPP) - mean(HIPP_AGE(indBS1_HIPP)));
        iHIPP_SEX_1 = term(HIPP_SEX(indBS1_HIPP));
        iHIPP_TBV_1 = term(HIPP_TBV(indBS1_HIPP) - mean(HIPP_TBV(indBS1_HIPP)));
        iHIPP_TSA_1 = term(HIPP_TSA(indBS1_HIPP) - mean(HIPP_TSA(indBS1_HIPP)));

        iAMYG_VOL_L_1 = zscore(AMYG_L(indBS1_AMYG));
        iAMYG_VOL_R_1 = zscore(AMYG_R(indBS1_AMYG));
        iAMYG_MOR_1 = zscore(AMYG_MOR(indBS1_AMYG,:));
        iAMYG_AGE_1 = term(AMYG_AGE(indBS1_AMYG) - mean(AMYG_AGE(indBS1_AMYG)));
        iAMYG_SEX_1 = term(AMYG_SEX(indBS1_AMYG));
        iAMYG_TBV_1 = term(AMYG_TBV(indBS1_AMYG) - mean(AMYG_TBV(indBS1_AMYG)));
        iAMYG_TSA_1 = term(AMYG_TSA(indBS1_AMYG) - mean(AMYG_TSA(indBS1_AMYG)));

        % Prepare the terms for second half sample
        iSTR_VOL_L_2 = zscore(STR_L(indBS2_STR));
        iSTR_VOL_R_2 = zscore(STR_R(indBS2_STR));
        iSTR_MOR_2 = zscore(STR_MOR(indBS2_STR,:));
        iSTR_AGE_2 = term(STR_AGE(indBS2_STR) - mean(STR_AGE(indBS2_STR)));
        iSTR_SEX_2 = term(STR_SEX(indBS2_STR));
        iSTR_TBV_2 = term(STR_TBV(indBS2_STR) - mean(STR_TBV(indBS2_STR)));
        iSTR_TSA_2 = term(STR_TSA(indBS2_STR) - mean(STR_TSA(indBS2_STR)));

        iGP_VOL_L_2 = zscore(GP_L(indBS2_GP));
        iGP_VOL_R_2 = zscore(GP_R(indBS2_GP));
        iGP_MOR_2 = zscore(GP_MOR(indBS2_GP,:));
        iGP_AGE_2 = term(GP_AGE(indBS2_GP) - mean(GP_AGE(indBS2_GP)));
        iGP_SEX_2 = term(GP_SEX(indBS2_GP));
        iGP_TBV_2 = term(GP_TBV(indBS2_GP) - mean(GP_TBV(indBS2_GP)));
        iGP_TSA_2 = term(GP_TSA(indBS2_GP) - mean(GP_TSA(indBS2_GP)));

        iTHAL_VOL_L_2 = zscore(THAL_L(indBS2_THAL));
        iTHAL_VOL_R_2 = zscore(THAL_R(indBS2_THAL));
        iTHAL_MOR_2 = zscore(THAL_MOR(indBS2_THAL,:));
        iTHAL_AGE_2 = term(THAL_AGE(indBS2_THAL) - mean(THAL_AGE(indBS2_THAL)));
        iTHAL_SEX_2 = term(THAL_SEX(indBS2_THAL));
        iTHAL_TBV_2 = term(THAL_TBV(indBS2_THAL) - mean(THAL_TBV(indBS2_THAL)));
        iTHAL_TSA_2 = term(THAL_TSA(indBS2_THAL) - mean(THAL_TSA(indBS2_THAL)));

        iHIPP_VOL_L_2 = zscore(HIPP_L(indBS2_HIPP));
        iHIPP_VOL_R_2 = zscore(HIPP_R(indBS2_HIPP));
        iHIPP_MOR_2 = zscore(HIPP_MOR(indBS2_HIPP,:));
        iHIPP_AGE_2 = term(HIPP_AGE(indBS2_HIPP) - mean(HIPP_AGE(indBS2_HIPP)));
        iHIPP_SEX_2 = term(HIPP_SEX(indBS2_HIPP));
        iHIPP_TBV_2 = term(HIPP_TBV(indBS2_HIPP) - mean(HIPP_TBV(indBS2_HIPP)));
        iHIPP_TSA_2 = term(HIPP_TSA(indBS2_HIPP) - mean(HIPP_TSA(indBS2_HIPP)));

        iAMYG_VOL_L_2 = zscore(AMYG_L(indBS2_AMYG));
        iAMYG_VOL_R_2 = zscore(AMYG_R(indBS2_AMYG));
        iAMYG_MOR_2 = zscore(AMYG_MOR(indBS2_AMYG,:));
        iAMYG_AGE_2 = term(AMYG_AGE(indBS2_AMYG) - mean(AMYG_AGE(indBS2_AMYG)));
        iAMYG_SEX_2 = term(AMYG_SEX(indBS2_AMYG));
        iAMYG_TBV_2 = term(AMYG_TBV(indBS2_AMYG) - mean(AMYG_TBV(indBS2_AMYG)));
        iAMYG_TSA_2 = term(AMYG_TSA(indBS2_AMYG) - mean(AMYG_TSA(indBS2_AMYG)));

        % PART 3
        % Compute stats for structural volume
        M_VOL_STR_1 = 1 + iSTR_SEX_1 + iSTR_AGE_1 + iSTR_TBV_1;
        M_VOL_GP_1 = 1 + iGP_SEX_1 + iGP_AGE_1 + iGP_TBV_1;
        M_VOL_THAL_1 = 1 + iTHAL_SEX_1 + iTHAL_AGE_1 + iTHAL_TBV_1;
        M_VOL_HIPP_1 = 1 + iHIPP_SEX_1 + iHIPP_AGE_1 + iHIPP_TBV_1;
        M_VOL_AMYG_1 = 1 + iAMYG_SEX_1 + iAMYG_AGE_1 + iAMYG_TBV_1;

        M_VOL_STR_2 = 1 + iSTR_SEX_2 + iSTR_AGE_2 + iSTR_TBV_2;
        M_VOL_GP_2 = 1 + iGP_SEX_2 + iGP_AGE_2 + iGP_TBV_2;
        M_VOL_THAL_2 = 1 + iTHAL_SEX_2 + iTHAL_AGE_2 + iTHAL_TBV_2;
        M_VOL_HIPP_2 = 1 + iHIPP_SEX_2 + iHIPP_AGE_2 + iHIPP_TBV_2;
        M_VOL_AMYG_2 = 1 + iAMYG_SEX_2 + iAMYG_AGE_2 + iAMYG_TBV_2;

        % Compute stats for regional surface area
        M_MOR_STR_1 = 1 + iSTR_SEX_1 + iSTR_AGE_1 + iSTR_TSA_1;
        M_MOR_GP_1 = 1 + iGP_SEX_1 + iGP_AGE_1 + iGP_TSA_1;
        M_MOR_THAL_1 = 1 + iTHAL_SEX_1 + iTHAL_AGE_1 + iTHAL_TSA_1;
        M_MOR_HIPP_1 = 1 + iHIPP_SEX_1 + iHIPP_AGE_1 + iHIPP_TSA_1;
        M_MOR_AMYG_1 = 1 + iAMYG_SEX_1 + iAMYG_AGE_1 + iAMYG_TSA_1;

        M_MOR_STR_2 = 1 + iSTR_SEX_2 + iSTR_AGE_2 + iSTR_TSA_2;
        M_MOR_GP_2 = 1 + iGP_SEX_2 + iGP_AGE_2 + iGP_TSA_2;
        M_MOR_THAL_2 = 1 + iTHAL_SEX_2 + iTHAL_AGE_2 + iTHAL_TSA_2;
        M_MOR_HIPP_2 = 1 + iHIPP_SEX_2 + iHIPP_AGE_2 + iHIPP_TSA_2;
        M_MOR_AMYG_2 = 1 + iAMYG_SEX_2 + iAMYG_AGE_2 + iAMYG_TSA_2;

        % PART 4A
        % Compute stats for structural volume (first half)
        SLM_STR_VOL_L_1 = SurfStatLinMod(iSTR_VOL_L_1, M_VOL_STR_1);
        SLM_STR_VOL_L_SEX_1 = SurfStatT(SLM_STR_VOL_L_1, iSTR_SEX_1);

        SLM_GP_VOL_L_1 = SurfStatLinMod(iGP_VOL_L_1, M_VOL_GP_1);
        SLM_GP_VOL_L_SEX_1 = SurfStatT(SLM_GP_VOL_L_1, iGP_SEX_1);

        SLM_THAL_VOL_L_1 = SurfStatLinMod(iTHAL_VOL_L_1, M_VOL_THAL_1);
        SLM_THAL_VOL_L_SEX_1 = SurfStatT(SLM_THAL_VOL_L_1, iTHAL_SEX_1);

        SLM_HIPP_VOL_L_1 = SurfStatLinMod(iHIPP_VOL_L_1, M_VOL_HIPP_1);
        SLM_HIPP_VOL_L_SEX_1 = SurfStatT(SLM_HIPP_VOL_L_1, iHIPP_SEX_1);

        SLM_AMYG_VOL_L_1 = SurfStatLinMod(iAMYG_VOL_L_1, M_VOL_AMYG_1);
        SLM_AMYG_VOL_L_SEX_1 = SurfStatT(SLM_AMYG_VOL_L_1, iAMYG_SEX_1);

        SLM_STR_VOL_R_1 = SurfStatLinMod(iSTR_VOL_R_1, M_VOL_STR_1);
        SLM_STR_VOL_R_SEX_1 = SurfStatT(SLM_STR_VOL_R_1, iSTR_SEX_1);

        SLM_GP_VOL_R_1 = SurfStatLinMod(iGP_VOL_R_1, M_VOL_GP_1);
        SLM_GP_VOL_R_SEX_1 = SurfStatT(SLM_GP_VOL_R_1, iGP_SEX_1);

        SLM_THAL_VOL_R_1 = SurfStatLinMod(iTHAL_VOL_R_1, M_VOL_THAL_1);
        SLM_THAL_VOL_R_SEX_1 = SurfStatT(SLM_THAL_VOL_R_1, iTHAL_SEX_1);

        SLM_HIPP_VOL_R_1 = SurfStatLinMod(iHIPP_VOL_R_1, M_VOL_HIPP_1);
        SLM_HIPP_VOL_R_SEX_1 = SurfStatT(SLM_HIPP_VOL_R_1, iHIPP_SEX_1);

        SLM_AMYG_VOL_R_1 = SurfStatLinMod(iAMYG_VOL_R_1, M_VOL_AMYG_1);
        SLM_AMYG_VOL_R_SEX_1 = SurfStatT(SLM_AMYG_VOL_R_1, iAMYG_SEX_1);

        % Compute stats for structural volume (second half)
        SLM_STR_VOL_L_2 = SurfStatLinMod(iSTR_VOL_L_2, M_VOL_STR_2);
        SLM_STR_VOL_L_SEX_2 = SurfStatT(SLM_STR_VOL_L_2, iSTR_SEX_2);

        SLM_GP_VOL_L_2 = SurfStatLinMod(iGP_VOL_L_2, M_VOL_GP_2);
        SLM_GP_VOL_L_SEX_2 = SurfStatT(SLM_GP_VOL_L_2, iGP_SEX_2);

        SLM_THAL_VOL_L_2 = SurfStatLinMod(iTHAL_VOL_L_2, M_VOL_THAL_2);
        SLM_THAL_VOL_L_SEX_2 = SurfStatT(SLM_THAL_VOL_L_2, iTHAL_SEX_2);

        SLM_HIPP_VOL_L_2 = SurfStatLinMod(iHIPP_VOL_L_2, M_VOL_HIPP_2);
        SLM_HIPP_VOL_L_SEX_2 = SurfStatT(SLM_HIPP_VOL_L_2, iHIPP_SEX_2);

        SLM_AMYG_VOL_L_2 = SurfStatLinMod(iAMYG_VOL_L_2, M_VOL_AMYG_2);
        SLM_AMYG_VOL_L_SEX_2 = SurfStatT(SLM_AMYG_VOL_L_2, iAMYG_SEX_2);

        SLM_STR_VOL_R_2 = SurfStatLinMod(iSTR_VOL_R_2, M_VOL_STR_2);
        SLM_STR_VOL_R_SEX_2 = SurfStatT(SLM_STR_VOL_R_2, iSTR_SEX_2);

        SLM_GP_VOL_R_2 = SurfStatLinMod(iGP_VOL_R_2, M_VOL_GP_2);
        SLM_GP_VOL_R_SEX_2 = SurfStatT(SLM_GP_VOL_R_2, iGP_SEX_2);

        SLM_THAL_VOL_R_2 = SurfStatLinMod(iTHAL_VOL_R_2, M_VOL_THAL_2);
        SLM_THAL_VOL_R_SEX_2 = SurfStatT(SLM_THAL_VOL_R_2, iTHAL_SEX_2);

        SLM_HIPP_VOL_R_2 = SurfStatLinMod(iHIPP_VOL_R_2, M_VOL_HIPP_2);
        SLM_HIPP_VOL_R_SEX_2 = SurfStatT(SLM_HIPP_VOL_R_2, iHIPP_SEX_2);

        SLM_AMYG_VOL_R_2 = SurfStatLinMod(iAMYG_VOL_R_2, M_VOL_AMYG_2);
        SLM_AMYG_VOL_R_SEX_2 = SurfStatT(SLM_AMYG_VOL_R_2, iAMYG_SEX_2);

        % PART 4B
        % Compute stats for regional surface area (first half)
        SLM_STR_MOR_1 = SurfStatLinMod(iSTR_MOR_1, M_MOR_STR_1);
        SLM_STR_MOR_SEX_1 = SurfStatT(SLM_STR_MOR_1, iSTR_SEX_1);

        SLM_GP_MOR_1 = SurfStatLinMod(iGP_MOR_1, M_MOR_GP_1);
        SLM_GP_MOR_SEX_1 = SurfStatT(SLM_GP_MOR_1, iGP_SEX_1);

        SLM_THAL_MOR_1 = SurfStatLinMod(iTHAL_MOR_1, M_MOR_THAL_1);
        SLM_THAL_MOR_SEX_1 = SurfStatT(SLM_THAL_MOR_1, iTHAL_SEX_1);

        SLM_HIPP_MOR_1 = SurfStatLinMod(iHIPP_MOR_1, M_MOR_HIPP_1);
        SLM_HIPP_MOR_SEX_1 = SurfStatT(SLM_HIPP_MOR_1, iHIPP_SEX_1);

        SLM_AMYG_MOR_1 = SurfStatLinMod(iAMYG_MOR_1, M_MOR_AMYG_1);
        SLM_AMYG_MOR_SEX_1 = SurfStatT(SLM_AMYG_MOR_1, iAMYG_SEX_1);

        % Compute stats for regional surface area (second half)
        SLM_STR_MOR_2 = SurfStatLinMod(iSTR_MOR_2, M_MOR_STR_2);
        SLM_STR_MOR_SEX_2 = SurfStatT(SLM_STR_MOR_2, iSTR_SEX_2);

        SLM_GP_MOR_2 = SurfStatLinMod(iGP_MOR_2, M_MOR_GP_2);
        SLM_GP_MOR_SEX_2 = SurfStatT(SLM_GP_MOR_2, iGP_SEX_2);

        SLM_THAL_MOR_2 = SurfStatLinMod(iTHAL_MOR_2, M_MOR_THAL_2);
        SLM_THAL_MOR_SEX_2 = SurfStatT(SLM_THAL_MOR_2, iTHAL_SEX_2);

        SLM_HIPP_MOR_2 = SurfStatLinMod(iHIPP_MOR_2, M_MOR_HIPP_2);
        SLM_HIPP_MOR_SEX_2 = SurfStatT(SLM_HIPP_MOR_2, iHIPP_SEX_2);

        SLM_AMYG_MOR_2 = SurfStatLinMod(iAMYG_MOR_2, M_MOR_AMYG_2);
        SLM_AMYG_MOR_SEX_2 = SurfStatT(SLM_AMYG_MOR_2, iAMYG_SEX_2);

        % PART 5
        % Extract stats (first half)
        STR_VOL_L_B_1(iBS,iS) = SLM_STR_VOL_L_SEX_1.coef(2);
        GP_VOL_L_B_1(iBS,iS) = SLM_GP_VOL_L_SEX_1.coef(2);
        THAL_VOL_L_B_1(iBS,iS) = SLM_THAL_VOL_L_SEX_1.coef(2);
        HIPP_VOL_L_B_1(iBS,iS) = SLM_HIPP_VOL_L_SEX_1.coef(2);
        AMYG_VOL_L_B_1(iBS,iS) = SLM_AMYG_VOL_L_SEX_1.coef(2);

        STR_VOL_R_B_1(iBS,iS) = SLM_STR_VOL_R_SEX_1.coef(2);
        GP_VOL_R_B_1(iBS,iS) = SLM_GP_VOL_R_SEX_1.coef(2);
        THAL_VOL_R_B_1(iBS,iS) = SLM_THAL_VOL_R_SEX_1.coef(2);
        HIPP_VOL_R_B_1(iBS,iS) = SLM_HIPP_VOL_R_SEX_1.coef(2);
        AMYG_VOL_R_B_1(iBS,iS) = SLM_AMYG_VOL_R_SEX_1.coef(2);

        STR_MOR_B_1(:,iBS,iS) = SLM_STR_MOR_SEX_1.coef(2,:);
        GP_MOR_B_1(:,iBS,iS) = SLM_GP_MOR_SEX_1.coef(2,:);
        THAL_MOR_B_1(:,iBS,iS) = SLM_THAL_MOR_SEX_1.coef(2,:);
        HIPP_MOR_B_1(:,iBS,iS) = SLM_HIPP_MOR_SEX_1.coef(2,:);
        AMYG_MOR_B_1(:,iBS,iS) = SLM_AMYG_MOR_SEX_1.coef(2,:);

        % Extract stats (second half)
        STR_VOL_L_B_2(iBS,iS) = SLM_STR_VOL_L_SEX_2.coef(2);
        GP_VOL_L_B_2(iBS,iS) = SLM_GP_VOL_L_SEX_2.coef(2);
        THAL_VOL_L_B_2(iBS,iS) = SLM_THAL_VOL_L_SEX_2.coef(2);
        HIPP_VOL_L_B_2(iBS,iS) = SLM_HIPP_VOL_L_SEX_2.coef(2);
        AMYG_VOL_L_B_2(iBS,iS) = SLM_AMYG_VOL_L_SEX_2.coef(2);

        STR_VOL_R_B_2(iBS,iS) = SLM_STR_VOL_R_SEX_2.coef(2);
        GP_VOL_R_B_2(iBS,iS) = SLM_GP_VOL_R_SEX_2.coef(2);
        THAL_VOL_R_B_2(iBS,iS) = SLM_THAL_VOL_R_SEX_2.coef(2);
        HIPP_VOL_R_B_2(iBS,iS) = SLM_HIPP_VOL_R_SEX_2.coef(2);
        AMYG_VOL_R_B_2(iBS,iS) = SLM_AMYG_VOL_R_SEX_2.coef(2);

        STR_MOR_B_2(:,iBS,iS) = SLM_STR_MOR_SEX_2.coef(2,:);
        GP_MOR_B_2(:,iBS,iS) = SLM_GP_MOR_SEX_2.coef(2,:);
        THAL_MOR_B_2(:,iBS,iS) = SLM_THAL_MOR_SEX_2.coef(2,:);
        HIPP_MOR_B_2(:,iBS,iS) = SLM_HIPP_MOR_SEX_2.coef(2,:);
        AMYG_MOR_B_2(:,iBS,iS) = SLM_AMYG_MOR_SEX_2.coef(2,:);
    end
end

%% Extract results

for iS = 1:length(SizeBS)

    disp(['[Extracting stats] Current sample size: ' num2str(SizeBS(iS))]);

    for iBS = 1:nBS
        % Compute spatial correlation for regional surface area
        HCP.REP.MOR.STR_SC(iBS,iS) = corr(STR_MOR_B_1(:,iBS,iS), STR_MOR_B_2(:,iBS,iS));
        HCP.REP.MOR.GP_SC(iBS,iS) = corr(GP_MOR_B_1(:,iBS,iS), GP_MOR_B_2(:,iBS,iS));
        HCP.REP.MOR.THAL_SC(iBS,iS) = corr(THAL_MOR_B_1(:,iBS,iS), THAL_MOR_B_2(:,iBS,iS));
        HCP.REP.MOR.HIPP_SC(iBS,iS) = corr(HIPP_MOR_B_1(:,iBS,iS), HIPP_MOR_B_2(:,iBS,iS));
        HCP.REP.MOR.AMYG_SC(iBS,iS) = corr(AMYG_MOR_B_1(:,iBS,iS), AMYG_MOR_B_2(:,iBS,iS));
    end
end

% Save split-half results
save(fullfile(dataDir, 'HCP_replication.mat'), 'HCP', '-v7.3');

%% Within-sample reproducibility

% To skip the subsampling above, load the saved results instead:
% load(fullfile(dataDir, 'HCP_replication.mat'));

close all;

% raincloud (morph)
clear rainData;
rainData{1} = HCP.REP.MOR.STR_SC;
rainData{2} = HCP.REP.MOR.GP_SC;
rainData{3} = HCP.REP.MOR.THAL_SC;
rainData{4} = HCP.REP.MOR.HIPP_SC;
rainData{5} = HCP.REP.MOR.AMYG_SC;
group_inx = [ones(1,nPerm), 2.*ones(1,nPerm) 3.*ones(1,nPerm) 4.*ones(1,nPerm) 5.*ones(1,nPerm) 6.*ones(1,nPerm) 7.*ones(1,nPerm) 8.*ones(1,nPerm) 9.*ones(1,nPerm)];
c =  [0.45, 0.80, 0.69;...
      0.98, 0.40, 0.35;...
      0.55, 0.60, 0.79;...
      0.90, 0.70, 0.30;...
      0.20, 0.20, 0.20];
group_names = {'Striatum', 'Globus pallidus', 'Thalamus', 'Hippocampus', 'Amygdala'};
condition_names = {'50', '100', '150', '200', '250', '300', '350', '400', '450'};
h = figure; d = daviolinplot(rainData, 'groups', group_inx(1:90), 'colors', c, 'box', 3,...
    'boxcolor', 'w', 'scatter', 2, 'jitter', 1, 'scattercolor', 'same', 'smoothing', 1e-2,...
    'scattersize', 10, 'scatteralpha', 0.4, 'linkline', 1, 'outliers', 0, ...
    'xtlabels', condition_names,...
    'legend', group_names);
ylabel('Spatial correlation');
set(d.sc,'MarkerEdgeColor','none');
set(gca,'FontSize',10);
xl = xlim; xlim([xl(1)-0.1, xl(2)+0.3]);
ylim([-0.7 1]); yticks(-0.7:0.1:1);
ax = gca;
ax.YGrid = 'on';
ax.TickDir = 'out';
h.Position = [0 0 2000 1000];
legend('off');
exportgraphics(h, fullfile(figDir, 'raincloud_R_mor.png'), 'Resolution', 300);

% Display median [IQR]
[median(HCP.REP.MOR.STR_SC(:,end)) prctile(HCP.REP.MOR.STR_SC(:,end), 25) prctile(HCP.REP.MOR.STR_SC(:,end), 75)]
[median(HCP.REP.MOR.GP_SC(:,end)) prctile(HCP.REP.MOR.GP_SC(:,end), 25) prctile(HCP.REP.MOR.GP_SC(:,end), 75)]
[median(HCP.REP.MOR.THAL_SC(:,end)) prctile(HCP.REP.MOR.THAL_SC(:,end), 25) prctile(HCP.REP.MOR.THAL_SC(:,end), 75)]
[median(HCP.REP.MOR.HIPP_SC(:,end)) prctile(HCP.REP.MOR.HIPP_SC(:,end), 25) prctile(HCP.REP.MOR.HIPP_SC(:,end), 75)]
[median(HCP.REP.MOR.AMYG_SC(:,end)) prctile(HCP.REP.MOR.AMYG_SC(:,end), 25) prctile(HCP.REP.MOR.AMYG_SC(:,end), 75)]

%% Between-sample reproducibility

% Conjunction map (morph)
UKB_STR_AD = zeros(1,26000);
HCP_STR_AD = zeros(1,26000);
STR_AD = zeros(1,26000);

UKB_STR_AD(UKB.STAT.MOR.TSA.STR_SEX_B > 0) =  1;
UKB_STR_AD(UKB.STAT.MOR.TSA.STR_SEX_B < 0) = -1;
UKB_STR_AD((UKB.STAT.MOR.TSA.STR_SEX_Q < thresh) & (UKB.STAT.MOR.TSA.STR_SEX_B > 0))  =  2;
UKB_STR_AD((UKB.STAT.MOR.TSA.STR_SEX_Q < thresh) & (UKB.STAT.MOR.TSA.STR_SEX_B < 0))  = -2;

HCP_STR_AD(HCP.STAT.MOR.TSA.STR_SEX_B > 0) =  1;
HCP_STR_AD(HCP.STAT.MOR.TSA.STR_SEX_B < 0) = -1;
HCP_STR_AD((HCP.STAT.MOR.TSA.STR_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.STR_SEX_B > 0))  =  2;
HCP_STR_AD((HCP.STAT.MOR.TSA.STR_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.STR_SEX_B < 0))  = -2;

STR_AD((UKB_STR_AD == 2) & (HCP_STR_AD == 2)) = -1.9;
STR_AD((UKB_STR_AD == 2) & (HCP_STR_AD == 1)) = -1.6;
STR_AD((UKB_STR_AD == 1) & (HCP_STR_AD == 2)) = -1.6;
STR_AD((UKB_STR_AD == 1) & (HCP_STR_AD == 1)) = -1.3;

STR_AD((UKB_STR_AD == 2) & (HCP_STR_AD == -2)) = -0.9;
STR_AD((UKB_STR_AD == 2) & (HCP_STR_AD == -1)) = -0.6;
STR_AD((UKB_STR_AD == 1) & (HCP_STR_AD == -2)) = -0.6;
STR_AD((UKB_STR_AD == 1) & (HCP_STR_AD == -1)) = -0.3;

STR_AD((UKB_STR_AD == -2) & (HCP_STR_AD == 2)) = 0.3;
STR_AD((UKB_STR_AD == -2) & (HCP_STR_AD == 1)) = 0.6;
STR_AD((UKB_STR_AD == -1) & (HCP_STR_AD == 2)) = 0.6;
STR_AD((UKB_STR_AD == -1) & (HCP_STR_AD == 1)) = 0.9;

STR_AD((UKB_STR_AD == -2) & (HCP_STR_AD == -2)) = 1.3;
STR_AD((UKB_STR_AD == -2) & (HCP_STR_AD == -1)) = 1.6;
STR_AD((UKB_STR_AD == -1) & (HCP_STR_AD == -2)) = 1.6;
STR_AD((UKB_STR_AD == -1) & (HCP_STR_AD == -1)) = 1.9;

h = figure('WindowState', 'Maximized'); HyoMagetViewFlat(STR_AD, [-2 2], cmap.tri_div);
exportgraphics(h, fullfile(figDir, 'conjunction_STR.png'), 'Resolution', 300);

STR_AD_Count = [sum(S_STR_NA(STR_AD == 1.3)) sum(S_STR_NA(STR_AD == 1.6)) sum(S_STR_NA(STR_AD == 1.9)) ...
               sum(S_STR_NA(STR_AD == -1.9)) sum(S_STR_NA(STR_AD == -1.6)) sum(S_STR_NA(STR_AD == -1.3)) ...
               sum(S_STR_NA((STR_AD == -0.9) | (STR_AD == 0.3))) sum(S_STR_NA((STR_AD == -0.6) | (STR_AD == 0.6))) sum(S_STR_NA((STR_AD == -0.3) | (STR_AD == 0.9)))];

h = figure;
b = bar(1, STR_AD_Count, 'stacked', 'EdgeColor', 'none');
for i = 1:9
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_tri_div(i,:);
end
box off
axis off
exportgraphics(h, fullfile(figDir, 'conjunction_STR_bar.png'), 'Resolution', 300);

UKB_GP_AD = zeros(1,6000);
HCP_GP_AD = zeros(1,6000);
GP_AD = zeros(1,6000);

UKB_GP_AD(UKB.STAT.MOR.TSA.GP_SEX_B > 0) =  1;
UKB_GP_AD(UKB.STAT.MOR.TSA.GP_SEX_B < 0) = -1;
UKB_GP_AD((UKB.STAT.MOR.TSA.GP_SEX_Q < thresh) & (UKB.STAT.MOR.TSA.GP_SEX_B > 0))  =  2;
UKB_GP_AD((UKB.STAT.MOR.TSA.GP_SEX_Q < thresh) & (UKB.STAT.MOR.TSA.GP_SEX_B < 0))  = -2;

HCP_GP_AD(HCP.STAT.MOR.TSA.GP_SEX_B > 0) =  1;
HCP_GP_AD(HCP.STAT.MOR.TSA.GP_SEX_B < 0) = -1;
HCP_GP_AD((HCP.STAT.MOR.TSA.GP_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.GP_SEX_B > 0))  =  2;
HCP_GP_AD((HCP.STAT.MOR.TSA.GP_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.GP_SEX_B < 0))  = -2;

GP_AD((UKB_GP_AD == 2) & (HCP_GP_AD == 2)) = -1.9;
GP_AD((UKB_GP_AD == 2) & (HCP_GP_AD == 1)) = -1.6;
GP_AD((UKB_GP_AD == 1) & (HCP_GP_AD == 2)) = -1.6;
GP_AD((UKB_GP_AD == 1) & (HCP_GP_AD == 1)) = -1.3;

GP_AD((UKB_GP_AD == 2) & (HCP_GP_AD == -2)) = -0.9;
GP_AD((UKB_GP_AD == 2) & (HCP_GP_AD == -1)) = -0.6;
GP_AD((UKB_GP_AD == 1) & (HCP_GP_AD == -2)) = -0.6;
GP_AD((UKB_GP_AD == 1) & (HCP_GP_AD == -1)) = -0.3;

GP_AD((UKB_GP_AD == -2) & (HCP_GP_AD == 2)) = 0.3;
GP_AD((UKB_GP_AD == -2) & (HCP_GP_AD == 1)) = 0.6;
GP_AD((UKB_GP_AD == -1) & (HCP_GP_AD == 2)) = 0.6;
GP_AD((UKB_GP_AD == -1) & (HCP_GP_AD == 1)) = 0.9;

GP_AD((UKB_GP_AD == -2) & (HCP_GP_AD == -2)) = 1.3;
GP_AD((UKB_GP_AD == -2) & (HCP_GP_AD == -1)) = 1.6;
GP_AD((UKB_GP_AD == -1) & (HCP_GP_AD == -2)) = 1.6;
GP_AD((UKB_GP_AD == -1) & (HCP_GP_AD == -1)) = 1.9;

h = figure('WindowState', 'Maximized'); HyoMagetViewFlat(GP_AD, [-2 2], cmap.tri_div);
exportgraphics(h, fullfile(figDir, 'conjunction_GP.png'), 'Resolution', 300);

GP_AD_Count = [sum(S_GP_NA(GP_AD == 1.3)) sum(S_GP_NA(GP_AD == 1.6)) sum(S_GP_NA(GP_AD == 1.9)) ...
               sum(S_GP_NA(GP_AD == -1.9)) sum(S_GP_NA(GP_AD == -1.6)) sum(S_GP_NA(GP_AD == -1.3)) ...
               sum(S_GP_NA((GP_AD == -0.9) | (GP_AD == 0.3))) sum(S_GP_NA((GP_AD == -0.6) | (GP_AD == 0.6))) sum(S_GP_NA((GP_AD == -0.3) | (GP_AD == 0.9)))];
h = figure;
b = bar(1, GP_AD_Count, 'stacked', 'EdgeColor', 'none');
for i = 1:9
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_tri_div(i,:);
end
box off
axis off
exportgraphics(h, fullfile(figDir, 'conjunction_GP_bar.png'), 'Resolution', 300);

UKB_THAL_AD = zeros(1,13000);
HCP_THAL_AD = zeros(1,13000);
THAL_AD = zeros(1,13000);

UKB_THAL_AD(UKB.STAT.MOR.TSA.THAL_SEX_B > 0) =  1;
UKB_THAL_AD(UKB.STAT.MOR.TSA.THAL_SEX_B < 0) = -1;
UKB_THAL_AD((UKB.STAT.MOR.TSA.THAL_SEX_Q < thresh) & (UKB.STAT.MOR.TSA.THAL_SEX_B > 0))  =  2;
UKB_THAL_AD((UKB.STAT.MOR.TSA.THAL_SEX_Q < thresh) & (UKB.STAT.MOR.TSA.THAL_SEX_B < 0))  = -2;

HCP_THAL_AD(HCP.STAT.MOR.TSA.THAL_SEX_B > 0) =  1;
HCP_THAL_AD(HCP.STAT.MOR.TSA.THAL_SEX_B < 0) = -1;
HCP_THAL_AD((HCP.STAT.MOR.TSA.THAL_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.THAL_SEX_B > 0))  =  2;
HCP_THAL_AD((HCP.STAT.MOR.TSA.THAL_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.THAL_SEX_B < 0))  = -2;

THAL_AD((UKB_THAL_AD == 2) & (HCP_THAL_AD == 2)) = -1.9;
THAL_AD((UKB_THAL_AD == 2) & (HCP_THAL_AD == 1)) = -1.6;
THAL_AD((UKB_THAL_AD == 1) & (HCP_THAL_AD == 2)) = -1.6;
THAL_AD((UKB_THAL_AD == 1) & (HCP_THAL_AD == 1)) = -1.3;

THAL_AD((UKB_THAL_AD == 2) & (HCP_THAL_AD == -2)) = -0.9;
THAL_AD((UKB_THAL_AD == 2) & (HCP_THAL_AD == -1)) = -0.6;
THAL_AD((UKB_THAL_AD == 1) & (HCP_THAL_AD == -2)) = -0.6;
THAL_AD((UKB_THAL_AD == 1) & (HCP_THAL_AD == -1)) = -0.3;

THAL_AD((UKB_THAL_AD == -2) & (HCP_THAL_AD == 2)) = 0.3;
THAL_AD((UKB_THAL_AD == -2) & (HCP_THAL_AD == 1)) = 0.6;
THAL_AD((UKB_THAL_AD == -1) & (HCP_THAL_AD == 2)) = 0.6;
THAL_AD((UKB_THAL_AD == -1) & (HCP_THAL_AD == 1)) = 0.9;

THAL_AD((UKB_THAL_AD == -2) & (HCP_THAL_AD == -2)) = 1.3;
THAL_AD((UKB_THAL_AD == -2) & (HCP_THAL_AD == -1)) = 1.6;
THAL_AD((UKB_THAL_AD == -1) & (HCP_THAL_AD == -2)) = 1.6;
THAL_AD((UKB_THAL_AD == -1) & (HCP_THAL_AD == -1)) = 1.9;

h = figure('WindowState', 'Maximized'); HyoMagetViewFlat(THAL_AD, [-2 2], cmap.tri_div);
exportgraphics(h, fullfile(figDir, 'conjunction_THAL.png'), 'Resolution', 300);

THAL_AD_Count = [sum(S_THAL_NA(THAL_AD == 1.3)) sum(S_THAL_NA(THAL_AD == 1.6)) sum(S_THAL_NA(THAL_AD == 1.9)) ...
               sum(S_THAL_NA(THAL_AD == -1.9)) sum(S_THAL_NA(THAL_AD == -1.6)) sum(S_THAL_NA(THAL_AD == -1.3)) ...
               sum(S_THAL_NA((THAL_AD == -0.9) | (THAL_AD == 0.3))) sum(S_THAL_NA((THAL_AD == -0.6) | (THAL_AD == 0.6))) sum(S_THAL_NA((THAL_AD == -0.3) | (THAL_AD == 0.9)))];

h = figure;
b = bar(1, THAL_AD_Count, 'stacked', 'EdgeColor', 'none');
for i = 1:9
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_tri_div(i,:);
end
box off
axis off
exportgraphics(h, fullfile(figDir, 'conjunction_THAL_bar.png'), 'Resolution', 300);

UKB_HIPP_AD = zeros(1,2367);
HCP_HIPP_AD = zeros(1,2367);
HIPP_AD = zeros(1,2367);

UKB_HIPP_AD(UKB.STAT.MOR.TSA.HIPP_SEX_B > 0) =  1;
UKB_HIPP_AD(UKB.STAT.MOR.TSA.HIPP_SEX_B < 0) = -1;
UKB_HIPP_AD((UKB.STAT.MOR.TSA.HIPP_SEX_Q < thresh) & (UKB.STAT.MOR.TSA.HIPP_SEX_B > 0))  =  2;
UKB_HIPP_AD((UKB.STAT.MOR.TSA.HIPP_SEX_Q < thresh) & (UKB.STAT.MOR.TSA.HIPP_SEX_B < 0))  = -2;

HCP_HIPP_AD(HCP.STAT.MOR.TSA.HIPP_SEX_B > 0) =  1;
HCP_HIPP_AD(HCP.STAT.MOR.TSA.HIPP_SEX_B < 0) = -1;
HCP_HIPP_AD((HCP.STAT.MOR.TSA.HIPP_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.HIPP_SEX_B > 0))  =  2;
HCP_HIPP_AD((HCP.STAT.MOR.TSA.HIPP_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.HIPP_SEX_B < 0))  = -2;

HIPP_AD((UKB_HIPP_AD == 2) & (HCP_HIPP_AD == 2)) = -1.9;
HIPP_AD((UKB_HIPP_AD == 2) & (HCP_HIPP_AD == 1)) = -1.6;
HIPP_AD((UKB_HIPP_AD == 1) & (HCP_HIPP_AD == 2)) = -1.6;
HIPP_AD((UKB_HIPP_AD == 1) & (HCP_HIPP_AD == 1)) = -1.3;

HIPP_AD((UKB_HIPP_AD == 2) & (HCP_HIPP_AD == -2)) = -0.9;
HIPP_AD((UKB_HIPP_AD == 2) & (HCP_HIPP_AD == -1)) = -0.6;
HIPP_AD((UKB_HIPP_AD == 1) & (HCP_HIPP_AD == -2)) = -0.6;
HIPP_AD((UKB_HIPP_AD == 1) & (HCP_HIPP_AD == -1)) = -0.3;

HIPP_AD((UKB_HIPP_AD == -2) & (HCP_HIPP_AD == 2)) = 0.3;
HIPP_AD((UKB_HIPP_AD == -2) & (HCP_HIPP_AD == 1)) = 0.6;
HIPP_AD((UKB_HIPP_AD == -1) & (HCP_HIPP_AD == 2)) = 0.6;
HIPP_AD((UKB_HIPP_AD == -1) & (HCP_HIPP_AD == 1)) = 0.9;

HIPP_AD((UKB_HIPP_AD == -2) & (HCP_HIPP_AD == -2)) = 1.3;
HIPP_AD((UKB_HIPP_AD == -2) & (HCP_HIPP_AD == -1)) = 1.6;
HIPP_AD((UKB_HIPP_AD == -1) & (HCP_HIPP_AD == -2)) = 1.6;
HIPP_AD((UKB_HIPP_AD == -1) & (HCP_HIPP_AD == -1)) = 1.9;

h = figure('WindowState', 'Maximized'); HyoMagetViewFlat(HIPP_AD, [-2 2], cmap.tri_div);
exportgraphics(h, fullfile(figDir, 'conjunction_HIPP.png'), 'Resolution', 300);

HIPP_AD_Count = [sum(S_HIPP_NA(HIPP_AD == 1.3)) sum(S_HIPP_NA(HIPP_AD == 1.6)) sum(S_HIPP_NA(HIPP_AD == 1.9)) ...
               sum(S_HIPP_NA(HIPP_AD == -1.9)) sum(S_HIPP_NA(HIPP_AD == -1.6)) sum(S_HIPP_NA(HIPP_AD == -1.3)) ...
               sum(S_HIPP_NA((HIPP_AD == -0.9) | (HIPP_AD == 0.3))) sum(S_HIPP_NA((HIPP_AD == -0.6) | (HIPP_AD == 0.6))) sum(S_HIPP_NA((HIPP_AD == -0.3) | (HIPP_AD == 0.9)))];

h = figure;
b = bar(1, HIPP_AD_Count, 'stacked', 'EdgeColor', 'none');
for i = 1:9
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_tri_div(i,:);
end
box off
axis off
exportgraphics(h, fullfile(figDir, 'conjunction_HIPP_bar.png'), 'Resolution', 300);

UKB_AMYG_AD = zeros(1,2878);
HCP_AMYG_AD = zeros(1,2878);
AMYG_AD = zeros(1,2878);

UKB_AMYG_AD(UKB.STAT.MOR.TSA.AMYG_SEX_B > 0) =  1;
UKB_AMYG_AD(UKB.STAT.MOR.TSA.AMYG_SEX_B < 0) = -1;
UKB_AMYG_AD((UKB.STAT.MOR.TSA.AMYG_SEX_Q < thresh) & (UKB.STAT.MOR.TSA.AMYG_SEX_B > 0))  =  2;
UKB_AMYG_AD((UKB.STAT.MOR.TSA.AMYG_SEX_Q < thresh) & (UKB.STAT.MOR.TSA.AMYG_SEX_B < 0))  = -2;

HCP_AMYG_AD(HCP.STAT.MOR.TSA.AMYG_SEX_B > 0) =  1;
HCP_AMYG_AD(HCP.STAT.MOR.TSA.AMYG_SEX_B < 0) = -1;
HCP_AMYG_AD((HCP.STAT.MOR.TSA.AMYG_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.AMYG_SEX_B > 0))  =  2;
HCP_AMYG_AD((HCP.STAT.MOR.TSA.AMYG_SEX_Q < thresh) & (HCP.STAT.MOR.TSA.AMYG_SEX_B < 0))  = -2;

AMYG_AD((UKB_AMYG_AD == 2) & (HCP_AMYG_AD == 2)) = -1.9;
AMYG_AD((UKB_AMYG_AD == 2) & (HCP_AMYG_AD == 1)) = -1.6;
AMYG_AD((UKB_AMYG_AD == 1) & (HCP_AMYG_AD == 2)) = -1.6;
AMYG_AD((UKB_AMYG_AD == 1) & (HCP_AMYG_AD == 1)) = -1.3;

AMYG_AD((UKB_AMYG_AD == 2) & (HCP_AMYG_AD == -2)) = -0.9;
AMYG_AD((UKB_AMYG_AD == 2) & (HCP_AMYG_AD == -1)) = -0.6;
AMYG_AD((UKB_AMYG_AD == 1) & (HCP_AMYG_AD == -2)) = -0.6;
AMYG_AD((UKB_AMYG_AD == 1) & (HCP_AMYG_AD == -1)) = -0.3;

AMYG_AD((UKB_AMYG_AD == -2) & (HCP_AMYG_AD == 2)) = 0.3;
AMYG_AD((UKB_AMYG_AD == -2) & (HCP_AMYG_AD == 1)) = 0.6;
AMYG_AD((UKB_AMYG_AD == -1) & (HCP_AMYG_AD == 2)) = 0.6;
AMYG_AD((UKB_AMYG_AD == -1) & (HCP_AMYG_AD == 1)) = 0.9;

AMYG_AD((UKB_AMYG_AD == -2) & (HCP_AMYG_AD == -2)) = 1.3;
AMYG_AD((UKB_AMYG_AD == -2) & (HCP_AMYG_AD == -1)) = 1.6;
AMYG_AD((UKB_AMYG_AD == -1) & (HCP_AMYG_AD == -2)) = 1.6;
AMYG_AD((UKB_AMYG_AD == -1) & (HCP_AMYG_AD == -1)) = 1.9;

h = figure('WindowState', 'Maximized'); HyoMagetViewFlat(AMYG_AD, [-2 2], cmap.tri_div);
exportgraphics(h, fullfile(figDir, 'conjunction_AMYG.png'), 'Resolution', 300);

AMYG_AD_Count = [sum(S_AMYG_NA(AMYG_AD == 1.3)) sum(S_AMYG_NA(AMYG_AD == 1.6)) sum(S_AMYG_NA(AMYG_AD == 1.9)) ...
               sum(S_AMYG_NA(AMYG_AD == -1.9)) sum(S_AMYG_NA(AMYG_AD == -1.6)) sum(S_AMYG_NA(AMYG_AD == -1.3)) ...
               sum(S_AMYG_NA((AMYG_AD == -0.9) | (AMYG_AD == 0.3))) sum(S_AMYG_NA((AMYG_AD == -0.6) | (AMYG_AD == 0.6))) sum(S_AMYG_NA((AMYG_AD == -0.3) | (AMYG_AD == 0.9)))];

h = figure;
b = bar(1, AMYG_AD_Count, 'stacked', 'EdgeColor', 'none');
for i = 1:9
    b(i).BarWidth = 0.1;
    b(i).FaceColor = cmap.bar_tri_div(i,:);
end
box off
axis off
exportgraphics(h, fullfile(figDir, 'conjunction_AMYG_bar.png'), 'Resolution', 300);

close all;

sum(STR_AD_Count(1:6)) ./ sum(STR_AD_Count)
sum(GP_AD_Count(1:6)) ./ sum(GP_AD_Count)
sum(THAL_AD_Count(1:6)) ./ sum(THAL_AD_Count)
sum(HIPP_AD_Count(1:6)) ./ sum(HIPP_AD_Count)
sum(AMYG_AD_Count(1:6)) ./ sum(AMYG_AD_Count)

%% Spatial permutation testing for HCP and UKB

% Categorize spatially permuted HCP sex-difference maps
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
PERM_STR_AD((UKB_STR_AD' == 2) & (PERM_HCP_STR_AD == 2)) = -1.9;
PERM_STR_AD((UKB_STR_AD' == 2) & (PERM_HCP_STR_AD == 1)) = -1.6;
PERM_STR_AD((UKB_STR_AD' == 1) & (PERM_HCP_STR_AD == 2)) = -1.6;
PERM_STR_AD((UKB_STR_AD' == 1) & (PERM_HCP_STR_AD == 1)) = -1.3;
PERM_STR_AD((UKB_STR_AD' == 2) & (PERM_HCP_STR_AD == -2)) = -0.9;
PERM_STR_AD((UKB_STR_AD' == 2) & (PERM_HCP_STR_AD == -1)) = -0.6;
PERM_STR_AD((UKB_STR_AD' == 1) & (PERM_HCP_STR_AD == -2)) = -0.6;
PERM_STR_AD((UKB_STR_AD' == 1) & (PERM_HCP_STR_AD == -1)) = -0.3;
PERM_STR_AD((UKB_STR_AD' == -2) & (PERM_HCP_STR_AD == 2)) = 0.3;
PERM_STR_AD((UKB_STR_AD' == -2) & (PERM_HCP_STR_AD == 1)) = 0.6;
PERM_STR_AD((UKB_STR_AD' == -1) & (PERM_HCP_STR_AD == 2)) = 0.6;
PERM_STR_AD((UKB_STR_AD' == -1) & (PERM_HCP_STR_AD == 1)) = 0.9;
PERM_STR_AD((UKB_STR_AD' == -2) & (PERM_HCP_STR_AD == -2)) = 1.3;
PERM_STR_AD((UKB_STR_AD' == -2) & (PERM_HCP_STR_AD == -1)) = 1.6;
PERM_STR_AD((UKB_STR_AD' == -1) & (PERM_HCP_STR_AD == -2)) = 1.6;
PERM_STR_AD((UKB_STR_AD' == -1) & (PERM_HCP_STR_AD == -1)) = 1.9;

PERM_GP_AD = zeros(6000, nPerm);
PERM_GP_AD((UKB_GP_AD' == 2) & (PERM_HCP_GP_AD == 2)) = -1.9;
PERM_GP_AD((UKB_GP_AD' == 2) & (PERM_HCP_GP_AD == 1)) = -1.6;
PERM_GP_AD((UKB_GP_AD' == 1) & (PERM_HCP_GP_AD == 2)) = -1.6;
PERM_GP_AD((UKB_GP_AD' == 1) & (PERM_HCP_GP_AD == 1)) = -1.3;
PERM_GP_AD((UKB_GP_AD' == 2) & (PERM_HCP_GP_AD == -2)) = -0.9;
PERM_GP_AD((UKB_GP_AD' == 2) & (PERM_HCP_GP_AD == -1)) = -0.6;
PERM_GP_AD((UKB_GP_AD' == 1) & (PERM_HCP_GP_AD == -2)) = -0.6;
PERM_GP_AD((UKB_GP_AD' == 1) & (PERM_HCP_GP_AD == -1)) = -0.3;
PERM_GP_AD((UKB_GP_AD' == -2) & (PERM_HCP_GP_AD == 2)) = 0.3;
PERM_GP_AD((UKB_GP_AD' == -2) & (PERM_HCP_GP_AD == 1)) = 0.6;
PERM_GP_AD((UKB_GP_AD' == -1) & (PERM_HCP_GP_AD == 2)) = 0.6;
PERM_GP_AD((UKB_GP_AD' == -1) & (PERM_HCP_GP_AD == 1)) = 0.9;
PERM_GP_AD((UKB_GP_AD' == -2) & (PERM_HCP_GP_AD == -2)) = 1.3;
PERM_GP_AD((UKB_GP_AD' == -2) & (PERM_HCP_GP_AD == -1)) = 1.6;
PERM_GP_AD((UKB_GP_AD' == -1) & (PERM_HCP_GP_AD == -2)) = 1.6;
PERM_GP_AD((UKB_GP_AD' == -1) & (PERM_HCP_GP_AD == -1)) = 1.9;

PERM_THAL_AD = zeros(13000, nPerm);
PERM_THAL_AD((UKB_THAL_AD' == 2) & (PERM_HCP_THAL_AD == 2)) = -1.9;
PERM_THAL_AD((UKB_THAL_AD' == 2) & (PERM_HCP_THAL_AD == 1)) = -1.6;
PERM_THAL_AD((UKB_THAL_AD' == 1) & (PERM_HCP_THAL_AD == 2)) = -1.6;
PERM_THAL_AD((UKB_THAL_AD' == 1) & (PERM_HCP_THAL_AD == 1)) = -1.3;
PERM_THAL_AD((UKB_THAL_AD' == 2) & (PERM_HCP_THAL_AD == -2)) = -0.9;
PERM_THAL_AD((UKB_THAL_AD' == 2) & (PERM_HCP_THAL_AD == -1)) = -0.6;
PERM_THAL_AD((UKB_THAL_AD' == 1) & (PERM_HCP_THAL_AD == -2)) = -0.6;
PERM_THAL_AD((UKB_THAL_AD' == 1) & (PERM_HCP_THAL_AD == -1)) = -0.3;
PERM_THAL_AD((UKB_THAL_AD' == -2) & (PERM_HCP_THAL_AD == 2)) = 0.3;
PERM_THAL_AD((UKB_THAL_AD' == -2) & (PERM_HCP_THAL_AD == 1)) = 0.6;
PERM_THAL_AD((UKB_THAL_AD' == -1) & (PERM_HCP_THAL_AD == 2)) = 0.6;
PERM_THAL_AD((UKB_THAL_AD' == -1) & (PERM_HCP_THAL_AD == 1)) = 0.9;
PERM_THAL_AD((UKB_THAL_AD' == -2) & (PERM_HCP_THAL_AD == -2)) = 1.3;
PERM_THAL_AD((UKB_THAL_AD' == -2) & (PERM_HCP_THAL_AD == -1)) = 1.6;
PERM_THAL_AD((UKB_THAL_AD' == -1) & (PERM_HCP_THAL_AD == -2)) = 1.6;
PERM_THAL_AD((UKB_THAL_AD' == -1) & (PERM_HCP_THAL_AD == -1)) = 1.9;

PERM_HIPP_AD = zeros(2367, nPerm);
PERM_HIPP_AD((UKB_HIPP_AD' == 2) & (PERM_HCP_HIPP_AD == 2)) = -1.9;
PERM_HIPP_AD((UKB_HIPP_AD' == 2) & (PERM_HCP_HIPP_AD == 1)) = -1.6;
PERM_HIPP_AD((UKB_HIPP_AD' == 1) & (PERM_HCP_HIPP_AD == 2)) = -1.6;
PERM_HIPP_AD((UKB_HIPP_AD' == 1) & (PERM_HCP_HIPP_AD == 1)) = -1.3;
PERM_HIPP_AD((UKB_HIPP_AD' == 2) & (PERM_HCP_HIPP_AD == -2)) = -0.9;
PERM_HIPP_AD((UKB_HIPP_AD' == 2) & (PERM_HCP_HIPP_AD == -1)) = -0.6;
PERM_HIPP_AD((UKB_HIPP_AD' == 1) & (PERM_HCP_HIPP_AD == -2)) = -0.6;
PERM_HIPP_AD((UKB_HIPP_AD' == 1) & (PERM_HCP_HIPP_AD == -1)) = -0.3;
PERM_HIPP_AD((UKB_HIPP_AD' == -2) & (PERM_HCP_HIPP_AD == 2)) = 0.3;
PERM_HIPP_AD((UKB_HIPP_AD' == -2) & (PERM_HCP_HIPP_AD == 1)) = 0.6;
PERM_HIPP_AD((UKB_HIPP_AD' == -1) & (PERM_HCP_HIPP_AD == 2)) = 0.6;
PERM_HIPP_AD((UKB_HIPP_AD' == -1) & (PERM_HCP_HIPP_AD == 1)) = 0.9;
PERM_HIPP_AD((UKB_HIPP_AD' == -2) & (PERM_HCP_HIPP_AD == -2)) = 1.3;
PERM_HIPP_AD((UKB_HIPP_AD' == -2) & (PERM_HCP_HIPP_AD == -1)) = 1.6;
PERM_HIPP_AD((UKB_HIPP_AD' == -1) & (PERM_HCP_HIPP_AD == -2)) = 1.6;
PERM_HIPP_AD((UKB_HIPP_AD' == -1) & (PERM_HCP_HIPP_AD == -1)) = 1.9;

PERM_AMYG_AD = zeros(2878, nPerm);
PERM_AMYG_AD((UKB_AMYG_AD' == 2) & (PERM_HCP_AMYG_AD == 2)) = -1.9;
PERM_AMYG_AD((UKB_AMYG_AD' == 2) & (PERM_HCP_AMYG_AD == 1)) = -1.6;
PERM_AMYG_AD((UKB_AMYG_AD' == 1) & (PERM_HCP_AMYG_AD == 2)) = -1.6;
PERM_AMYG_AD((UKB_AMYG_AD' == 1) & (PERM_HCP_AMYG_AD == 1)) = -1.3;
PERM_AMYG_AD((UKB_AMYG_AD' == 2) & (PERM_HCP_AMYG_AD == -2)) = -0.9;
PERM_AMYG_AD((UKB_AMYG_AD' == 2) & (PERM_HCP_AMYG_AD == -1)) = -0.6;
PERM_AMYG_AD((UKB_AMYG_AD' == 1) & (PERM_HCP_AMYG_AD == -2)) = -0.6;
PERM_AMYG_AD((UKB_AMYG_AD' == 1) & (PERM_HCP_AMYG_AD == -1)) = -0.3;
PERM_AMYG_AD((UKB_AMYG_AD' == -2) & (PERM_HCP_AMYG_AD == 2)) = 0.3;
PERM_AMYG_AD((UKB_AMYG_AD' == -2) & (PERM_HCP_AMYG_AD == 1)) = 0.6;
PERM_AMYG_AD((UKB_AMYG_AD' == -1) & (PERM_HCP_AMYG_AD == 2)) = 0.6;
PERM_AMYG_AD((UKB_AMYG_AD' == -1) & (PERM_HCP_AMYG_AD == 1)) = 0.9;
PERM_AMYG_AD((UKB_AMYG_AD' == -2) & (PERM_HCP_AMYG_AD == -2)) = 1.3;
PERM_AMYG_AD((UKB_AMYG_AD' == -2) & (PERM_HCP_AMYG_AD == -1)) = 1.6;
PERM_AMYG_AD((UKB_AMYG_AD' == -1) & (PERM_HCP_AMYG_AD == -2)) = 1.6;
PERM_AMYG_AD((UKB_AMYG_AD' == -1) & (PERM_HCP_AMYG_AD == -1)) = 1.9;

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

%% Spatial correlation

% Calculate correlations
[R_STR ~] = corr(HCP.STAT.MOR.TSA.STR_SEX_B', UKB.STAT.MOR.TSA.STR_SEX_B');
[R_GP ~] = corr(HCP.STAT.MOR.TSA.GP_SEX_B', UKB.STAT.MOR.TSA.GP_SEX_B');
[R_THAL ~] = corr(HCP.STAT.MOR.TSA.THAL_SEX_B', UKB.STAT.MOR.TSA.THAL_SEX_B');
[R_HIPP ~] = corr(HCP.STAT.MOR.TSA.HIPP_SEX_B', UKB.STAT.MOR.TSA.HIPP_SEX_B');
[R_AMYG ~] = corr(HCP.STAT.MOR.TSA.AMYG_SEX_B', UKB.STAT.MOR.TSA.AMYG_SEX_B');

[R_STR R_GP R_THAL R_HIPP R_AMYG]

% Null distributions for spatial correlations
for i = 1:nPerm
    PERM_R_STR(i) = corr(HCP.STAT.MOR.MSR.STR_SEX_B(:,i), UKB.STAT.MOR.TSA.STR_SEX_B');
    PERM_R_GP(i) = corr(HCP.STAT.MOR.MSR.GP_SEX_B(:,i), UKB.STAT.MOR.TSA.GP_SEX_B');
    PERM_R_THAL(i) = corr(HCP.STAT.MOR.MSR.THAL_SEX_B(:,i), UKB.STAT.MOR.TSA.THAL_SEX_B');
    PERM_R_HIPP(i) = corr(HCP.STAT.MOR.MSR.HIPP_SEX_B(:,i), UKB.STAT.MOR.TSA.HIPP_SEX_B');
    PERM_R_AMYG(i) = corr(HCP.STAT.MOR.MSR.AMYG_SEX_B(:,i), UKB.STAT.MOR.TSA.AMYG_SEX_B');
end

disp('statistics')
sum(R_STR < PERM_R_STR) / nPerm
sum(R_GP < PERM_R_GP) / nPerm
sum(R_THAL < PERM_R_THAL) / nPerm
sum(R_HIPP < PERM_R_HIPP) / nPerm
sum(R_AMYG < PERM_R_AMYG) / nPerm

% Vectors
CMP_STR = [HCP.STAT.MOR.TSA.STR_SEX_B; UKB.STAT.MOR.TSA.STR_SEX_B];
CMP_GP = [HCP.STAT.MOR.TSA.GP_SEX_B; UKB.STAT.MOR.TSA.GP_SEX_B];
CMP_THAL = [HCP.STAT.MOR.TSA.THAL_SEX_B; UKB.STAT.MOR.TSA.THAL_SEX_B];
CMP_HIPP = [HCP.STAT.MOR.TSA.HIPP_SEX_B; UKB.STAT.MOR.TSA.HIPP_SEX_B];
CMP_AMYG = [HCP.STAT.MOR.TSA.AMYG_SEX_B; UKB.STAT.MOR.TSA.AMYG_SEX_B];
IDL = [-1:1e-5:1; -1:1e-5:1];

% Calculate residuals from the identity line
for i = 1:size(CMP_STR,2); ID_CMP_STR(i) = min(sqrt(sum((CMP_STR(:,i) - IDL).^2))); end
for i = 1:size(CMP_GP,2); ID_CMP_GP(i) = min(sqrt(sum((CMP_GP(:,i) - IDL).^2))); end
for i = 1:size(CMP_THAL,2); ID_CMP_THAL(i) = min(sqrt(sum((CMP_THAL(:,i) - IDL).^2))); end
for i = 1:size(CMP_HIPP,2); ID_CMP_HIPP(i) = min(sqrt(sum((CMP_HIPP(:,i) - IDL).^2))); end
for i = 1:size(CMP_AMYG,2); ID_CMP_AMYG(i) = min(sqrt(sum((CMP_AMYG(:,i) - IDL).^2))); end

% Test the significance of slope - "1"
H0 = [0 1];  c = 1;
mdl_STR = fitlm(HCP.STAT.MOR.TSA.STR_SEX_B, UKB.STAT.MOR.TSA.STR_SEX_B);
mdl_GP = fitlm(HCP.STAT.MOR.TSA.GP_SEX_B, UKB.STAT.MOR.TSA.GP_SEX_B);
mdl_THAL = fitlm(HCP.STAT.MOR.TSA.THAL_SEX_B, UKB.STAT.MOR.TSA.THAL_SEX_B);
mdl_HIPP = fitlm(HCP.STAT.MOR.TSA.HIPP_SEX_B, UKB.STAT.MOR.TSA.HIPP_SEX_B);
mdl_AMYG = fitlm(HCP.STAT.MOR.TSA.AMYG_SEX_B, UKB.STAT.MOR.TSA.AMYG_SEX_B);

disp('Statistics for slope vs 1')
[mdl_STR.Coefficients.Estimate(2) mdl_GP.Coefficients.Estimate(2) mdl_THAL.Coefficients.Estimate(2) mdl_HIPP.Coefficients.Estimate(2) mdl_AMYG.Coefficients.Estimate(2)]
5*[coefTest(mdl_STR, H0, c) coefTest(mdl_GP, H0, c) coefTest(mdl_THAL, H0, c) coefTest(mdl_HIPP, H0, c) coefTest(mdl_AMYG, H0, c)]

% Hexaplots
h = figure('WindowState', 'maximized'); hold on;
hexscatter(HCP.STAT.MOR.TSA.STR_SEX_B, UKB.STAT.MOR.TSA.STR_SEX_B, 'res', 50, 'xlim', [-1 1], 'ylim', [-1 1]);
scatter(HCP.STAT.MOR.TSA.STR_SEX_B, UKB.STAT.MOR.TSA.STR_SEX_B, 5, [0 0.5 0.5], 'filled', 'MarkerEdgeColor', [0 0.5 0.5], 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0);
xlim([-1 1]); xticks(-1:0.25:1); ylim([-1 1]); yticks(-1:0.25:1); set(gca,'TickDir','out');
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--'; b.LineWidth = 2;
x = xline(0); x.LineStyle = '--'; x.LineWidth = 2; y = yline(0); y.LineStyle = '--'; y.LineWidth = 2;
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xlim([-1 1]); ylim([-1 1]);
axis square;
exportgraphics(h, fullfile(figDir, 'spatialcorr_str.png'), 'Resolution', 300);

h = figure('WindowState', 'maximized'); hold on;
hexscatter(HCP.STAT.MOR.TSA.GP_SEX_B, UKB.STAT.MOR.TSA.GP_SEX_B, 'res', 50, 'xlim', [-1 1], 'ylim', [-1 1]);
scatter(HCP.STAT.MOR.TSA.GP_SEX_B, UKB.STAT.MOR.TSA.GP_SEX_B, 5, [0 0.5 0.5], 'filled', 'MarkerEdgeColor', [0 0.5 0.5], 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0);
xlim([-1 1]); xticks(-1:0.25:1); ylim([-1 1]); yticks(-1:0.25:1); set(gca,'TickDir','out');
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--'; b.LineWidth = 2;
x = xline(0); x.LineStyle = '--'; x.LineWidth = 2; y = yline(0); y.LineStyle = '--'; y.LineWidth = 2;
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xlim([-1 1]); ylim([-1 1]);
axis square;
exportgraphics(h, fullfile(figDir, 'spatialcorr_gp.png'), 'Resolution', 300);

h = figure('WindowState', 'maximized'); hold on;
hexscatter(HCP.STAT.MOR.TSA.THAL_SEX_B, UKB.STAT.MOR.TSA.THAL_SEX_B, 'res', 50, 'xlim', [-1 1], 'ylim', [-1 1]);
scatter(HCP.STAT.MOR.TSA.THAL_SEX_B, UKB.STAT.MOR.TSA.THAL_SEX_B, 5, [0 0.5 0.5], 'filled', 'MarkerEdgeColor', [0 0.5 0.5], 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0);
xlim([-1 1]); xticks(-1:0.25:1); ylim([-1 1]); yticks(-1:0.25:1); set(gca,'TickDir','out');
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--'; b.LineWidth = 2;
x = xline(0); x.LineStyle = '--'; x.LineWidth = 2; y = yline(0); y.LineStyle = '--'; y.LineWidth = 2;
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xlim([-1 1]); ylim([-1 1]);
axis square;
exportgraphics(h, fullfile(figDir, 'spatialcorr_thal.png'), 'Resolution', 300);

h = figure('WindowState', 'maximized'); hold on;
hexscatter(HCP.STAT.MOR.TSA.HIPP_SEX_B, UKB.STAT.MOR.TSA.HIPP_SEX_B, 'res', 50, 'xlim', [-1 1], 'ylim', [-1 1]);
scatter(HCP.STAT.MOR.TSA.HIPP_SEX_B, UKB.STAT.MOR.TSA.HIPP_SEX_B, 5, [0 0.5 0.5], 'filled', 'MarkerEdgeColor', [0 0.5 0.5], 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0);
xlim([-1 1]); xticks(-1:0.25:1); ylim([-1 1]); yticks(-1:0.25:1); set(gca,'TickDir','out');
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--'; b.LineWidth = 2;
x = xline(0); x.LineStyle = '--'; x.LineWidth = 2; y = yline(0); y.LineStyle = '--'; y.LineWidth = 2;
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xlim([-1 1]); ylim([-1 1]);
axis square;
exportgraphics(h, fullfile(figDir, 'spatialcorr_hipp.png'), 'Resolution', 300);

h = figure('WindowState', 'maximized'); hold on;
hexscatter(HCP.STAT.MOR.TSA.AMYG_SEX_B, UKB.STAT.MOR.TSA.AMYG_SEX_B, 'res', 50, 'xlim', [-1 1], 'ylim', [-1 1]);
scatter(HCP.STAT.MOR.TSA.AMYG_SEX_B, UKB.STAT.MOR.TSA.AMYG_SEX_B, 5, [0 0.5 0.5], 'filled', 'MarkerEdgeColor', [0 0.5 0.5], 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 0);
xlim([-1 1]); xticks(-1:0.25:1); ylim([-1 1]); yticks(-1:0.25:1); set(gca,'TickDir','out');
b = refline(1,0); b.Color = [0 0 0]; b.LineStyle = '--'; b.LineWidth = 2;
x = xline(0); x.LineStyle = '--'; x.LineWidth = 2; y = yline(0); y.LineStyle = '--'; y.LineWidth = 2;
l = lsline; l.Color = [0 0 0]; l.LineWidth = 2;
xlim([-1 1]); ylim([-1 1]);
axis square;
exportgraphics(h, fullfile(figDir, 'spatialcorr_amyg.png'), 'Resolution', 300);

close all;
