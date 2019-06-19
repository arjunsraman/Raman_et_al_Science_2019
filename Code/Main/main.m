%%%%%%%%%%%
% Identification of 'ecogroup' through taxon-taxon co-variation as appplied
% to a healthy cohort of Bangladeshi children sampled from postnatal months
% 1 to 60. 


% Authors: Arjun Raman, Jeffrey Gordon, Washington University School of Medicine
% Date: 07/12/2019
% Manuscript Title: A sparse co-varying unit of the human gut microbiota
% that describes healthy and impaired community development
% 
% Brief statement of problem: The function of biological systems is emergent, arising from the collective action of 
% variables that comprise the system. Longitudinally sampled cohorts offer a unique opportunity to 
% potentially identify which taxa encode important information regarding important microbial community properties. 
% The workflow below identifies taxa that co-vary with each other ('ecogroup') in a temporally reproducible fashion, i.e. temporally conserved co-varying taxa.

% Input for 'needed_workspace.mat': 
% Data matrix, 3-D matrix where
%     Rows (i): Fecal samples
%     Columns (j): Taxa
%     Z-dimension (k): Time
%     Element: Fractional abundance of taxa (j) in fecal sample (i) at time (k).

% OTU_ID: Cell array where each entry is the ID of the OTUs in the same order as the columns in Data Matrix
% OTU_names: Cell array where each entry is the taxonomic name in the same order as the columns in Data Matrix

% Functions needed:
% serial_eig.m (For Step 1)
% calc_Cij_bin_t (For Step 2)

% 
% Steps
%     Step 1: Microbiota organizational dynamics over time
%         
%     Step 2: Creation and eigendecomposition of binarized taxon-taxon covariance matrix averaged over time
%     
%     Step 3: 'Ecogroup' identification
%     
%     
%     
%     
%%%%%%%%%%%
%% STEP 1: Microbiota organizational dynamics over time


prompt = 'What is the timeframe of interest?';
timeframe = input(prompt); %User defines timeframe of interest
    %For Raman et al., answer is [1 36] and [1 60]
    
data_trm = data(:,:,timeframe(1,1):timeframe(1,2)); %Isolate data for specified timeframe
v_t = serial_eig(data_trm,1,length(data_trm(1,1,:))); %Eigenspectra as a function of time

%Plot first eigenvector through time
figure; plot(fliplr(v_t(1,:)));

clear timeframe
clear prompt
clear data_trm

%% STEP 2: Creation and eigendecomposition of binarized taxon-taxon covariance matrix averaged over time
        

warning('off','all');
prompt = 'Correlation(1) or Covariance (2)?';
decision = input(prompt)
if decision == 1; %Correlation
    clear prompt
    prompt = 'What is desired pseudo_wt (0.01 for Raman et al)?'
    pseudo_wt = input(prompt);
    [Cij_bin_t,taxa_names,OTU_ID_trm,y_frac,ev_wts] = calc_Cij_bin_t_corr(data,OTU_ID,OTU_names,pseudo_wt);
    clear prompt
elseif decision == 2; %Covariance
    [Cij_bin_t,taxa_names,OTU_ID_trm,y_frac,ev_wts] = calc_Cij_bin_t_cov(data,OTU_ID,OTU_names);
end;
clear prompt
    
%I/O Explanation
% Input
% data: 3D matrix of dimensions l x n x t timeseries data where rows are fecal samples
% (dimension l), columns are taxa (dimension n), 3rd dimension is time
% (dimension t)
% OTU_ID: cell array (nx1 dimensions) where each entry is the OTU ID in the same order as columns for 'data'
% OTU_names: cell array (nx1 dimensions) where each entry is the OTU name in the same order as columsn for 'data'
% 
% Output
% Cij_bin_t: 2D matrix (mxm dimensions, where m<n) where each row and column is a taxa and each element is the temporally conserved covariance between taxa i and j
% taxa_names: cell array (mx1 dimensions) where each entry is the taxa name in the same order as Cij_bin_t
% OTU_ID_trm: cell array (mx1 dimensions) where each entry is the taxa OTU ID in the same order as Cij_bin_t
% y_frac: mx1 matrix where each entry is the fractional data variance carried by eigenvector m
% ev_wts: mxk matrix. rows specify taxa, k specify eigenvectors (as specified by user). Each entry in matrix specifies weighting of taxa m on eigenvector k

%Prompt information
%       User is prompted to answer three questions
% Q1: 'What is the timeframe of interest for ecogroup identification?'
%     This is based on Step 1 of this script. For Raman et al, we choose months 20 to 60 as month 20 is when the microbiota achieves a stable configuration as based by the dynamics of eigenvector 1 through time
%     When prompted type in [20 60].
% 
% Q2: 'At what threshold is desired cutoff for taxon inclusion?'
%     For each month a hierarchically clustered covariance matrix is computed. Due to the sparsity of covariance secondary to the abundance of absent or lowly abundant taxa, particularly in the early months
%     there exist many 0 entries in this matrix. In order to ease distribution fitting (see next prompt), a conservative threshold for taxon inclusion is implemented
%     thereby reducing the number of taxa considered each month. For Raman et al, this value was set at 20. This threshold takes the 20 clustered taxa that occupy the bottom right diagonal
%     of the hierarchically clustered covariance matrix per month. Over all months, this equates to 80 out of a possible of 118 taxa. 
%     
% Q3: 'What percentage cutoff is desired?'
%     After monthly normalization, each monthly covariance matrix is fit to a t-location scale distribution as the monthly taxon-taxon covariance distribution
%     exhibits significant positive and negative tails. For each month, a cutoff is applied for binarization. In Raman et al, this cutoff is 0.1, equating to the top and bottom 10%
%     of the CDF as determined by a t-location scale distribution. As a goal of this analysis is to average the most coupled taxa over time, this binarization
%     step is necessary so as to not average temporally de-correlated positive and negative taxon-taxon covariance that might articially produce an average covariance of close to zero.
%     For Raman et al, this cutoff was set to 0.1 and adjusted to 0.05, 0.15, and 0.075 for sensitivity analysis. 

%The output ev_wts provide the necessary projections along eigenvectors to
%define the ecogroup. However, ecogroup identification requires empiric fitting of
%distribution and cutoff applied to the CDF of projections along a single, or multiple eigenctors. 

%% STEP 3: 'Ecogroup' Identification
%This is dependent on a) how many eigenvectors are deemed significant and
%b) the threshold of projection significance applied to each eigenvector.
%For the case in Raman et al, we take the first eigenvector to be
%significant as it holds ~80% of the variance, fit the taxa projections
%along ev1 to a generalized extreme value distribution ('gev') and identify
%significantly co-varying taxa as those taxa that project onto the extreme values of ev1. 
%For 'ecogroup' definition we use the cutoff of ev1_projection < -0.15 which
%corresponds to projection within the top 20%.  

pd = fitdist(ev_wts(:,1),'gev');

cutoff_ecogroup(1,1) = icdf(pd,0.15);
cutoff_ecogroup(1,2) = icdf(pd,0.175);
cutoff_ecogroup(1,3) = icdf(pd,0.2);











