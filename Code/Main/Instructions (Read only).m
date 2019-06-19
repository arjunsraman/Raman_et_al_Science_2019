%Run Instructions prompt by prompt for generating 'ecogroup' as defined in
%Raman et al.

%Open main.m (this is the main matlab script used for calculations that
%will operate on 'needed_workspace.mat')

%Prompt 1 (for Step 1): What is the timeframe of interest? 
%Answer: [1 36]
%Matlab output: (1)Plot shown in fig. S4B. (2) v_t array where each row is
%an eigenvector and each column corresponds to a timepoint with column 1
%corresponding to timepoint 1. Each cell is the value of eigenvector n at
%timepoint m.

%Prompt 2 (for Step 2): Correlation(1) or Covariance (2)?
%Answer: 2

%Prompt 3 (for Step 2): What is the timeframe of interest for ecogroup
%identification?
%Answer: [20 60]

%Prompt 4 (for Step 2): At what threshold is desired cutoff for taxon
%inclusion?
%Answer: 20
%Matlab output: Dendrogram shown in fig. S5C

%Prompt 5 (for Step 2): What percentage cutoff is desired?
%Answer: 0.1; 
%Matlab output: Eigenspectrum of temporally conserved covariance matrix
%shown in Fig. 1C inset

%Prompt 6 (for Step 2): Which eigenvectors to use for ecogroup
%identification?
%Answer: 1
%Matlab output: 
    %Cij_bin_t: Temporally conserved covariance matrix (80x80,
    %shown in Fig. 1B)
    
    %decision: 2 (Indicative of choice of covariance instead of
    %correlation)
    
    %ev_wts: Eigenprojection of the 80 taxa onto eigenvector 1 of the
    %temporally conserved covariance matrix (shown in Fig. 1C as histogram)
    
    %OTU_ID: 118 OTU ID designations
    %OTU_ID_trm: 80 OTU that define 80x80 Cij_bin_t (in same order)
    %taxa_names: Names of 80 OTUs
    %y_frac: Eigenspectrum
    
    
    
