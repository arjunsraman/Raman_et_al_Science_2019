%%To calculate time average covariance matrix, names and OTU ID of taxa
%%contributing to this matrix, fractional variance carried by eigenvectors
%%of this matrix, and taxa projections onto eigenvectors of choice.

function [Cij_bin_t,taxa_names, OTU_ID_trm,y_frac,ev_wts] = calc_Cij_bin_t(data,OTU_ID,OTU_names);

prompt = 'What is the timeframe of interest for ecogroup identification?';
timepoint = input(prompt);
    %For Raman et al, answer is [20 60]
data_trm = data(:,:,timepoint(1,1):timepoint(1,2));
for i=1:length(data_trm(1,1,:));
    cov_data(:,:,i) = cov(data_trm(:,find(sum(data_trm(:,:,i) > 0)),i)'); %Creates a 118x118 matrix
    Z = linkage(cov_data(:,:,i)); [H,T,reordered_indices(i,:)] = dendrogram(Z,118); 
    cov_data_reordered(:,:,i) = cov_data(reordered_indices(i,:),reordered_indices(i,:));
    cov_data_reordered(:,:,i) = cov_data_reordered(:,:,i)./max(max(cov_data_reordered(:,:,i))); %Creates 118x118 matrix that is reordered according to Hierarchical clustering and normalized relative to max monthly covariance
    %As example, dendrogram for month 60 will show
end;
clear prompt

prompt = 'At what threshold is desired cutoff for taxon inclusion?';
thrsh = input(prompt);
    %For Raman et al, threshold was top 20 of hierarchically clustered taxa for each month
monthly_relevant_taxa_covariance = reordered_indices(:,end-thrsh:end);
monthly_relevant_taxa_unique = unique(monthly_relevant_taxa_covariance);
cov_data_uniq = cov_data(monthly_relevant_taxa_unique,monthly_relevant_taxa_unique,:); %For the case of Raman et al, this creates a 76x76x41 matrix
cov_data_uniq_norm = cov_data_uniq./max(max(cov_data_uniq));
taxa_names = OTU_names(monthly_relevant_taxa_unique);
OTU_ID_trm = OTU_ID(monthly_relevant_taxa_unique);

%Now apply cutoffs
clear prompt
prompt = 'What percentage cutoff is desired?';
    %For Raman et al, 0.1 was applied
percentage_cutoff = input(prompt);
for i=1:length(cov_data_reordered(1,1,:)); 
    test = cov_data_reordered(end-thrsh:end,end-thrsh:end,i);
    pd = fitdist(test(:),'tLocationScale'); %Type of distribution can be altered depending on nature of data. 
    %For Raman et al, heavy-tailed distributions naturally lent to
    %t-location scale distribution choice
    
    cutoff(i,1) = icdf(pd,percentage_cutoff);
    cutoff(i,2) = icdf(pd,1-(percentage_cutoff)); 
    
    tmp = cov_data_uniq_norm(:,:,i);
    tmp(tmp <= cutoff(i,1) | tmp >= cutoff(i,2)) = 1;
    tmp(tmp > cutoff(i,1) & tmp < cutoff(i,2)) = 0;
    
    
    cov_data_uniq_norm_bin(:,:,i) = tmp;
end;
clear prompt
Cij_bin_t = mean(cov_data_uniq_norm_bin,3);

%Now decompose Cij_bin_t to find ecogroup
[x,y] = eig(cov(Cij_bin_t)); x = fliplr(x); y = diag(y); y = flipud(y); 
y_frac = y/sum(y);
figure; histogram(y_frac,50); %Eigenspectrum of Cij_bin_t
title('Eigenspectrum of Temporally Conserved Taxon-Taxon Covariance Matrix');
xlabel('Eigenvalue');
ylabel('Number');


prompt = 'Which eigenvectors to use for ecogroup identification?';
ev_relevant = input(prompt);
    %For Raman et al, answer is 1

clear tmp
ev_wts = x(:,ev_relevant);

disp('Projection along eigenvectors of choice calculated.')







    
    
    


    

    
    
    


