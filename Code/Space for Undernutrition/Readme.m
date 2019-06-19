%%Instructions for plotting Fig. 3A

%Variables included in .mat file
    %data_char
        %SAM: Data pertaining to SAM trial
            %Sample: Sample number (months)
            %PID: Patient ID
            %Age: Age of sample
            %Treatment Arm: Digitized treatment arm
            %Treatment Arm lookup: Mapping to treatment arm
            %Ecogroup columns: Column number of ecogroup taxa
            %Sample mapping: Mapping index to cohort characteristic
            %OTU_list_total: List of all OTUs
        
        %MAM: Data pertaining to MAM trial
            %SampleID: Sample ID
            %PID: Patient ID
            %Study Arm: Digitized version of MDCF trial treatment arm
                    %1: RUSF
                    %2: MDCF-1
                    %3: MDCF-2
                    %4: MDCF-3
            %Sample number: Week of sampling for given patient
            %Sample number lookup: Mapping between sample number and
                %clinical cohort
            
            %Ecogroup columns: Columns for ecogroup taxa
            %Full data: 531x945 matrix where rows are feecal samples and
                %columns are OTUs
            
            %OTU_list_full: List of OTUs
            %OTU_labels: Taxonomic labels for OTUs
            %Structures
                %Four included structures for RUSF, MDCF-1, MDCF-2, MDCF-3
                %treatment arms. Each structure contains 'pretreatment
                %data' (wk2 of sampling), 'posttreatment data' (wk9 of
                %sampling) where each row is a fecal sample and each column
                %is an ecogroup taxon in the same order as 'Ecogroup
                %columns'.  'mean_std_pretreatment' is a 2x15 array where
                %the 1st row is the mean fractional abundance of the
                %ecogroup taxa (in the same order as 'Ecogroup columns) and
                %the 2nd row is the standard deviation of the fractional
                %abundance of corresponding ecogroup taxa. 
        %ind: An index that digitizes the cohort characteristic for all 766
        %samples used to plot Fig. 3A. plot_space.m describes the digitized
        %barcode in the colormap portion of the script titled 'Digitization
        %of cohort'
                
    %ecogroup_data and full_data
        %These are both structures that have the fractional abundance for
        %all samples named accordingly within each structure for the
        %ecogroup (15 OTU) and the full_data (944 OTU for SAM, 945 OTU for
        %MAM). The 'ecogroup_data' structure contains a matrix termed 'full_matrix' which
        %is the concatenated of all patient cohort.
      
%Generating plot in Fig. 3A
[x,y] = plot_space(ecogroup_data.full_matrix,data_char.ind);
%The second figure of the above command is the plot shown in Fig. 3A where
%x and y pertain to the projections along the computed eigenvectors and the
%eigenspectrum respectively for the full_matrix within the ecogroup_data
%structure.
    
        
        
            
            
            
         
            
        