%%This script is to color the SAM/MAM space (Fig 8B of Raman et al) using
%%different input data (full v ecogroup v RF v whatever else will be used)

function [x,y] = plot_space(data,ind)

%This function will create a 3-D plot of the centroids of the cohorts that
%are delineated as the input for 'data'. 

%Inputs: data, a nxm matrix where n is the number of fecal samples and m is
%the number of OTUs that are included

%ind, a n x 1 matrix that is a mapping between each fecal sample and the
%metadata of choice

%ind_names, a nx1 matrix indicating the metadata category that will be used
%to color cohort

%Note! We can only use so many colors, thus there will be a cycling of
%colors according to
%https://www.mathworks.com/help/matlab/ref/colorspec.html with the
%additional use of [0.5 0.5 0.5] or grey instead of white

%Sometimes you will need many dimensions of metadata...this will be
%reflected in the ind matrix. The ind matrix needs to be made before
%placing it in this function. 


[x,y] = eig(cov(data')); x = fliplr(x); y = diag(y); y = flipud(y);


%Digitization of cohort 
colorMap = zeros(8,3);
colorMap(1,:) = [0, 0.4470, 0.7410]; %Navy blue--> Healthy, Young
colorMap(2,:) = [0.8500, 0.3250, 0.0980]; %Dark orange--> Healthy, Middle age
colorMap(3,:) = [0.9290, 0.6940, 0.1250]; %Light orange --> Healthy, Old age
colorMap(4,:) = [0.4940, 0.1840, 0.5560]; %Purple--> SAM, untreated
colorMap(5,:) = [0.4660, 0.6740, 0.1880]; %Green --> SAM, at discharge
colorMap(6,:) = [0.3010, 0.7450, 0.9330]; %Sky blue--> SAM, 1 month post discharge
colorMap(7,:) = [0.6350, 0.0780, 0.1840]; %Maroon--> SAM, 6 months post discharge
colorMap(8,:) = [0, 0, 1]; %blue--> SAM, 12 months post discharge
colorMap(9,:) = [0, 0.5, 0]; %Forest green --> MAM, pretreatment
colorMap(10,:) = [1, 0, 0]; %Red --> MAM treated with RUSF
colorMap(11,:) = [0, 0.75, 0.75]; %Cyan --> MAM treated with MDCF1
colorMap(12,:) = [0.75, 0.75, 0]; %Dark yellow --> MAM treated with MDCF2
colorMap(13,:) = [0,0,0]; %Black --> MAM treated with MDCF3

%Color space
colorMap_total = zeros(length(data),3);
for i = 1:length(colorMap(:,1)); %Should be cycling through 13 colorMap features
    tmp_ind = find(ind(:,1) == i); %Isolate rows that correspond to metadata 'i'
    colorMap_total(tmp_ind,:) = repmat(colorMap(i,:),length(tmp_ind),1);
    tmp = x(find(ind == i),1:3);
    %Calculate centroid
    centroid(i,1:3) = mean(tmp);
    
end;

%For visualization of all fecal samples at all timepoints
figure;scatter3(x(:,1),x(:,2),x(:,3),80,colorMap_total(:,1:3),'filled');

%For visualization of centroids of populations described in the colormap
%designation

figure; scatter3(centroid(:,1),centroid(:,2),centroid(:,3),80,colorMap(:,1:3),'filled');

%Bargraph of PC1 projections
figure; bar(centroid(:,1));