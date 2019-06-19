%Conduct PCA on serially concatenated data
%This script will conduct PCA on serially concatenated matrices and only
%give the eigenvalue spectrum of each matrix 

function [v_t] = serial_eig(data,timerange1,timerange2);
data = flipdim(data,3);
for i=timerange1:timerange2;
    if i == 1;
        
        tmp = data(:,find(sum(data(:,:,i)) > 0), i)';
        tmp_cov = cov(tmp);
        [x,v] = eig(tmp_cov);
        x = fliplr(x);
        v = diag(v);
        v = flipud(v);
        v_t(:,i) = v;
        %disp(tmp)
        
    elseif i>1;
        tmp1 = data(:,find(sum(data(:,:,i)) > 0), i)';
        tmp = cat(1,tmp,tmp1);
        %disp(tmp);
        tmp_cov = cov(tmp);
        [x,v] = eig(tmp_cov);
        x = fliplr(x);
        v = diag(v);
        v = flipud(v);
        v_t(:,i) = v;
        
        tmp1 = [];
    end;
end;
        
    
    
