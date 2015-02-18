clc
clear
%% Training Porgramme Solution by Solomon Vimal 
% Tutorials 2 - 5 - Problem Sheets


%% Problem 1: Stream Flow data

StreamImport = importdata('Somasila_monthly_inflows.txt.txt');
StreamData = StreamImport.data;
FlowValues = StreamData(:,3); % Pic the stream flow vector

% Time 
Time_Flow= StreamData(:,1:2); %Pic the time : Index --> 1-199
[MaxStreamFlow, MaxIndex] = max(FlowValues);

       
% Plot Stream Flow
figure(1)
plot(FlowValues)
text(MaxIndex,MaxStreamFlow,[' \leftarrow maximum monthly flow is ' num2str(MaxStreamFlow)]'','FontSize',10)
ylabel('stream flow');
xlabel('time');
title('Time Series of Stream Flow');

%% 2. Precipitation Data
[data,delimiter_out,headerlinesout] = importdata('Pennar_IMD_monthly_precipitation_at_three_locations.txt.txt');
precip = data.data;

% Loop to store precip data of 3 locations as separate vectors 
for i = 3:5
    eval(['precip', int2str(i-2), '= precip(:,i);'])
end
clear i;

% Pracipitation Data 
precip_data = precip(:,3:5);
N=length(precip1);

% Precip Time 6-1989 to 12-2005 ------------>> Index: 222:420
precip_time = precip(:,1:2);

% Mean Precipitaion

mean_precip = mean(precip_data);

% Calculate Covariance Matrix
cov_precip = cov(precip_data);
X=cov_precip;

% Standardize
X=X-mean(X(:));
X=X/std(X(:));
std_precip = X;
clear X;

%% Problem 3
% Correlation coefficient for precipitation at 3 locations with Stream Flow

r1 = corrcoef(precip1(222:420,1),FlowValues(1:199,1));
Corr_precip1_Flow = r1(2,1);


r2 = corrcoef(precip2(222:420,1),FlowValues(1:199,1));
Corr_precip2_Flow = r2(2,1);


r3 = corrcoef(precip3(222:420,1),FlowValues(1:199,1));
Corr_precip3_Flow = r3(2,1);

display(['Maximum correlation is with precipitatoin at location no. 3; Corresponding correlation coeff is:  ', num2str(Corr_precip3_Flow)])

%% Linear Regression 

x = precip3(222:420,1);
y = FlowValues(1:199,1);

p = polyfit(x,y,1); % linear regression that predicts y from x
m = p(1); % slope
b = p(2); % y intercept

% Call polyval to use p to predict y, calling the result yfit:
yfit = polyval(p,x);
% Using polyval saves you from typing the fit equation yourself, which in this case looks like:
% yfit =  p(1) * x + p(2);

% Compute the residual values as a vector signed numbers:
yresid = y - yfit;
% Square the residuals and total them obtain the residual sum of squares:
SSresid = sum(yresid.^2);
% Compute the total sum of squares of y by multiplying the variance of y by the number of observations minus 1:
SStotal = (length(y)-1) * var(y);
% Compute R2 using the formula
rsq = 1 - SSresid/SStotal;
rsq;

%% Problem 4
NSAT_NCEP_file = importdata('Pennar_NCEP_monthly_near_surface_air_temp.txt.txt');
MSLP_NCEP_file = importdata('Pennar_NCEP_monthly_mean_sea_level_pr.txt.txt');
HUS_NCEP_file = importdata('Pennar_NCEP_monthly_specific_humidity_at_850hPa.txt.txt');
UAS_NCEP_file = importdata('Pennar_NCEP_monthly_surface_u_wind.txt.txt');
VAS_NCEP_file = importdata('Pennar_NCEP_monthly_surface_v_wind.txt.txt');
GEOPOT_NCEP_file = importdata('Pennar_NCEP_monthly_geopotential_height_at_500hPa.txt.txt');

NSAT_NCEP = NSAT_NCEP_file.data; MSLP_NCEP = MSLP_NCEP_file.data;
HUS_NCEP = HUS_NCEP_file.data; UAS_NCEP = UAS_NCEP_file.data;
VAS_NCEP = VAS_NCEP_file.data; GEOPOT_NCEP = GEOPOT_NCEP_file.data;

% Common range with Preciitation over time at location A ->index 222:420

NSAT_NCEP_A = NSAT_NCEP(222:420,3);
MSLP_NCEP_A = MSLP_NCEP(222:420,3);
HUS_NCEP_A = HUS_NCEP(222:420,3);
UAS_NCEP_A = UAS_NCEP(222:420,3);
VAS_NCEP_A = VAS_NCEP(222:420,3);
GEOPOT_NCEP_A = GEOPOT_NCEP(222:420,3);

%% Correlation of Climate Variables at location A with Precipitation at Location 1
NSAT_cor = corr(NSAT_NCEP_A,x);
MSLP_cor = corr(MSLP_NCEP_A,x);
HUS_cor = corr(HUS_NCEP_A,x);
UAS_cor = corr(UAS_NCEP_A,x);
VAS_cor = corr(VAS_NCEP_A,x);
GEOPOT_cor = corr(GEOPOT_NCEP_A,x);

if abs(NSAT_cor)>0.2
    display('NSAT is significantly correlated with Precipitation')
end
if abs(MSLP_cor)>0.2
    display('MSLP is significantly correlated with Precipitation')
end
if abs(HUS_cor)>0.2
    display('HUS is significantly correlated with Precipitation')
end
if abs(UAS_cor)>0.2
    display('UAS is significantly correlated with Precipitation')
end
if abs(VAS_cor)>0.2
    display('VAS is significantly correlated with Precipitation')
end
if abs(GEOPOT_cor)>0.2
    display('GEOPOT is significantly correlated with Precipitation')
end

% Predictor Dataset
P = [NSAT_NCEP_A, MSLP_NCEP_A, HUS_NCEP_A, UAS_NCEP_A, VAS_NCEP_A, GEOPOT_NCEP_A];

% Standardize with respect to Climatology of 1971-1990
ClimaVar_70_90_raw = [NSAT_NCEP,MSLP_NCEP,HUS_NCEP,UAS_NCEP,VAS_NCEP,GEOPOT_NCEP] ;
ClimaVar_70_90 = ClimaVar_70_90_raw([1:240],[3 7 11 15 19 23]); % extracting from double_index over 1971 - 12,1990

% Mean and Std-dev over 1971-1990
Mean_Clima_70_90 = mean(ClimaVar_70_90);
StdDev_clima_70_90 = std(ClimaVar_70_90); 

P_1 = P-repmat(Mean_Clima_70_90,length(P) ,1);
P_Std = P_1./repmat(StdDev_clima_70_90,length(P),1); % Standardized matrix

% Normalization of Predictor Dataset
P_norm = normc(P_Std); % Column-wise normalization 

%% Tutorial 3: Data Analysis:

%% Principal Component Analysis:  (modified from PCA_example - UMD)

[N,K] = size(P);
Cov_P = cov(P_norm);
S = Cov_P;

[eigenvectors, D] = eig(S);

eigenvalues = diag(D);

% sort the eigenvalues in descending order

% so that the first EOF explains most of the variance


[eigenvalues_sorted,sort_order] = sort(eigenvalues,1,'descend'); % keep a tab on the sorting order 

eigenvectors_sorted = eigenvectors(sort_order,sort_order);

% eigenvectors = pattern vectors = empirical orthogonal functions (EOFs)

% eigenvalues = principal components (PCs) = amplitudes

eigenvectors_sorted;
eigenvalues_sorted;

% explain the total variance

total_variance = sum(eigenvalues);

% explain the variance of each EOF (Principal Component)

for i=1:min(N,K)

    cumulative_variance = sum(eigenvalues_sorted(1:i)/total_variance)*100;

    disp(['Principal Component ' num2str(i) ' explains ' num2str(eigenvalues_sorted(i)/total_variance*100) '% of the variance in the dataset.']);

end
%% Tutorial 4 : PCA, MLR, Fuzzy C-means 
sample_data = [2.0 1.0 3.0;4.0 2.0 3.0; 4.0 1.0 0.0; 2.0 3.0 3.0; 5.0 1.0 9.0] ;
std_samp_data = zscore(sample_data);
cov_samp_data = cov(std_samp_data);
[N_samp K_samp] = size(cov_samp_data);

[eig_vect_samp eig_values_samp] = eigs(cov_samp_data);

eig_val_samp = diag(eig_values_samp);

% sort the eigenvalues in descending order
% so that the first EOF explains most of the variance

[eigenvalues_sorted_samp,sort_order_samp] = sort(eig_val_samp,1,'descend'); % keep a tab on the sorting order 

eigenvectors_sorted_samp = eig_vect_samp(sort_order_samp,sort_order_samp);

% eigenvectors = pattern vectors = empirical orthogonal functions (EOFs)

% eigenvalues = principal components (PCs) = amplitudes

eigenvectors_sorted_samp;
eigenvalues_sorted_samp;

% explain the total variance

total_variance_samp = sum(eig_val_samp);

% explain the variance of each EOF (Principal Component)

for i=1:min(N_samp,K_samp)

    cumulative_variance_samp = sum(eigenvalues_sorted_samp(1:i)/total_variance_samp)*100;

    disp(['Principal Component ' num2str(i) ' explains ' num2str(eigenvalues_sorted_samp(i)/total_variance_samp*100) '% of the variance in the dataset.']);
end

% Matrix of loadings from highest to lowest eigen values
loadings = eigenvectors_sorted_samp;

% Transformed data matrix: 
Transformed_data_matrix = zscore(sample_data)*loadings; % using 3 steps -> standardize-> get eig_vectors -> multiply

% Using Prnincop function -> Loadings and PCs 
[coeff,score] = princomp(zscore(sample_data)); % Using direct function on standardized data

% Plot the variability of PCs using scatter and pareto functions:
figure(2)
scatter(score(:,1),score(:,2))
xlabel('PC 1');
ylabel('PC 2');
title('Scatter Plot about first two Principal Directions')
figure(3)
pareto(eig_val_samp, {'PC1'  'PC2' 'PC3' 'PC4' 'PC5' 'PC6' 'PC7'})
title('Percentage of Variability Explained by PCs') 

%% MLR from Precipitation at 3 locations with Streamflow
% Multiple Linear Regression:
precip_all3 = precip_data(222:420, :);
X = precip_all3;
Y = FlowValues(1:199,1);

% alpha = 0.05 by default;
stats=regstats(Y,X);
Beta = stats.beta;
rmse = sqrt(stats.mse);
rsquare =stats.rsquare;
% Residuals (errors):
residuals = stats.r;
standards_residuals = stats.standres; 
figure(4)
hist(standards_residuals)
std_res = std(standards_residuals);

% Significance of Co-efficients: 
significance_t_stat_coeff = stats.tstat.beta; %%% how can I show from t-stat that coeff is significant??

%%% RMSE = fprintf('%f', rmse) Weirdest thing !
% Multiple Linear Fit 
display('Y = Beta(1) + Beta(2)*X1 + Beta(3)*X2 + Beta(4)*X3')
X1 = X(:,1);
X2 = X(:,2);
X3 = X(:,3);
Yfit=Beta(1) + Beta(2)*X1 + Beta(3)*X2 + Beta(3)*X3;

figure(5)
plot(yfit)
hold on;
plot(y,'r')
title('Actual and Fitted Stream Flow')
legend('Fitted Stream Flow','Actual Stream Flow')

%% Fuzzy C-means Clustering
data=P_Std;
% m = 1.5 (constant) and C = 3,4,5,6 
for c=3:6
[center,U,obj_fcn] = fcm(data,c,1.5); %% Cluster Center: How can it be plotted?
figure(6)
subplot(2,2,c-2)
plot(data(:,1), data(:,2),'o');
hold on;
title(['Fuzzy Clustering; m = 1.5, C = ',num2str(c)])
maxU = max(U);
index1 = find(U(1,:) == maxU);
index2 = find(U(2, :) == maxU);
index3 = find(U(3, :) == maxU);
line(data(index1,1),data(index1, 2),'linestyle','none',...
     'marker','*','color','g');
line(data(index2,1),data(index2, 2),'linestyle','none',...
     'marker', '*','color','r');
line(data(index3,1),data(index3, 2),'linestyle','none',...
     'marker', '*','color','y');
end
clear m c;
%% Fuzzy Clustering with C = 3, and m = 0.5, 1.5, 2.0

C=3; 
for m = [0.5,1.5,2.0]
    for i=1:3
    list = [0.5 1.5 2.0];
            if m==list(i)
            x = i;
            end
    end
        if (m>1)
[center,U,obj_fcn] = fcm(data,3,m); %% Cluster Center: How can it be plotted?
figure(7)
subplot(3,1,x)
plot(data(:,1), data(:,2),'o');
hold on;
title(['Fuzzy Clustering; c = 3, m = ',num2str(m)])
maxU = max(U);
index1 = find(U(1,:) == maxU);
index2 = find(U(2, :) == maxU);
index3 = find(U(3, :) == maxU);
line(data(index1,1),data(index1, 2),'linestyle','none',...
     'marker','*','color','g');
line(data(index2,1),data(index2, 2),'linestyle','none',...
     'marker', '*','color','r');
line(data(index3,1),data(index3, 2),'linestyle','none',...
     'marker', '*','color','y');
        else
            figure(7)
            subplot(3,1,x)
            plot(data(:,1), data(:,2),'x');
            title(['C = 3, m = ',num2str(m),' : m value should be greater than 1'])
        end
end
%% Tutorial 5: Downscaling \97 Pennar Case Study:

% Plot PCs
[pc,score,latent,t_square] = princomp(P_norm);
variance_pcs = cumsum(latent)./sum(latent); % latent is the eigen values sorted in descending order
for i = 1:length(variance_pcs)
    if variance_pcs(i)>0.95 && variance_pcs(i-1)<=0.95;
        num_pcs = i;
    end
end

display(['First ', num2str(num_pcs),' Principal Components contribute to over 95% variance, so it is sufficient to retain them'])
figure(8)
biplot(pc(:,1:2),'Scores',score(:,1:2),'VarLabels',...
  {'X1' 'X2' 'X3' 'X4' 'X5' 'X6'})

% Transformed Data matrix

T = P_norm*pc(:,1:num_pcs);

%% Fuzzy Clustering on the Transformed Data matrix:

data = T;
c = 3;
m = 1.5;

[center,U,obj_fcn] = fcm(data,c,m); %% Cluster Center: How can it be plotted?
figure(9)
plot(data(:,1), data(:,2),'o');
hold on;
title(['Fuzzy Clustering; m =', num2str(m), ' ; C = ',num2str(c)])
maxU = max(U);
index1 = find(U(1,:) == maxU);
index2 = find(U(2, :) == maxU);
index3 = find(U(3, :) == maxU);
line(data(index1,1),data(index1, 2),'linestyle','none',...
     'marker','*','color','g');
line(data(index2,1),data(index2, 2),'linestyle','none',...
     'marker', '*','color','r');
line(data(index3,1),data(index3, 2),'linestyle','none',...
     'marker', '*','color','y');
clear m c data;

%% ANN Model
% No. of hidden layers = 10
IMD_precip = precip1(222:420,1);  % time frame from 6-1989 to 12-2005
% U - membership function matrix
% T - transformed data matrix
GCM_dataset = xlsread('GCM_predictors.xlsx');

% pc from NCEP data set (6 variables) normalized with climatology of 1971-1990 

GCM = GCM_dataset(:, [3 5 7 9 11]);

% size of pcs and GCM dataset: 
% size(GCM_dataset(:, [3 5 7 9 11]))
% size(pc)

input_dataset = GCM*pc(:,1:num_pcs); 
size(GCM_dataset);









