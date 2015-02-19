
   %% ###############  Ensemble averaging for Projected Rainfall ####################
    
ensemblevariable  = 'Rainfall_';
clear data ensemble_rain;
filenames = who;
n = strmatch('sim_rain_',filenames); % dont use strcmp!! & change search word here
names = filenames(n);
size(n);
[~, nlocations] = size(eval([cell2mat([names(1)])]));

for k = 1:nlocations
    
   
        for ii = 1:length(n) % Looping over Models
            data_multilocations = eval([cell2mat([names(ii)])]); % send data to separate columns iteratively
            data_multilocations = data_multilocations(1:1140,:);
            data(:,ii) = data_multilocations(:,k);
        end
            
    ensemble = [mean(data')]' ;
    eval([strcat('ensemble_', ensemblevariable,'_at_location_', num2str(k)) '= ensemble;' ])
    ensemble_rain(:,k) = ensemble;
    clear data;
    
end 

% Plot projected ensemble and observed time series
figure(1)
plot(ensemble_rain);
title('Monthly Multimodal Projected Rainfall: 2006-2100')
figure(2)
plot(observed_rain);
title('Monthly Observed Rainfall: 1971-2005')

% To get time slabs like 2010-30, 40-60, 70-90:

% Create month column
u = [1:12]; u=u';
month = repmat(u,95,1);

% Create year column
Y = repmat(2006,12,1);
for iii = 2007:2100
V = repmat(iii,12,1);
Y = vertcat (Y,V);
end
clear iii;
Y = Y(1:1140, :); % truncate the excess years. Error! Check later!
size(month);
size(Y);
year_month = horzcat(Y, month);
% Get the index of years.
for i = 1:length(Y)
    
    if (month(i)==1 && Y(i)==2010)
        Jan_2010_index = i;
    end
    if (month(i)==12 && Y(i)==2030)
        Dec_2030_index = i;
    end
    if (month(i)==1 && Y(i)==2040)
        Jan_2040_index = i;
    end
    if (month(i)==12 && Y(i)==2060)
        Dec_2060_index = i;
    end
    if (month(i)==1 && Y(i)==2070)
        Jan_2070_index = i;
    end
    if (month(i)==12 && Y(i)==2090)
        Dec_2090_index = i;
    end
end
clear i;

    
%% ########  Box plot for wet and dry spells, for 2010-30, 2040-60 & 2070-90 ##########

[w_ensemble y_ensemble d_ensemble] = wyd(ensemble_rain);
[w_observed y_observed d_observed] = wyd(observed_rain);

data_for_box = y_ensemble;   % ###### change 1 ######
data_observed = y_observed;  % ###### change 2 ######

% File name for saving
ensemblevariable  = 'Yearly Rainfall'; % ###### change 3 ######

[~, n]=size(data_for_box);
for k=1:n
future_temp = data_for_box(:,k); % change names of RHS
observed_temp = data_observed(:,k);    %(1:252,k); % ##### Is this ok? - truncate the observed data to first 20 years (similar to projected data.


x = future_temp; % MONTHLY DATA FOR 2006-2100

x10_30 =  x(5:25);
% x(Jan_2010_index:Dec_2030_index,:); % 2010-30, at all stations
x40_60 = x(35:55);
% x(Jan_2040_index:Dec_2060_index,:); % 2040-60
x70_90 = x(75:95); % check this once
% x(Jan_2070_index:Dec_2090_index,:); % 2070-90

figure(k+3)
% 
% boxplot([observed_rain_temp x10_30 x40_60 x70_90])
% xlabel(gca,'xticklabel',{'1971-2005 (Observed)' '2010-2030' '2040-2060' '2070-2090'})

subplot(1,4,1); p1 = boxplot(observed_temp, 'labels', '1971-2005 (Observed)');
subplot(1,4,2); p1 = boxplot(x10_30, 'labels', '2010-2030');
subplot(1,4,3); p1 = boxplot(x40_60, 'labels', '2040-2060');
subplot(1,4,4); p1 = boxplot(x70_90, 'labels', '2070-2090');

% Global Title
set(gcf,'NextPlot','add');
axes;
%            ###### change 4:last ###### below line:
h = title(strcat('Yearly Rainfall for different Time Slices at Location  ', num2str(k)));
set(gca,'Visible','off');
set(h,'Visible','on');
clear future_rain_temp observed_rain_temp;
filename = strcat('Box_plot_', ensemblevariable, '_at location_', num2str(k));
saveas(gcf,filename, 'tiffn');
saveas(gcf,filename, 'fig');
end


%% ################    Monthly Temperature - RCPs separated   ########################### 
% genralized for copying - changes specified

ensemblevariable  = 'Rainfall'; %    ####### Change 1 ###########

% Load Projected data: List files from workspace

filenames = who;
n = strmatch('sim_rain',filenames); % dont use strcmp!! ####### Change 2 ###########
names = filenames(n);

for_size = eval([cell2mat([names(1)])]);
[~, nlocations] = size(for_size);

for k = 1:nlocations
clear data data_multilocations;
for ii = 1:length(n) % Looping over Models
            data_multilocations = eval([cell2mat([names(ii)])]); % send data to separate columns iteratively
            data_multilocations = [data_multilocations(1:1140,:)];
            data(:,ii) = data_multilocations(:,k);
end
monthly_rain = data; %  ####### Change 3 ###########
clear data

x = monthly_rain'; %  ####### Change 4 ###########

clear data data_RCPs i;
for i=1:4
   data = [x(i,:); x(i+4,:) ; x(i+8,:); x(i+12,:); x(i+16,:); x(i+20,:); x(i+24,:); x(i+28,:)];
   dataT = data';
   
   % RCPs Separated
              %  ####### Change 5 ########### below here
   eval([strcat('rain_monthly_RCP_', num2str(suffix(i)),'_at_location_',num2str(k)) '= dataT;']) % RCPs separated 
   data_RCPs(i,:) = mean(data);
end
   data_RCPs = data_RCPs';
eval([strcat('monthly_rain_across_RCPs_at_location_', num2str(k)) '= data_RCPs']);

clear monthly_rain data data_RCPs dataT data_multilocations
end
%%  ############################### CDF ###################################

% For Monthly Tmean

% Get no. of locations:

filename = who;
n = strmatch('rain_monthly_RCP_26', filename) %  ####### change 1 #########
names = filename(n);
nlocations = length(n); 
L = nlocations;


 
for k=1:L
    
% List Files

rcp_26 = eval([strcat('rain_monthly_RCP_26_at_location_',num2str(k))]);
rcp_45 = eval([strcat('rain_monthly_RCP_45_at_location_',num2str(k))]);
rcp_60 = eval([strcat('rain_monthly_RCP_60_at_location_',num2str(k))]);
rcp_85 = eval([strcat('rain_monthly_RCP_85_at_location_',num2str(k))]); %stream_monthly_RCP_85;


% Color Map
cmap = hsv(9);


for i=1:8
    figure(k)
    hold on
    p2 = subplot(2,2,1);
    F26(i) = cdfplot(rcp_26(:,i)); %ylim([0.7 1]); xlim([0 0.8*10^8]);
    set(F26(i), 'color', cmap(i,:))
    title('CDF for RCP 2.6')
    hold on
    
    p3 = subplot(2,2,2);
    F45(i) = cdfplot(rcp_45(:,i)); %ylim([0.7 1]); xlim([0 0.8*10^8]);
    set(F45(i), 'color', cmap(i,:))
    title('CDF for RCP 4.5')
    hold on
    
    p4 = subplot(2,2,3);
    F60(i) = cdfplot(rcp_60(:,i)); %ylim([0.7 1]); xlim([0 0.8*10^8]);
    set(F60(i), 'color', cmap(i,:))
    title('CDF for RCP 6.0')
    hold on
    
    p5 = subplot(2,2,4);
    F85(i) = cdfplot(rcp_85(:,i)); %ylim([0.7 1]); xlim([0 0.8*10^8]);
    set(F85(i), 'color', cmap(i,:))
    title('CDF for RCP 8.5')
    hold on
    
    modelnames = {'BCC-CSM-1' 'GFDL-CM3' 'GFDL-ESM-2M' 'GISS-E2-H' 'GISS-E2-R' 'IPSL-CM5A-LR' 'MRI-CGCM3' 'NorESM1-M' 'Observed'};
    modelname = cell2mat(modelnames(i));
    
    s{i}=sprintf('%s%g',modelname);
    legend(s)
    
end

% Add Observed CDFs

observed_data = observed_rain(:,k) ; % ####### change 3 #########

for i=1:4
    hold on
    subplot(2,2,i); F = cdfplot(observed_data); set(F, 'color', [0 0 0], 'LineWidth', [2])
    %ylim([0.7 1]); xlim([0 0.8*10^8])
    title(strcat('RCP', num2str(suffix(i))));xlabel('x = Rainfall (mm)');
end

% Legend for Observed
modelname = cell2mat(modelnames(9)); s{9}=sprintf('%s%g',modelname); legend(s)


% Global Title
set(gcf,'NextPlot','add');
axes;                                            % ####### change 4; last ######### below here
h = title(strcat('CDFs for different RCPs: Observed and Projected Rainfall at location- ', num2str(k)));
set(gca,'Visible','off');
set(h,'Visible','on');

 %                    ####### change 5; last ######### below here
filename = strcat('CDF_for_Rainfall_at_location_', num2str(k));
saveas(gcf,filename, 'tiffn'); 
saveas(gcf,filename, 'fig');

end

clear observed_data data;
