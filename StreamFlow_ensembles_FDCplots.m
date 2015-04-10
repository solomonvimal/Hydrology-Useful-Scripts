%% ######################           STREAMFLOW          ######################
% This script takes a matlab file [array] as an input and makes box plots, CDFs and FDC curves for the time series
% Take in the file names 
clc; 

% Work Space names

prefix = 'ws_for_final_plots_RCP';
suffix = [26 45 60 85];

% corrected data - outliers removed
load('C:\Users\Solo\Desktop\Pennar\Matlab\proj_streamflow_corrected_for_outliers.mat');
%% ##########  Ensemble averaging for streamflow - Across all models X all RCPS  ##############

ensemblevariable  = 'StreamFlow';

% List and extract filenames 
clear data n data_multilocations ensemble
filenames = who;
n = strmatch('proj_streamflow',filenames); % dont use strcmp!! Doesn't work.
names = filenames(n);

for ii = 1:length(n) % Looping over Models
    data_multilocations = eval([cell2mat([names(ii)])]); % send data to separate columns iteratively
    data_multilocations = data_multilocations(1:1140,:);
    data(:,ii) = data_multilocations;
end
data=data';
ensemble = [mean(data)]' ;

% Multi-model Stream Flow of 32 Combinations - 4 RCPs x 8 models = 32 columns  
eval([strcat('ensemble_', ensemblevariable) '= ensemble;'])
clear data n data_multilocations ensemble;

%% ########## ############# Separate Models - Across all RCPs ############# ######################

% List and extract

filenames = who;
n = strmatch('proj_streamflow',filenames); % dont use strcmp!!
names = filenames(n);

% Concatenate all the 32 combinations

for ii=1:length(n)
    data = eval([cell2mat([names(ii)])]);
    data_all_models(:,ii) = data(1:1140,:);
end

% Transpose to get 1140 x 32 matrix

data = data_all_models';
[n, ~] = size(data);


% Average to get RCP-wise model ensembles

j=1;
for ii = 1:4:n-3
    data_to_avg = data(ii:ii+3,:);
    data_avgd = mean(data_to_avg);
    data_matrix_avgd(j,:) = data_avgd;
    j=j+1;
end

clear data data_matri_summed data_summed data_to_sum n j ii i;

Stream_across_RCPs_Monthly = data_matrix_avgd';


%% ############################ Yearly Stream Flow - RCPs separated ##############
% 
% ensemblevariable  = 'StreamFlow';
% 
% % Load Observed Data 
% 
% load 'somasila_inflow_1989_6_2004_12.mat';
% inflow = somasila_inflow_1989_6_2004_12';
% inflow = inflow';
% inflow = inflow(8:end); % trucate data to start from January
% [~, y_inflow, ~] = wyd(inflow);
% 
% % Load Projected data: List files from workspace
% 
% filenames = who;
% n = strmatch('yearly_streamflow',filenames); % dont use strcmp!!
% names = filenames(n);
% 
% clear data data_multilocations ii;
% for ii = 1:length(n) % Looping over Models
%             data_multilocations = eval([cell2mat([names(ii)])]); % send data to separate columns iteratively
%             data_multilocations = data_multilocations';
%             data_multilocations = [data_multilocations(1:95,:)];
%             data(:,ii) = data_multilocations;
% end
% yearly_StreamFlow = data';
% clear data
% 
% x = yearly_StreamFlow;
% 
% clear data data_RCPs;
% for i=1:4
%    data = [x(i,:); x(i+4,:) ; x(i+8,:); x(i+12,:); x(i+16,:); x(i+20,:); x(i+24,:); x(i+28,:)];
%    dataT = data';
%    
%    % RCPs Separated
%    
%    eval([strcat('stream_annual_RCP_', num2str(suffix(i))), '= dataT;']) % RCPs separated 
%    data_RCPs(i,:) = mean(data);
% end
%    
% yearly_data_across_RCPs = data_RCPs';

%% ################    Monthly Stream FLow - RCPs separated   ########################### 

ensemblevariable  = 'StreamFlow';

% Load Projected data: List files from workspace

filenames = who;
n = strmatch('proj_streamflow',filenames); % dont use strcmp!!
names = filenames(n);

clear data data_multilocations;
for ii = 1:length(n) % Looping over Models
            data_multilocations = eval([cell2mat([names(ii)])]); % send data to separate columns iteratively
            data_multilocations = [data_multilocations(1:1140,:)];
            data(:,ii) = data_multilocations;
end
monthly_StreamFlow = data;
clear data

x = monthly_StreamFlow';

clear data data_RCPs i;
for i=1:4
   data = [x(i,:); x(i+4,:) ; x(i+8,:); x(i+12,:); x(i+16,:); x(i+20,:); x(i+24,:); x(i+28,:)];
   dataT = data';
   
   % RCPs Separated
   
   eval([strcat('stream_monthly_RCP_', num2str(suffix(i))), '= dataT;']) % RCPs separated 
   data_RCPs(i,:) = mean(data);
end
   
monthly_streamflow_across_RCPs = data_RCPs';

clear data data_RCPs dataT data_multilocations i

%% ################## Annual Stream Flow Series Plot ###################
% 
% clear n;
% filenames = who;
% n = strmatch('stream_annual_RCP_', filenames);
% names = filenames(n);
% 
% figure(1) % 4 subplots for each RCP -> Observed + 8 models
% for ii = 1:length(n)
%     data = eval([cell2mat([names(ii)])]);
% subplot(2,2,ii); plot([2006:2100],data)
% hold on 
% plot([1990:2004],y_inflow, 'color', [0 0 0], 'LineWidth', [2])
% xlim([1990 2100]); 
% title(strcat('RCP ', num2str(suffix(ii)))); 
% ylabel('Annual Stream Flow (Mcft)'); xlabel('Years');
% end
% 
% legend('BCC-CSM-1', 'GFDL-CM3', 'GFDL-ESM-2M', 'GISS-E2-H', 'GISS-E2-R', 'IPSL-CM5A-LR', 'MRI-CGCM3', 'NorESM1-M', 'Observed'); 
% 
% % Global Title
% set(gcf,'NextPlot','add');
% axes;
% h = title(strcat('Annual Stream Flow for different RCPs: Observed and Projected'));
% set(gca,'Visible','off');
% set(h,'Visible','on');
% filename = 'Stream_Series_4plots_8models';
% saveas(gcf,filename, 'tiffn'); 
% saveas(gcf,filename, 'fig');

%% ################       StreamFLow across RCPs        ####################
% % Transpose to get 32 x 95 matrix
% clear data data_matrix_avgd data_avgd data_to_avg n ii i;
% data = yearly_StreamFlow;
% [n, ~] = size(data);
% 
% % Sum to get RCP-wise model ensembles
% 
% j=1;
% for ii = 1:4:n-3
%     data_to_avg = data(ii:ii+3,:);
%     data_avgd = mean(data_to_avg);
%     data_matrix_avgd(j,:) = data_avgd;
%     j=j+1;
% end
% 
% Yearly_StreamFlow_8_models = data_matrix_avgd';
% clear data data_matrix_avgd data_avgd data_to_avg n
% 
% 
% % Load inflow data
% 
% load 'somasila_inflow_1989_6_2004_12.mat';
% inflow = somasila_inflow_1989_6_2004_12;
% 
% % Truncate to get from Jan to December
% inflow = inflow(8:end); 
% 
% % Get yearly inflow using wyd function
% [~, y_inflow, ~] = wyd(inflow);
% 
% % plot(data); hold on %% 10^8 values???
% 
% figure(2)
% plot([1990:2004],y_inflow,'color',[0 0 0],'LineWidth', 3)
% xlabel('Years')
% ylabel('Stream Flow (MCft)')
% title('Yearly Inflow and yearly ensemble Streamflow')
% 
% [~, y_ensemble_streamflow, ~] = wyd(ensemble_StreamFlow);
% 
% hold on
% plot([2006:2100],y_ensemble_streamflow);
% 
% filename = 'Stream_Yearly_Series_inflow_ensemble';
% saveas(gcf,filename, 'tiffn'); 
% saveas(gcf,filename, 'fig');

%% #################  To get time slabs like 2010-30, 40-60, 70-90: ###############

% Create month column

   u = [1:12]; u=u'; 
   month = repmat(u,95,1);
      
% Create year column

   Y = repmat(2006,12,1);
       for iii = 2007:2100
   V = repmat(iii,12,1);
   Y = vertcat (Y,V);
       end
   
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

%% #########################   BOX PLOTS  - STREMFLOW - not needed  ##############################################

ensemblevariable  = 'Stream_Flow_';
load 'somasila_inflow_1989_6_2004_12.mat';
inflow = somasila_inflow_1989_6_2004_12;


[~, n]=size(ensemble_StreamFlow);
for k=1:n
future_stream_temp = ensemble_StreamFlow(:,k); % change names of RHS
observed_stream_temp = inflow(:,k);    %(1:252,k); % ##### Is this ok? - truncate the observed data to first 20 years (similar to projected data.
% figure(3)
% plot(ensemble_StreamFlow, 'color', [1 0 0]); hold on
% plot(inflow)

x = future_stream_temp; % MONTHLY DATA FOR 2006-2100

x10_30 = x(Jan_2010_index:Dec_2030_index,:); % 2010-30, at all stations
x40_60 = x(Jan_2040_index:Dec_2060_index,:); % 2040-60
x70_90 = x(Jan_2070_index:Dec_2090_index,:); % 2070-90


x10_30 = x10_30(1:187,:);
x40_60 = x40_60(1:187,:);
x70_90 = x70_90(1:187,:);

figure(4)
% BOX PLOT NOT NEEDED 

subplot(1,4,1); p1 = boxplot(observed_stream_temp, 'labels', '1971-2005 (Observed)');
subplot(1,4,2); p1 = boxplot(x10_30, 'labels', '2010-2030');
subplot(1,4,3); p1 = boxplot(x40_60, 'labels', '2040-2060');
subplot(1,4,4); p1 = boxplot(x70_90, 'labels', '2070-2090');


% Global Title
set(gcf,'NextPlot','add');
axes;
h = title(strcat('Stream Flow for different Time Slices - Observed and Multi=model Average'));
set(gca,'Visible','off');
set(h,'Visible','on');
clear future_rain_temp observed_rain_temp;
filename = strcat('Boxplot of Monthly Stream Flow');
saveas(gcf,filename, 'tiffn'); 
end

% ##########################  Color Space - not used  #################################### 

% x1 = ensemble_StreamFlow;
% data_to_plot = x1;
% 
% [m nplots] = size(data_to_plot);
% cmap = hsv(nplots); % colors evenly spread across RGB space
% for i=1:nplots
%     x=data_to_plot(:,i);
% % figure(3)    
% % subplot(nplots,1,i); % dimension of subplot
% plotname = plot(x);
% set(plotname, 'color',cmap(i,:))
% hold on
% end
% hold off
% clear nplots i;

%%  ############################### CDF ###################################

% For Monthly Stream Flow
rcp_26 = stream_monthly_RCP_26;
rcp_45 = stream_monthly_RCP_45;
rcp_60 = stream_monthly_RCP_60;
rcp_85 = stream_monthly_RCP_85;

% Color Map
cmap = hsv(9);


for i=1:8
    figure(5)
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
for i=1:4
    hold on
    subplot(2,2,i); F = cdfplot(inflow); set(F, 'color', [0 0 0], 'LineWidth', [2])
%     ylim([0.7 1]); xlim([0 0.8*10^8])
    title(strcat('RCP', num2str(suffix(i))))
    xlabel('x = Stream Flow (Mcft)')
end

% Legend for Observed
modelname = cell2mat(modelnames(9)); s{9}=sprintf('%s%g',modelname); legend(s)


% Global Title
set(gcf,'NextPlot','add');
axes;
h = title(strcat('CDFs for different RCPs: Observed and Projected Monthly Stream Flow'));
set(gca,'Visible','off');
set(h,'Visible','on');


filename = 'CDF_Streamflow';
saveas(gcf,filename, 'tiffn'); 
saveas(gcf,filename, 'fig');

%% #################### Flow Duration Curves ################################
 clear i data X P s; 
% Import data

x = ensemble_StreamFlow; % MONTHLY DATA FOR 2006-2100

% Split Data

x10_30 = x(Jan_2010_index:Dec_2030_index,:); % 2010-30, at all stations
x40_60 = x(Jan_2040_index:Dec_2060_index,:); % 2040-60
x70_90 = x(Jan_2070_index:Dec_2090_index,:); % 2070-90

% Create filename cell
N = {'x10_30' 'x40_60' 'x70_90' 'inflow'};


for ii=1:4
    data = eval([N{ii}]); 
    
    % Make colum vector
    
    X = data(:);
    
    % Sorted in ascending to plot easily
    
    X = sort(X,'ascend'); 
    for j = 1:length(X(1,:))
        [~, order]= sort(X(:,j));
        m(order,j) = 1:length(X(:,j)); % rank is 'm'
    end
    
    
    % Get Exceedance Probability
    
    n =length(X);
    for i=1:length(X)
        P(i)=m(i)./(n+1);
        Pc_exceeded(i,:)=(1-P(i))*100;
    end
    
    clear i 
    
    % Create axis -> is this needed ? 
    
    %###################### Curve Fitting ###################### -> not used in the plots
    
    % Spline interpolation
    
    x = Pc_exceeded';
    y = X;
    %x1 = x(1); x2 = x(end);
    % xx= [x1:-0.1:x2,1];
    % yy = spline(x,y,xx);
    % figure(1)
    % hold on
    % plot(xx,yy, 'LineWidth', 3);
    
    figure(6)
    hold on
    % scatter(Pc_exceeded,X); ylim([0 2*10^8])
    plot(x,y, 'color', cmap(ii,:)); ylim([0 4*10^7])
    
    %########################### FDC Plot Annotation ##############################
    
    title('Flow Duration Curve'); %xlim([0 20])
    xlabel('Percentage of time equalled or exceeded');
    ylabel('Flow (MCft)');
    
    % Legends
    
    modelnames = {'Projected 2010-2030' 'Projected 2040-2060' 'Projected 2070-2090' 'Observed 1990-2004'};
    modelname = cell2mat(modelnames(ii));
    s{ii}=sprintf('%s%g',modelname);
    legend(s)
    clear Pc_exceeded X P n
end


filename = 'Flow_Duration_Curve';
saveas(gcf,filename, 'tiffn'); 
saveas(gcf,filename, 'fig');

% Reference for Flow Duration Curve - http://crk.iri.columbia.edu/water/course/view.php?id=14
