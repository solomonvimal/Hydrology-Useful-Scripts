
% This code creates specific file names, searches for them in a few workspaces, does some operations on the files, and then write an excel sheet for each combination.

%% Rainfall at 19 locations - .xls sheets - Time Series:

clc;

% Work Space names
prefix = 'ws_for_final_plots_RCP'
suffix = [26 45 60 85]
for k = 1:19
    
    for i=1:4 % Looping over RCPs
      filename = [prefix, num2str(suffix(i))];
      filenames = who('-file',filename);
      n = strmatch('sim_rain',filenames);  %%% Change search word here
   
        for ii = 1:8 % Looping over Models
        j=n(ii); 
        
      names = filenames(n);  
      RHS = cell2mat([names(ii)]);
      
      % data(:,ii) = eval([RHS])
      % Assign Convenient names to files in the workspace for different GCMs 
      RHS = eval([RHS]);
      RHS = RHS(1:1140,:) % reduce the dimension of the larger time series
      
      % Change the below line
      eval ([strcat('Rainfall_RCP_', num2str(suffix(i)), '_Model_' , num2str(ii)) '= RHS;']);
      
      % GCM data Time Series
      data(:,ii) = RHS(:,k);  %%% data for each location
      
      eval([strcat('RCP_', num2str(suffix(i)), '_Rain_data') '=' 'data'])
      
      title_row={'Year','Month','BCC-CSM-1','GFDL-CM3','GFDL-ESM-2M','GISS-E2-H','GISS-E2-R','IPSL-CM5A-LR','MRI-CGCM3','NorESM1-M'};
       
      % Create month column
      u = [1:12]; u=u'; 
      month = repmat(u,95,1);
      
      % Create year column
      Y = repmat(2006,12,1);
        for iii = 2007:2100
      V = repmat(iii,12,1);
      Y = vertcat (Y,V);
        end
   
      Y = Y(1:1140, :) % Filter out the excess years. Error! Check later!
      size(month);
      size(Y);
      year_month = horzcat(Y, month)
      
     % Create Excel Sheets and write data
        ExcelFileName = ['Rainfall Time Series_at_Location_', num2str(k)]
        sheetname = strcat('RCP_', num2str(suffix(i)))
        
        xlswrite(ExcelFileName, title_row, sheetname, 'A1')
        xlswrite (ExcelFileName, year_month, sheetname, 'A2')
        xlswrite(ExcelFileName, data, sheetname, 'C2')
        
        end
    end
    clear data RHS eval V u Y i ii iii ans j n ExcelFileName title ; % important
    
end
