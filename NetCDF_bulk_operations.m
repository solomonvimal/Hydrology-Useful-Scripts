%% NetCDF bulk operations - modified and extended from Arpita MOndal's original code

%{
This code can be run directly from the folders that contain the NetCDF files. It
pulls out all the .nc files and truncates them to the required study
region: by default for the Pennar region and saves the files as .mat files with the same filename.

- Solomon Vimal
%}

%% Loading .nc files automatically from the folder and subsetting..
clear all; clc;
list = dir('*.nc'); % lists all .nc files
[len_list, ~] = size(list); % get length of list

filenames = ({ }); % pre-alocate cell for speed
for i=1:len_list
    filenames{i} = list(i,1).name;  % Assign file names to cell
end
clear i

% filenamesv = filenames'; l1 = length(filenames); l2 = length(filenamesv);

%%

for j=1:length(filenames) % loop over all files
    %% open & inquire (filenames{j})
    ncid = netcdf.open(filenames{j}, 'NC_NOWRITE');
    nc_1 = netcdf(filenames{j},'nowrite');
    
    
    [numdim , ~, ~] = netcdf.inq(ncid);
    
    %% loop over all variables in each file to pick out lat-long variables
    
    % Explore the Contents
    
    [numdims,nvars,~] = netcdf.inq(ncid);
    
    for ii=1:nvars
        [name,~,~,~] = netcdf.inqVar(ncid,ii-1);
        blah=nc_1{name};
        
        % Set Lat-Lon Extents
        if (strcmp(name,'lat')==1)
            blahblah=blah(:,:,:,:);
            for i=1:length(blahblah)
                if blahblah(i)<12.5 && blahblah(i+1)>12.5
                    lat_start_index=i;
                end
                if blahblah(i)>17.5 && blahblah(i-1)<17.5
                    lat_end_index=i;
                end
            end
        end
        
        if (strcmp(name,'lon')==1)
            blahblah=blah(:,:,:,:);
            for i=1:length(blahblah)
                if blahblah(i)<75 && blahblah(i+1)>75
                    lon_start_index=i;
                end
                if blahblah(i)>82.5 && blahblah(i-1)<82.5
                    lon_end_index=i;
                end
            end
        end
        % Set plev at 85000
        if (strcmp(name,'plev')==1)
            blahblah=blah(:,:,:,:);
            
            for i=2:length(blahblah)
                if blahblah(i)==85000
                    plev_index=i;
                end
            end
        end
        
        
        
        
    end
    temp_var =(blah(:,lat_start_index:lat_end_index,lon_start_index:lon_end_index));
    eval(['savefile = ' '''' filenames{j} '.mat'''])
    savefile = savefile(1:end-7); % to remove all the .nc and .mat extensions in the filename to be saved
    % save syntax: save(name_tobe_savedas, 'variable')
    save(savefile , 'temp_var');
    
    clear blah blahblah name aaa bbb natts
    
end




%{

Below are the other variants of this code that I tried. It would be useful
to use the below code if numdims, numatts, etc. of the files are important
in the analysis.
&& (strcmp(name,'hus')==1)
           hus_temp=(blah(:,plev_index,lat_start_index:lat_end_index,lon_start_index:lon_end_index));
        elseif (strcmp(name,'tas')==1)
           tas_temp=(blah(:,lat_start_index:lat_end_index,lon_start_index:lon_end_index));
        elseif (strcmp(name,'psl')==1)
           psl_temp=(blah(:,lat_start_index:lat_end_index,lon_start_index:lon_end_index));
        elseif (strcmp(name,'uas')==1)
           uas_temp=(blah(:,lat_start_index:lat_end_index,lon_start_index:lon_end_index));
        elseif (strcmp(name,'vas')==1)
           vas_temp=(blah(:,lat_start_index:lat_end_index,lon_start_index:lon_end_index));
      filename
%}




%% This is a variant that I tried. But got confused and gave up after a few days :)
% Too much use of eval function is kind of confusing. This code can help us check the variables, dimensions and other details within the netcdf file.

%{
eval(['ncid' num2str(i) '= netcdf.open(' 'filenames{i}' ',''NC_NOWRITE'' );']);
    eval(['nc_' num2str(i) '= netcdf(' 'filenames{1}' ',''nowrite'' );']);
    a = num2str(i);
eval(['[numdim_' a, ' nvars_' a,' natts_' a ']' '= netcdf.inq(ncid' num2str(i) ');']);

 
    % ncid = ['ncid', num2str(j)]


Creating variables for each files separately

for ii = 1: eval(['nvars_' num2str(i)])
        eval(['[name', ' xtype', ' dimids' ' natts]' ' = netcdf.inqVar(ncid', a ' , ii-1)']);
        eval([ 'blah_' a ' = nc_' a '{name' '}'])
        
        % Set Lat-Lon extents
                
        eval(['blahblah' a '= blah_' a, '(:,:,:,:);']);
        
%% Extract Latitutde
        
        if eval(['(strcmp(name, ',     '''lat'')== 1)'])  % problem in introducing ' -> solved
            eval(['blahblah =' 'blah_' a '(:,:,:,:)'])
            for j = 2: eval(['length(blahblah' a ')'])  % loop over all latitudes  %% j = 2 for running the loop without error
                eval(['blahblah = blahblah' a])
                if blahblah(j) < 12.5 && blahblah(j-1) > 12.5
                    lat_start_index = j;
                else
                    fprintf('There is no variable! \n')
                end
                if blahblah(j) > 17.5 && blahblah(j-1) < 17.5
                    lat_end_index = j;
                else
                    fprintf('There is no variable! \n')
                end
            end
        end
        
        
   
            
%% Extract Longitude
            if eval(['(strcmp(name, ',     '''lon'')== 1)']) % problem in introducing '
                eval(['blahblah =' 'blah_' a '(:,:,:,:)'])
                
                for j = 2:length(blahblah) % loop over all longitudes
                    if blahblah(j) < 75 && blahblah(j-1) > 75
                        lon_start_index = j;
                    else
                        fprintf('There is no variable! \n')
                    end
                    if blahblah(j) > 82.5 && blahblah(j-1)< 82.5
                        lon_end_index = j;
                    else
                        fprintf('there is no variable! \n')
                    end
                    
                end
            end
            
 %% Assign Lat-Lon indices as limits and save file:
 
 
 
            
    end


%}



%% Other functions of use to us in concatenation:

% can also use cellfun

% genvarname(filenames{i}) - > cannot be used because of . (dots)

% Eval function in nameing with loop



% vertcat

% s=whos; s={s.name}; g=s(strmatch('matrix',s)); z=cellfun(@(x)
% evalin('base',x),g,'uni',0); vertcat(z{:})
%
% clear s g z clear; clc for i=1:length(filenames)
%      eval(['filename_' num2str(i) '= i'])
% end


% Different approach:

% assign variables names of each file to the same common variable names,
% and then after hyperslabbing, save them..

% Specify model name first -> incorporate model names into the programme..
% model_name = 'bcc_csm'
