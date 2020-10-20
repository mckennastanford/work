%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Israel Silber
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[Flist_struct] = give_me_files_and_subfolders(Search_String, Path)
%[Flist_struct] = give_me_files_and_subfolders(Search_String, Path, Filename_struct)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% file description: The function finds all files in the given Path and
% subfolder, whom name matches the 'Search_String' variable (for loading).
% In case that all the files are needed, make Search_String=''; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This improved version receives a cell type search string, or a simple 
% string meaning it can search for several strings in a single run, thus
% highly improving the function's performance, when working on several
% subfolders.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Flist_struct] = give_me_files_and_subfolders(Search_String, Path, varargin)

if ispc
    slash_backslash = '\';
else
    addpath(genpath('/chinook/meteo/ixs34/Matlab/IS_home_made_functions/'));
    slash_backslash = '/';
end

% if no filename struct was given as an input, creating a new one.
if nargin  == 2 
   Flist_struct = [];
else
   Flist_struct = varargin{1};
end

if ~iscell(Search_String)
    TMP = Search_String; %if the search stirng is nto a cell array, then turning it into one.
    clear Search_String %nead to clear the string in order to call a cell array with the same name.
    Search_String{1} = TMP; 
    clear TMP
end

for jj = 1 : length(Search_String)
    
    tmp_list = dir([Path,'*',Search_String{jj},'*']); %searching files with a matching string in the given path.    
    tmp_size = size(Flist_struct,1);
    
    if size(tmp_list,1) > 0        
        for ii = 1 : size(tmp_list)            
            Flist_struct(tmp_size+ii,1).name = tmp_list(ii).name; %adding matching files to the Flist struct.
            Flist_struct(tmp_size+ii,1).path = Path; %adding matching files to the Flist struct.     
            Flist_struct(tmp_size+ii,1).KB = tmp_list(ii).bytes ./ 1024; %adding matching file's size in KB (the dir function gives the file size in bytes).     
        end        
    end
    
end

%finding subfolders.
tmp_list = dir(Path);
tmp_list = tmp_list([tmp_list.isdir] == 1); %the rectangle brackets let you receive the values of a non scalar struct values as a double array.
if size(tmp_list,1) > 0
    
    for ii=1:size(tmp_list)
        
        if ~strcmp(tmp_list(ii).name,'..') && ~strcmp(tmp_list(ii).name,'.')
            New_Path = [Path,tmp_list(ii).name,slash_backslash] %no ';' to print on screen the current folder
            [Flist_struct] = give_me_files_and_subfolders(Search_String,New_Path,Flist_struct); %recursion for subfolders
            
        end
    end
    
end

return
