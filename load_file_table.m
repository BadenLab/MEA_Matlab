function [name_loaded_files,Epochs] = load_file_table (Mode,existing_table)
%Mode 0 means that the table has not been created before, Mode 1 means it
%does already exist
[Saved_Epoch,Epoch_path] = uigetfile('*.mat','Select stimulus file');
     cd (Epoch_path);

   
     %Test if a file containing a Epochs structure has been loaded
     
     Epochs = matfile(Saved_Epoch);
     test = exist("Epochs",'var');
     %Throw error if not
     if test == 0 
         error("No data found, check loaded file")
     end
     
     
     
if Mode == 1
    name_loaded_files = existing_table;
end

if Mode == 0
%create the template of the table 

name_loaded_files = cell(1,3);
%The first entries are Name, binsize and Kernels (0,1)
name_loaded_files{1,1} = "Name";
name_loaded_files{1,2} = "Binsize";
name_loaded_files{1,3} = "Kernels";
%Change Mode so that next part of the Code will be run and error will be
%thrown in case the Mode has the wrong value
Mode = 1;
end

if Mode == 1
    [row_table,column_table] = size(name_loaded_files);
    
    %Make sure to start the entry after the last entry in the table
    rr = row_table+1;
    
        for cc = 1:column_table
                       
                if cc == 1
                    
                    name_loaded_files{rr,cc} = strcat(Epoch_path,Saved_Epoch);
                    
                elseif cc == 2
                    
                    name_loaded_files{rr,cc} = num2str(Epochs.binsize);
                
                elseif cc == 3
                    %Check if Kernels have been calculated and return false
                    %or true boolean
                    t = (who('-file', Saved_Epoch, 'Kernels'));
                    test  = isempty(t);                 
                    name_loaded_files{rr,cc} = num2str(logical(test == 0));
                end
                    
        end
else
    error("wrong Mode selected")
end
            
end
    
    
    






