function [medfilt_data] = med_filter_calcium(input)
%input should be the full path to the file, including the file name ending
%with .tif. The output is a saved tiff file after applying a median filter
%using a 5-frame window. MEDFILT is added to the file name at the end. 


data = imread_tifflib(input); 
data = double(data); 
medfilt_data = medfilt1(data,5,[],3); 

C = strsplit(input,'.'); 
output = [C{1},'_MEDFILT.',C{2}]; 

medfilt_data = uint16(medfilt_data); 
save_tiff(medfilt_data,output)