function [] = concatenate_mukamel(N,input,output)
% To apply before applying the Mukamel PCA-ICA bits.  February 2019, Carey
% Lab - leonard.dupont@ens.fr
%
%  INPUT
%
%     N             the number of tiff files (trials) to be concatenated
%                   into one
%     input         the input path
%     output        the output path
%
%  OUTPUT
%
%     The function saves a tiff file from the generated concatenation matrix. 

%% About directories and existing tiff files

if (isempty(input))
    input = uigetdir();
end


%output = set_param_if_not_empty(output,input);

%finding number of .tif files in the input path which we can concatenate

filePattern = fullfile(input, '*Ready.tif');

if (isempty(filePattern))
    error('No registered tif files in this directory')
end

V = length(dir(filePattern));
videofiles = dir(filePattern); 
original_output = output;

%% Dealing with concatenation groupsize 
if (N > V)
    warning(['Number of registered tif files smaller than concatenation group size, using number of tif files (',num2str(V),') instead.']); 
    headcounts = V;
else
    aside = rem(V,N);
    if aside < 4
        headcounts = N*linspace(1,1,floor(V/N)); %we create a list with the number of concatenated videos per group
        headcounts(end) = headcounts(end) + aside; %last group gets the remainer of the division
    else
       headcounts = N*linspace(1,1,floor(V/N)); 
       headcounts(end+1) = aside;
    end
    
end

%% Concatenation bit + saving in uint16 tiff format 
c=1;
for konkat=1:length(headcounts) 
    clear filename
    clear output
    clear concatenated
    clear concatenation
    
    output = original_output;
    
    tic
    disp(['Starting to concatenate tif files in group ', num2str(konkat), ' out of ', num2str(length(headcounts)), '.'])
   
    thename = videofiles(c).name;
    underscore = find(thename == '_');
    
    if ispc
        filename = ['\',thename(1:underscore(8)),'_concatenated_trials'];
    else
        filename = ['/',thename(1:underscore(8)),'_concatenated_trials']; %building base filename
    end
    
   
    samples = headcounts(konkat);
    konkat_t = 0 ;
    
    for h=1:samples
        baseFileName = videofiles(c).name;
        thename = baseFileName;
        thefile = imread_tifflib(thename); 
        [x,y,t] = size(thefile);
        konkat_t = konkat_t + t;
        concatenation.(['v_',num2str(h)]).data = thefile;
        concatenation.(['v_',num2str(h)]).time = t;
        underscore = find(thename == '_');
        trial_nb = thename(underscore(end-2)+1:underscore(end-1)-1); %we get the reference number of the trial
        filename = strcat(filename,['_',num2str(trial_nb)]); %adding trial reference number to filename
        c=c+1;
    end
    
    concatenated = zeros(x,y,konkat_t); %the final matrix
    clear konkat_t
    konkat_t = 1;
 
    
    for h=1:samples
        t = concatenation.(['v_',num2str(h)]).time;
        concatenated(:,:, konkat_t:(konkat_t-1+t)) = concatenation.(['v_',num2str(h)]).data; 
        konkat_t = konkat_t + t;
    end
    
    concatenated = uint16(concatenated); %back to 16 bits format
    
    output = strcat(output,filename); %concatenate with the output path
    output = strcat(output,'.tif');
    save_tiff(concatenated,output); %limiting timestep, clearly
    
    elapsed = toc;
    disp(['Concatenation done. Elapsed time was ', num2str(elapsed) , ' s.'])
    
end



