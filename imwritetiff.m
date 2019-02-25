function [ varargout ] = imwritetiff( varargin )
%IMWRITETIFF Write tiff files on disk using Matlab's Tiff class
%   This function can be seen as an extension to imwrite, to be used when
%   writing multistack tiff images. This function uses, when possible, the
%   compiled mex file that allows for faster tiff writing. When used this
%   way, the tiff file is left open by default and the user is responsible to close 
%   the data stream once he doen not need to write any more pages 
%   to the tiff file, by calling imwritetiff('close').
%   Only if the parameter 'close' is set to true the file will be closed automatically.
%
%SYNTAX:
%   [ varargout ] = imwritetiff( data, filename, varargin )
%     data and filename are the only mandatory parameters. If no other
%     parameters are specified, the function uses the default ones.
%   [ varargout ] = imwritetiff( data, map, filename, varargin ) 
%     alternative where the user specifies a colormap in the second
%     argument.
%   imwritetiff('close') closes the file on which the user is currently
%     writing
%
%INPUT:
%   data the stack of images that will go into the tiff file. It is assumed
%     to be of size (rows, cols, numImages) for grayscale images and of size
%     (rows, cols, 3, numImages) for RGB images
%   filename the filename used to save data on disk
%
%OPTIONAL INPUT:
%   numImages total number of images that should be written to file
%   close (true/false): The file will only be closed if close  is set to true
%     By default close is false and the file is left open for further write
%     operations. In this case the user has to call util.imwritetiff('close') 
%     once he is done with writing to this file.
%   compression compression used when saved images. default is 'lzw',
%     other values allowed are 'none', 'lzw', 'deflate', 'packbits'
%   writeMode file opening mode. Default is 'create' (to create a new file).
%     Other accepted values are 'append', to append to a file previously
%     closed, or 'write' to add data to an already opened file. PLEASE NOTE
%     that calling 'append' on a file which was already opened will close
%     the file and then re-open it, losing the advantages of fast tiff
%     writing. If the file is already opoened the correct behaviour is to
%     use the 'write' mode.
%   isRGB explicitly specifies if the data should be saved as an RGB color
%     image. Default is false.
%   logging specifies whether to output some logging info such as time
%     information. It can be either 'off' (default) or 'on'
%   checkExistence (true/false) specifies if checks should be performed for the existance
%     of the file to a) warn if it is overwritten and b) create it if 
%     writtemode write is used on a nonexisting file. 
%     Can be very timeconsuming on large folers and can therefore be turned
%     off. Default true.
%   isBig boolean specifying whether the final files will be bigger than 4
%     Gb. Default false
%   resolution A two-element vector containing the XResolution and YResolution,
%     or a scalar indicating both resolutions; the default value is 72
%   resolutionUnit Specifies the unit used to specify the resolution
%     parameter. Can be 'inch', 'centimeter', 'cm', 'millimeter', 'mm', 'micrometer',
%     'um', or 'unknown' (the default)
%OUTPUT:
%   success integer value that equals zero if everything went smooth, -1
%   otherwise
%
%NOTE: IMWRITETIFF creates the object "Keep_me_I_survey_open_tiff_files" in the 
%   workspace of the calling function to close any open tiff file even if 
%   the calling function errors. Do not delete this variable.
%
%EXAMPLES:
%   data = uint8(ones(1024, 512, 50));
%   util.imwritetiff( data, 'test.tiff' ) writes all the 50 images of data on
%     test.tiff
%   util.imwritetiff( data, 'test.tiff', 'numImages', 20 ) writes the first 20
%     images of data on test.tiff
%   util.imwritetiff( data, 'test.tiff', 'numImages', 70 ) writes all the 50
%    images of data on test.tiff and issue a warning, because the
%    specified number of images is greater than the number of images in
%    data
%   util.imwritetiff( data, 'test.tiff', 'resolution', [48 48], 'resolutionUnit', 'millimeter' )
%     writes the imagess using a resolution of 48x48 dots per millimeter
%   util.imwritetiff( data, 'test.tiff', 'compression', 'none', 'writeMode', 'a')
%     appends data to test.tiff without using compression
%   util.imwritetiff('close')  closes the file previously opened and returns
%   For fast writing of multipage tiff on the same file:  
%     util.imwritetiff( data, 'test.tiff' );
%     util.imwritetiff( newdata, 'test.tiff', 'writemode', 'write' )
%     util.imwritetiff( otherdata, 'test.tiff', 'writemode', 'write' )
%     util.imwritetiff('close');
%
%AUTHORS: Stefano.Masneri@brain.mpg.de, Marcel.Lauterbach@brain.mpg.de


%% check if we just want to close the current image we are writing
if ischar( varargin{1} ) && strcmp('close', varargin{1})
  multitiff('close');
  if nargout > 0
  	varargout{1} = 0;
  end
  %delete the guard object, which is not needed any more
  evalin('caller','clear Keep_me_I_survey_open_tiff_files');
  return;
end

%% parse input
[data, map, filename, paramPairs] = parse_inputs(varargin{:});

p = inputParser;
p.addRequired('data', @(x) isnumeric(x) || islogical(x));
p.addRequired('filename', @ischar);
p.addParameter('numImages', 0, @(x) isinteger(uint32(x) ) );
p.addParameter('writeMode', 'create', @(x) any(strcmp(x, {'write', 'create', 'append', 'w', 'c', 'a'} ) ) );
p.addParameter('compression', 'lzw', @(x) any(strcmp(x, {'none', 'lzw', 'deflate', 'packbits'} ) ) );
p.addParameter('isRGB', false, @islogical);
p.addParameter('close', false, @islogical);
p.addParameter('checkExistence', true, @islogical);
p.addParameter('logging', 'off', @(x) any(strcmp(x, {'on', 'off'} ) ) );
p.addParameter('isBig', false, @islogical);
p.addParameter('resolution', 72, @(x) isinteger(uint32(x) ) );
p.addParameter('resolutionUnit', 'unknown', @(x) any(strcmp(x, {'inch', 'centimeter', ...
  'cm', 'millimeter', 'mm', 'micrometer', 'um', 'unknown'} ) ) );
p.parse( data, filename, paramPairs{:});

numImages = p.Results.numImages;
isBig = p.Results.isBig;
colMap = map;
resolution = uint32(p.Results.resolution);
if isscalar(resolution)
  resolution = [resolution resolution];
end

[~, ~, ext] = fileparts(filename);
if ~(strcmpi(ext, '.tif') || strcmpi(ext, '.tiff'))
  warning('imwritetiff: Using incorrect extension! Appending ''.tif'' to filename')
  filename = [filename '.tif'];
end

%% Specify resolution according to resolution unit
if strcmp(p.Results.resolutionUnit, 'unknown')
  resUnit = int16(1); %RESUNIT_INCH in Tiff standard
elseif strcmp(p.Results.resolutionUnit, 'inch')
  resUnit = int16(2); %RESUNIT_INCH in Tiff standard
else
  resUnit = int16(3); %%RESUNIT_CENTIMETER in Tiff standard
  if any(strcmp(p.Results.resolutionUnit, {'millimeter', 'mm'}));
    resolution = 10 * resolution;
  elseif any(strcmp(p.Results.resolutionUnit, {'micrometer', 'um'}));
    resolution = 10000 * resolution;
  end % do nothing is unit is centimeter
end

%% Check compression
switch p.Results.compression
  case 'none'
    compression = Tiff.Compression.None;
  case 'lzw'
    compression = Tiff.Compression.LZW;
  case 'deflate'
    compression = Tiff.Compression.Deflate;
  otherwise %paranoia
    compression = Tiff.Compression.PackBits;
end
compression = uint16(compression);

%% Check for map 
% Tiff library wants colormap to be a vector of length 256 (for each color
% channel). So we interpolate if the length is shorter than that
if ~isempty(colMap) && length(colMap) ~= 256
  len = size(colMap, 1);
  tempCM = zeros(256, 3);
  for k = 1:3
    tempCM(:,k) = interp1(1:len, double(colMap(:,k)), linspace(1, len, 256));
  end
  colMap = tempCM;
end
% Matlab always uses colormaps between 0 and 1, type double. TIFF
% specification instead wants uint16. So if the data is double and less
% than one we rescale. If it is another type we convert to uint16
if ~isempty(colMap)
  if isa(colMap, 'double') && max(colMap(:)) <= 1
    colMap = uint16(65536 * colMap);
  elseif ~isa(colMap, 'uint16')
    colMap = uint16(colMap);
  end
end

writeMode = p.Results.writeMode;
if strcmp( p.Results.writeMode, 'w')
  writeMode = 'write';
elseif strcmp( p.Results.writeMode, 'c')
  writeMode = 'create';
elseif strcmp( p.Results.writeMode, 'a')
  writeMode = 'append';
end

numDims = ndims(data);
dataSize = size(data);

%% number of images to write
if numDims == 2
  numDirs = 1;
elseif numDims == 3 && p.Results.isRGB && dataSize(3) == 3
	numDirs = 1;
else
  tiffImgs = size(data, numDims);
  if numImages > 0
    if tiffImgs < numImages
      numDirs = tiffImgs;
      warning('Trying to write more images than available!'); 
    else
      numDirs = numImages;
    end
  else
    numDirs = tiffImgs;
  end
end

% check datatype
switch class(data)
  case {'uint8', 'int8'}
    bitsPerSample = 8;
  case {'uint16', 'int16'}
    bitsPerSample = 16;
  case {'uint32', 'int32', 'single'}
    bitsPerSample = 32;
  case {'double'}
    bitsPerSample = 64;
  case {'logical'}
    bitsPerSample = 1;
  otherwise
    error('ERROR: Unsupported data type')
end

%% issue a warning if file is open in write mode and a file with same name
%already exists
fileExistsAlready = true; % a priori assume that the file exists (necessary for writemode write)
if p.Results.checkExistence
    checktimeStart = tic;
    fileExistsAlready = exist(filename,'file'); %if file exists (we use here filename and not fname because filename contains the path)
    checkTime = toc(checktimeStart);
    if checkTime > 0.2
        warning('imwritetiff:SlowExistenceCheck','imwritetiff:SlowExistenceCheck: Checking for exisence of file is very slow here (you have probably many files in this folder). Condider setting the option checkExistence to false.')
    end
    if  strcmp( p.Results.logging, 'on') && strcmp('create', writeMode) %
        if fileExistsAlready
            warning('File already exists, you are overwriting it!')
        end
    end
end
%% Check if we can use fast method
if ~isBig % if false, recheck!
  isBig = (bitsPerSample/8 * numDirs * dataSize(1) * dataSize(2)) > 4294967295;
end

if ~islogical(data) && ~isfloat(data) && ~isBig %doesn't work with logicals / floating points / > 4GB
  try
    if strcmp(writeMode, 'create') || strcmp(writeMode, 'append')
      checkExist = evalin('caller', 'exist(''Keep_me_I_survey_open_tiff_files'', ''var'')');
      if checkExist && strcmp(writeMode, 'append')
        warning('imwritetiff:automatic_file_closing', ...
                'Closing tiff file automatically by the imwritetiffGuard')
        disp('This is probably unexpected behaviour!')
        disp('If you are in ''append'' mode please note that:')
        disp('The writeMode ''append'' is intended to be used on closed files')
        disp('The file is then opened and will start writing after the end of previous data')
        disp('-----')
        disp('If the file is already open the same behaviour is achieved calling')
        disp('util.imwritetiff( ..., ''writeMode'', ''write'')')
        disp('Please check the help of imwritetiff for more details')
        disp('-----')
        disp('Now closing the file and returning to imwritetiff')
      end
      evalin('caller','Keep_me_I_survey_open_tiff_files = util.imwritetiffGuard;');
      multitiff(writeMode, filename); %initialize
    elseif strcmp(writeMode, 'write') && ~fileExistsAlready;
      % file does not exist, must create it!
      evalin('caller','Keep_me_I_survey_open_tiff_files = util.imwritetiffGuard;');
      multitiff('create', filename);
    end
    %single page grayscale
    if numDims == 2
      if isempty(colMap)
        multitiff('write', data', [], compression, resolution, resUnit);
      else
        multitiff('write', data', colMap, compression, resolution, resUnit);
      end
      if nargout > 0
          varargout{1} = 0;
      end

    elseif numDims == 3 && p.Results.isRGB && dataSize(3) == 3
      dataToWrite = permute(data,[3 2 1]); %different ordering compared to Matlab
      if ~isempty(colMap)
        warning('WARNING: Colormaps are used only with single channel data!')
      end
      multitiff('write', dataToWrite, [], compression, resolution, resUnit);
      if nargout > 0
        varargout{1} = 0;
      end

    %multipage
    else
      for k = 1:numDirs
        if numDims == 3
          if isempty(colMap)
            multitiff('write', squeeze(data(:,:,k))', [], compression, resolution, resUnit);
          else
            multitiff('write', squeeze(data(:,:,k))', colMap, compression, resolution, resUnit);
          end
        else
          dataToWrite = permute( squeeze(data(:,:,:,k)), [3 2 1] );
          if ~isempty(colMap)
            warning('WARNING: Colormaps are used only with single channel data!')
          end
          multitiff('write', dataToWrite, [], compression, resolution, resUnit);
        end
      end
    end  
    if p.Results.close %close file if wanted by user
        multitiff('close');
        evalin('caller','clear Keep_me_I_survey_open_tiff_files'); %delete the guard object, which is not needed any more
    end
    if nargout > 0
      varargout{1} = 0;
    end
  catch ME %catch any error during execution and try to close any remaining open files.
      if strcmpi(ME.identifier,'MATLAB:invalidMEXFile');
          warning('You probably do not have the Microsoft C++ runtime installed. Install it from https://www.microsoft.com/en-us/download/details.aspx?id=48145')
      end
      
      try %onother try for closing any potentially open files (if this fails, probably just no file was open)
        fileWasClosedMessage = 'Any open tiff file was closed now.';% a priory, assume that we will close any open file.
        multitiff('close');
        %after closing the file issue the error now
        error('imwritetiff:ErrorCaught', 'util.imwritetiff: %s Could not save file. %s',ME.message, fileWasClosedMessage);
      catch MEmultitiff
          if strcmpi(MEmultitiff.message ,'No file was open, first create one.') 
              fileWasClosedMessage = 'No file needed to be closed.';
          else
              warning('imwritetiff:CouldNotCloseOnError','imwritetiff:CouldNotCloseOnError: File could not be closed automatically')
          end
          
      end%of inner try (multitiff('close'))
    if nargout > 0
    end
    error('imwritetiff:ErrorCaught', 'util.imwritetiff: %s Could not save file. %s', ME.message, fileWasClosedMessage);
  end%of outer try
%% old way for logicals  and floating point data
else 
  %% we must setup the Tiff tags
  % RGB or greyscale?
  if numDims == 2
    samplesPerPixel = 1;
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
  elseif numDims == 3
    if dataSize(3) == 3 && p.Results.isRGB % it's a single RGB
      samplesPerPixel = 3;
      tagstruct.Photometric = Tiff.Photometric.RGB;
    else %it is a stack of 3 grayscale images
      samplesPerPixel = 1;
      tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    end
  elseif numDims == 4
    samplesPerPixel = 3;
    tagstruct.Photometric = Tiff.Photometric.RGB;
  else
    error('ERROR: data must be either 2D (single images), 3D (for greyscale stacks) or 4D (for RGB stacks)')
  end

  % setup tiff tag
  tagstruct.ImageLength = size(data,1);
  tagstruct.ImageWidth = size(data,2);
  tagstruct.BitsPerSample = bitsPerSample;
  tagstruct.SamplesPerPixel = samplesPerPixel;
  tagstruct.RowsPerStrip = size(data,2);
  tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
  tagstruct.Software = 'MATLAB';
  tagstruct.XResolution = double(resolution(1));
  tagstruct.YResolution = double(resolution(2));
  switch p.Results.compression
    case 'none'
      tagstruct.Compression = Tiff.Compression.None;
    case 'lzw'
      tagstruct.Compression = Tiff.Compression.LZW;
    case 'jpeg'
      tagstruct.Compression = Tiff.Compression.JPEG;
    case 'adobe'
      tagstruct.Compression = Tiff.Compression.AdobeDeflate;
    otherwise %paranoia
      tagstruct.Compression = Tiff.Compression.LZW;
  end
  % for binary images, compression is either none or PackBits
  if islogical(data) && ~strcmp('none', p.Results.compression)
    tagstruct.Compression = Tiff.Compression.PackBits;
    if strcmp( p.Results.logging, 'on')
      warning('Compression for binary images set to PackBits');
    end
  end

  switch class(data)
    case {'uint8', 'uint16', 'uint32', 'logical'}
      tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
    case {'int8', 'int16', 'int32'}
      tagstruct.SampleFormat = Tiff.SampleFormat.Int;
    case {'single', 'double'}
      tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
    otherwise
      error('ERROR: unsupported data type');
  end
  
  %from create, append, write to w / a
  if strcmp(writeMode, 'append')
    wrtMode = 'a';
  elseif isBig
    wrtMode = 'w8';
  else
    wrtMode = 'w';
  end
  
  try
    t = Tiff(filename, wrtMode);
    % Actually apply this tag to the first image
    t.setTag(tagstruct);
    % Write the first image
    if numDims == 2 || ( numDims == 3 && p.Results.isRGB && dataSize(3) == 3)
      t.write(data);
    elseif numDims == 3
      t.write(data(:,:,1));
    else %already checked that numDim is either 3 or 4
      t.write(data(:,:,:,1));
    end
    
    % For all further images
    for k = 2:numDirs
      % Create a new image inside the file
      t.writeDirectory();
      % Apply the tag for the new sub image
      t.setTag(tagstruct);
      % Write the image
      if numDims == 3
        write(t,data(:,:,k));
      else %already checked that numDim is either 3 or 4
        write(t,data(:,:,:,k));
      end
    end
    % Close the Tiff file
    t.close()
    if nargout > 0
      varargout{1} = 0;
    end
  catch ME
      if strcmpi(ME.identifier,'MATLAB:invalidMEXFile');
          warning('You probably do not have the Microsoft C++ runtime installed. Install it from https://www.microsoft.com/en-us/download/details.aspx?id=48145')
      end
      
    t.close()
    error('imwritetiff:ErrorCaught', 'util.imwritetiff: %s Could not save file. Any open tiff file was closed now.', ME.message);
  end
  
end
end


%%%
%%% Function parse_inputs
%%%
function [data, map, filename, paramPairs] = parse_inputs(varargin)

  data = [];
  map = [];
  filename = '';
  paramPairs = {};

  if (nargin < 2)
    error('Imwritetiff input parameter error: wrong number of inputs');
  end

  % find the first string -> the output filename
  firstString = [];
  for k = 1:length(varargin)
    if (ischar(varargin{k}))
      firstString = k;
      break;
    end
  end

  if (isempty(firstString))
    error('Imwritetiff input parameter error: missing filename');
  end

  switch firstString
  case 1
    error('Imwritetiff input parameter error: first arg must be image data');
  case 2
    % imwritetiff(data, filename, ...)
    data = varargin{1};
    filename = varargin{2};
  case 3
    % imwrite(data, map, filename, ...)
    data = varargin{1};
    map = varargin{2};
    filename = varargin{3};
    if (size(map,2) ~= 3)
      error('Imwritetiff input parameter error: invalid colormap');
    end
    validateattributes(map,{'numeric'},{'>=',0},'','COLORMAP');
  otherwise
    error('Imwritetiff input parameter error: bad filename argument position');
  end

  if (length(varargin) > firstString)
    % There are additional arguments after the filename.
    if (~ischar(varargin{firstString + 1}))
      error('Imwritetiff input parameter error: invalid arguments');
    end
    % Do some validity checking on param-value pairs
    if (rem(length(paramPairs), 2) ~= 0)
      error('Imwritetiff input parameter error: invalid syntax');
    end
    %copy all the remaining parameters in param_pairs
    paramPairs = varargin((firstString + 1):end);
  end
end

