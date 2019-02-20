function FinalImage = imread_tifflib(filename, frametoread)
%function FinalImage = imread_tifflib(filename, frametoread)
%read tiff image, using tifflib (faster).
%
%input:
%filename: the file that is read
%frametoread: frame of the file that is read(optional)
%if frametoread is not specified all frames will be read
%output: 
%the tiff file as an array

[~, ~, ext ] = fileparts(filename);

%check if file is a tiff file
if isempty(strfind(ext, 'tif'))
    error('imread_tifflib:Tiffonly','imread_tifflib can only read tiff imges. This filename has not a tiff extension and is presumably not a tiff file. Aborting.')
end


warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning')
%get the tiff info
InfoImage=imfinfo(filename);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;

switch InfoImage(1).ColorType
    case {'grayscale','indexed'}%imreads greyscaled files
        if ~exist('frametoread','var')%read all frames if not specified otherwise
            
            NumberImages=length(InfoImage);
            toread = 1:NumberImages;
            
            if InfoImage(1).BitDepth == 8%check if 8 bit or 16 bit
                FinalImage=zeros(nImage,mImage,NumberImages,'uint8');
            else
                FinalImage=zeros(nImage,mImage,NumberImages,'uint16');
            end
        else
            if InfoImage(1).BitDepth == 8%check if 8 bit or 16 bit
                FinalImage=zeros(nImage,mImage, 1,'uint8');
            else
                FinalImage=zeros(nImage,mImage, 1,'uint16');
            end
            toread = frametoread;
        end
    case 'truecolor'%imreads rgb tiffs
        if ~exist('frametoread','var')%read all frames if not specified otherwise
            
            NumberImages=length(InfoImage);
            toread = 1:NumberImages;
            
            if InfoImage(1).BitDepth == 8%check if 8 bit or 16 bit
                FinalImage=zeros(nImage,mImage, 4,NumberImages,'uint8');
            else
                FinalImage=zeros(nImage,mImage, 4,NumberImages,'uint16');
            end
        else
            if InfoImage(1).BitDepth == 8%check if 8 bit or 16 bit
                FinalImage=zeros(nImage,mImage, 4, 1,'uint8');
            else
                FinalImage=zeros(nImage,mImage, 4, 1,'uint16');
            end
            toread = frametoread;
        end
    case 'CMYK' % imreads cmyk files (insnt working yet)
        if ~exist('frametoread','var')%read all frames if not specified otherwise
            
            NumberImages=length(InfoImage);
            %toread = 1:NumberImages;
            toread = 1;
            if InfoImage(1).BitDepth == 8%check if 8 bit or 16 bit
                FinalImage=zeros(nImage,mImage, 4,NumberImages,'uint8');
            else
                FinalImage=zeros(nImage,mImage, 4,NumberImages,'uint16');
            end
        else
            if InfoImage(1).BitDepth == 8%check if 8 bit or 16 bit
                FinalImage=zeros(nImage,mImage, 4, 1,'uint8');
            else
                FinalImage=zeros(nImage,mImage, 4, 1,'uint16');
            end
            toread = frametoread;
        end
    otherwise
        warning('File is not a RGB Binary Image');
end

TifLink = Tiff(filename, 'r');
%read the image and output it
for i = toread
    TifLink.setDirectory(i);
    switch InfoImage(1).ColorType
        case 'truecolor'
            
            FinalImage(:,:,:,i)=TifLink.read();
        case {'grayscale','indexed'}
            FinalImage(:,:,i)=TifLink.read();
        case 'CMYK'
             FinalImage(:,:,:,:,i)=TifLink.read();
    end
end
%close tiff file
TifLink.close();
