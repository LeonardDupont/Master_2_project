classdef imwritetiffGuard < handle
%classdef imwritetiffGuard < handle
%a class to survey that files that are opened with util.imwritetiff are also closed again.
%This class is called by util.imwritetiff. You should not use it manually for not interfering with util.imwritetiffGuard
%When the function that called util.imwritetiff finishes (with or without error) an eventually open tiff file will be closed automatically by the destructor of the util.imwritetiffGuard object
%
%Properties: keepopen: If set to true, no tiff file will be closed automatically (eg for further write operations by another function). (default false)
%            verbose:  if set to true more output, mainly for debugging. (defautl false)
%
%Methods: None except for constructor and destructor
%
%usage example: 
%  util.imwritetiff(data,'dummyfile.tiff','numImages',10,'writemode','create');
%  util.imwritetiff(data,'dummyfile.tiff','numImages',10,'writemode','write');
%  %some calcualtion that might lead to errors. guard will than take care of closing dummyfile.tiff
%  util.imwritetiff('close') %regular closing of the dummyfile.tiff 
%
%AUTHOR: Marcel.Lauterbach@brain.mpg.de

  properties
    keepopen = false;
    verbose = false;
    filename;
  end
    
  methods
    function obj = imwritetiffGuard(varargin)
      %constructor.
      %Save the filename. At the moment this info is not at all used,
      %but it might become necessary if we introduce later handling of two open tiff files in parallel.
      %To avoid later changes of the interface, we implement this already here.
      if 2 == nargin 
        filename = varargin{1};
        obj.keepopen = varargin{2};
      elseif 1 == nargin
        filename = varargin{1};
      else
        filename = '';
      end
      %for uniquely identifying the file make sure that it contains the path:
      if ~isempty(filename) %we really have a filename
          [fpath, fname, fext] = fileparts(filename);
          if isempty(fpath) %filename was given without path
              fpath = pwd;%take current directory as path
          end
          filename = fullfile(fpath, fname, fext);
      end
      obj.filename = filename;
      %If the user gives a relative path, this might not be handled correctly. 
      %Add this funcitonality if needed. 
      %one solution might be: 
      % oldpath=pwd;
      % cd('mydir')
      % fullpath=pwd;
      % cd(oldpath)
      %To handle also pathologic cases this seems a serious problem. 
      %There is a GetFullPath by Jan Simon on Matlab File Exchange, 
      %but it seems overkill for the moment, especially since the 
      %filename is not used at all for now. 
    end

    function delete(obj) %destructor. Will take care of closing any open tiff file from imwritetiff
      if obj.verbose
        disp('Destroying imwritetiffGuard object')
      end

      if ~obj.keepopen %if the user did not ask explicitely to keep the tiff file open
        try
          util.imwritetiff('close');
        catch ME
          if strcmpi(ME.message, 'No file was open, first create one.')
            %most likely cause for util.imwritetiff('close') to fail: 
            %everything ok, the file is already closed (as it
            %should)
          else
            warning('imwritetiffGuard:ErrorAtDelete', ...
                     'Warning: The following error was caught while executing util.imwritetiffGuard class destructor: %s', ...
                     ME.message);
          end%of if strcpmi..
        end%of try
      end%of if ~obj.keepopen
    end%of destructor

  end%of methods
    
end%of classdef
