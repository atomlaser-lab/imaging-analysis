classdef RawImageData < handle
    %RAWIMAGEDATA Defines a class for handling loading raw image data
    properties
        directory   %Directory to load data from
        files       %Files to load
        images      %Loaded image data
    end
    
    properties(SetAccess = protected)
        status     %Status information about raw image data
    end

    properties(Constant, Hidden=true)
        DEFAULT_SIZE = [1036,1384];
        DEFAULT_DIRECTORY = 'E:\RawImages\2019\01Jan\top';
        DEFAULT_NUM_IMAGES = 2;
        DEFAULT_BINARY_TYPE = 'mono16';
        FILE_EXT = '.raw';
    end


    methods
        function self = RawImageData(varargin)
            %RAWIMAGEDATA Creates a RAWIMAGEDATA object
            %
            %   RAW = RAWIMAGEDATA() creates a blank object RAW
            %
            %   RAW = RAWIMAGEDATA(NAME,VALUE,...) creates an object RAW
            %   and loads data based on the NAME/VALUE pairs given in the
            %   function call
            self.status = ImageAnalysisErrorHandler;
            if nargin == 0
                return
            elseif nargin > 0 && mod(numel(varargin),2) ~= 0
                error('Must specify input arguments as name/value pairs');
            else
                self.load(varargin{:});
            end
        end

        function self = copy(self,obj)
            %COPY Copies the properties of an input object to the current
            %object
            %
            %   RAW = COPY(OBJ) copies the properties of OBJ to the current
            %   object RAW
            %
            p = properties(self);
            for nn = 1:numel(p)
                self.(p{nn}) = obj.(p{nn});
            end
        end
        
        function self = load(self,varargin)
            %LOAD loads raw image data
            %
            %   RAW = RAW.LOAD(NAME,VALUE,...) loads data based on
            %   NAME/VALUE pairs.  Options for NAME are 'directory' (from
            %   which to load data), 'filenames' for the file names to
            %   load, 'length' for the number of images in each image set,
            %   'index' for the number of image sets ago to load, 'dims'
            %   for the dimensions of the image, and 'datatype' for the
            %   data type to use
            %
            %   'filenames' should be given as a cell array of file names.
            %   Alternatively, it can be an array of file structures in the
            %   format returned by the DIR function and stored as the
            %   normal property RAWIMAGEDATA.FILES
            %
            %   'datatype' can be 'mono8', 'mono16', or 'raw8'
            if mod(numel(varargin),2) == 0
                self.directory = self.DEFAULT_DIRECTORY;
                filenames = 'last';
                len = self.DEFAULT_NUM_IMAGES;
                idx = 1;
                dims = self.DEFAULT_SIZE;
                dataType = self.DEFAULT_BINARY_TYPE;
                for nn = 1:2:numel(varargin)
                    cmd = lower(varargin{nn});
                    switch cmd
                        case 'directory'
                            self.directory = varargin{nn+1};
                        case {'filenames','files'}
                            filenames = varargin{nn+1};
                        case {'length','len'}
                            len = varargin{nn+1};
                        case {'index','idx'}
                            idx = varargin{nn+1};
                        case 'dims'
                            dims = varargin{nn+1};
                        case 'datatype'
                            dataType = varargin{nn+1};
                        otherwise
                            error('Option ''%s'' unsupported',cmd);
                    end
                end
                
                if ~iscell(filenames) && strcmpi(filenames,'last')
                    [self.files,self.status] = RawImageData.getLastFilenames(self.directory,len,idx);
                elseif iscell(filenames)
                    self.files = RawImageData.getFileInfo(self.directory,filenames);
                elseif isstruct(filenames)
                    self.files = filenames;
                end
                self.readImages(dataType,dims);
            end
        end

        function self = readImages(self,dataType,dims)
            %READIMAGES reads images from file
            %
            %   RAW = RAW.READIMAGES() reads images based on values stored
            %   in RAW.FILES using a default data type of 'mono16' and
            %   assuming the image is the same size as DEFAULT_SIZE.
            %
            %   RAW = RAW.READIMAGES(DATATYPE) uses the data type given by
            %   DATATYPE
            %
            %   RAW = RAW.READIMAGES(__,DIMS) uses the image dimensions
            %   given by DIMS
            %
            if nargin < 2
                dataType = self.DEFAULT_BINARY_TYPE;
            end
            if nargin < 3
                dims = self.DEFAULT_SIZE;
            end
            %
            % Change binary data type and duplication level (dupl). 'raw8'
            % has 3 colour channels but they are all the same value for
            % each pixel
            %
            switch lower(dataType)
                case 'mono8'
                    binaryType = 'uint8';
                    dupl = 1;
                case 'mono16'
                    binaryType = 'uint16';
                    dupl = 1;
                case 'raw8'
                    binaryType = 'uint8';
                    dupl = 3;
                otherwise
                    error('Unsupported data type');
            end
            %
            % Pre-allocate arrays
            %
            self.images = zeros([dims,numel(self.files)]);
            for nn = 1:numel(self.files)
                fid = fopen(fullfile(self.directory,self.files(nn).name),'r');          %Open file
                tmp = fread(fid,dupl*prod(dims),binaryType);                            %Read file
                fclose(fid);                                                            %Close file
                self.images(:,:,nn) = double(reshape(tmp(1:dupl:end),flip(dims))).';    %Reshape the data to be the appropriate size
            end
        end

        function imgNum = getImageNumbers(self)
            %GETIMAGENUMBERS Gets the image numbers from the file names
            %
            %   NUMS = RAW.GETIMAGENUMBERS() returns the numbers associated
            %   with the image files in RAW
            %
            for nn = 1:numel(self.files)
                tmp = regexp(self.files(nn).name,'[0-9]+(?=\.raw)','match');
                imgNum(nn) = str2double(tmp{1});
            end
        end
        
        function s = struct(self)
            %STRUCT Converts object instance into structure
            %
            %   S = RAW.STRUCT() converts object RAW into structure S
            s.directory = self.directory;
            s.files = self.files;
            s.images = self.images;
        end
        
        function s = saveobj(self)
            %SAVEOBJ Saves a simpler struct representation of the class
            %
            %   S = RAW.SAVEOBJ() creates simpler structure S for saving
            %   object RAW
            s = self.struct;
            s.images = uint16(s.images);
        end
    end

    methods(Static)
        function b = loadobj(a)
            %LOADOBJ Converts saved structure into object instance
            %
            %   RAW = LOADOBJ(A) converts saved structure A representing an
            %   instance of RAWIMAGEDATA into object instance RAW
            b = RawImageData;
            b.directory = a.directory;
            b.files = a.files;
            b.images = double(a.images);
        end
        
        function [f,msg] = getLastFilenames(directory,len,idx)
            %GETLASTFILENAMES Gets the last filenames from the directory
            %
            %   F = GETLASTFILENAMES(DIRECTORY,LEN,IDX) looks at directory
            %   DIRECTORY and grabs LEN files starting IDX ago.  So if IDX
            %   = 1 then it grabs files from END - (LEN - 1) : END.
            %
            %   If IDX is not specified, it assumes that IDX = 1.
            %
            %   [F,MSG] = GETLASTFILENAMES(__) returns a warning/error
            %   message MSG if something goes wrong while loading a file
            
            %
            % Get all file names, sort by date
            %
            files = dir(fullfile(directory,['*' RawImageData.FILE_EXT])); %#ok<*PROP>
            [~,k] = sortrows(datevec([files.datenum]));
            files(:) = files(k');
            %
            % Get LEN file names starting at END - (IDX*LEN - 1)
            %
            if nargin <= 2
                f = files(end-(len-1):end);
            else
                f = files(end-(idx*len-1) + (0:(len-1)));
            end
            
            %
            % This checks that the files were taken together by comparing
            % the times at which they were taken
            %
            if any(abs(diff([f.datenum]))*3600*24 > 3)
                str = sprintf('Set of images may not have been taken together. Time difference of %.3f s',max(abs(diff([f.datenum]))));
                warning(str); %#ok<SPWRN>
                msg = ImageAnalysisErrorHandler(ImageAnalysisErrorHandler.STATUS_WARNING,str);
            else
                msg = ImageAnalysisErrorHandler(ImageAnalysisErrorHandler.STATUS_OK);
            end
        end

        function f = getFileInfo(directory,filenames)
            %GETFILEINFO Returns information about user-supplied files
            %
            %   F = GETFILEINFO(DIRECTORY,FILENAMES) returns the file
            %   information for FILENAMES in DIRECTORY
            %
            for nn = 1:numel(filenames)
                files(nn) = dir(fullfile(directory,filenames{nn})); %#ok<*AGROW>
            end
            [~,k] = sortrows(datevec([files.datenum]));
            f = files(k);
        end
        
        function raw = loadImageSets(varargin)
            %LOADIMAGESETS Loads (possibly) multiple images and creates
            %RAWIMAGEDATA objects for each image set
            %
            %   RAW = LOADIMAGESETS() Loads the last image set using default
            %   values in RAWIMAGEDATA
            %
            %   RAW = LOADIMAGESETS('last',IDX) loads the last image sets
            %   given by IDX, where IDX is a vector of positive integers
            %   specifying the image set relative to the last image set.
            %   So IDX = 1 is the last image set, IDX = 2 is the
            %   second-to-last image set, etc.
            %
            %   RAW = LOADIMAGESETS(NAME,VALUE,...) Loads image sets given
            %   by NAME/VALUE pairs.  Valid NAME parameters are
            %   'filenames', 'directory', 'index', 'length', 'dims', and
            %   'datatype'
            %
            %   'filenames' should be a cell array of file names
            %   corresponding to one image set.  If more than one image
            %   set, it should be a cell array of cell arrays.
            %
            %   'filenames' can also be a cell array of struct arrays (1
            %   cell element for each image set) or an array of structs
            %   that containing file information
            %
            %   'length' is the length of each image set.
            %
            %   'directory' is the directory from which to load images
            %
            %   'index' is the same as IDX above
            %
            %   'dims' is the image size, and 'datatype' is the data type
            %   of the images, either 'mono8', 'mono16', or 'raw8'

            if mod(nargin,2) ~= 0
                error('Arguments must occur as name/value pairs!');
            end
            %
            % Set default values
            %
            filenames = 'last';
            index = 1;
            len = RawImageData.DEFAULT_NUM_IMAGES;
            directory = RawImageData.DEFAULT_DIRECTORY;
            dims = RawImageData.DEFAULT_SIZE;
            dataType = RawImageData.DEFAULT_BINARY_TYPE;
            %
            % Parse inputs
            %
            for nn = 1:2:nargin
                v = varargin{nn+1};
                switch lower(varargin{nn})
                    case {'filenames','files'}
                        filenames = v;
                    case 'directory'
                        directory = v;
                    case {'index','idx'}
                        index = v;
                    case {'length','len'}
                        len = v;
                    case 'dims'
                        dims = v;
                    case 'datatype'
                        dataType = v;
                end
            end
            
            if ischar(filenames) && strcmpi(filenames,'last')
                %
                % Load the last images according to the index
                %
                numImages = numel(index);
                raw(numImages,1) = RawImageData;
                for mm = 1:numImages
                    raw(mm).load('filenames','last','directory',directory,'index',index(mm),'len',len,'dims',dims,'datatype',dataType);
                end
            elseif iscell(filenames)
                %
                % Load files as a cell array
                %
                numImages = numel(filenames);
                raw(numImages,1) = RawImageData;
                for mm = 1:numImages
                    raw(mm).load('filenames',filenames{mm},'directory',directory,'len',len,'dims',dims,'datatype',dataType);
                end
            elseif isstruct(filenames)
                %
                % If files are as a array of structures and the structures
                % have a name, then load as array.  First check to make
                % sure the "name" field exists
                %
                if ~isfield(filenames(1),'name')
                    error('File structures must have the ''name'' field');
                end
                numImages = floor(numel(filenames)/len);
                raw(numImages,1) = RawImageData;
                for mm = 1:numImages
                    fileidx = ((mm-1)*len+1):(mm*len);
                    raw(mm).load('filenames',filenames(fileidx),'directory',directory,'len',len,'dims',dims,'datatype',dataType);
                end
            else
                %
                % Throw an error since this is an unsupported mode
                %
                error('Unsupported filename input type!');
            end
        end
    end

end