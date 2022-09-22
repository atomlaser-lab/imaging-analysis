classdef BinaryImageData < RawImageData
    %BINARYIMAGEDATA Defines a class for handling loading raw image data

    methods
        function self = BinaryImageData(varargin)
            %BINARYIMAGEDATA Creates a BINARYIMAGEDATA object
            %
            %   RAW = BINARYIMAGEDATA() creates a blank object RAW
            %
            %   RAW = BINARYIMAGEDATA(NAME,VALUE,...) creates an object RAW
            %   and loads data based on the NAME/VALUE pairs given in the
            %   function call
            self@RawImageData(varargin{:});
        end

        function self = copy(self,obj)
            %COPY Copies the properties of an input object to the current
            %object
            %
            %   RAW = COPY(OBJ) copies the properties of OBJ to the current
            %   object RAW
            %
            p = {'directory','files','images'};
            for nn = 1:numel(p)
                self.(p{nn}) = obj.(p{nn});
            end
            self.status.copy(obj.status);
        end

        function r = is_multi_camera(self)
            %IS_MULTI_CAMERA Returns true if image files contains images
            %from multiple cameras

            r = iscell(self.images);
        end
        
        function self = load(self,varargin)
            %LOAD loads raw image data
            %
            %   RAW = RAW.LOAD(NAME,VALUE,...) loads data based on
            %   NAME/VALUE pairs.  Options for NAME are 'directory' (from
            %   which to load data), 'filenames' for the file names to
            %   load and'index' for the number of image sets ago to load
            %
            %   'filenames' should be given as a numeric or cell array
            %   of image numbers to read. Alternatively, it can be an array
            %   of file structures in the format returned by the DIR
            %   function and stored as the normal property
            %   RAWIMAGEDATA.FILES
            %

            if mod(numel(varargin),2) == 0
                self.directory = self.DEFAULT_DIRECTORY;
                filenames = 'last';
                idx = 1;
                rotation = 0;
                for nn = 1:2:numel(varargin)
                    cmd = lower(varargin{nn});
                    switch cmd
                        case 'directory'
                            self.directory = varargin{nn+1};
                        case {'filenames','files'}
                            filenames = varargin{nn+1};
                        case {'index','idx'}
                            idx = varargin{nn+1};
                        case 'rotation'
                            rotation = varargin{nn+1};
                        otherwise
                            error('Option ''%s'' unsupported',cmd);
                    end
                end
                
                if ~iscell(filenames) && strcmpi(filenames,'last')
                    [self.files,self.status] = BinaryImageData.getLastFilenames(self.directory,idx);
                elseif iscell(filenames) || isstring(filenames) || isnumeric(filenames)
                    self.files = BinaryImageData.getFileInfo(self.directory,filenames);
                elseif isstruct(filenames)
                    self.files = filenames;
                end
                self.readImages(rotation);
            end
        end

        function self = readImages(self,rotation)
            %READIMAGES reads images from file
            %
            %   RAW = RAW.READIMAGES() reads images based on files stored
            %   in RAW.FILES
            %
            %   RAW = RAW.READIMAGES(__,ROTATION) rotates the image
            %   dimensions by ROTATION, which is either 0, -90, 90, or 180
            %
            if nargin < 2
                rotation = 0;
            end
            binaryType = 'uint16';
            %
            % Read and parse images
            %
            fid = fopen(fullfile(self.directory,self.files.name),'r');          %Open file
            tmp = uint16(fread(fid,self.files.bytes/2,binaryType));             %Read file
            fclose(fid);                                                        %Close file
            self.images = BinaryImageData.parseBinaryImageData(tmp);            %Parse data
            if ~self.is_multi_camera
                num_sets = 1;
                self.images = {self.images};
            else
                num_sets = numel(self.images);
            end
            for mm = 1:num_sets
                %
                % Apply rotations
                %
                if rotation == 90
                    self.images{mm} = pagetranspose(self.images{mm});
                elseif rotation == 180
                    self.images{mm} = flipud(self.images{mm});
                elseif rotation == -90
                    self.images{mm} = flipud(pagetranspose(self.images{mm}));
                elseif rotation ~= 0
                    error('Only rotations supported at 0, 90, 180, and -90');
                end
                %
                % Check for saturation of the camera
                %
                switch binaryType
                    case 'uint8'
                        max_value = 240;
    
                    case 'uint16'
                        max_value = 65000;
                end
                for nn = 1:size(self.images{mm},3)
                    tmp = self.images{mm}(:,:,nn);
                    num_saturated_pixels = sum(tmp(:) >= max_value);
                    if num_saturated_pixels > 0.1*numel(tmp)
                        self.status.status = ImageAnalysisErrorHandler.STATUS_WARNING;
                        self.status.message = 'Image is saturated';
                        warning(self.status.message);
                        break
                    end
                end
            end
            %
            % Fix image structure for single camera image sets
            %
            if num_sets == 1
                self.images = self.images{1};
            end
        end

        function imgNum = getImageNumbers(self)
            %GETIMAGENUMBERS Gets the image numbers from the file names
            %
            %   NUMS = RAW.GETIMAGENUMBERS() returns the numbers associated
            %   with the image files in RAW
            %
            for nn = 1:numel(self.files)
                tmp = regexpi(self.files(nn).name,'[0-9]+(?=\.bin)','match');
                imgNum(nn) = str2double(tmp{1});
            end
        end
        
        function s = struct(self)
            %STRUCT Converts object instance into structure
            %
            %   S = RAW.STRUCT() converts object RAW into structure S
            if numel(self) > 1
                for nn = 1:numel(self)
                    s(nn,1) = self(nn).struct;
                end
            else
                s.directory = self.directory;
                s.files = self.files;
                s.images = uint16(self.images);
            end
        end
        
        function s = saveobj(self)
            %SAVEOBJ Saves a simpler struct representation of the class
            %
            %   S = RAW.SAVEOBJ() creates simpler structure S for saving
            %   object RAW
            s = self.struct;
        end
    end

    methods(Static)
        function raw = loadobj(a)
            %LOADOBJ Converts saved structure into object instance
            %
            %   RAW = LOADOBJ(A) converts saved structure A representing an
            %   instance of RAWIMAGEDATA into object instance RAW
            raw = RawImageData.empty;
            for nn = 1:numel(a)
                raw(nn,1).directory = a(nn).directory;
                raw(nn,1).files = a(nn).files;
                raw(nn,1).images = double(a(nn).images);
            end
        end
        
        function [files,msg] = getLastFilenames(directory,idx)
            %GETLASTFILENAMES Gets the last filenames from the directory
            %
            %   F = GETLASTFILENAMES(DIRECTORY,IDX) looks at directory
            %   DIRECTORY and grabs all files starting IDX ago.
            %
            %   If IDX is not specified, it assumes that IDX = 1.
            %
            %   [F,MSG] = GETLASTFILENAMES(__) returns a warning/error
            %   message MSG if something goes wrong while loading a file
            
            %
            % Get the last image file written
            %
            fid = fopen(fullfile(directory,'last-image.txt'),'r');
            last_image = str2double(fgetl(fid));
            fclose(fid);
            image_numbers = last_image - idx + 1;
            mm = 1;
            for nn = image_numbers
                files(mm,1) = dir(fullfile(directory,sprintf('bec%d.bin',nn)));
                mm = mm + 1;
            end
            msg = ImageAnalysisErrorHandler;
        end

        function f = getFileInfo(directory,filenames)
            %GETFILEINFO Returns information about user-supplied files
            %
            %   F = GETFILEINFO(DIRECTORY,FILENAMES) returns the file
            %   information for FILENAMES in DIRECTORY
            %
            for nn = 1:numel(filenames)
                if iscell(filenames)
                    fname = filenames{nn};
                    if isnumeric(fname)
                        fname = sprintf('bec%d.bin',fname(nn));
                    end
                elseif isstring(filenames)
                    fname = filenames(nn);
                elseif isnumeric(filenames)
                    fname = sprintf('bec%d.bin',filenames(nn));
                end
                files(nn,1) = dir(fullfile(directory,fname)); %#ok<*AGROW>
            end
            [~,k] = sortrows(datevec([files.datenum]));
            f = files(k);
        end
        
        function raw = loadImageSets(varargin)
            %LOADIMAGESETS Loads (possibly) multiple images and creates
            %BINARYIMAGEDATA objects for each image set
            %
            %   RAW = LOADIMAGESETS() Loads the last image set using default
            %   values in BINARYIMAGEDATA
            %
            %   RAW = LOADIMAGESETS('last',IDX) loads the last image sets
            %   given by IDX, where IDX is a vector of positive integers
            %   specifying the image set relative to the last image set.
            %   So IDX = 1 is the last image set, IDX = 2 is the
            %   second-to-last image set, etc.
            %
            %   RAW = LOADIMAGESETS(NAME,VALUE,...) Loads image sets given
            %   by NAME/VALUE pairs.  Valid NAME parameters are
            %   'filenames', 'directory', 'index', and 'rotation'
            %
            %   'filenames' should be a cell array of file names
            %   corresponding to one image set.  If more than one image
            %   set (N sets, for example), it should be a cell array of
            %   length N of cell arrays of length len
            %
            %   'filenames' can also be a cell array of struct arrays (1
            %   cell element for each image set) or an array of structs
            %   that containing file information
            %
            %   'filenames' can also be an array of structs which MUST have
            %   a "name" field. For N image sets, it should be N*len long
            %
            %   'filenames' can also be an array of strings. For N image
            %   sets, each set of length len, the array should be N*len
            %   long
            %
            %   'filenames' can also just be a numeric array of image
            %   numbers
            %
            %   'directory' is the directory from which to load images
            %
            %   'index' is the same as IDX above
            %
            %   'rotation' is the angle to rotate the image, either 0, 90,
            %   -90, or 180

            if mod(nargin,2) ~= 0
                error('Arguments must occur as name/value pairs!');
            end
            %
            % Set default values
            %
            filenames = 'last';
            index = 1;
            directory = BinaryImageData.DEFAULT_DIRECTORY;
            rotation = 0;
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
                    case 'rotation'
                        rotation = v;
                end
            end
            
            if (ischar(filenames) || (isstring(filenames) && numel(filenames) == 1)) && strcmpi(filenames,'last')
                %
                % Load the last images according to the index
                %
                numImages = numel(index);
                raw(numImages,1) = BinaryImageData;
                for mm = 1:numImages
                    raw(mm).load('filenames','last','directory',directory,'index',index(mm),'rotation',rotation);
                end
            elseif iscell(filenames)
                %
                % Load files as a cell array. For N image sets, this should
                % be a cell array of length N with each cell being a cell
                % array containing the file names of each raw image file
                % (so it should be len long). If one image set, it can be
                % just the one cell array of length len
                %
                numImages = numel(filenames);
                raw(numImages,1) = BinaryImageData;
                for mm = 1:numImages
                    raw(mm).load('filenames',filenames{mm},'directory',directory,'rotation',rotation);
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
                numImages = numel(filenames);
                raw(numImages,1) = BinaryImageData;
                for mm = 1:numImages
                    raw(mm).load('filenames',filenames(mm),'directory',directory,'rotation',rotation);
                end
            elseif isstring(filenames)
                %
                % If files are as an array of strings, then assume that it
                % is arranged as
                % [set1_file1,set1_file2,set2_file1,set2_file2,...]
                %
                numImages = floor(numel(filenames)/len);
                raw(numImages,1) = BinaryImageData;
                for mm = 1:numImages
                    fileidx = ((mm-1)*len+1):(mm*len);
                    raw(mm).load('filenames',filenames(fileidx),'directory',directory,'rotation',rotation);
                end
            elseif isnumeric(filenames)
                %
                % If filenames are just numeric values, assume they
                % correspond to image numbers
                %
                numImages = numel(filenames);
                raw(numImages,1) = BinaryImageData;
                for mm = 1:numImages
                    raw(mm).load('filenames',filenames(mm),'directory',directory,'rotation',rotation);
                end
            else
                %
                % Throw an error since this is an unsupported mode
                %
                error('Unsupported filename input type!');
            end
        end

        function imgs = parseBinaryImageData(data)
            %PARSEBINARYIMAGEDATA Parses a 1D array of binary data into a
            %set of 2D images.  The first 2-bytes of the image data must be
            %the image format version number, and this function should be modified
            %to handle updates to the image format.
            %
            %   IMGS = PARSEBINARYIMAGEDATA(DATA) parses 1D array DATA into
            %   3D array IMGS
            %
            img_version = data(1);
            switch img_version
                case 0
                    img_width = double(data(2));
                    img_height = double(data(3));
                    num_images = data(4);
                    data = data(5:end);
                    imgs = zeros(img_width,img_height,num_images);
                    start_idx = 1;
                    final_idx = start_idx + img_height*img_width - 1;
                    for nn = 1:num_images
                        imgs(:,:,nn) = reshape(double(data(start_idx:final_idx)),img_width,img_height);
                        start_idx = final_idx + 1;
                        final_idx = start_idx + img_height*img_width - 1;
                    end

                case 1
                    img_width(1) = double(data(2));
                    img_height(1) = double(data(3));
                    num_images(1) = data(4);

                    img_width(2) = double(data(5));
                    img_height(2) = double(data(6));
                    num_images(2) = data(7);

                    start_idx = 8;
                    final_idx = start_idx + img_height(1)*img_width(1) - 1;
                    for mm = 1:numel(img_width)
                        if mm > 1
                            final_idx = start_idx + img_height(mm)*img_width(mm) - 1;
                        end
                        imgs{mm} = zeros(img_width(mm),img_height(mm),num_images(mm));
                        for nn = 1:num_images(mm)
                            imgs{mm}(:,:,nn) = reshape(double(data(start_idx:final_idx)),img_width(mm),img_height(mm));
                            start_idx = final_idx + 1;
                            final_idx = start_idx + img_height(mm)*img_width(mm) - 1;
                        end
                    end

                otherwise
                    error('Unknown image version %d!',img_version);
            end
        end
    end

end