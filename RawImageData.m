classdef RawImageData < handle

    properties
        directory
        files
        images
    end

    properties(Constant, Hidden=true)
        SIZE = [1036,1384];
        DEFAULT_DIRECTORY = 'E:\RawImages\2019\01Jan\top';
        DEFAULT_NUM_IMAGES = 2;
        FILE_EXT = '.raw';
    end


    methods
        function self = RawImageData(varargin)
            if nargin == 0
                return
            elseif nargin > 0 && mod(numel(varargin),2) ~= 0
                error('Must specify input arguments as name/value pairs');
            else
                self.load(varargin{:});
            end
        end

        function self = copy(self,obj)
            p = properties(self);
            for nn = 1:numel(p)
                self.(p{nn}) = obj.(p{nn});
            end
        end
        
        function self = load(self,varargin)
            if mod(numel(varargin),2) == 0
                self.directory = RawImageData.DEFAULT_DIRECTORY;
                filenames = 'last';
                len = RawImageData.DEFAULT_NUM_IMAGES;
                idx = 1;
                for nn = 1:2:numel(varargin)
                    cmd = lower(varargin{nn});
                    switch cmd
                        case 'directory'
                            self.directory = varargin{nn+1};
                        case 'filenames'
                            filenames = varargin{nn+1};
                        case 'length'
                            len = varargin{nn+1};
                        case 'index'
                            idx = varargin{nn+1};
                        otherwise
                            error('Option unsupported');
                    end
                end
                
                if ~iscell(filenames) && strcmpi(filenames,'last')
                    self.files = RawImageData.getLastFilenames(self.directory,len,idx);
                elseif iscell(filenames)
                    self.files = RawImageData.getFileInfo(self.directory,varargin{2});
                end
                self.readImages;
            end
        end

        function readImages(self,camType)
            if nargin < 2
                binaryType = 'uint16';
            else
                switch camType
                    case 'new'
                        binaryType = 'uint8';
                    case 'old'
                        binaryType = 'uint16';
                    otherwise
                        error('Camera type not supported');
                end
            end

            self.images = zeros([RawImageData.SIZE,numel(self.files)]);
            for nn = 1:numel(self.files)
                fid = fopen(fullfile(self.directory,self.files(nn).name),'r');
                tmp = fread(fid,flip(RawImageData.SIZE),binaryType);
                fclose(fid);
                self.images(:,:,nn) = double(tmp.');
            end
        end

        function imgNum = getImageNumbers(self)
            for nn = 1:numel(self.files)
                tmp = regexp(self.files(nn).name,'[0-9]+(?=\.raw)','match');
                imgNum(nn) = str2double(tmp{1});
            end
        end
    end

    methods(Static)
        function f = getLastFilenames(directory,len,idx)
            files = dir(fullfile(directory,['*' RawImageData.FILE_EXT])); %#ok<*PROP>
            [~,k] = sortrows(datevec([files.datenum]));
            files(:) = files(k');
            if nargin <= 2
                f = files(end-(len-1):end);
            else
                f = files(end-(idx*len-1) + (0:(len-1)));
            end
            
            dates = zeros(numel(f),1);
            for nn = 1:numel(f)
                dates(nn,1) = datenum(f(nn).date);
            end
            if any(abs(diff(dates)) > 3)
                warning('Set of images may not have been taken together. Time difference of %.3f s',max(abs(diff(dates))));
            end
        end

        function f = getFileInfo(directory,filenames)
            for nn = 1:numel(filenames)
                files(nn) = dir(fullfile(directory,filenames{nn})); %#ok<*AGROW>
            end
            [~,k] = sortrows(datevec([files.datenum]));
            f = files(k);
        end
    end

end