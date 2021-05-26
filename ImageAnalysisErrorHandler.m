classdef ImageAnalysisErrorHandler < handle
    %IMAGEANALYSISERRORHANDLER Defines a class that allows for warnings and
    %other errors to be passed to the user
    
    properties
        status      %The error status
        message     %A message describing the error
    end
    
    properties(Constant,Hidden=true)
        STATUS_OK = 'ok';           %Indicates that no error has occurred
        STATUS_WARNING = 'warning'; %Indicates that a warning has occurred
        STATUS_ERROR = 'error';     %Indicates that an error has occurred
    end
    
    methods
        function self = ImageAnalysisErrorHandler(status,message)
            %IMAGEANALYSISERRORHANDLER Creates an instance of the class
            %
            %   ERR = IMAGEANALYSISERRORHANDLER() Creates a blank instance
            %
            %   ERR = IMAGEANALYSISERRORHANDLER(S) Creates an instance with
            %   status S
            %
            %   ERR = IMAGEANALYSISERRORHANDLER(S,M) Creates an instance
            %   with status S and message M
            if nargin == 0
                self.status = self.STATUS_OK;
                self.message = '';
            elseif nargin == 1
                self.status = status;
                self.message = '';
            else
                self.status = status;
                self.message = message;
            end
        end
        
        function r = ok(self)
            %OK Returns true if no error or warning, false otherwise
            %
            %   R = OK() Returns true if no error or warning, false
            %   otherwise
            if strcmpi(self.status,self.STATUS_OK)
                r = true;
            else
                r = false;
            end
        end
        
        function r = warning(self)
            %OK Returns true if warning, false otherwise
            %
            %   R = WARNING() Returns true if warning, false
            %   otherwise
            if strcmpi(self.status,self.STATUS_WARNING)
                r = true;
            else
                r = false;
            end
        end
        
        function r = error(self)
            %OK Returns true if error, false otherwise
            %
            %   R = ERROR() Returns true if error, false
            %   otherwise
            if strcmpi(self.status,self.STATUS_ERROR)
                r = true;
            else
                r = false;
            end
        end
    end
end