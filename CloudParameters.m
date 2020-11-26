classdef CloudParameters

    properties
        gaussAmp
        becAmp
        offset
        pos
        gaussWidth
        becWidth
        cloudAngle
    end

    methods
        function self = CloudParameters(varargin)
            self.offset = 0;
            self.pos = 0;
            self.gaussAmp = 0;
            self.gaussWidth = 0;
            self.becAmp = 0;
            self.becWidth = 0;
            self.cloudAngle = 0;
            if isa(varargin{1},'CloudParameters') && isa(varargin{2},'CloudParameters')
                p1 = varargin{1};p2 = varargin{2};
                self.offset = [p1.offset,p2.offset];
                self.pos = [p1.pos,p2.pos];
                self.gaussAmp = [p1.gaussAmp,p2.gaussAmp];
                self.gaussWidth = [p1.gaussWidth,p2.gaussWidth];
                self.becAmp = [p1.becAmp,p2.becAmp];
                self.becWidth = [p1.becWidth,p2.becWidth];
            else
                self.offset = varargin{1};
                self.pos = varargin{2};
                self.gaussAmp = varargin{3};
                self.gaussWidth = varargin{4};
                if nargin > 4
                    self.becAmp = varargin{5};
                    self.becWidth = varargin{6};
                end
                if nargin > 6
                    self.cloudAngle = varargin{7};
                end
            end
        end
    end



end