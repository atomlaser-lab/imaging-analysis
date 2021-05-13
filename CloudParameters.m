classdef CloudParameters
    %CLOUDPARAMETERS Defines a simple value class with properties that
    %describe a cloud of cold atoms

    properties
        gaussAmp    %Amplitude(s) of Gaussian fits
        becAmp      %Amplitude(s) of Thomas-Fermi fits
        offset      %Offset(s) on data
        lin         %Linear background variations in X and Y
        pos         %X and Y position offsets
        gaussWidth  %Gaussian widths
        becWidth    %Thomas-Fermi widths
        cloudAngle  %Rotation of Gaussian cloud when using 2D fits
    end

    methods
        function self = CloudParameters(varargin)
            %CLOUDPARAMETERS Creates a CLOUDPARAMETERS object with
            %user-supplied properties
            %
            %   C = CLOUDPARAMETERS() creates CLOUDPARAMETERS object C with
            %   default parameters (all 0)
            %
            %   C = CLOUDPARAMETERS(C1,C2) creates a new object C from
            %   input CLOUDPARAMETERS objects C1 and C2 assuming that C1
            %   describes the X distribution and C2 the Y distribution
            %
            %   C = CLOUDPARAMETERS(NAME,VALUE,...) creates object C using
            %   NAME/VALUE pairs
            %
            
            %
            % Set default values
            %
            self.offset = 0;
            self.pos = [0,0];
            self.gaussAmp = 0;
            self.gaussWidth = [0,0];
            self.becAmp = 0;
            self.becWidth = [0,0];
            self.cloudAngle = 0;
            self.lin = [0,0];
            if nargin == 0
                return
            elseif numel(varargin) == 2 && isa(varargin{1},'CloudParameters') && isa(varargin{2},'CloudParameters')
                %
                % If there are two inputs and they are both CloudParameter
                % objects then create new object assuming that the first
                % object is for the X distribution and the second is for
                % the Y distribution
                %
                p1 = varargin{1};p2 = varargin{2};
                self.offset = [p1.offset(1),p2.offset(1)];
                self.pos = [p1.pos,p2.pos];
                self.gaussAmp = [p1.gaussAmp(1),p2.gaussAmp(1)];
                self.gaussWidth = [p1.gaussWidth(1),p2.gaussWidth(1)];
                self.becAmp = [p1.becAmp(1),p2.becAmp(1)];
                self.becWidth = [p1.becWidth(1),p2.becWidth(1)];
            elseif mod(numel(varargin),2) == 0
                %
                % If arguments are in name/value pairs then parse the
                % arguments
                %
                for nn = 1:2:numel(varargin)
                    v = varargin{nn+1};
                    switch lower(varargin{nn})
                        case 'offset'
                            self.offset = v;
                        case 'pos'
                            self.pos = v;
                        case 'gaussamp'
                            self.gaussAmp = v;
                        case 'gausswidth'
                            self.gaussWidth = v;
                        case 'becamp'
                            self.becAmp = v;
                        case 'becwidth'
                            self.becWidth = v;
                        case 'cloudangle'
                            self.cloudAngle = v;
                        case {'lin','linear'}
                            self.lin = v;
                        otherwise
                            error('Input option %s not supported!',varargin{nn});
                    end
                end
            else
                error('Inputs must either be no inputs, two CloudParameters objects, or name/value pairs');
            end
        end
    end



end