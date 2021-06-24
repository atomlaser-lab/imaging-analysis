classdef CloudParameters < handle
    %CLOUDPARAMETERS Defines a  class with properties that parameters used
    %to fit to a distribution of atoms

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
            %   C = CLOUDPARAMETERS(D) creates a CLOUDPARAMETERS
            %   object C with parameters all set to value D
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
            if nargin == 0
                return;
            elseif nargin == 1
                %
                % Set all parameters to input
                %
                if isempty(varargin{1})   
                    p = properties(self);
                    for nn = 1:numel(p)
                        self.(p{nn}) = [];
                    end
                else
                    v = varargin{1};
                    self.offset = v;
                    self.pos = [v,v];
                    self.gaussAmp = v;
                    self.gaussWidth = [v,v];
                    self.becAmp = v;
                    self.becWidth = [v,v];
                    self.cloudAngle = v;
                    self.lin = [v,v];
                end
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
                self.set(varargin{:});
            else
                error('Inputs must either be no inputs, the single numeric input, two CloudParameters objects, or name/value pairs');
            end
        end
        
        function self = set(self,varargin)
            %SET Sets properties based on name/value pairs
            %
            %   SELF = SELF.SET(NAME,VALUE,...) sets properties indicated
            %   by case-insensitive NAME to VALUE
            if mod(numel(varargin),2) ~= 0
                error('Arguments must appear as name/value pairs!');
            else
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
            end
        end
        
        function self = compare(self,cp)
            %COMPARE Compares current object with input object, and for
            %every non-empty property in input object sets current object
            %property to that value
            %
            %   SELF = SELF.COMPARE(CP) Uses input CLOUDPARAMETERS object
            %   CP for comparison.
            if ~isa(cp,'CloudParameters')
                error('Input must be of type CloudParameters');
            end
            p = properties(self);
            for nn = 1:numel(p)
                v2 = cp.(p{nn});
                if ~isempty(v2)
                    self.(p{nn}) = v2;
                end
            end
        end
        
        function v = convert2array(self,varargin)
            %CONVERT2ARRAY Converts object to an array with properties
            %ordered according to input arguments
            %
            %   V = SELF.CONVERT2ARRAY(VARARGIN) Orders object properties
            %   in output array V according to input variable argument list
            %   VARARGIN. So inputs 'offset','pos','gaussamp', returns an
            %   array where the first element is the offset, the second is
            %   the position (or second and third, if pos is 2 elements),
            %   and so on
            v = [];
            for nn = 1:numel(varargin)
                switch lower(varargin{nn})
                    case 'offset'
                        v = [v,self.offset(:)']; %#ok<*AGROW>
                    case 'pos'
                        v = [v,self.pos(:)'];
                    case 'posx'
                        v = [v,self.pos(1)];
                    case 'posy'
                        v = [v,self.pos(2)];
                    case 'gaussamp'
                        v = [v,self.gaussAmp(:)'];
                    case 'gausswidth'
                        v = [v,self.gaussWidth(:)'];
                    case 'gausswidthx'
                        v = [v,self.gaussWidth(1)];
                    case 'gausswidthy'
                        v = [v,self.gaussWidth(2)];
                    case 'becamp'
                        v = [v,self.becAmp(:)'];
                    case 'becwidth'
                        v = [v,self.becWidth(:)'];
                    case 'becwidthx'
                        v = [v,self.becWidth(1)];
                    case 'becwidthy'
                        v = [v,self.becWidth(2)];
                    case {'lin','linear'}
                        v = [v,self.lin];
                    case {'linx','linearx'}
                        v = [v,self.lin(1)];
                    case {'liny','lineary'}
                        v = [v,self.lin(2)];
                    case 'cloudangle'
                        v = [v,self.cloudAngle(:)'];
                end
            end
        end
        
        function self = setFromArray(self,order,a)
            %SETFROMARRAY Sets properties in current object from an array
            %according to an input ordering scheme
            %
            %   SELF = SELF.SETFROMARRAY(ORDER,ARRAY) uses cell array ORDER
            %   to determine which elements in ARRAY are assigned to which
            %   properties of the current CLOUDPARAMETERS object SELF
            for nn = 1:numel(order)
                switch lower(order{nn})
                    case 'offset'
                        self.offset = a(nn);
                    case 'pos'
                        self.pos = a(nn);
                    case 'posx'
                        self.pos(1) = a(nn);
                    case 'posy'
                        self.pos(2) = a(nn);
                    case 'gaussamp'
                        self.gaussAmp = a(nn);
                    case 'gausswidth'
                        self.gaussWidth = a(nn);
                    case 'gausswidthx'
                        self.gaussWidth(1) = a(nn);
                    case 'gausswidthy'
                        self.gaussWidth(2) = a(nn);
                    case 'becamp'
                        self.becAmp = a(nn);
                    case 'becwidth'
                        self.becWidth = a(nn);
                    case 'becwidthx'
                        self.becWidth(1) = a(nn);
                    case 'becwidthy'
                        self.becWidth(2) = a(nn);
                    case {'lin','linear'}
                        self.lin = a(nn);
                    case {'linx','linearx'}
                        self.lin(1) = a(nn);
                    case {'liny','lineary'}
                        self.lin(2) = a(nn);
                    case 'cloudangle'
                        self.cloudAngle = a(nn);
                end
            end
        end
        
        function s = struct(self)
            %STRUCT Returns a structure version of the class instance
            %
            %   S = CP.STRUCT() Returns a struct version of CLOUDPARAMETERS
            %   instance CP
            p = properties(self);
            for nn = 1:numel(p)
                s.(p{nn}) = self.(p{nn});
            end
        end
        
        function s = saveobj(self)
            %SAVEOBJ Returns a version of the object that is easier to save
            %
            %   S = CP.SAVEOBJ() returns a structure S that represents
            %   instance CP
            s = self.struct;
        end
    end
    
    methods(Static)
        function b = loadobj(a)
            %LOADOBJ Converts saved structure A into class instance B
            %
            %   CONST = LOADOBJ(A) Converts saved structure A representing
            %   CLOUDPARAMTERS object into proper instance CONST
            
            b = CloudParameters;
            p = fields(a);
            for nn = 1:numel(p)
                b.(p{nn}) = a.(p{nn});
            end
        end
        
        function cp = convertFromArray(order,a)
            %CONVERTFROMARRAY Converts an input array into a
            %CLOUDPARAMETERS object based on input labelling scheme
            %
            %   CP = CONVERTFROMARRAY(ORDER,ARRAY) uses cell array ORDER to
            %   determine which elements in ARRAY are assigned to which
            %   properties of new CLOUDPARAMETERS object CP
            cp = CloudParameters;
            for nn = 1:numel(order)
                switch lower(order{nn})
                    case 'offset'
                        cp.offset = a(nn);
                    case 'pos'
                        cp.pos = a(nn);
                    case 'posx'
                        cp.pos(1) = a(nn);
                    case 'posy'
                        cp.pos(2) = a(nn);
                    case 'gaussamp'
                        cp.gaussAmp = a(nn);
                    case 'gausswidth'
                        cp.gaussWidth = a(nn);
                    case 'gausswidthx'
                        cp.gaussWidth(1) = a(nn);
                    case 'gausswidthy'
                        cp.gaussWidth(2) = a(nn);
                    case 'becamp'
                        cp.becAmp = a(nn);
                    case 'becwidth'
                        cp.becWidth = a(nn);
                    case 'becwidthx'
                        cp.becWidth(1) = a(nn);
                    case 'becwidthy'
                        cp.becWidth(2) = a(nn);
                    case {'lin','linear'}
                        cp.lin = a(nn);
                    case {'linx','linearx'}
                        cp.lin(1) = a(nn);
                    case {'liny','lineary'}
                        cp.lin(2) = a(nn);
                    case 'cloudangle'
                        cp.cloudAngle = a(nn);
                end
            end
        end
    end
end