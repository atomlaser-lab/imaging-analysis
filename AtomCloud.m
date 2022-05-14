classdef AtomCloud < handle

    properties
        %
        % Atomic sample properties
        %
        N               %Number of atoms extracted using a fit method in FITDATA
        Nsum            %Number of atoms extracted by summing of the region of interest.
        pos             %Position of the cloud as [x0,y0]
        gaussWidth      %Width of the Gaussian fits as [xwidth,ywidth]
        T               %Temperature of the cloud as [Tx,Ty]
        peakOD          %Peak optical depth
        PSD             %Phase-space density
        cloudAngle      %Angle of the cloud
        becFrac         %Fraction of atoms that are condensed
        becWidth        %Width of the condensed fraction as [xwidth,ywidth]
    end

    properties(SetAccess = protected)
        fitdata         %Fit data as instance of AtomCloudFit
        constants       %Constants associated with atom cloud as instance of AtomImageConstants
    end

    methods
        function self = AtomCloud(constants)
            %ATOMCLOUD Creates an instance of ATOMCLOUD
            %
            %   SELF = ATOMCLOUD() Creates an instance of ATOMCLOUD with
            %   default values
            %
            %   SELF = ATOMCLOUD(CONSTANTS) Uses the input
            %   ATOMIMAGECONSTANTS object as the current object's version.
            %   Does not copy - simply copies the reference
            self.fitdata = AtomCloudFit;
            if nargin > 0
                if isa(constants,'AtomImageConstants')
                    self.constants = constants;
                else
                    error('Input must of type AtomImageConstants');
                end
            else
                self.constants = AtomImageConstants('Rb87');
            end
        end
        
        function self = setConstants(self,constantsIn)
            %SETCONSTANTS Sets the internal CONSTANTS property to be the
            %same reference as the input argument
            %
            %   SELF = SELF.SETCONSTANTS(IN) Sets the internal CONSTANTS
            %   property to the same reference as IN
            self.constants = constantsIn;
        end
        
        function N = sum(self,bg)
            %SUM Computes the number of atoms by summing over the ROI
            %
            %   N = C.SUM() Returns the number of atoms as calculated by
            %   summing over the ROI
            %
            %   N = C.SUM(OFFSET) Subtracts the offset OFFSET from the
            %   image before summing.
            if nargin == 1
                bg = 0;
            end
            %
            % Shorten names
            %
            f = self.fitdata;
            c = self.constants;
            %
            % Get pixel area
            %
            dx = (c.pixelSize/c.magnification)*f.roiStep(1);
            dy = (c.pixelSize/c.magnification)*f.roiStep(2);
            %
            % Correct for user-supplied OD offset and compute sum using the
            % corrected OD map.  Clamp the computed number at a minimum of
            % 0 - no negative atom numbers allowed!
            %
            img = self.fitdata.image - bg;
            N = sum(sum(img))*dx*dy./c.absorptionCrossSection.*(1+4*(c.detuning/c.gamma).^2);
            N = max(N,0);
        end
        
        function self = fit(self,varargin)
            %FIT Performs a fit to the image and extracts information about
            %the sample
            %
            %   SELF = SELF.FIT() performs the fit for ATOMCLOUD object
            %   SELF using default parameters -- namely, the fit type is that
            %   set in SELF.FITDATA and the method for calculating the number
            %   of atoms when using 1D fits is to use the 'xy' marginal
            %   distributions
            %
            %   SELF = SELF.FIT('fittype',FITTYPE,'ex',EX,'method',METHOD)
            %   performs the fit and gets parameters based on the fit type
            %   FITTYPE, 1D extraction method METHOD, and OD exclusion
            %   value EX.  Name/value pairs can be excluded or in any
            %   order.  FITTYPE can be any value accepted by class
            %   ATOMCLOUDFIT. EX is the OD above which data is excluded for
            %   2D fits. METHOD is either 'y', 'x', or 'xy' for calculating
            %   the number of atoms from 1D fits
            %
            
            %
            % Parse input arguments
            %
            if mod(numel(varargin),2) ~= 0
                error('Arguments must appear in name/value pairs!');
            else
                %
                % Default values are to use the fit type specified in
                % FITDATA, exclude no points based on OD, and to calculate
                % the number of atoms using the xy marginal distributions
                %
                fittype = [];
                ex = [];
                calcMethod = 'xy';
                dofit = true;
                for nn = 1:2:numel(varargin)
                    v = varargin{nn+1};
                    switch lower(varargin{nn})
                        case {'fittype','type'}
                            fittype = v;
                        case 'ex'
                            ex = v;
                        case 'method'
                            calcMethod = v;
                        case 'dofit'
                            dofit = v;
                    end
                end
            end
            %
            % Shorten names and fit
            %
            c = self.constants;
            f = self.fitdata;
            if dofit
                f.fit(fittype,ex);
            end
            %
            % Compute effective pixel area, extract fit parameters
            %
            dx = diff(f.x(1:2));dy = diff(f.y(1:2));
            p = f.params;
            %
            % Check that the fit is good
            %
            if ~strcmpi(f.fittype,'sum')
                if f.is1D()
                    if p.gaussAmp(1) < 1.5*std(f.residuals.x)
                        p.gaussAmp(1) = 0;
                    end
                    if p.gaussAmp(2) < 1.5*std(f.residuals.y)
                        p.gaussAmp(2) = 0;
                    end
                else
                    if p.gaussAmp < 1.5*std(f.residuals(:))
                        p.gaussAmp = 0;
                    end
                end
            end
            %
            % Calculate number of atoms
            %
            if strcmpi(f.fittype,'sum')
                %
                % If fittype is 'sum', just use the summed value for the
                % number of atoms
                %
                self.N = self.sum();
                self.Nsum = self.N;
                self.becFrac = NaN;
                self.PSD = NaN;
                self.gaussWidth = nan(1,2);
                self.becWidth = nan(1,2);
                self.pos = nan(1,2);
                self.T = nan(1,2);
                self.peakOD = max(max(f.image));
            else
                %
                % Copy over parameters that don't need extra processing
                %
                self.gaussWidth = p.gaussWidth;
                self.pos = p.pos;
                self.becWidth = p.becWidth;
                self.cloudAngle = p.cloudAngle;
                self.peakOD = max(max(f.image));
                if f.is1D()
                    %
                    % If the fit method is a 1D type, calculate the number
                    % of atoms based on the given calcMethod
                    %
                    switch lower(calcMethod)
                        case 'x'
                            Nth = sqrt(2*pi)*dy*p.gaussAmp(1).*p.gaussWidth(1);
                            Nbec = 16/15*p.becAmp(1).*p.becWidth(1)*dy;
                        case 'y'
                            Nth = sqrt(2*pi)*dx*p.gaussAmp(2).*p.gaussWidth(2);
                            Nbec = 16/15*p.becAmp(2).*p.becWidth(2)*dx;
                        case {'xy','yx'}
                            Nth = sqrt(2*pi*dx*dy)*sqrt(prod(p.gaussAmp.*p.gaussWidth));
                            Nbec = 16/15*sqrt(prod(p.becAmp.*p.becWidth))*sqrt(dx*dy);
                        otherwise
                            error('Only allowed calculation methods for number of atoms are ''x'', ''y'', and ''xy''');
                    end
                else
                    %
                    % If the fit method is not a 1D method (i.e. it is 2D)
                    % then calculate the number of atoms from the widths
                    % and amplitudes directly. This also re-defines the
                    % peak OD based on the fit amplitudes
                    %
                    Nbec = 2*pi/5*p.becAmp*prod(p.becWidth);
                    Nth = p.gaussAmp*(2*pi*prod(p.gaussWidth));
                    self.peakOD = p.becAmp +p.gaussAmp;
                end
                %
                % This converts from the raw number to an actual number
                % based on the absorption cross section and detuning
                %
                self.N = (Nth + Nbec)./c.absorptionCrossSection.*(1+4*(c.detuning/c.gamma).^2);
                %
                % We can calculate the number of atoms by summing and
                % subtracting off the background.  This works best for 2D
                % fits
                %
                if ~f.is1D()
                    self.Nsum = self.sum(self.fitdata.bg);
                else
                    self.Nsum = self.sum();
                end
                %
                % Calculate the BEC fraction and the phase-space density
                %
                self.T = c.calcTemperature(self.gaussWidth);
                self.becFrac = Nbec./(Nth+Nbec);
                self.PSD = self.calcPSD;
            end
        end
        
        function PSD = calcPSD(self,N,T)
            %CALCPSD Calculates the phase-space density based on either the
            %object parameters or supplied parameters
            %
            %   PSD = C.CALCPSD() calculates the PSD based on object C's
            %   number of atoms and temperature
            %
            %   PSD = C.CALCPSD(N,T) calculates the PSD based on the
            %   supplied number N and temperature T along with object C's
            %   trap frequencies and atom mass
            %
            if nargin == 1
                N = self.N;
                T = sqrt(prod(self.T));
                F = self.becFrac;
            else
                F = 0;
            end
            c = self.constants;
            deBroglie = sqrt(2*pi*const.hbar^2./(c.mass.*const.kb.*T));
            freqMean = prod(c.freqs)^(1/3);
            estGaussWidths = sqrt(const.kb*T./(c.mass*freqMean.^2));

            nGauss = (1-F)*N./((2*pi)^1.5*estGaussWidths.^3);
            nBEC = (15*F*N/(8*pi)).^(2/5).*(c.mass*freqMean^2/2).^(3/5);
            n0 = nGauss + nBEC;
            PSD = n0.*deBroglie.^3;
            if nargin == 1
                self.PSD = PSD;
            end
        end

        function varargout = get(self,varargin)
            %GET Takes a vector of AbsorptionImage objects and returns the
            %properties as vectors.
            %
            %   VARARGOUT = C.GET('NAME1','NAME2',...) takes a vector of
            %   AbsorptionImage objects C and returns in VARARGOUT the vectors
            %   of properties 'NAME1', 'NAME2', etc.  Valid values for 'NAME1',
            %   etc correspond to property names in the AbsorptionImage class.
            for nn = 1:numel(self)
                for mm = 1:numel(varargin)
                    varargout{mm}(nn,:) = self(nn).(varargin{mm}); %#ok<*AGROW>
                end
            end
        end
        
        %% Saving and loading functions
        function s = struct(self)
            %STRUCT Creates a struct from the ATOMCLOUD object
            %
            %   S = C.STRUCT() Creates a struct from the ATOMCLOUD
            %   object C
            if numel(self) > 1
                for nn = 1:numel(self)
                    s(nn,1) = struct(self(nn)); %#ok<*AGROW>
                end
            else
                s.N = self.N;
                s.Nsum = self.Nsum;
                s.pos = self.pos;
                s.gaussWidth = self.gaussWidth;
                s.T = self.T;
                s.peakOD = self.peakOD;
                s.PSD = self.PSD;
                s.cloudAngle = self.cloudAngle;
                s.becFrac = self.becFrac;
                s.becWidth = self.becWidth;
                s.constants = struct(self.constants);
                s.fitdata = struct(self.fitdata);
            end
        end
        
        function s = saveobj(self)
            if numel(self) > 1
                for nn = 1:numel(self)
                    s(nn,1) = saveobj(self(nn));
                end
            else
                s = self.struct;
            end
        end
        
    end
    
    methods(Static)
        function self = loadobj(a)
            %LOADOBJ Converts simpler object into ATOMCLOUD
            %
            %   C = LOADOBJ(A) converts simpler object A into instance C
            if numel(a) > 1
                self = AtomCloud.empty;
                for nn = 1:numel(a)
                    self(nn,1) = AtomCloud.loadobj(a(nn));
                end
            else
                self = AtomCloud;
                c = AtomImageConstants.loadobj(a.constants);
                f = AtomCloudFit.loadobj(a.fitdata);
                self.constants.copy(c);
                self.fitdata.copy(f);
                %
                % Properties
                %
                self.N = a.N;
                self.Nsum = a.Nsum;
                self.pos = a.pos;
                self.gaussWidth = a.gaussWidth;
                self.T = a.T;
                self.peakOD = a.peakOD;
                self.PSD = a.PSD;
                self.cloudAngle = a.cloudAngle;
                self.becFrac = a.becFrac;
                self.becWidth = a.becWidth;
            end
        end 
    end


end