classdef AtomImageConstants < handle
    %ATOMIMAGECONSTANTS Defines an object that describes constants that
    %can be used to describe an absorption image
    properties
        atomType                    %Atom type, currently only 'Rb87'
        mass                        %Mass of the atom
        freqs                       %Trapping frequencies in rad/s, 3-element vector

        pixelSize                   %Size of pixels on camera
        magnification               %Magnification of the imaging system
        tof                         %Time-of-flight of the sample

        exposureTime                %Exposure time of light on camera
        photonsPerCount             %Number of photons per "count" on the camera
        detuning                    %Detuning of the imaging beam in MHz
        gamma                       %FWHM of the resonance, in MHz

        satOD                       %Saturation OD of the imaging system
        Isat                        %Saturation intensity in W/m^2 of the transition
        wavelength                  %Wavelength of the transition
        absorptionCrossSection      %On-resonance, stretched-state absorption cross section 3*wavelength^2/(2*pi)
        polarizationCorrection      %Polarization correction  s.t. cross-section -> cross-section/polarizationCorrection
    end

    methods

        function self = AtomImageConstants(atomType,varargin)
            %ATOMIMAGECONSTANTS Creates an object with constants describing
            %the imaging process
            %
            %   C = ATOMIMAGECONSTANTS(ATOMTYPE) Creates object C using
            %   atom type ATOMTYPE.  Currently it can only be 'Rb87'.
            %
            %   C = ATOMIMAGECONSTANTS(__,NAME,VALUE,...) NAME/VALUE pairs
            %   describing different properties. See SETUP()
            %
            if nargin > 0
                self.useDefaults(atomType);
                self.setup(atomType,varargin{:});
            end
        end

        function self = copy(self,obj)
            %COPY Copies the object properties from input
            %
            %   C = C.COPY(OBJ) copies the object properties from input
            %   ATOMIMAGECONSTANTS object OBJ
            %
            p = properties(self);
            for nn = 1:numel(p)
                self.(p{nn}) = obj.(p{nn});
            end
        end

        function self = setup(self,atomType,varargin)
            %SETUP Sets up ATOMIMAGECONSTANTS based on input arguments
            %
            %   C = C.SETUP(ATOMTYPE,NAME,VALUE,...) uses ATOMTYPE ('Rb87')
            %   and NAME/VALUE pairs to define properties. Allowed values
            %   for NAME are any of the object property names (case
            %   insensitive).
            %
            self.atomType = atomType;
            if mod(numel(varargin),2) ~= 0
                error('Variable argument list must in name/value pairs!');
            else
                for nn = 1:2:numel(varargin)
                    v = varargin{nn+1};
                    switch lower(varargin{nn})
                        case 'freqs'
                            self.freqs = v;
                        case 'pixelsize'
                            self.pixelSize = v;
                        case 'magnification'
                            self.magnification = v;
                        case 'tof'
                            self.tof = v;
                        case 'exposuretime'
                            self.exposureTime = v;
                        case 'photonspercount'
                            self.photonsPerCount = v;
                        case 'detuning'
                            self.detuning = v;
                        case 'gamma'
                            self.gamma = v;
                        case 'satod'
                            self.satOD = v;
                        case 'isat'
                            self.Isat = v;
                        case 'wavelength'
                            self.wavelength = v;
                        case 'absorptioncrosssection'
                            self.absorptionCrossSection = v;
                        case 'polarizationcorrection'
                            self.polarizationCorrection = v;
                        otherwise
                            error('Option %s not recognized',varargin{nn});
                    end
                end
            end
        end

        function self = useDefaults(self,atomType)
            %USEDEFAULTS Sets the default values for the object
            %
            %   C = C.USEDEFAULTS(ATOMTYPE) sets the default values for
            %   object C given atom type ATOMTYPE (only 'Rb87' is
            %   supported)
            self.pixelSize = 6.45e-6;
            self.magnification = 0.97;
            self.tof = 35e-3;
            self.exposureTime = 30e-6;
            self.photonsPerCount = 17e3/(2^16*0.174);
            self.detuning = 0;
            self.gamma = 6.065;

            if strcmpi(atomType,'Rb87')
                self.atomType = 'Rb87';
                self.mass = const.mRb;
                self.freqs = 2*pi*[220,230,100];
                self.satOD = Inf;
                self.Isat = 1.665e-3*1e4;
                self.wavelength = 780.241e-9;
                self.absorptionCrossSection = 3*self.wavelength^2/(2*pi);
                self.polarizationCorrection = 15/8;
            else
                error('Atom type ''%s'' not supported!',atomType);
            end
        end

        function T = calcTemperature(self,width,tof)
            %CALCTEMPERATURE calculates the temperature of a sample given a
            %width and a time of flight
            %
            %   T = C.CALCTEMPERATURE(WIDTH) Calculates the temperature T
            %   using the supplied WIDTH in meters and the property C.tof.5
            %
            %   T = C.CALCTEMPERATURE(__,TOF) uses the supplied time of
            %   flight TOF for calculating the temperature
            %
            if nargin < 3
                tof = self.tof;
            end
            T = (self.freqs(1:2).*width).^2./(1+(self.freqs(1:2)*tof).^2).*self.mass./const.kb;
        end
        
        function Nsat = satN(self)
            %SATN Computes the number of camera "counts" that corresponds
            %to the saturation intensity, after correcting for polarization
            %and detuning
            %
            %   NSAT = C.SATN() returns the saturation number of counts
            %   NSAT
            Nsat = self.Isat.*(self.pixelSize/self.magnification)^2*self.exposureTime/(const.h*const.c/self.wavelength)/self.photonsPerCount*(1+4*(self.detuning/self.gamma).^2);
        end

    end

end
