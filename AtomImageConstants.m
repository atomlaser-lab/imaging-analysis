classdef AtomImageConstants < handle

    properties
        atomType
        mass
        freqs

        pixelSize
        magnification

        exposureTime
        photonsPerCount

        satOD
        Isat
        wavelength
        absorptionCrossSection
        polarizationCorrection
    end

    methods

        function self = AtomImageConstants(atomType,varargin)
            if nargin > 0
                self.useDefaults;
                self.setup(varargin);
            end
        end

        function self = copy(self,obj)
            for p = properties(self)
                self.(p) = obj.(p);
            end
        end

        function self = setup(self,atomType,varargin)
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
                        case 'exposuretime'
                            self.exposureTime = v;
                        case 'photonspercount'
                            self.photonsPerCount = v;
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
            self.pixelSize = 5.67e-6;
            self.magnification = 0.97;
            self.exposureTime = 30e-6;
            self.photonsPerCount = 1;

            if strcmpi(atomType,'Rb87')
                self.atomType = 'Rb87';
                self.mass = const.mRb;
                self.freqs = 2*pi*[220,230,100];
                self.satOD = Inf;
                self.Isat = 1.665e-3*1e4;
                self.wavelength = 780.241e-9;
                self.absorptionCrossSection = 3*self.wavelength^2/(2*pi);
                self.polarizationCorrection = 15/8;
                error('Atom type not supported!');
            end
        end

        function T = calcTemperature(self,width,tof)
            T = width.^2./(1+(self.freqs(1:2)*tof).^2).*self.mass./const.kb;
        end

    end

end
