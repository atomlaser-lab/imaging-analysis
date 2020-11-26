classdef AbsorptionImage < handle
    
    properties
        x
        y
        image

        tof
        N
        pos
        gaussWidth
        T
        peakOD
        PSD
        cloudAngle
        becFrac
        becWidth
    end

    properties(SetAccess = immutable)
        raw
        constants
        fitdata
    end

    methods
        function self = AbsorptionImage(varargin)
            if nargin < 1
                self.raw = RawImageData;
            else
                self.raw = varargin{1};
            end
            if nargin < 2
                self.constants = AtomImageConstants;
            else
                self.constants = varargin{2};
            end
            if nargin < 3
                self.fitdata = AtomCloudFit;
            else
                self.fitdata = varargin{3};
            end
        end

        function self = makeImage(self)
            c = self.constants;
            Nsat = c.Isat.*(c.pixelSize/c.magnification)^2*c.exposureTime/(const.h*const.c/c.wavelength)/c.photonsPerCount;
            if numel(self.raw.images) == 3
                imgWithAtoms = rawImages.images(:,:,1) - rawImages.images(:,:,3);
                imgWithoutAtoms = rawImages.images(:,:,2) - rawImages.images(:,:,3);
            elseif numel(self.raw.images) == 2
                imgWithAtoms = rawImages.images(:,:,1);
                imgWithoutAtoms = rawImages.images(:,:,2);
            end

            self.image = -c.polarizationCorrection*log(imgWithAtoms./imgWithoutAtoms) + (imgWithoutAtoms - imgWithAtoms)./Nsat;
            if ~isinf(c.satOD)
                self.image = log((1-exp(-c.satOD))./(exp(-self.image)-exp(-c.satOD)));
            end
            self.image = self.image./c.polarizationCorrection;
            self.x = (self.pixelSize/self.magnification)*(1:size(self.image,2));
            self.y = (self.pixelSize/self.magnification)*(1:size(self.image,1));
        end

        

        function self = fit(self,fittype,tof,calcmethod)
            f = self.fitdata.fit(fittype);
            c = self.constants;
            A = (c.pixelSize./c.magnification).^2;
            self.tof = tof;
            self.peakOD = max(max(f.image));
            self.gaussWidth = f.gaussWidth;
            self.pos = f.pos;
            self.becWidth = f.becWidth;
            self.cloudAngle = f.cloudAngle;
            self.T = c.calcTemperature(self.gaussWidth,tof);

            if f.is1D()
                switch lower(calcmethod)
                    case 'x'
                        Nth = sqrt(2*pi*A)*f.gaussAmp(1).*f.gaussWidth(1)./c.absorptionCrossSection*c.polarizationCorrection;
                        Nbec = 8/15*pi*f.becAmp(1).*f.becWidth(1)*sqrt(A)./c.absorptionCrossSection*c.polarizationCorrection;
                    case 'y'
                        Nth = sqrt(2*pi*A)*f.gaussAmp(2).*f.gaussWidth(2)./c.absorptionCrossSection*c.polarizationCorrection;
                        Nbec = 8/15*pi*f.becAmp(2).*f.becWidth(2)*sqrt(A)./c.absorptionCrossSection*c.polarizationCorrection;
                    case 'xy'
                        Nth = sqrt(2*pi*A)*sqrt(prod(f.gaussAmp.*f.gaussWidth))./c.absorptionCrossSection*c.polarizationCorrection;
                        Nbec = 8/15*pi*sqrt(prod(f.becAmp.*f.becWidth))*sqrt(A)./c.absorptionCrossSection*c.polarizationCorrection;
                    otherwise
                        error('Only allowed calculation methods for number of atoms are ''x'', ''y'', and ''xy''');
                end
            else
                Nbec = 0;
                Nth = f.gaussAmp*(2*pi*prod(f.gaussWidth));
            end

            self.N = Nth + Nbec;
            self.becFrac = Nbec./self.N;
            self.PSD = self.calcPSD;
        end

        function PSD = calcPSD(self)
            c = self.constants;
            deBroglie = sqrt(2*pi*const.hbar^2./(f.mass.*const.kb.*sqrt(self.T));
            estGaussWidths = sqrt(const.kb*sqrt(prod(self.T))./(c.mass*c.freqs.^2));

            nGauss = (1-self.becFrac)*self.N./((2*pi)^1.5*prod(estGaussWidths));
            nBEC = (15*self.becFrac*self.N/(8*pi))^(2/5)*(c.mass*prod(c.freqs)^(2/3)/2).^(3/5);
            n0 = nGauss + nBEC;
            self.PSD = n0.*deBroglie^3;
       end


        %% Plotting functions
        function self = plotROI(self)
            row = self.fitdata.roiRow;col = self.fitdata.roiCol;
            plot([col(1),col(end),col(end),col(1),col(1)],[row(1),row(1),row(end),row(end),row(1)],'r--');
        end

        function str = makeImageSummary(self)
            str = sprintf('N = %1.3g (%d%%)    OD_{peak} = %1.3g    T_{y} = %3.2f uK',...
                self.N,round(self.becFrac*100),self.peakOD,sqrt(prod(self.T))*1e6);
        end

        function self = plotYData(self,col1,col2)
            if nargin < 2
                col1 = 'b.';
                col2 = 'r-';
            elseif nargin < 3
                col2 = 'r-';
            end
            f = self.fitdata;
            plot(f.ydata,-f.y,col1);
            hold on
            plot(f.yfit,-f.y,col2);
            str{1} = sprintf('\sigma_{y} = %3.1f um',self.gaussWidth(2)*1e6);
            str{2} = sprintf('R_{y} = %3.1f um',self.becWidth(2)*1e6);
            hold off;
            xlabel(str,'fontsize',8);
        end

        function self = plotXData(self,col1,col2)
            if nargin < 2
                col1 = 'b.';
                col2 = 'r-';
            elseif nargin < 3
                col2 = 'r-';
            end
            f = self.fitdata;
            plot(f.x,f.xdata,col1);
            hold on
            plot(f.x,f.xfit,col2);
            str{1} = sprintf('\sigma_{x} = %3.1f um',self.gaussWidth(1)*1e6);
            str{2} = sprintf('R_{x} = %3.1f um',self.becWidth(1)*1e6);
            hold off;
            xlabel(str,'fontsize',8);
        end

        function self = plotAllData(self,maxOD,col1,col2)
            % subplot(6,6,[2:5 8:11 14:17 20:23]);
            axes('position',[0.3,0.3,0.6,0.7]);
            imagesc(self.image,[0,maxOD]);
            axis equal;
            axis tight;
            colorbar;
            colormap(jet);
            imgNums = self.raw.getImageNumbers;
            strTitle = sprintf('Image: %d',imgNums(1));
            title(strTitle,'fontsize',14);
            xlabel(self.makeImageSummary,'fontsize',8);

            hold on;
            self.plotROI;

            %             subplot(6,6,[1 7 13 19]);
            axes('position',[0.075,0.35,0.15,0.6]);
            self.plotYData(col1,col2);

            %             subplot(6,6, 31:36)
            axes('position',[0.1,0.1,0.8,0.15]);
            self.plotZData(col1,col2);
        end  %End plotAllData

        function self = plotAbsData(self,maxOD)
            imagesc(self.image,[0,maxOD]);
            axis equal;
            axis tight;
            colorbar;
            colormap(jet);
            imgNums = self.raw.getImageNumbers;
            strTitle = sprintf('Image: %d',imgNums(1));
            title(strTitle,'fontsize',14);
        end

        %% Labelling functions
        function [labelStr,numberStrTotal] = labelOneROI(self)
            labelCell = {'Image','y width/um','x width/um','Natoms','PeakOD','T/nk','PSD'};
            formatCell = {'% 5d','%0.3e','%0.3e','%0.2e','%0.2e','%0.2e','%0.2e'};
            imgNum = str2double(self.raw.getImageNumbers);
            numberCell = {imgNum(1),self.yWidth(1)*1e6,self.zWidth(1)*1e6,self.numAtoms(1),self.peakOD(1),self.Ty(1)*1e9,self.PSD};
            [labelStr,numberStrTotal] = self.formatLabel(labelCell,formatCell,numberCell);
        end

    end

    methods(Static)
        function [LabelStr,NumberStrTotal] = formatLabel(LabelCell,FormatCell,NumberCell)
            N=numel(LabelCell);

            LabelStr='';
            for nn=1:N
                LabelMidPoint(nn)=numel(LabelStr)+round(numel(LabelCell{nn})/2);    %#ok
                if nn==N
                    LabelStr=[LabelStr,LabelCell{nn}];  %#ok
                else
                    LabelStr=[LabelStr,LabelCell{nn},'  |  '];  %#ok
                end

            end

            NumberStrTotal=repmat(' ',1,numel(LabelStr));
            for nn=1:N
                NumberStr=num2str(NumberCell{nn},FormatCell{nn});
                NumberLength=numel(NumberStr);
                NumberMidPoint=floor(NumberLength/2);
                NumberStrTotal((LabelMidPoint(nn)-NumberMidPoint+1):(LabelMidPoint(nn)+(NumberLength-NumberMidPoint)))=NumberStr;
            end
            NumberStrTotal(LabelStr=='|')='|';
        end
    end

end