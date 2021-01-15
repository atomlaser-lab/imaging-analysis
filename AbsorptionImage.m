classdef AbsorptionImage < handle
    
    properties
        x
        y
        image
        imageCorr

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
    
%     properties(SetAccess = protected)
%         
%     end

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
            r = self.raw;
            Nsat = c.Isat.*(c.pixelSize/c.magnification)^2*c.exposureTime/(const.h*const.c/c.wavelength)/c.photonsPerCount*c.polarizationCorrection*(1+4*(c.detuning/c.gamma).^2);
            if size(r.images,3) == 3
                imgWithAtoms = r.images(:,:,1) - r.images(:,:,3);
                imgWithoutAtoms = r.images(:,:,2) - r.images(:,:,3);
            elseif size(r.images,3) == 2
                imgWithAtoms = r.images(:,:,1);
                imgWithoutAtoms = r.images(:,:,2);
            else
                error('Image sets with %d images are unsupported',size(r.images,3));
            end

            ODraw = real(-log(imgWithAtoms./imgWithoutAtoms));
            self.image = ODraw;
            self.peakOD = max(max(ODraw));
            
            if ~isinf(c.satOD)
                ODmod = real(log((1-exp(-c.satOD))./(exp(-ODraw)-exp(-c.satOD))));
            else
                ODmod = ODraw;
            end
            
            self.imageCorr = c.polarizationCorrection*ODmod + (1 - exp(-ODmod)).*imgWithoutAtoms./Nsat;
            self.x = (c.pixelSize/c.magnification)*(1:size(self.image,2));
            self.y = (c.pixelSize/c.magnification)*(1:size(self.image,1));
        end

        function self = fit(self,fittype,tof,calcmethod,ex)
            c = self.constants;
            f = self.fitdata;
            
            f.makeFitObjects(self.x,self.y,self.imageCorr);
            if nargin < 5
                ex = [];
            end
            f.fit(fittype,ex);
            
            dx = diff(f.x(1:2));dy = diff(f.y(1:2));
            self.tof = tof;
            p = f.params;
            self.gaussWidth = p.gaussWidth;
            self.pos = p.pos;
            self.becWidth = p.becWidth;
            self.cloudAngle = p.cloudAngle;
            self.T = c.calcTemperature(self.gaussWidth,tof);
            self.peakOD = max(max(f.image));

            if f.is1D()
                switch lower(calcmethod)
                    case 'x'
                        Nth = sqrt(2*pi)*dy*p.gaussAmp(1).*p.gaussWidth(1);
%                         Nbec = 8/15*pi*p.becAmp(1).*p.becWidth(1)*dy;
                        Nbec = 16/15*p.becAmp(1).*p.becWidth(1)*dy;
                    case 'y'
                        Nth = sqrt(2*pi)*dx*p.gaussAmp(2).*p.gaussWidth(2);
%                         Nbec = 8/15*pi*p.becAmp(2).*p.becWidth(2)*dx;
                        Nbec = 16/15*p.becAmp(2).*p.becWidth(2)*dx;
                    case 'xy'
                        Nth = sqrt(2*pi*dx*dy)*sqrt(prod(p.gaussAmp.*p.gaussWidth));
%                         Nbec = 8/15*pi*sqrt(prod(p.becAmp.*p.becWidth))*sqrt(dx*dy);
                        Nbec = 16/15*sqrt(prod(p.becAmp.*p.becWidth))*sqrt(dx*dy);
                    otherwise
                        error('Only allowed calculation methods for number of atoms are ''x'', ''y'', and ''xy''');
                end
            else
                Nbec = 2*pi/5*p.becAmp*prod(p.becWidth);
                Nth = p.gaussAmp*(2*pi*prod(p.gaussWidth));
                self.peakOD = p.becAmp +p.gaussAmp;
            end

%             self.N = (Nth + Nbec)./c.absorptionCrossSection*c.polarizationCorrection.*(1+4*(c.detuning/c.gamma).^2);
            self.N = (Nth + Nbec)./c.absorptionCrossSection.*(1+4*(c.detuning/c.gamma).^2);
            self.becFrac = Nbec./(Nth+Nbec);
            self.PSD = self.calcPSD;
        end

        function PSD = calcPSD(self)
            c = self.constants;
            deBroglie = sqrt(2*pi*const.hbar^2./(c.mass.*const.kb.*sqrt(prod(self.T))));
            estGaussWidths = sqrt(const.kb*sqrt(prod(self.T))./(c.mass*c.freqs.^2));

            nGauss = (1-self.becFrac)*self.N./((2*pi)^1.5*prod(estGaussWidths));
            nBEC = (15*self.becFrac*self.N/(8*pi))^(2/5)*(c.mass*prod(c.freqs)^(2/3)/2).^(3/5);
            n0 = nGauss + nBEC;
            self.PSD = n0.*deBroglie^3;
            PSD = self.PSD;
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
            str{1} = sprintf('Gauss_{y} = %3.1f um',self.gaussWidth(2)*1e6);
            str{2} = sprintf('TF_{y} = %3.1f um',self.becWidth(2)*1e6);
            hold off;
%             xlim([min(f.ydata),max(f.ydata)]);
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
            str{1} = sprintf('Gauss_{x} = %3.1f um',self.gaussWidth(1)*1e6);
            str{2} = sprintf('TF_{x} = %3.1f um',self.becWidth(1)*1e6);
            hold off;
            xlabel(str,'fontsize',8);
        end

        function self = plotAllData(self,maxOD,col1,col2)
%             subplot(6,6,[2:5 8:11 14:17 20:23]);
            axes('position',[0.3,0.3,0.6,0.65]);
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
            axes('position',[0.1,0.075,0.8,0.15]);
            self.plotXData(col1,col2);
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
            labelCell = {'Image','x width/um','y width/um','Natoms','BEC %','PeakOD','T/nk','PSD'};
            formatCell = {'% 5d','%0.3e','%0.3e','%0.2e','%0.2f','%0.2e','%0.2e','%0.2e'};
            imgNum = self.raw.getImageNumbers;
            numberCell = {imgNum(1),self.gaussWidth(1)*1e6,self.gaussWidth(2)*1e6,self.N,self.becFrac*1e2,self.peakOD,sqrt(prod(self.T))*1e9,self.PSD};
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