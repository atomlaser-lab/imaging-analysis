classdef AbsorptionImage < handle
    
    properties
        x
        y
        image
        imageCorr

        tof
        N
        Nsum
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

        function self = makeImage(self,imgIndexes)
            c = self.constants;
            r = self.raw;
            Nsat = c.satN;
            if nargin < 2
                imgIndexes = [1,2];
            end
            imgWithAtoms = r.images(:,:,imgIndexes(1));     %Subtract dark here if desired
            imgWithoutAtoms = r.images(:,:,imgIndexes(2));  %Subtract dark here if desired

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
        
        function self = removeBackground(self,exRegion)
            options = optimset('Display','off', 'MaxFunEvals',1000, 'TolFun', 1e-9, 'TolX', 1e-9);
            
            x = 1:size(self.imageCorr,2); %#ok<*PROP>
            y = 1:size(self.imageCorr,1);
            [X,Y] = meshgrid(x,y);
            
            lb = [0,min(x),0,min(y),0,-1e6,-1e6,-10];
            ub = [10,max(x),10*range(x),max(y),10*range(y),1e6,1e6,10];
            
            idx = true(size(self.imageCorr));
            for nn = 1:size(exRegion,1)
                row = exRegion{nn,1};
                col = exRegion{nn,2};
                idx(row,col) = false;
            end
            X = X(idx);Y = Y(idx);img = self.imageCorr(idx);
            s = 20;
            X = X(1:s:end,1:s:end);
            Y = Y(1:s:end,1:s:end);
            img = img(1:s:end);
            position(:,:,1) = X;
            position(:,:,2) = Y;
            guess = [0.1,mean(x),range(x)/2,mean(y),range(y)/2,0,0,0];
            params = lsqcurvefit(@(c,x) AtomCloudFit.gauss2D(c,x),guess,position,img,lb,ub,options);
            
            [X,Y] = meshgrid(x,y);
            p(:,:,1) = X;
            p(:,:,2) = Y;
            bg = AtomCloudFit.gauss2D(params,p);
            self.image = self.image - bg;
            self.imageCorr = self.imageCorr - bg;
            
        end
        
        function N = sum(self,offset)
            %SUM Computes the number of atoms by summing over the ROI
            %
            %   N = C.SUM() Returns the number of atoms as calculated by
            %   summing over the ROI
            %
            %   N = C.SUM(OFFSET) Subtracts the offset OFFSET from the
            %   image before summing
            if nargin == 1
                offset = 0;
            end
            f = self.fitdata;
            c = self.constants;
            row = f.roiRow(1):f.roiRow(2);
            col = f.roiCol(1):f.roiCol(2);
            dx = self.x(2) - self.x(1);
            dy = self.y(2) - self.y(1);
            img = self.imageCorr(row,col) - offset;
            N = sum(sum(img))*dx*dy./c.absorptionCrossSection.*(1+4*(c.detuning/c.gamma).^2);
            N = max(N,0);
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
            
            if strcmpi(f.fittype,'sum')
                self.N = self.sum();
                self.Nsum = self.N;
            else
                if f.is1D()
                    switch lower(calcmethod)
                        case 'x'
                            Nth = sqrt(2*pi)*dy*p.gaussAmp(1).*p.gaussWidth(1);
                            Nbec = 16/15*p.becAmp(1).*p.becWidth(1)*dy;
                        case 'y'
                            Nth = sqrt(2*pi)*dx*p.gaussAmp(2).*p.gaussWidth(2);
                            Nbec = 16/15*p.becAmp(2).*p.becWidth(2)*dx;
                        case 'xy'
                            Nth = sqrt(2*pi*dx*dy)*sqrt(prod(p.gaussAmp.*p.gaussWidth));
                            Nbec = 16/15*sqrt(prod(p.becAmp.*p.becWidth))*sqrt(dx*dy);
                        otherwise
                            error('Only allowed calculation methods for number of atoms are ''x'', ''y'', and ''xy''');
                    end
                else
                    Nbec = 2*pi/5*p.becAmp*prod(p.becWidth);
                    Nth = p.gaussAmp*(2*pi*prod(p.gaussWidth));
                    self.peakOD = p.becAmp +p.gaussAmp;
                end

                self.N = (Nth + Nbec)./c.absorptionCrossSection.*(1+4*(c.detuning/c.gamma).^2);
                if ~f.is1D()
                    self.Nsum = self.sum(f.params.offset);
                else
                    self.Nsum = self.sum;
                end
                self.becFrac = Nbec./(Nth+Nbec);
                self.PSD = self.calcPSD;
            end
        end

        function PSD = calcPSD(self,N,T)
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


        %% Plotting functions
        function self = plotROI(self,row,col)
            if nargin == 1
                row = self.fitdata.roiRow;col = self.fitdata.roiCol;
            end
            plot([col(1),col(end),col(end),col(1),col(1)],[row(1),row(1),row(end),row(end),row(1)],'r--');
        end

        function str = makeImageSummary(self)
            str{1} = sprintf('N = %1.3g (%d%%)    OD_{peak} = %1.3g    T_{y} = %3.2f uK',...
                self.N,round(self.becFrac*100),self.peakOD,sqrt(prod(self.T))*1e6);
            str{2} = sprintf('Nsum = %1.3g',self.Nsum);
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

        function self = plotAllData(self,dispOD,col1,col2,plotROI)
            if nargin < 5
                plotROI = false;
            end
%             subplot(6,6,[2:5 8:11 14:17 20:23]);
            axes('position',[0.3,0.3,0.6,0.65]);
            self.plotAbsData(dispOD,plotROI);
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

        function self = plotAbsData(self,dispOD,plotROI)
            if nargin < 3
                plotROI = false;
            end
            imagesc(self.image,dispOD);
            axis equal;
            axis tight;
            colorbar;
            colormap(jet);
            if ~iscell(plotROI) && plotROI
                xlim(self.fitdata.roiCol);
                ylim(self.fitdata.roiRow);
            elseif iscell(plotROI)
                xlim(plotROI{2});
                ylim(plotROI{1});
            end
            imgNums = self.raw.getImageNumbers;
            strTitle = sprintf('Image: %d',imgNums(1));
            title(strTitle,'fontsize',14);
        end

        %% Labelling functions
        function [labelStr,numberStrTotal] = labelOneROI(self)
            imgNum = self.raw.getImageNumbers;
            labelCell = {'Image','x width/um','y width/um','Nsum' ,'Nfit' ,'BEC %','PeakOD','T/nk' ,'PSD'};
            formatCell = {'% 5d','%0.3e'     ,'%0.3e'     ,'%0.2e','%0.2e','%0.2f','%0.2e' ,'%0.2e','%0.2e'};
            numberCell = {imgNum(1),self.gaussWidth(1)*1e6,self.gaussWidth(2)*1e6,self.Nsum,self.N,self.becFrac*1e2,self.peakOD,sqrt(prod(self.T))*1e9,self.PSD};
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
        
        function imgOut = match(img1,img2,exRegion)
            Filt = ones(size(img1));
            for nn = 1:size(exRegion,1)
                row = exRegion{nn,1};
                col = exRegion{nn,2};
                Filt(row,col) = 0;
            end
            s = fminbnd(@(x) sum(sum(((img1-x*img2).*Filt)).^2),0.1,10);
            imgOut = img2*s;
        end
        
    end

end