classdef AbsorptionImage < handle
    %ABSORPTIONIMAGE Defines a class that describes an absorption image.
    %
    %   Properties of the absorption image are things like the image
    %   itself, the x and y coordinates, the number of atoms, widths,
    %   temperatures, etc.
    %
    %   Three "constant" handle properties are used - RAW, CONSTANTS, and
    %   FITDATA of type RAWIMAGEDATA, ATOMIMAGECONSTANTS, and ATOMCLOUDFIT,
    %   respectively.  These describe the raw image files, the constants
    %   used for extracted useful information about the sample (like image
    %   detuning, pixel size, etc), and the fit function, data, and
    %   methods, respectively
    
    properties
        %
        % Image properties
        %
        x               %x position
        y               %y position
        image           %raw absorption image: -log(imgWithAtoms/imgWithoutAtoms)
        imageCorr       %corrected absorption image
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

    properties(SetAccess = immutable)
        raw             %RAWIMAGEDATA object describing the raw image data
        constants       %ATOMIMAGECONSTANTS object describing the conditions under which the measurement was done
        fitdata         %ATOMCLOUDFIT object describing the fitting
    end

    methods
        function self = AbsorptionImage(varargin)
            %ABSORPTIONIMAGE Creates and ABSORPTIONIMAGE object
            %
            %   C = ABSORPTIONIMAGE() Creates object C with constant
            %   properties set to default values
            %
            %   C = ABSORPTIONIMAGE(RAW) creates object C with raw image
            %   data RAW
            %
            %   C = ABSORPTIONIMAGE(RAW,CONSTANTS) creates object C with
            %   raw image data RAW and constants CONSTANTS
            %
            %   C = ABSORPTIONIMAGE(RAW,CONSTANTS,FITDATA) creates object C
            %   with raw image data RAW, constants CONSTANTS, and fit data
            %   FITDATA
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
            %MAKEIMAGE Creates the uncorrected and corrected absorption
            %images
            %
            %   C = C.MAKEIMAGE() creates the images for object C based on
            %   constant properties
            %
            %   C = C.MAKEIMAGE(IMGINDEXES) creates the images for object C
            %   assuming that the image with atoms is the raw image with
            %   index IMGINDEXES(1) and that the image without atoms is the
            %   raw image with index IMGINDEXES(2)
            
            % Shorten names of constants and raw data
            c = self.constants;
            r = self.raw;
            % Compute the saturation intensity in camera counts
            Nsat = c.satN;
            %
            % Select "image with atoms" and "image without atoms". Either
            % use the default (first and second raw images) or use the
            % user-supplied indices.
            %
            if nargin < 2
                imgIndexes = [1,2];
            end
            imgWithAtoms = r.images(:,:,imgIndexes(1));     %Subtract dark here if desired
            imgWithoutAtoms = r.images(:,:,imgIndexes(2));  %Subtract dark here if desired
            %
            % Create the uncorrected optical depth map and get peak OD.
            % Set NaNs and Infs to zero
            %
            ODraw = real(-log(imgWithAtoms./imgWithoutAtoms));
            ODraw(isnan(ODraw) | isinf(abs(ODraw))) = 0;    %Clean up data so that there are no NaNs or Infs
            self.image = ODraw;
            self.peakOD = max(max(ODraw));
            %
            % Correct for finite saturation OD
            %
            if ~isinf(c.satOD)
                ODmod = real(log((1-exp(-c.satOD))./(exp(-ODraw)-exp(-c.satOD))));
            else
                ODmod = ODraw;
            end
            %
            % Correct for polarization and intensity saturation - this is
            % the corrected optical depth map
            %
            self.imageCorr = c.polarizationCorrection*ODmod + (1 - exp(-ODmod)).*imgWithoutAtoms./Nsat;
            %
            % Create x and y vectors based on pixel size and magnification
            %
            self.x = (c.pixelSize/c.magnification)*(1:size(self.image,2));
            self.y = (c.pixelSize/c.magnification)*(1:size(self.image,1));
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
            %
            % Shorten names
            %
            f = self.fitdata;
            c = self.constants;
            %
            % Extract the region of interest
            %
            row = f.roiRow(1):f.roiRow(2);
            col = f.roiCol(1):f.roiCol(2);
            %
            % Get pixel area
            %
            Apx = (c.pixelSize/c.magnification)^2;
            %
            % Correct for user-supplied OD offset and compute sum using the
            % corrected OD map.  Clamp the computed number at a minimum of
            % 0 - no negative atom numbers allowed!
            %
            img = self.imageCorr(row,col) - offset;
            N = sum(sum(img))*Apx./c.absorptionCrossSection.*(1+4*(c.detuning/c.gamma).^2);
            N = max(N,0);
        end

        function self = fit(self,varargin)
            %FIT Performs a fit to the image and extracts information about
            %the sample
            %
            %   C = C.FIT() performs the fit for ABSORPTIONIMAGE object C
            %   using default parameters -- namely, the fit type is that
            %   set in C.FITDATA and the method for calculating the number
            %   of atoms when using 1D fits is to use the 'y' marginal
            %   distribution
            %
            %   C = C.FIT('fittype',FITTYPE,'ex',EX,'method',METHOD)
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
                % the number of atoms using the y marginal distribution
                %
                fittype = [];
                ex = [];
                calcMethod = 'y';
                for nn = 1:2:numel(varargin)
                    v = varargin{nn+1};
                    switch lower(varargin{nn})
                        case {'fittype','type'}
                            fittype = v;
                        case 'ex'
                            ex = v;
                        case 'method'
                            calcMethod = v;
                    end
                end
            end
            %Shorten names
            c = self.constants;
            f = self.fitdata;
            %
            % Create fit objects and then fit
            %
            f.makeFitObjects(self.x,self.y,self.imageCorr);
            f.fit(fittype,ex);
            %
            % Compute effective pixel area, extract fit parameters
            %
            dx = diff(f.x(1:2));dy = diff(f.y(1:2));
            p = f.params;
            self.gaussWidth = p.gaussWidth;
            self.pos = p.pos;
            self.becWidth = p.becWidth;
            self.cloudAngle = p.cloudAngle;
            self.T = c.calcTemperature(self.gaussWidth);
            self.peakOD = max(max(f.image));
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
            else
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
                        case 'xy'
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
                % If the fit is not a 1D fit, we can use the offset from
                % the 2D fit to correct the number of atoms calculated
                % through summing over the region of interest
                %
                if ~f.is1D()
                    self.Nsum = self.sum(f.params.offset);
                else
                    self.Nsum = self.sum;
                end
                %
                % Calculate the BEC fraction and the phase-space density
                %
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


        %% Plotting functions
        function self = plotROI(self,row,col)
            %PLOTROI Plots the region-of-interest (ROI) on the current axes
            %
            %   C = C.PLOTROI() plots the ROI stored in ABSORPTIONIMAGE
            %   object C
            %
            %   C = C.PLOTROI(ROW,COL) uses the ROI given by ROW and COL to
            %   plot the ROI
            %
            if nargin == 1
                row = self.fitdata.roiRow;col = self.fitdata.roiCol;
            end
            plot([col(1),col(end),col(end),col(1),col(1)],[row(1),row(1),row(end),row(end),row(1)],'r--');
        end

        function str = makeImageSummary(self)
            %MAKEIMAGESUMMARY Creates a string used for summarizing a given
            %absorption image
            %
            %   STR = C.MAKEIMAGESUMMARY() Creates a string STR with the
            %   number of atoms from the fit, the peak OD, the temperature,
            %   and the number of atoms from summing over the ROI
            %
            str{1} = sprintf('N = %1.3g (%d%%)    OD_{peak} = %1.3g    T_{y} = %3.2f uK',...
                self.N,round(self.becFrac*100),self.peakOD,sqrt(prod(self.T))*1e6);
            str{2} = sprintf('Nsum = %1.3g',self.Nsum);
        end

        function self = plotYData(self,col1,col2)
            %PLOTYDATA Plots the marginal y distribution as a vertically
            %oriented plot
            %
            %   C = C.PLOTYDATA plots the data for object C with default
            %   line specs 'b.' for the data and 'r-' for the fit
            %
            %   C = C.PLOTYDATA(COL1) uses line spec COL1 for the data
            %
            %   C = C.PLOTYDATA(__,COL2) uses line spec COL2 for the fit
            %
            if nargin < 2
                col1 = 'b.';
                col2 = 'r-';
            elseif nargin < 3
                col2 = 'r-';
            end
            f = self.fitdata;
            %
            % Plots the marginal Y distribution as a vertically-oriented
            % plot
            %
            plot(f.ydata,-f.y,col1);
            hold on
            plot(f.yfit,-f.y,col2);
            %
            % Creates a summary of the Y distribution
            %
            str{1} = sprintf('Gauss_{y} = %3.1f um',self.gaussWidth(2)*1e6);
            str{2} = sprintf('TF_{y} = %3.1f um',self.becWidth(2)*1e6);
            hold off;
            xlabel(str,'fontsize',8);
        end

        function self = plotXData(self,col1,col2)
            %PLOTXDATA Plots the marginal X distribution as a horizontally
            %oriented plot
            %
            %   C = C.PLOTXDATA plots the data for object C with default
            %   line specs 'b.' for the data and 'r-' for the fit
            %
            %   C = C.PLOTXDATA(COL1) uses line spec COL1 for the data
            %
            %   C = C.PLOTXDATA(__,COL2) uses line spec COL2 for the fit
            %
            if nargin < 2
                col1 = 'b.';
                col2 = 'r-';
            elseif nargin < 3
                col2 = 'r-';
            end
            f = self.fitdata;
            %
            % Plot the marginal X distribution and fit
            %
            plot(f.x,f.xdata,col1);
            hold on
            plot(f.x,f.xfit,col2);
            %
            % Creates a summary of the X distribution
            %
            str{1} = sprintf('Gauss_{x} = %3.1f um',self.gaussWidth(1)*1e6);
            str{2} = sprintf('TF_{x} = %3.1f um',self.becWidth(1)*1e6);
            hold off;
            xlabel(str,'fontsize',8);
        end

        function self = plotAllData(self,dispOD,col1,col2,plotROI)
            %PLOTALLDATA Plots the absorption image and the X and Y
            %marginal distributions on the current figure
            %
            %   C = C.PLOTALLDATA(DISPOD,COL1,COL2) Plots, for object C,
            %   the absorption image and the X and Y distributions on the
            %   current figure. DISPOD is the OD range to display on the
            %   image (using caxis), and COL1 and COL2 are the line specs
            %   to use for the marginal data and fits
            %
            %   C = C.PLOTALLDATA(__,PLOTROI) restricts the image to
            %   PLOTROI.  If PLOTROI is not a cell array, it is interpreted
            %   as a boolean value and restricts the image to the ROI
            %   defined in FITDATA if true.  If it is a cell array, then
            %   this is interpreted as the {ROW,COL} to use for limiting
            %   the image to a given ROI
            %
            if nargin < 5
                plotROI = false;
            end
            %
            % Create the main axes for the image and plot the absorption
            % data
            %
%             subplot(6,6,[2:5 8:11 14:17 20:23]);
            axes('position',[0.3,0.3,0.6,0.65]);
            self.plotAbsData(dispOD,plotROI);
            xlabel(self.makeImageSummary,'fontsize',8);
            %
            % Draw the ROI
            %
            hold on;
            self.plotROI;
            %
            % Plot the Y distribution
            %
%             subplot(6,6,[1 7 13 19]);
            axes('position',[0.075,0.35,0.15,0.6]);
            self.plotYData(col1,col2);
            %
            % Plot the X distribution
            %
%             subplot(6,6, 31:36)
            axes('position',[0.1,0.075,0.8,0.15]);
            self.plotXData(col1,col2);
        end

        function self = plotAbsData(self,dispOD,plotROI)
            %PLOTABSDATA Plots the absorption image data on the current
            %axes
            %
            %   C = C.PLOTABSDATA(DISPOD) Plots the absorption image data
            %   using DISPOD as the caxis limits
            %
            %   C = C.PLOTABSDATA(__,PLOTROI) restricts the image to
            %   PLOTROI.  If PLOTROI is not a cell array, it is interpreted
            %   as a boolean value and restricts the image to the ROI
            %   defined in FITDATA if true.  If it is a cell array, then
            %   this is interpreted as the {ROW,COL} to use for limiting
            %   the image to a given ROI
            %
            if nargin < 3
                plotROI = false;
            end
            %
            % Create the image make it look nice
            %
            imagesc(self.image,dispOD);
            axis equal;
            axis tight;
            colorbar;
            colormap(jet);
            %
            % Set x and y limits to the stored ROI or the user-specified
            % ROI
            %
            if ~iscell(plotROI) && plotROI
                xlim(self.fitdata.roiCol);
                ylim(self.fitdata.roiRow);
            elseif iscell(plotROI)
                xlim(plotROI{2});
                ylim(plotROI{1});
            end
            %
            % Make the title the image number
            %
            imgNums = self.raw.getImageNumbers;
            strTitle = sprintf('Image: %d',imgNums(1));
            title(strTitle,'fontsize',14);
        end

        %% Labelling functions
        function [labelStr,numberStrTotal] = labelOneROI(self)
            %LABELONEROI Creates labels that can be printed on the command
            %line so that the user can see image properties
            %
            %   [LABEL,NUMBER] = C.LABELONEROI() creates a label string
            %   LABEL and a number string with parameters NUMBER
            %
            imgNum = self.raw.getImageNumbers;
            labelCell = {'Image','x width/um','y width/um','Nsum' ,'Nfit' ,'BEC %','PeakOD','T/nk' ,'PSD'};
            formatCell = {'% 5d','%0.3e'     ,'%0.3e'     ,'%0.2e','%0.2e','%0.2f','%0.2e' ,'%0.2e','%0.2e'};
            numberCell = {imgNum(1),self.gaussWidth(1)*1e6,self.gaussWidth(2)*1e6,self.Nsum,self.N,self.becFrac*1e2,self.peakOD,sqrt(prod(self.T))*1e9,self.PSD};
            [labelStr,numberStrTotal] = self.formatLabel(labelCell,formatCell,numberCell);
        end

    end

    methods(Static)
        function [LabelStr,NumberStrTotal] = formatLabel(LabelCell,FormatCell,NumberCell)
            %FORMATLABEL Uses input labels, formats, and values to create a
            %nice label that can be printed to the command line
            %
            %   [LABELSTR,NUMBERSTR] =
            %   ABSORPTIONIMAGE.FORMATLABEL(LABELCELL,FORMATCELL,NUMBERCELL)
            %   creates a label string LABELSTR and number string NUMBERSTR
            %   using cell arrays LABELCELL, FORMATCELL, and NUMBERCELL
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
            %MATCH Attempts to rescale a background image to more closely
            %match an image with atoms while excluding the region where the
            %atoms are
            %
            %   IMGOUT = MATCH(IMG1,IMG2,EXREGION) re-scales IMG2 to match
            %   IMG1 while excluding region EXREGION
            %
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