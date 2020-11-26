classdef AtomCloudFit < handle

    properties
        roiRow
        roiCol
        roiStep
        image

        x
        xdata
        xfit
        y
        ydata
        yfit
        
        fittype
        fitfunc
        params
    end

    methods
        function self = AtomCloudFit(varargin)
            if nargin > 0
                self.setup(varargin{:});
            end
        end

        function self = copy(self,obj)
            for p = properties(self)
                self.(p) = obj.(p);
            end
        end

        function self = setup(self,varargin)
            if mod(numel(varargin),2) ~= 0
                error('Arguments must be as name/value pairs');
            else
                for nn = 1:2:numel(varargin)
                    p = lower(varargin{nn});
                    v = varargin{nn+1};
                    switch p
                        case 'roirow'
                            self.roiRow = v;
                        case 'roicol'
                            self.roiCol = v;
                        case 'roiStep'
                            self.roiStep = v;
                        case 'fittype'
                            self.fittype = v;
                        otherwise
                            error('Unsupported parameter %s',p);
                    end
                end
            end
        end

        function set.fittype(self,v)
            v = lower(v);
            switch v
                case 'gauss1d'
                    self.fittype = 'gauss1d';
                case {'twocomp1d','2comp1d'}
                    self.fittype = 'twocomp1d';
                case {'tf1d','bec1d'}
                    self.fittype = 'tf1d';
                case 'gauss2d'
                    self.fittype = 'gauss2d';
                otherwise
                    error('Fit type ''%s'' not supported',v);
            end
        end

        function set.roiRow(self,v)
            if numel(v) ~= 2
                error('Must supply a 2-element vector');
            else
                self.roiRow = v;
            end
        end

        function set.roiCol(self,v)
            if numel(v) ~= 2
                error('Must supply a 2-element vector');
            else
                self.roiCol = v;
            end
        end

        function set.roiStep(self,v)
            if numel(v) == 1
                self.roiStep = [v,v];
            elseif numel(v) == 2
                self.roiStep = v;
            else
                error('Must supply a 1 or 2 element vector');
            end
        end

        function r = is1D(self)
            r = ~strcmpi(self.fittype,'gauss2D');
        end

        function [row,col] = makeROIVectors(self)
            if numel(self.roiStep) == 1
                row = self.roiRow(1):self.roiStep:self.roiRow(2);
                col = self.roiCol(1):self.roiStep:self.roiCol(2);
            elseif numel(self.roiStep) == 2
                row = self.roiRow(1):self.roiStep(1):self.roiRow(2);
                col = self.roiCol(1):self.roiStep(2):self.roiCol(2);
            else
                error('ROI step should be either 1 or 2 elements!');
            end
        end

        function self = makeFitObjects(self,x,y,image)          
            [row,col] = self.makeROIVectors;
            self.image = image(row,col);
            self.x = x(col);
            self.xdata = sum(image(row,col),1);
            self.xfit = [];
            self.y = y(row);
            self.ydata = sum(image(row,col),2);
            self.yfit = [];
        end

        function self = fit(self,fittype)
            if nargin > 0 || isempty(fittype)
                self.fittype = fittype;
            end

            switch self.fittype
                case 'gauss1d'
                    [px,self.xfit] = self.fitGauss1D(self.x,self.xdata);
                    [py,self.yfit] = self.fitGauss1D(self.x,self.ydata);
                    self.params = CloudParameters(px,py);
                case 'twocomp1d'
                    [px,self.xfit] = self.twoComp1D(self.x,self.xdata);
                    [py,self.yfit] = self.twoComp1D(self.x,self.ydata);
                    self.params = CloudParameters(px,py);
                case 'tf1d'
                    [px,self.xfit] = self.fitTF1D(self.x,self.xdata);
                    [py,self.yfit] = self.fitTF1D(self.x,self.ydata);
                    self.params = CloudParameters(px,py);
                case 'gauss2d'
                    [row,col] = self.makeROIVectors;
                    [self.params,f] = self.fitGauss2D(self.x,self.y,self.image,false);
                    self.xfit = sum(f,1);self.yfit = sum(f,2);
                case 'gauss2dangle'
                    [row,col] = self.makeROIVectors;
                    [self.params,f] = self.fitGauss2D(self.x,self.y,self.image,true);
                    self.xfit = sum(f,1);self.yfit = sum(f,2);
            end
        end

    end

    methods(Static)
        function y = gauss1D(c,x)
            y=c(1)+c(2)*exp(-(x-c(3)).^2./(2*c(4).^2))+c(5).*x;
        end

        function y = bec1D(c,x)
            y = c(1) + pi/2*c(2).*(1-((x-c(3))/c(4)).^2).*(abs(x) <= c(4));
        end

        function y = twoComp1D(c,x)
            y0 = c(1);
            ampGauss = c(2);
            x0 = c(3);
            wG = c(4);
            linG = c(5);
            ampBEC = c(6);
            wBEC = c(7);
            y = AbsorptionImage.gauss1D([y0,ampGauss,x0,wG,linG],x) + AbsorptionImage.bec1D([y0,ampBEC,x0,wBEC],x);
        end

        function F = gauss2DAngle(c,Z)
            x = Z(:,:,1);
            y = Z(:,:,2);
            
            F = c(1)*exp(-((x-c(2))*cos(c(9))+(y-c(4))*sin(c(9))).^2./(2*c(3).^2)-(-(x-c(2))*sin(c(9))+(y-c(4))*cos(c(9))).^2./(2*c(5).^2))+c(6)*x+c(7)*y+c(8);
        end

        function F = gauss2D(c,Z)
            x = Z(:,:,1);
            y = Z(:,:,2);
            
            F = c(1)*exp(-(x-c(2)).^2./(2*c(3).^2)-(y-c(4)).^2./(2*c(5).^2))+c(6)*x+c(7)*y+c(8);
        end

        function guess = guessGaussParams(x,y)
            % Estimate of fitting params...
            Az = min(y);
            [Max,I] = max(y);
            Bz = Max - Az;
            Cz = x(I);
            tmp=y-Az;
            tmp=tmp/max(tmp);
            width_point = exp(-0.5);

            I1 = find(smooth(tmp,3) >= width_point, 2, 'first');
            I2 = find(smooth(tmp,3) >= width_point, 2, 'last');
            Dz1 = mean(x(I1));
            Dz2 = mean(x(I2));
            Dztrial = Dz2 - Dz1;

            if length(I1) < 1 || Dztrial < 10e-6,
                Dz = 10e-6;
            else
                Dz = Dztrial;
            end

            guess = CloudParameters(Az,Cz,Bz,Dz,0,0);
        end

        function guess = guessTFParams(x,y)
            Az = min(y);
            [B,I] = max(y);
            Bz = 2/pi*(B-guess.y0);
            Cz = x(I);

            tmp = y - guess.y0;
            tmp = tmp./max(tmp);
            wp = 0.5;

            I1 = find(smooth(tmp,3) >= wp, 2, 'first');
            I2 = find(smooth(tmp,3) >= wp, 2, 'last');
            Dz1 = mean(x(I1));
            Dz2 = mean(x(I2));
            Dztrial = Dz2 - Dz1;

            if length(I1) < 1 || Dztrial < 5e-6,
                Dz = 5e-6;
            else
                Dz = Dztrial;
            end

            guess = CloudParameters(Az,Cz,0,0,Bz,Dz);
        end

        function params = attemptFit(f,guess,x,y,lb,ub,options)
            try
                params = lsqcurvefit(f, guess, x, y, lb, ub, options);
            catch err
                if (strcmp(err.identifier,'optimlib:snls:UsrObjUndefAtX0'))
                    disp('No atoms found in ROI for fit');
                    params = zeros(size(params));
                else
                    rethrow(err);
                end
            end
        end

        function [p,f] = fitGauss1D(x,y)
            x = x(:);y = y(:);
            options = optimset('Display','off', 'MaxFunEvals',100000, 'TolFun', 1e-9, 'TolX', 1e-9);
            ub = [1e6,1e8,10,10,1e6];lb = [-1e6,0,-10,0,-1e6];

            g = AtomCloudFit.guessGaussParams(x,y);
            guess = [g.offset, g.gaussAmp, g.pos, g.gaussWidth, 0];
            func = @(c,x) AtomCloudFit.gauss1D(c,x);
            params = AtomCloudFit.attemptFit(func,guess,x,y,lb,ub,options);
            p = CloudParameters(params(1),params(3),params(2),params(4),0,0);
            f = func(params,x);
         end %End fitGauss1D
         
        function [p,f] = fit2Comp1D(x,y)
            x = x(:);y = y(:);
            options = optimset('Display','off', 'MaxFunEvals',100000, 'TolFun', 1e-9, 'TolX', 1e-9);
            ub = [1e6,1e8,10,10,1e6,50,1];lb=[-1e6,0,-10,0,-1e6,0,0];

            g = AtomCloudFit.guessGaussParams(x,y);
            guess = [g.offset/10,g.gaussAmp/2,g.pos,g.gaussWidth,0,g.gaussAmp,g.gaussWidth];
            func = @(c,x) AtomCloudFit.twoComp1D(c,x);
            params = AtomCloudFit.attemptFit(func,guess,x,y,lb,ub,options);
            p = CloudParameters(params(1),params(3),params(2),params(4),params(6),params(7));
            f = func(params,x);
        end

        function [p,f] = fitTF1D(x,y)
            x = x(:);y = y(:);
            options = optimset('Display','off', 'MaxFunEvals',100000, 'TolFun', 1e-9, 'TolX', 1e-9);
            ub = [1e6,1e6,1,1];lb = [-1e6,-1e6,-1,-1];
            g = AtomCloudFit.guessTFParams(x,y);
            guess = [g.offset,g.becAmp,g.pos,g.becWidth];
            func = @(c,x) AtomCloudFit.bec1D(c,x);
            params = AtomCloudFit.attemptFit(func,guess,x,y,lb,ub,options);
            p = CloudParameters(params(1),params(3),0,0,params(2),params(4));
            f = func(params,x);
        end
         
        function [p,f] = fitGauss2D(x,y,z,includeRotation)
            if nargin < 4
                includeRotation = false;
            end
            options = optimset('Display','off', 'MaxFunEvals',100000, 'TolFun', 1e-9, 'TolX', 1e-9);
            lb = [0,0,0,0,0,-1e6,-1e6,-10];ub = [10,1,1,1,1,1e6,1e6,10];
            gx = AtomCloudFit.guessGaussParams(x,sum(z,1));
            gx = AtomCloudFit.guessGaussParams(y,sum(z,2));
            amp = max(z(:));
            z0 = min(z(:));
  
            % Fitting
            [X,Y]=meshgrid(x,y);
            Position(:,:,1)=X;
            Position(:,:,2)=Y;
            guess = [Amp gx.pos gx.gaussWidth gy.pos gy.gaussWidth 0 0 z0];
            if includeRotation
                lb(end+1) = -pi/6;
                ub(end+1) = pi/6;
                guess(end+1) = 0;
                func = @(c,x) AtomCloudFit.gauss2DAngle(c,x);
            else
                func = @(c,x) = AtomCloudFit.gauss2D(c,x);
            end
            params = AtomCloudFit.attemptFit(func,guess,Position,z,lb,ub,options);
            p = CloudParameters(params(8),params([2,4]),params(1),params([3,5]),0,0);
            f = func(params,Position);
        end  %End fitGauss2D
    end




end