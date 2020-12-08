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
            p = properties(self);
            for nn = 1:numel(p)
                self.(p{nn}) = obj.(p{nn});
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
                        case 'roistep'
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
                case 'none'
                    self.fittype = 'none';
                case 'gauss1d'
                    self.fittype = 'gauss1d';
                case {'twocomp1d','2comp1d'}
                    self.fittype = 'twocomp1d';
                case {'tf1d','bec1d'}
                    self.fittype = 'tf1d';
                case 'gauss2d'
                    self.fittype = 'gauss2d';
                case {'bec2d','tf2d'}
                    self.fittype = 'tf2d';
                case {'twocomp2d','2comp2d'}
                    self.fittype = 'twocomp2d';
                otherwise
                    error('Fit type ''%s'' not supported',v);
            end
        end

        function set.roiRow(self,v)
            if numel(v) ~= 2
                error('Must supply a 2-element vector');
            else
                v(1) = max(v(1),1);
                v(2) = min(v(2),RawImageData.SIZE(1));
                self.roiRow = v;
            end
        end

        function set.roiCol(self,v)
            if numel(v) ~= 2
                error('Must supply a 2-element vector');
            else
                v(1) = max(v(1),1);
                v(2) = min(v(2),RawImageData.SIZE(2));
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
            r = ~strcmpi(self.fittype,'gauss2D') && ~strcmpi(self.fittype,'twocomp2d') && ~strcmpi(self.fittype,'tf2d') && ~strcmpi(self.fittype,'none');
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
            self.y = y(row);
            self.ydata = sum(image(row,col),2);
            
            
            idx = ~(isnan(self.xdata) | (imag(self.xdata) ~= 0) | isinf(self.xdata));
            self.xdata = self.xdata(idx);
            self.x = self.x(idx);
            self.xfit = zeros(size(self.x));
            idx = ~(isnan(self.ydata) | (imag(self.ydata) ~= 0) | isinf(self.ydata));
            self.ydata = self.ydata(idx);
            self.y = self.y(idx);
            self.yfit = zeros(size(self.y));
        end

        function self = fit(self,fittype)
            if nargin > 0 && ~isempty(fittype)
                self.fittype = fittype;
            end

            switch self.fittype
                case 'none'
                    self.params = CloudParameters();
                case 'gauss1d'
                    [px,self.xfit] = self.fitGauss1D(self.x,self.xdata);
                    [py,self.yfit] = self.fitGauss1D(self.y,self.ydata);
                    self.params = CloudParameters(px,py);
                case 'twocomp1d'
                    [px,self.xfit] = self.fit2Comp1D(self.x,self.xdata);
                    [py,self.yfit] = self.fit2Comp1D(self.y,self.ydata);
                    self.params = CloudParameters(px,py);
                case 'tf1d'
                    [px,self.xfit] = self.fitTF1D(self.x,self.xdata);
                    [py,self.yfit] = self.fitTF1D(self.y,self.ydata);
                    self.params = CloudParameters(px,py);
                case 'gauss2d'
                    [self.params,f] = self.fitGauss2D(self.x,self.y,self.image,false);
                    self.xfit = sum(f,1);self.yfit = sum(f,2);
                case 'gauss2dangle'
                    [self.params,f] = self.fitGauss2D(self.x,self.y,self.image,true);
                    self.xfit = sum(f,1);self.yfit = sum(f,2);
                case 'twocomp2d'
                    [self.params,f] = self.fitTwoComp2D(self.x,self.y,self.image);
                    self.xfit = sum(f,1);self.yfit = sum(f,2);
                case 'tf2d'
                    [self.params,f] = self.fitTF2D(self.x,self.y,self.image);
                    self.xfit = sum(f,1);self.yfit = sum(f,2);
                otherwise
                    error('Fit type %s not supported',self.fittype);
            end
        end

    end

    methods(Static)
        function y = gauss1D(c,x)
            y=c(1)+c(2)*exp(-(x-c(3)).^2./(2*c(4).^2))+c(5).*x;
        end

        function y = bec1D(c,x)
            y = c(1) + c(2).*(1-((x-c(3))/c(4)).^2).^2.*(abs(x-c(3)) <= c(4))+c(5).*x;
        end

        function y = twoComp1D(c,x)
            y0 = c(1);
            ampGauss = c(2);
            x0 = c(3);
            wG = c(4);
            lin = c(5);
            ampBEC = c(6);
            wBEC = c(7);
            y = AtomCloudFit.gauss1D([y0,ampGauss,x0,wG,lin],x) + AtomCloudFit.bec1D([y0,ampBEC,x0,wBEC,lin],x);
        end

        function F = gauss2DAngle(c,Z)
            x = Z(:,:,1); %#ok<*PROP>
            y = Z(:,:,2);
            
            F = c(1)*exp(-((x-c(2))*cos(c(9))+(y-c(4))*sin(c(9))).^2./(2*c(3).^2)-(-(x-c(2))*sin(c(9))+(y-c(4))*cos(c(9))).^2./(2*c(5).^2))+c(6)*x+c(7)*y+c(8);
        end

        function F = gauss2D(c,Z)
            x = Z(:,:,1);
            y = Z(:,:,2);
            
            F = c(1)*exp(-(x-c(2)).^2./(2*c(3).^2)-(y-c(4)).^2./(2*c(5).^2))+c(6)*x+c(7)*y+c(8);
        end
        
        function F = bec2D(c,Z)
            n0 = c(1);
            x0 = c(2);
            xw = c(3);
            y0 = c(4);
            yw = c(5);
            linx = c(6);
            liny = c(7);
            offset = c(8);
            
            x = (Z(:,:,1)-x0)/xw;
            y = (Z(:,:,2)-y0)/yw;
            s2 = x.^2+y.^2;
            F = n0*((1-s2).*(s2<=1)).^1.5+linx.*x+liny.*y+offset;
        end
        
        function F = twoComp2D(c,Z)
            pg = c(1:8);
            pbec = [c(9),c(2),c(10),c(4),c(11),0,0,0];
            F = AtomCloudFit.gauss2D(pg,Z) + AtomCloudFit.bec2D(pbec,Z);
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
            Bz = 2/pi*(B-Az);
            Cz = x(I);

            tmp = y - Az;
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
                    params = zeros(size(guess));
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
            guess = [g.offset/10,g.gaussAmp/2,g.pos,g.gaussWidth,0,g.gaussAmp/2,g.gaussWidth/2];
            func = @(c,x) AtomCloudFit.twoComp1D(c,x);
            params = AtomCloudFit.attemptFit(func,guess,x,y,lb,ub,options);
            p = CloudParameters(params(1),params(3),params(2),params(4),params(6),params(7));
            f = func(params,x);
        end

        function [p,f] = fitTF1D(x,y)
            x = x(:);y = y(:);
            options = optimset('Display','off', 'MaxFunEvals',100000, 'TolFun', 1e-9, 'TolX', 1e-9);
            ub = [1e6,1e6,1,1,1e6];lb = [-1e6,-1e6,-1,-1,-1e6];
            g = AtomCloudFit.guessTFParams(x,y);
            guess = [g.offset,g.becAmp,g.pos,g.becWidth,0];
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
            gy = AtomCloudFit.guessGaussParams(y,sum(z,2));
            amp = max(z(:));
            z0 = min(z(:));
  
            % Fitting
            [X,Y]=meshgrid(x,y);
            Position(:,:,1)=X;
            Position(:,:,2)=Y;
            guess = [amp gx.pos gx.gaussWidth gy.pos gy.gaussWidth 0 0 z0];
            if includeRotation
                lb(end+1) = -pi/6;
                ub(end+1) = pi/6;
                guess(end+1) = 0;
                func = @(c,x) AtomCloudFit.gauss2DAngle(c,x);
            else
                func = @(c,x) AtomCloudFit.gauss2D(c,x);
            end
            params = AtomCloudFit.attemptFit(func,guess,Position,z,lb,ub,options);
            p = CloudParameters(params(8),params([2,4]),params(1),params([3,5]),0,[0,0]);
            f = func(params,Position);
        end
        
        function [p,f] = fitTwoComp2D(x,y,z)
            options = optimset('Display','off', 'MaxFunEvals',100000, 'TolFun', 1e-9, 'TolX', 1e-9);
            lb = [0,0,0,0,0,-1e6,-1e6,-10,0,0,0];ub = [10,1,1,1,1,1e6,1e6,10,10,1,1];
            gx = AtomCloudFit.guessGaussParams(x,sum(z,1));
            gy = AtomCloudFit.guessGaussParams(y,sum(z,2));
            amp = max(z(:));
            z0 = min(z(:));
  
            % Fitting
            [X,Y]=meshgrid(x,y);
            Position(:,:,1)=X;
            Position(:,:,2)=Y;
            guess = [amp/5 gx.pos gx.gaussWidth gy.pos gy.gaussWidth 0 0 z0 amp gx.gaussWidth/2 gy.gaussWidth/2];
            func = @(c,x) AtomCloudFit.twoComp2D(c,x);
            params = AtomCloudFit.attemptFit(func,guess,Position,z,lb,ub,options);
            p = CloudParameters(params(8),params([2,4]),params(1),params([3,5]),params(9),params(10:11));
            f = func(params,Position);
        end
        
        function [p,f] = fitTF2D(x,y,z)
            options = optimset('Display','off', 'MaxFunEvals',100000, 'TolFun', 1e-9, 'TolX', 1e-9);
            lb = [0,0,0,0,0,-1e6,-1e6,-10];ub = [10,1,1,1,1,1e6,1e6,10];
            gx = AtomCloudFit.guessTFParams(x,sum(z,1));
            gy = AtomCloudFit.guessTFParams(y,sum(z,2));
            amp = max(z(:));
            z0 = min(z(:));
  
            % Fitting
            [X,Y]=meshgrid(x,y);
            Position(:,:,1)=X;
            Position(:,:,2)=Y;
            guess = [amp gx.pos gx.becWidth gy.pos gy.becWidth 0 0 z0];
            func = @(c,x) AtomCloudFit.bec2D(c,x);
            params = AtomCloudFit.attemptFit(func,guess,Position,z,lb,ub,options);
            p = CloudParameters(params(8),params([2,4]),0,[0,0],params(1),params([3,5]));
            f = func(params,Position);
        end
    end




end