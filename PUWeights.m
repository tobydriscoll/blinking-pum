classdef PUWeights
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        weights
        diffweights
        diff2weights
        chebweights
    end
    
    properties (Constant)
        overlap = 0.1;
    end
    
    methods
        function obj = PUWeights()
            
            t = obj.overlap;
            
            bump = chebfun( {0,@(x) exp(1-1./(1-x.^2)),0},[-20 -1 1 20]);
            
            R = 1 + obj.overlap;
            bl = @(x) bump((x+1)/R);
            br = @(x) bump((x-1)/R);
            
            %These could be precalculated. This is what I will do.
            weight = chebfun(@(x) bl(x)./(bl(x)+br(x)),'vectorize');
            
            obj.chebweights = [weight 1-weight];
            
            obj.weights = {@(x) 1./(1+exp(4*(1+t)^2*x./((t-x).*(2+t-x).*(t+x).*(2+t+x)))), ...
                           @(x) 1./(1+exp(-4*(1+t)^2*x./((t-x).*(2+t-x).*(t+x).*(2+t+x))))};
            

            %obj.diffweights = @(x) (1+t).^2.*(t+(-1).*x).^(-2).*(2+t+(-1).*x).^(-2).*(t+x).^(-2).*(2+t+x).^(-2).*...
            %                        (t.^2.*(2+t).^2+2.*(2+t.*(2+t)).*x.^2+(-3).*x.^4).*...
            %                        sech(2.*(1+t).^2.*(t+(-1).*x).^(-1).*(2+t+(-1).*x).^(-1).*x.*(t+x).^(-1).*(2+t+x).^(-1)).^2;
            
            obj.diffweights = diff(obj.chebweights);
            obj.diff2weights = diff(obj.chebweights,2);
            
            %obj.diff2weights = @(x) (1+t).^2.*(t+(-1).*x).^(-3).*(2+t+(-1).*x).^(-3).*(t+x).^(-3).*(2+t+x).^(-3).*...
            %                        sech(2.*(1+t).^2.*(t+(-1).*x).^(-1).*(2+t+(-1).*x).^(-1).*x.*(t+x).^(-1).*(2+t+x).^(-1)).^2.*...
            %                        (4.*x.*(3.*t.^2.*(2+t).^2.*(2+t.*(2+t))+(8+(-1).*t.*(2+t).*((-8)+3.*t.*(2+t))).*x.^2+(-3).*(2+t.*(2+t)).*x.^4+3.*x.^6)...
            %                        +(-4).*(1+t).^2.*(t+(-1).*x).^(-1).*(2+t+(-1).*x).^(-1).*(t+x).^(-1).*(2+t+x).^(-1).*(t.^2.*(2+t).^2+2.*(2+t.*(2+t)).*x.^2+(-3).*x.^4).^2.*...
            %                        tanh(2.*(1+t).^2.*(t+(-1).*x).^(-1).*(2+t+(-1).*x).^(-1).*x.*(t+x).^(-1).*(2+t+x).^(-1)));
        end
        
        function ef = evalf(obj,X,MIDINV,k,split_dim,diff_j)
            
            t = obj.overlap;
            
            if nargin == 4
                split_dim =1;
                diff_j =0;
            end
            
            if nargin == 5
                diff_j=0;
            end
            
           % h = @(x) 2/(MIDINV(2)-MIDINV(1))*x-(MIDINV(2)+MIDINV(1))/(MIDINV(2)-MIDINV(1));
           % SCALE = 2/(MIDINV(2)-MIDINV(1));
            h = @(x) 2*t/(MIDINV(2)-MIDINV(1))*x-t*(MIDINV(2)+MIDINV(1))/(MIDINV(2)-MIDINV(1));
            SCALE = 2*t/(MIDINV(2)-MIDINV(1));
            %collect points along the splitting dimension.
            x = X(:,split_dim);
            
            ef = zeros(length(X),1);
            
            X_CENTER = h(x(x>=MIDINV(1) & x<=MIDINV(2)));
            
            switch diff_j
                case 0
                    ef(x<MIDINV(1)) = (k==1);
                    
                    ef(x>=MIDINV(1) & x<=MIDINV(2)) = obj.weights{k}(X_CENTER);
                    
                    ef(x>MIDINV(2)) = (k==2);
                case 1
                    ef(x>=MIDINV(1) & x<=MIDINV(2)) = SCALE*feval(obj.diffweights(:,k),X_CENTER);
                case 2
                    ef(x>=MIDINV(1) & x<=MIDINV(2)) = SCALE^2*feval(obj.diffweights2(:,k),X_CENTER);
                otherwise
                    ef(x>=MIDINV(1) & x<=MIDINV(2)) = SCALE^diff_j*feval(diff(obj.weights(:,k),diff_j),X_CENTER);
            end
        end
    end
end

