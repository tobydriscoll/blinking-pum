classdef PUWeights
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        weights
        diffweights
        diff2weights
    end
    
    properties (Constant)
        overlap = 0.08;
    end
    
    methods
        function obj = PUWeights()
            
            bump = chebfun( {0,@(x) exp(1-1./(1-x.^2)),0},[-20 -1 1 20]);
            
            R = 1 + obj.overlap;
            bl = @(x) bump((obj.overlap*x+1)/R);
            br = @(x) bump((obj.overlap*x-1)/R);
            
            %These could be precalculated
            weight = chebfun(@(x) bl(x)./(bl(x)+br(x)),'vectorize');
            
            obj.weights = [weight 1-weight];
            
            obj.diffweights = diff(obj.weights);
            
            obj.diff2weights = diff(obj.weights,2);
        end
        
        function ef = evalf(obj,X,MIDINV,k,split_dim,diff_j)
            if nargin == 3
                diff_j=0;
            end
            
            h = @(x) 2/(MIDINV(2)-MIDINV(1))*x-(MIDINV(2)+MIDINV(1))/(MIDINV(2)-MIDINV(1));
            SCALE = 2/(MIDINV(2)-MIDINV(1));
            
            %collect points along the splitting dimension.
            x = X(:,split_dim);
            
            ef = zeros(length(X),1);
            
            X_CENTER = h(x(x>=MIDINV(1) & x<=MIDINV(2)));
            
            switch diff_j
                case 0
                    ef(x<MIDINV(1)) = (k==1);
                    ef(x>=MIDINV(1) & x<=MIDINV(2)) = feval(obj.weights(:,k),X_CENTER);
                    ef(x>MIDINV(2)) = (k==2);
                case 1
                    ef(x>=MIDINV(1) & x<=MIDINV(2)) = SCALE*feval(obj.diffweights(:,k),X_CENTER);
                case 2
                    ef(x>=MIDINV(1) & x<=MIDINV(2)) = SCALE^2*feval(obj.diff2weights(:,k),X_CENTER);
                otherwise
                    ef(x>=MIDINV(1) & x<=MIDINV(2)) = SCALE^diff_j*feval(diff(obj.weights(:,k),diff_j),X_CENTER);
            end
        end
    end
end

