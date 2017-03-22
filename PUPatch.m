classdef PUPatch<Patch
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    % To do. figure out way to deal with weights. Maybe have seperate
    % object for weights? IDK.
    properties
        children
        splitting_dim
        overlap
    end
    
    properties (Constant)
        weights = PUWeights();
    end
    
    methods
        
        function obj = PUPatch(domain,cheb_length,children,overlap,splitting_dim)
            obj.domain = domain;
            obj.cheb_length = cheb_length;
            obj.children = children;
            obj.overlap = overlap;
            obj.splitting_dim = splitting_dim;
            
        end
        
        function ln = length(obj)
            ln = obj.cheb_length;
        end
        
        %Implement this later.
        function pts = points(obj)
            pts = [];
        end
        
        function is_refined = PUsplit(obj)
            
            is_refined = true;
            
            for k=1:2
                if obj.children{k}.is_leaf && obj.children{k}.is_refined
                    obj.children{k} = obj.children{k}.splitleaf();
                    is_refined = is_refined && obj.children{k}.is_refined;
                else
                    is_refined = is_refined && obj.children{k}.PUsplit();
                end
            end
        end
        
        function sample(obj,f)
            for k=1:2
                obj.children{k}.sample(f);
            end
        end
        
        function vals = evalf(obj,X,dim,order)
            
            vals = zeros(length(X),order+1);
            
            ind = zeros(length(X),2);
            
            %Figure out indicies of the points in the left and right
            %domains
            for k=1:2
                ind(:,k) = obj.children{k}.InDomain(x);
            end
            
            child_vals{k} = cell(2,1);
            
            %calculate values for the children
            for k=1:2
                child_vals{k} = obj.child{k}.evalf(X(ind(:,k)),dim,order);
            end
            
            for j=0:order
                %Generalized product rule
                for i=0:j
                    
                    %Weights depend only on one variable; don't need to
                    %calculate anything if the splitting dim does not match
                    %the dim we are differentiating in (i.e. the derivative
                    %in these dimensions will be zero).
                    if i>0 && obj.splitting_dim == dim
                        for k=1:2
                            if any(ind(:,k))
                                vals(ind(:,k),j+1) = vals(ind(:,k),j+1) + ...
                                    nchoosek(j,i)*obj.weights.evalf(X(ind(:,k)),obj.overlap,k,obj.splitting_dim,i).*...
                                    child_vals{k}(:,order+1-i);
                            end
                        end
                    end
                end
            end
        end
    end
end

