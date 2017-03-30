classdef PUPatch<Patch
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    % To do. figure out way to deal with weights. Maybe have seperate
    % object for weights? IDK.
    properties
        children
        splitting_dim
        overlap_in
    end
    
    properties (Constant)
        weights = PUWeights();
    end
    
    methods
        
        function obj = PUPatch(domain,overlap_in,cheb_length,children,splitting_dim)
            obj.domain = domain;
            obj.overlap_in = overlap_in;
            obj.cheb_length = cheb_length;
            obj.children = children;
            obj.splitting_dim = splitting_dim;
            obj.is_leaf = false;
            obj.is_refined = false;
            
        end
        
        function ln = length(obj)
            ln = obj.cheb_length;
        end
        
        %Implement this later.
        function pts = points(obj)
            pts = [obj.children{1}.points();obj.children{2}.points()];
        end
        
        function is_refined = PUsplit(obj)
            
            is_refined = true;
            
            for k=1:2
                
                if ~obj.children{k}.is_refined
                    if obj.children{k}.is_leaf
                        obj.children{k} = obj.children{k}.splitleaf();
                        is_refined = is_refined && obj.children{k}.is_refined;
                    else
                        is_refined = is_refined && obj.children{k}.PUsplit();
                    end
                end
            end
            obj.cheb_length = obj.children{1}.cheb_length+obj.children{2}.cheb_length;
            obj.is_refined = is_refined;
        end
        
        function sample(obj,f)
            obj.children{1}.sample(f(1:length(obj.children{1})));
            obj.children{2}.sample(f(length(obj.children{1})+1:end));
        end
        
        function vals = evalf(obj,X,dim,order)
            
            vals = zeros(length(X),order+1);
            
            ind = false(length(X),2);
            
            %Figure out indicies of the points in the left and right
            %domains
            for k=1:2
                ind(:,k) = obj.children{k}.InDomain(X);
            end
            
            child_vals{k} = cell(2,1);
            
            %calculate values for the children
            for k=1:2
                if any(ind(:,k))
                    child_vals{k} = obj.children{k}.evalf(X(ind(:,k),:),dim,order);
                end
            end
            
            for j=0:order
                
                %Generalized product rule
                for ord_i=0:j
                    
                    %Weights depend only on one variable; don't need to
                    %calculate anything if the splitting dim does not match
                    %the dim we are differentiating in (i.e. the derivative
                    %in these dimensions will be zero).
                    if (j-ord_i)==0 || obj.splitting_dim == dim
                        for k=1:2
                            if any(ind(:,k))
                                vals(ind(:,k),j+1) = vals(ind(:,k),j+1) + ...
                                    nchoosek(j,ord_i)*obj.weights.evalf(X(ind(:,k),:),obj.overlap_in,k,obj.splitting_dim,j-ord_i).*...
                                    child_vals{k}(:,ord_i+1);
                            end
                        end
                    end
                end
            end
        end
    end
end

