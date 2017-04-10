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
            [obj.dim,~] = size(obj.domain);
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
            if ~isnumeric(f)
                obj.children{1}.sample(f);
                obj.children{2}.sample(f);
            else
                obj.children{1}.sample(f(1:length(obj.children{1})));
                obj.children{2}.sample(f(length(obj.children{1})+1:end));
            end
        end
        
        function vals = evalfGrid(obj,X,dim,order)
            
            grid_lengths = cellfun(@(x)length(x),X);
            
            [n,~] = size(grid_lengths);
            
            if n>1
                grid_lengths = grid_lengths';
            end
            
            %Put the order at the end; this makes things
            %easier to deal with in matlab
            vals = zeros([grid_lengths order+1]);
            
            sub_grids = cell(2,1);
            child_vals = cell(2,1);
            
            inds = cell(2,1);
            
            for k=1:2
                [sub_grids{k},inds{k}] = obj.children{k}.InDomainGrid(X,obj.splitting_dim);
            end
            

            %calculate values for the children
            for k=1:2
                if ~isempty(sub_grids{k})
                    child_vals{k} = obj.children{k}.evalfGrid(sub_grids{k},dim,order);
                end
            end
            
            
            %Go ahead and shift the splitting diminsion
            %Here I ciculate the dimensions before the order
            dim_permute = circshift(1:obj.dim,[0 -obj.splitting_dim+1]);
            
            if order>0
                vals = permute(vals,[dim_permute obj.dim+1]);
                lengths = size(vals);
                vals = reshape(vals,lengths(1),prod(lengths(2:end-1)),order+1);
                
                for k=1:2
                    child_vals{k} = permute(child_vals{k},[dim_permute obj.dim+1]);
                    c_lengths = size(child_vals{k});
                    child_vals{k} = reshape(child_vals{k},c_lengths(1),prod(c_lengths(2:end-1)),order+1);
                end
            else
                vals = permute(vals,dim_permute);
                lengths = size(vals);
                vals = reshape(vals,lengths(1),prod(lengths(2:end)));
                
                for k=1:2
                    child_vals{k} = permute(child_vals,[dim_permute obj.dim+1]);
                    c_lengths = size(child_vals{k});
                    child_vals{k} = reshape(child_vals,c_lengths(1),prod(c_lengths(2:end)));
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
                            if ~isempty(child_vals{k})
                                weight_vals =  nchoosek(j,ord_i)*obj.weights.evalf(X{obj.splitting_dim}(inds{k}),obj.overlap_in,k,1,j-ord_i);
                                weight_vals = repmat(weight_vals,1,length(child_vals{k}(:,:,ord_i+1)));
                                vals(inds{k},:,j+1) = vals(inds{k},:,j+1) + weight_vals.*child_vals{k}(:,:,ord_i+1);
                            end
                        end
                    end
                end
            end
            
            if order>0
                vals = reshape(vals,lengths);
                vals = ipermute(vals,[dim_permute obj.dim+1]);
            else
                vals = reshape(vals,lengths);
                vals = ipermute(vals,dim_permute);
            end
        end
        
        function vals = evalf(obj,X,dim,order)
            
            
            vals = zeros(length(X),order+1);
            
            ind = false(length(X),2);
            
            %Figure out indicies of the points in the left and right
            %domains
            for k=1:2
                ind(:,k) = obj.children{k}.InDomain(X);
            end
            
            child_vals = cell(2,1);
            
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
        
        function str = toString(obj)
            str = strvcat(strcat('1',obj.children{1}.toString()),strcat('2',obj.children{2}.toString()));
        end
    end
   
end

