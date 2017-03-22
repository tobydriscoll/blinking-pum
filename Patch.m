classdef (Abstract) Patch < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        domain
        boundaryIndex
        cheb_length
        is_leaf
        is_refined
        dim
    end
    
    methods (Abstract)
        points = points(obj)
        ef = evalf(obj,x,dim,order)
        ln = length(obj)
        sample(obj,f)
    end
    
    methods
        %Right now just assume the domain is just the square.
        function domain_ind = InDomain(obj,x)
            domain_ind = ones(length(x),1);
            for i=1:obj.dim
                domain_ind = domain_ind && ...
                    ( x(:,i)>= obj.domain(i,1) && x(:,i)<= obj.domain(i,2));
            end
        end
    end
    
end

