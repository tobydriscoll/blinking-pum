classdef LeafPatch<Patch
    % This is the abstract class for a leaf object. This is used with the
    % PUPatch object.
 
% LeafPatch(varargin) serves as the base constructor for a leafPatch
% object. Here PUFun(f,'perf1',perf1,'pref2',pref2,..) is called, with
% options:
% 
%
% *The domain used for the Chebyshev polynomial: 'domain', [a,b;c,d]
%
% *The zone (non overlapping part from partition) used: 'zone', [a,b;c,d]
%
% *The domain of the root of the tree: 'outerbox', [a,b;c,d]
%
% LeafPatch(struct) will construct an approximation with a structure
% struct. Here struct must contain the following fields:
% outerbox : domain of the outerbox
% zone : domain of the zone
% domain : square domain of the polynomial
    properties
        index = [];
        chebweights = [];
        cheb_bump = chebfun( {0,@(x) exp(1-1./(1-x.^2)),0},[-20 -1 1 20]);
        %cheb_bump = @(x) exp(1-1./(1-x.^2));
        bump
    end
    
    methods
        function obj = LeafPatch(varargin)
            
            if length(varargin)==1
                varargin = varargin{:};
            end
            
            if isstruct(varargin)
                
                obj.domain = varargin.domain;
                obj.outerbox = varargin.outerbox;
                obj.zone = varargin.zone;
                
            else
                
                args = varargin;
                while (~isempty(args))
                    if strcmpi(args{1}, 'domain')
                        obj.domain = args{2};
                        [obj.dim,~] = size(obj.domain);
                    elseif strcmpi(args{1}, 'zone')
                        obj.zone = args{2};
                    elseif strcmpi(args{1}, 'outerbox')
                        obj.outerbox = args{2};
                    end
                    args(1:2) = [];
                end
            end
            
            if isempty(obj.domain)
                error('You need a domain!');
            end
            
            if isempty(obj.zone)
                obj.zone = obj.domain;
            end
            
            if isempty(obj.outerbox)
                obj.outerbox = obj.domain;
            end
            
            obj.is_geometric_refined = true; %square is always geometrically refined
            [obj.dim,~] = size(obj.domain);
            obj.is_leaf = true;
            
            obj.bump = cell(3,1);
            
            for k=1:obj.dim
               % if isequal(obj.domain(k,:),obj.outerbox(k,:))
               %     w = @(x) ones(size(x));
               if obj.domain(k,1) == obj.outerbox(k,1)
                    w = @(x) obj.cheb_bump((obj.invf(x,obj.domain(k,:))+1)/2);
                elseif obj.domain(k,2) == obj.outerbox(k,2)
                    w = @(x) obj.cheb_bump((obj.invf(x,obj.domain(k,:))-1)/2);
                else
                    w = @(x) obj.cheb_bump(obj.invf(x,obj.domain(k,:)));
                end
                obj.bump{k} = w;
            end
        end
        
        
        % Evaluates the bump function of a patch.
        %
        %  Input:
        %      X: set of points to evaluate at
        %
        % Output:
        %     ef: array of length(X) of the patch's weight evaluated at X.
        function ef = evalfBump(obj,X)
            
            ef = ones(size(X,1),1);
            for k=1:obj.dim
                ef = ef.*obj.bump{k}(X(:,k));
            end
            
        end
        
        % Evaluates the bump function of a patch.
        %
        %  Input:
        %      X: cellarray of grids to be evaluated on.
        %
        % Output:
        %     ef: array of dim(X) of the patch's weight evaluated at X.
        function ef = evalfGridBump(obj,X)
            
            W = cell(3,1);
            
            for i=1:obj.dim
                W{i} = obj.bump{i}(X{i});
            end
            
            if obj.dim==2
                ef = W{1}*W{2}.';
            else
                ef = reshape(W{3},1,1,length(W{3})).*(W{2}'.*W{1});
            end
            
        end
        
        % Plots the overlaping domain of the tree. Works only
        % in 2D and 3D.
        %
        %     Input:
        %     color: color of marker for the center of the domain.
        %
        function plotdomain(obj,color)
            
            if nargin==1
                color = 'black';
            end
            
            if obj.dim==2
                hold on;
                lengths = [diff(obj.domain(1,:));diff(obj.domain(2,:))];
                rectangle('position',[obj.domain(:,1)' lengths'],'LineWidth',2,'EdgeColor',color);
                plot(mean(obj.domain(1,:)),mean(obj.domain(2,:)),'.','MarkerSize',10,'Color',color);
                hold off;
            elseif obj.dim==3
                hold on;
                lengths = [diff(obj.domain(1,:));diff(obj.domain(2,:));diff(obj.domain(3,:))];
                center = sum(obj.domain,2)/2;
                %Vertices for Line Cube. Order matters
                X = [-1 -1 1 1 -1 -1 1 1 1 1 1 1 -1 -1 -1 -1 -1]';
                Y = [-1 1 1 -1 -1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1]';
                Z = [-1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 1 -1 1 1 -1]';
                %Example two cube matrix. Unit cube and one scaled/translated cube
                X1 = X*lengths(1)/2+center(1);
                Y1 = Y*lengths(2)/2+center(2);
                Z1 = Z*lengths(3)/2+center(3);
                %Single plot command for all 'cube lines'
                plot3(X1,Y1,Z1,'color','black');
                plot3(mean(obj.domain(1,:)),mean(obj.domain(2,:)),mean(obj.domain(3,:)),'.','MarkerSize',10,'Color','black');
                hold off;
            end
        end
        
        function List = collectLeaves(obj)
            List = {obj};
        end
        
        % Plots the zones of the tree. Works only
        % in 2D and 3D.
        %
        function plotzone(obj)
            
            if obj.dim==2
                hold on;
                lengths = [diff(obj.zone(1,:));diff(obj.zone(2,:))];
                rectangle('position',[obj.zone(:,1)' lengths'],'LineWidth',2);
                hold off;
            elseif obj.dim==3
                hold on;
                lengths = [diff(obj.zone(1,:));diff(obj.zone(2,:));diff(obj.zone(3,:))];
                center = sum(obj.zone,2)/2;
                %Vertices for Line Cube. Order matters
                X = [-1 -1 1 1 -1 -1 1 1 1 1 1 1 -1 -1 -1 -1 -1]';
                Y = [-1 1 1 -1 -1 -1 -1 -1 -1 1 1 1 1 1 1 -1 -1]';
                Z = [-1 -1 -1 -1 -1 1 1 -1 1 1 -1 1 1 -1 1 1 -1]';
                %Example two cube matrix. Unit cube and one scaled/translated cube
                X1 = X*lengths(1)/2+center(1);
                Y1 = Y*lengths(2)/2+center(2);
                Z1 = Z*lengths(3)/2+center(3);
                %Single plot command for all 'cube lines'
                plot3(X1,Y1,Z1,'color','black');
                hold off;
            end
        end
        
    end
    methods (Abstract)
        %This method will split the child, creating a new PUPatch. If the
        %obj does not need to split, the method returns obj.
        Child = splitleaf(obj);
        
        ef = evalf(obj,X)
        
        ef = evalfGrid(obj,X)
    end
    
end
