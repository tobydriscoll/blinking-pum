classdef blinkonelog < blinkone
    
	% variant of blinkone that solves for log(h) in place of h, by
	% interposing the timederiv and jacobian functions
	    
    methods %(Access=private)
              
        function [w_t,H_t] = timederiv(r,t,w)
            % Compute f in M*u_t = f(t,u)
            [V,P] = unpack(r,w);
            V(r.disc.boundary.loc_outer) = log(r.boundaryH);
            H = exp(V);
            u_t = r.timederiv@blinkone(t,r.pack(H,P));
            [H_t,P_t] = unpack(r,u_t);
			H_t(r.disc.boundary.loc_outer) = 0;
            V_t = H_t./H;
            w_t = r.pack(V_t,P_t);
        end
        
        function J = jac(r,t,w)
            idxh = 1:r.disc.num.h; 
            idxp = r.disc.num.h+1:length(w);
            v = w(idxh);
            h = exp(v);
            u = [h;w(idxp)];
            J = r.jac@blinkone(t,u);
            [~,f] = r.timederiv(t,w);
            J(:,idxh) = J(:,idxh).*h';
            J(idxh,:) = (1./h).*J(idxh,:);
            J(idxh,idxh) = J(idxh,idxh) - diag(f(~r.disc.boundary.loc_outer)./h);
        end
                
    end
end