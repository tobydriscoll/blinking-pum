classdef blinkmultilog < blinkmulti
    
	% variant of blinkmulti that solves for log(h) in place of h
	
    methods
		function b = initialize(b,lid,flux,state)			
			Hleaf = b.H.leafArray;
			for i = 1:length(Hleaf)
				b.region{i} = blinkonelog(b.model,Hleaf{i},b.map,lid,flux);
			end			

			% use constant initial state if none given
			if isempty(state)
				state.H = chebfun2(b.model.h_boundary);
				state.P = chebfun2(0);
				state.dH = chebfun2(0);
				state.dP = chebfun2(0);
            end	
            
			state.dH = state.dH ./ state.H;
			state.H = log(state.H);
			b.initstate = state;
		end

        function H = evalH(r,t)
            V = r.evalH@blinkmulti(t);
            H = exp( V );
        end
        
    end
    
 end