classdef blinkmultilog < blinkmulti
    
	% variant of blinkmulti that solves for log(h) in place of h
	
    methods
		function b = initialize(b,lid,flux,state)
			b = b.initialize@blinkmulti(lid,flux,state);
			
			Hleaf = b.H.leafArray;
			for i = 1:length(Hleaf)
				b.region{i} = blinkonelog(b.model,Hleaf{i},b.map,lid,flux);
			end			

			b.initstate.dH = b.initstate.dH ./ b.initstate.H;
			b.initstate.H = log(b.initstate.H);
		end

        function H = evalH(r,t)
            V = r.evalH@blinkmulti(t);
            H = exp( V );
        end
        
    end
    
 end