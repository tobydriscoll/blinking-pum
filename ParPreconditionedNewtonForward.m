function [z,J] = ParPreconditionedNewtonForward(PUApprox,sol,evalF,num_sols)

%PUApprox.sample(sol);

step = zeros(length(PUApprox.leafArray),1);

for k=2:length(PUApprox.leafArray)
    step(k) = step(k-1) + num_sols*length(PUApprox.leafArray{k-1});
end

for k=1:length(PUApprox.leafArray)
    
    degs = PUApprox.leafArray{k}.degs;
    
    sol_loc{k} = [];
    diff{k} = [];
    
    [~,~,in_border{k},~] = FindBoundaryIndex2DSides(degs,PUApprox.leafArray{k}.domain,PUApprox.leafArray{k}.outerbox);
    
    loc_sol_step = 0;
    sol_step = 0;
    
    for j=1:num_sols
        tmp = sol(loc_sol_step+step(k)+(1:prod(degs)));
        sol_loc{k} = [sol_loc{k};tmp];
        diff{k} = [diff{k};PUApprox.leafArray{k}.Binterp*sol(sol_step+(1:length(PUApprox)))];
        loc_sol_step = loc_sol_step + prod(degs);
        sol_step = sol_step + length(PUApprox);
    end
    
end

%parallel step
for k=1:length(PUApprox.leafArray)
    
    [z{k},J{k}] = local_inverse(PUApprox.leafArray{k},sol_loc{k},in_border{k},diff{k},evalF,num_sols,length(PUApprox.leafArray{k}));
    
end

z = cell2mat(z');

end

function [c,Jk] = local_inverse(approx,sol_k,border_k,diff_k,evalF,num_sols,sol_length)

    function [F,J] = residual(z)
        
        z_orig = z;
        z = z+sol_k;
        
        [F,J] = evalF(approx,z);
        
        sol_step = 0;
        
        diff_len = length(diff_k)/num_sols;
        diff_step = 0;
        for i=1:num_sols
            F_s{i} = F(sol_step+(1:sol_length));
            F_s{i}(border_k) = F_s{i}(border_k) - diff_k(diff_step+(1:diff_len));
            sol_step = sol_step+sol_length;
            diff_step = diff_step + diff_len;
        end
        
        F = cell2mat(F_s');
        
    end

    options = optimoptions(@fsolve,'SpecifyObjectiveGradient',true,'MaxIterations',10000,'FunctionTolerance',1e-4);

    [s,~,~,~,Jk] = fsolve(@residual,zeros(size(sol_k)),options);
    
    c = s(:,end);
end


