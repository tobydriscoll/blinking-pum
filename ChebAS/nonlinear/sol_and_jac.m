function [ sol,J ] = sol_and_jac( f,jac,u )
sol = f(u);

if nargout==2
    J = jac(u);
end

end

