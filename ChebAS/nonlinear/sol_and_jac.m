function [ sol,J ] = sol_and_jac( f,jac,u )
sol = f(u);
J = jac(u);
end

