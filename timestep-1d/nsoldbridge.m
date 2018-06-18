function fun = nsoldbridge(res,jac)

fun = @nfun;

    function [r,J] = nfun(x)
        r = res(x);  J = jac(x);
    end

end
