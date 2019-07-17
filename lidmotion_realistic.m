function [rho,drho_dt,dtbc,dtco,dto,dtoc] = lidmotion_realistic(pcl)

% pcl = percent closed at minimum opening (t=0)

rho = @X;
drho_dt = @Xt;

% adjust for 5 mm/sec for Longfei vs 10 cm/sec for Alfa
s = 20;

dtco= 3.52/s; dto = 100/s; dtoc = 1.64/s;
dtbc = dtco + dto + dtoc;  % period of motion
L = 1;  pcl = 1-pcl;

    function x = X(t)
        
        x = zeros(size(t));
        tt = rem(t,dtbc);
        i1 = (tt >= 0) & (tt < dtco);
        x(i1) = -L*( 1-2*pcl-2*(1-pcl)*(tt(i1)/dtco).^2.*exp(1-(tt(i1)/dtco).^2) );
        i2 = (tt >= dtco) & (tt < dtco + dto);
        x(i2) = L;
        i3 = (tt >= dtco + dto) & (tt <= dtbc);
        x(i3) = -L*( -1+2*(1-pcl)*((tt(i3)-dtco-dto)/dtoc).^2.*exp(1-((tt(i3)-dtco-dto)/dtoc).^2) );
    end


    function xt = Xt(t)
        
        xt = zeros(size(t));
        tt = rem(t,dtbc);
        i1 = (tt >= 0) & (tt < dtco);
        xt(i1) = 2*L*(1-pcl)*(2/dtco)*(tt(i1)/dtco).*exp(1-(tt(i1)/dtco).^2).*(1-(tt(i1)/dtco).^2);
        i2 = (tt >= dtco) & (tt < dtco + dto);
        xt(i2) = 0;
        i3 = (tt >= dtco + dto) & (tt <= dtbc);
        a = (tt(i3)-dtco-dto)/dtoc;
        xt(i3) = -2*L*(1-pcl)*(2/dtoc)*a.*exp(1-a.^2).*(1-a.^2);
        
    end

end