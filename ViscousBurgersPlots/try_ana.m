D=0.1;
length_domain=2*pi;

Ntime=20;
Nx0=101;
dx0=length_domain/(Nx0-1);
dt=0.1;

x0_valueana=zeros(Nx0,1);
solution_ana=zeros(Nx0,Ntime);


for c=1:Ntime
    t_value=c*dt;
    for r=1:Nx0
        x0_value=(r-1)*dx0;
        x0_valueana(r,1)=x0_value;
        top_integrand=@(x0,x1,t) (x0-x1).*exp(-((x0-x1).*(x0-x1))/(4*D*t)-(1/(2*D)).*(1-cos(x1)));
        bottom_integrand=@(x0,x1,t) t.*exp(-((x0-x1).*(x0-x1))/(4*D*t)-(1/(2*D)).*(1-cos(x1)));

        top=integral(@(x1)top_integrand(x0_value,x1,t_value),-inf,inf);
        bottom=integral(@(x1)bottom_integrand(x0_value,x1,t_value),-inf,inf);
        solution_ana(r,c)=top/bottom;
    end
end

clearvars bottom bottom_integrand top top_integrand c D dt length_domain Ntime Nx0 r t_value x x0_value dx0