%SGE potential in 2d: sigma.B;
%evolution in time
function rashba2dbailey
global x c X Y waveiniup waveinidn Rho x0 y0 alpha Psiup Psidn r
dt= 3.0;          Lt =3.0;         Nt = Lt/dt;
dx =0.2;          Lx =12.0;        Nx = Lx/dx;
dy =0.2;          Ly =12.0;        Ny = Lx/dy;
alpha=1;     m=1;       hbar=1; 
nr=50000;                ampup=1/sqrt(2);               ampdn=1/sqrt(2);
%w packet's 1width

wx=2;
wy=2;
x=(-Nx/2:(Nx/2))*dx;
y=(-Ny/2:(Ny/2))*dy; 
[X, Y] = meshgrid(x,y);
Psiup = exp(-(X.^2+Y.^2));
Rho=(abs(Psiup).^2);
plotrho
pause(0.1)

for j=1:Nt
x0=-6.0+12.0*rand(nr,1);
y0=-6.0+12.0*rand(nr,1);
    c=[];
    d=[];
for k=1:Nx+1
      for m=1:Ny+1
         r=sqrt((X(k,m)-x0).^2+(Y(k,m)-y0).^2);
         waveiniup=1/sqrt(pi*wx*wy)*ampup*exp(-((x0.^2/wx^2)+(y0.^2/wy^2)));
         waveinidn=1/sqrt(pi*wx*wy)*ampdn*exp(-((x0.^2/wx^2)+(y0.^2/wy^2)));
          c=[c,(144/nr)*sum(1/(2*pi*1i*dt*j).*exp(-((X(k,m)-x0).^2.+(Y(k,m)-y0).^2-alpha^2*(dt)^2*j)/(2*1i*dt*j)).*(cos(alpha*r).*waveiniup+1i*sin(alpha*r).*waveinidn.*(1i*(X(k,m)-x0)+(Y(k,m)-y0))./r))];
          d=[d,(144/nr)*sum(1/(2*pi*1i*dt*j).*exp(-((X(k,m)-x0).^2.+(Y(k,m)-y0).^2-alpha^2*(dt)^2*j)/(2*1i*dt*j)).*(cos(alpha*r).*waveinidn+1i*sin(alpha*r).*waveiniup.*(-1i*(X(k,m)-x0)+(Y(k,m)-y0))./r))];
      end
end
Psiup=transpose(reshape(transpose(c),Nx+1,Nx+1));
Psidn=transpose(reshape(transpose(d),Nx+1,Nx+1));
%Rho = conj(Psiup).*Psiup-conj(Psidn).*Psidn;
 %Rho = (conj(Psiup).*Psidn + conj(Psidn).*Psiup);
 Rho= 1i*(-conj(Psiup).*Psidn+conj(Psidn).*Psiup);
plotrho
pause(0.1)
end

%Rho =(abs(Psidn+Psiup)).^2;
 %A = moviein(Nx/4);
%A(k) = getframe;
%movie(A);
%legend('Evolved W Pkt', 2); %'Location', Northwest);
function plotrho
global X Y Rho  
surf(X,Y,Rho)
shading interp;
colormap jet;
set(gca,'FontSize',40)
colorbar('FontSize',28)
%axis ([-5 5 -5 5 0 10]);

axis ([-6 6 -6 6]);

%axis square;
xlabel('');  ylabel('');
topline = sprintf('(d)');
title(topline);
 %1/(2*pi*1i*dt*j).*exp(-((X(k,m)-x0).^2.+(Y(k,m)-y0).^2)/(2*1i*dt*j)-((-X(k,m)-x0+Y(k,m)+y0)*dt*j/(2*1i))).*
% 1/(2*pi*1i*dt*j).*exp(-((X(k,m)-x0).^2.+(Y(k,m)-y0).^2)/(2*1i*dt*j)-((-X(k,m)-x0+Y(k,m)+y0)*dt*j/(2*1i))).*