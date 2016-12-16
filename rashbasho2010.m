%SGE potential in 2d: sigma.B;
%evolution in time
function rashba2dbailey
global x c X Y Rho x0 y0 alpha Psiup Psidn 
dt=3;      Lt=3;       Nt = Lt/dt;
dx =1/5;      Lx =16.0;        Nx = Lx/dx;
dy =1/5;      Ly =16.0;        Ny = Ly/dx;
alpha=0.5;    nr=1000;        
ampup=1;  
ampdn=0; 
w=0.5;
wx=2.0;wy=1.0;
sigma =1.0;           %w packet's 1width
x=(-Nx/2:(Nx/2))*dx;
y=(-Ny/2:(Ny/2))*dy;
[X, Y] = meshgrid(x,y);
Psiup = exp(-(X.^2+Y.^2));
Rho=(abs(Psiup).^2);
plotrho
pause(0.1)
for j=1:Nt
x0=-8+16*rand(nr,1);
y0=-8+16*rand(nr,1);
    c=[];
    d=[];
for k=1:Nx+1
      for m=1:Ny+1
         c2=w*j*dt*(cot(w*j*dt)+csc(w*j*dt)).*(X(k,m)-x0);
         c1=c2-2*(X(k,m)-x0)+2*(Y(k,m)-y0);    
         r=sqrt(c1.^2+c2.^2);
         waveup=1/sqrt(pi*wx*wy)*ampup*exp(-(x0.^2/wx^2+y0.^2/wy^2));
         wavedn=1/sqrt(pi*wx*wy)*ampdn*exp(-(x0.^2/wx^2+y0.^2/wy^2));
         c=[c,(256/nr)*sum(sqrt(w)/(2*pi*1i*sqrt(j*dt*sin(w*j*dt))).*exp(0.5*1i*((X(k,m).^2+x0.^2+alpha^2.*(j*dt).^2).*w.*cot(w*j*dt)-(2*X(k,m).*x0+alpha^2.*(j*dt)^2).*w.*csc(w*(j*dt))+alpha^2*(j*dt)^2)+0.5*1i/(j*dt)*((Y(k,m)-y0).^2+alpha^2*(j*dt)^2)).*(cos(0.5*alpha.*r).*waveup+1i*(c1+1i*c2)./r.*sin(0.5*alpha.*r).*wavedn))];        
         d=[d,(256/nr)*sum(sqrt(w)/(2*pi*1i*sqrt(j*dt*sin(w*j*dt))).*exp(0.5*1i*((X(k,m).^2+x0.^2+alpha^2.*(j*dt).^2).*w.*cot(w*j*dt)-(2*X(k,m).*x0+alpha^2.*(j*dt)^2).*w.*csc(w*(j*dt))+alpha^2*(j*dt)^2)+0.5*1i/(j*dt)*((Y(k,m)-y0).^2+alpha^2*(j*dt)^2)).*(cos(0.5*alpha.*r).*wavedn+1i*(c1-1i*c2)./r.*sin(0.5*alpha.*r).*waveup))];        
      end
end
Psiup=transpose(reshape(transpose(c),Nx+1,Nx+1));
Psidn=transpose(reshape(transpose(d),Nx+1,Nx+1));
 Rho = conj(Psiup).*Psiup-conj(Psidn).*Psidn;
 %Rho = (conj(Psiup).*Psidn + conj(Psidn).*Psiup);
 %Rho= 1i*(-conj(Psiup).*Psidn+conj(Psidn).*Psiup);
 
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
axis ([-8.0 8.0 -8.0 8.0]);

%axis square;
xlabel('');  ylabel('');
topline = sprintf('(b)');
title(topline);
  