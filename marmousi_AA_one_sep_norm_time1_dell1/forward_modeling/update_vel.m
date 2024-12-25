function [u,w]=update_vel(u,w,p,dx,dz,rho,nx,nz,npd,npd1,dt,mm,c,lee)
dtx=dt/dx;
dtz=dt/dz;
for j=(1-npd+mm+npd1):(nz+npd-mm+npd1)
    for i=(1-npd+mm+npd1):(nx+npd-mm+npd1)
        dpx=0.0;
        dpz=0.0;
         for ii=1:mm
              dpx=dpx+c(mm,ii)*(p(i+ii,j)-p(i-ii+1,j));
	          dpz=dpz+c(mm,ii)*(p(i,j+ii)-p(i,j-ii+1));
         end
           u(i,j)=u(i,j)+dtx*dpx*rho(i,j)*lee;
           w(i,j)=w(i,j)+dtz*dpz*rho(i,j)*lee;
    end
end
end