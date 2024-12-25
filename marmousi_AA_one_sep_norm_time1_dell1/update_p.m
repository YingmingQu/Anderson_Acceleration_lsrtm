function [p]=update_p(dt,npd,npd1,nx,nz,u,w,p,vp,dx,dz,mm,c,lee)
dtx=dt/dx;
dtz=dt/dz;
for j=(1-npd+mm+npd1):(nz+npd-mm+npd1)
    for i=(1-npd+mm+npd1):(nx+npd-mm+npd1)
        dux=0.0;
        dwz=0.0;
        for ii=1:mm
            dux=dux+c(mm,ii)*(u(i+ii-1,j)-u(i-ii,j));
            dwz=dwz+c(mm,ii)*(w(i,j+ii-1)-w(i,j-ii));
        end
        p(i,j)=p(i,j)+(dtx*dux+dtz*dwz)*vp(i,j)*lee; 
    end
end
end