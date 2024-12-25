function [u,w,p]=abs_bc(u,w,p,nx,nz,npd,npd1)
 qp=-0.15;
 %top boundary
 for j=0:-1:-npd+1
     for i=1:nx
           u(i+npd1,j+npd1)=(qp*((1-j)/(1.0*npd))^2+1)*u(i+npd1,j+npd1);
           w(i+npd1,j+npd1)=(qp*((1-j)/(1.0*npd))^2+1)*w(i+npd1,j+npd1);
           p(i+npd1,j+npd1)=(qp*((1-j)/(1.0*npd))^2+1)*p(i+npd1,j+npd1);
     end
 end
 %bottom boundary
   for j=nz+1:nz+npd
           for i=1:nx
              u(i+npd1,j+npd1)=(qp*((j-nz)/(1.0*npd))^2+1)*u(i+npd1,j+npd1);
              w(i+npd1,j+npd1)=(qp*((j-nz)/(1.0*npd))^2+1)*w(i+npd1,j+npd1);
              p(i+npd1,j+npd1)=(qp*((j-nz)/(1.0*npd))^2+1)*p(i+npd1,j+npd1);
           end
   end
 %left boundary
for i=0:-1:-npd+1
           for j=-npd+1:nz+npd
              u(i+npd1,j+npd1)=(qp*((1-i)/(1.0*npd))^2+1)*u(i+npd1,j+npd1);
              w(i+npd1,j+npd1)=(qp*((1-i)/(1.0*npd))^2+1)*w(i+npd1,j+npd1);
              p(i+npd1,j+npd1)=(qp*((1-i)/(1.0*npd))^2+1)*p(i+npd1,j+npd1);

           end
end
 %right boundary
    for i=nx+1:nx+npd
          for j=-npd+1:nz+npd
              u(i+npd1,j+npd1)=(qp*((i-nx)/(1.0*npd))^2+1)*u(i+npd1,j+npd1);
              w(i+npd1,j+npd1)=(qp*((i-nx)/(1.0*npd))^2+1)*w(i+npd1,j+npd1);
              p(i+npd1,j+npd1)=(qp*((i-nx)/(1.0*npd))^2+1)*p(i+npd1,j+npd1);

          end
    end
end