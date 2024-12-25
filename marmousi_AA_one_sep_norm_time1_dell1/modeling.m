function  [left_u,right_u,top_w,bottom_w,top_p,bottom_p,left_p,right_p,p_end,u_end,w_end] = modeling(vp,rho,mm,c,source,xsn,zs_hcp,zrec,dt,tmax,npd1,npd,nx,nz,dx,dz,nt,ndtt)
   u1=zeros(nx+npd+npd1,nz+npd+npd1);
   w1=zeros(nx+npd+npd1,nz+npd+npd1);
   p1=zeros(nx+npd+npd1,nz+npd+npd1);
   left_u=zeros(nz*5,nt);
   right_u=zeros(nz*5,nt);
   top_w=zeros(nx*5,nt);
   bottom_w=zeros(nx*5,nt);
   left_p=zeros(nz*5,nt);
   right_p=zeros(nz*5,nt);
   top_p=zeros(nx*5,nt);
   bottom_p=zeros(nx*5,nt);
    p_end=zeros(nx+npd+npd1,nz+npd+npd1);
    u_end=zeros(nx+npd+npd1,nz+npd+npd1);
    w_end=zeros(nx+npd+npd1,nz+npd+npd1);
 
   for it=1:nt
       p1=update_p(dt,npd,npd1,nx,nz,u1,w1,p1,vp,dx,dz,mm,c,1);
%           fprintf('it=%d,p0=%f\n',it,p1(1+npd1,1+npd1));
       
           
           
           p1(xsn+npd1,zs_hcp+npd1)=p1(xsn+npd1,zs_hcp+npd1)+source(it)*vp(xsn+npd1,zs_hcp+npd1)*dt; 
         
      
%          fprintf('it=%d,p1=%f\n',it,p1(1+npd1,1+npd1));
       [u1,w1]=update_vel(u1,w1,p1,dx,dz,rho,nx,nz,npd,npd1,dt,mm,c,1);
%          fprintf('it=%d,u1=%f\n',it,u1(1+npd1,1+npd1));
       [u1,w1,p1]=abs_bc(u1,w1,p1,nx,nz,npd,npd1);
     
       jj=0;
       for i=(1+npd1):(nx+npd1)
        top_p(jj+1,it)=p1(i,1+npd1);
		top_p(jj+2,it)=p1(i,2+npd1);
		top_p(jj+3,it)=p1(i,3+npd1);
		top_p(jj+4,it)=p1(i,4+npd1);
		top_p(jj+5,it)=p1(i,5+npd1);
	    top_w(jj+1,it)=w1(i,1+npd1);
		top_w(jj+2,it)=w1(i,2+npd1);
		top_w(jj+3,it)=w1(i,3+npd1);
		top_w(jj+4,it)=w1(i,4+npd1);
		top_w(jj+5,it)=w1(i,5+npd1);
		jj=jj+5;
       end
       jj=0;
	 for i=(1+npd1):(nx+npd1)
	    bottom_p(jj+1,it)=p1(i,nz-0+npd1);
		bottom_p(jj+2,it)=p1(i,nz-1+npd1);
		bottom_p(jj+3,it)=p1(i,nz-2+npd1);
		bottom_p(jj+4,it)=p1(i,nz-3+npd1);
		bottom_p(jj+5,it)=p1(i,nz-4+npd1);
	    bottom_w(jj+1,it)=w1(i,nz-1+npd1);
		bottom_w(jj+2,it)=w1(i,nz-2+npd1);
		bottom_w(jj+3,it)=w1(i,nz-3+npd1);
		bottom_w(jj+4,it)=w1(i,nz-4+npd1);
		bottom_w(jj+5,it)=w1(i,nz-5+npd1);
		jj=jj+5;
     end
	  jj=0;
	  for j=(1+npd1):(nz+npd1)
	    left_p(jj+1,it)=p1(1+npd1,j);
		left_p(jj+2,it)=p1(2+npd1,j);
		left_p(jj+3,it)=p1(3+npd1,j);
		left_p(jj+4,it)=p1(4+npd1,j);
		left_p(jj+5,it)=p1(5+npd1,j);
	    left_u(jj+1,it)=u1(1+npd1,j);
		left_u(jj+2,it)=u1(2+npd1,j);
		left_u(jj+3,it)=u1(3+npd1,j);
		left_u(jj+4,it)=u1(4+npd1,j);
		left_u(jj+5,it)=u1(5+npd1,j);
		jj=jj+5;
      end
	  jj=0;
	  for j=(1+npd1):(nz+npd1)
	    right_p(jj+1,it)=p1(nx-0+npd1,j);
		right_p(jj+2,it)=p1(nx-1+npd1,j);
		right_p(jj+3,it)=p1(nx-2+npd1,j);
		right_p(jj+4,it)=p1(nx-3+npd1,j);
		right_p(jj+5,it)=p1(nx-4+npd1,j);
	    right_u(jj+1,it)=u1(nx-1+npd1,j);
		right_u(jj+2,it)=u1(nx-2+npd1,j);
		right_u(jj+3,it)=u1(nx-3+npd1,j);
		right_u(jj+4,it)=u1(nx-4+npd1,j);
		right_u(jj+5,it)=u1(nx-5+npd1,j);
		jj=jj+5;
      end
 
      
	  if(it==nt)
	    p_end(:,:)=p1(:,:);
	    u_end(:,:)=u1(:,:);
	    w_end(:,:)=w1(:,:); 
      end
     
   end
 end