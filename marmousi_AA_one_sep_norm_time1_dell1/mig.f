      program make inclined fracture velocity field
	parameter(nx=368,nz=6001)
!	parameter(pi=3.1415926)
	dimension vp(1:nx,1:nz)
	dimension vp1(1:nx,1:nz)

      
       integer iter 
	character*256 fn1,fn2

	fn1='shot_born.dat'
	fn2='shot_born10.dat'
	
      

ccccccccccccccwriteccccccccccccccccccccccccccccccccccccc
	open(1,file=fn1,access='direct',recl=4*nz)
	open(2,file=fn2,access='direct',recl=4*nz)
       iter=10
        
	do i=1,nx
        
	   read(1,rec=(iter-1)*nx+i)(vp(i,j),j=1,nz)
	enddo
    
      
    
    
	


	do i=1,nx
	   write(2,rec=i)(vp(i,j),j=1,nz)
	enddo

	  close(1)
	  close(2)
	  end
	     


