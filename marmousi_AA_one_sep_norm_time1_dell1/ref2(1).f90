!**************************************************************************************
program ref
implicit none
real,allocatable:: vp0(:,:),vs0(:,:),ref0(:,:) 
integer i,j,vnx,nz
character*256 FN1,FN2,FN3
FN1='aoxian100_100.dat'
FN2='vpsmooth.dat'
FN3='ref.dat' 

vnx=100
nz=100	 
allocate(vp0(1:vnx,1:nz))
allocate(vs0(1:vnx,1:nz))
allocate(ref0(1:vnx,1:nz))	
vp0=0
vs0=0
ref0=0 
call read_file(FN1,vnx,nz,vp0)
call read_file(FN2,vnx,nz,vs0)
do i=1,vnx
  do j=1,nz
    ref0(i,j)=1.0/vp0(i,j)/vp0(i,j)-1.0/vs0(i,j)/vs0(i,j)
  enddo
enddo 	
open(4,file=FN3,access='direct',recl=4*nz)
do i=1,vnx
  write(4,rec=i)(ref0(i,j),j=1,nz)
enddo
close(4)     
end
!******************************************************************************************************	
      subroutine read_file(FN1,nx,nz,vv)

        integer nx,nz,i,j
        dimension vv(1:nx,1:nz)
        character*256 FN1  
        open(2,file=FN1,access='direct',recl=4*nz)
        do i=1,nx
           read(2,rec=i) (vv(i,j),j=1,nz)
        enddo	
        close(2)
        return
      end
!*********************************************************


      
