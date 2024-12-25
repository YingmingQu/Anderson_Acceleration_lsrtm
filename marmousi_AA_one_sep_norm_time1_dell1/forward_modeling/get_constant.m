function [dt,dtout,ndtt,nt,ntout]=get_constant(vp0,dx,dz,vnx,nz,frequency,dt,tmax,dtout)
vpmax=max(max(vp0));
vpmin=min(min(vp0));
% determine mininum spatial sampling interval
H_min=min(dx,dz);
% determine time sampling interval to ensure stability
dt_max=0.5*H_min/vpmax;
dx_max=vpmin/frequency*0.2;
dz_max=dx_max;
if(dx_max<dx) 
    disp('YOU NEED HAVE TO REDEFINE DX!');
end
if(dz_max<dz) 
    disp('YOU NEED HAVE TO REDEFINE DZ!');
end
 dt=dt/1000.0;  
 dtout=dtout/1000.0;	   
 ndtt=dtout/dt;	 
 nt=floor(tmax/dt+0.5)+1; 
 ntout=(nt-1)/floor(ndtt+0.5)+1; 
end