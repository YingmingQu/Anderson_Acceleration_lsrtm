function [xsn,vp,rho,ref]=current_shot(vp0,rho0,ref0,nx,nz,npd,npd1,vnx,fs_hcp,ds_hcp,is,stype)
if(stype==1)
  ivstart=1+(is-1)*ds_hcp;
   ivend=nx+(is-1)*ds_hcp;
   xsn=fs_hcp;
else
    ivstart=1;
	ivend=nx;
    xsn=fs_hcp+(is-1)*ds_hcp;
end
if(ivstart<=0)
    disp('ivstart less than zero');
end
if(ivend>vnx)
    disp('ivend great than vnx');
end
for ix=ivstart:ivend
    for iz=1:nz
      vp(ix-ivstart+1+npd1,iz+npd1)=vp0(ix,iz);
      rho(ix-ivstart+1+npd1,iz+npd1)=rho0(ix,iz);
      ref(ix-ivstart+1+npd1,iz+npd1)=ref0(ix,iz);
    end
end 
end