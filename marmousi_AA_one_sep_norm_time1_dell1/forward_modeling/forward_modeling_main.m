% ����                  
% �Ƴ����д��2021��4��20��,�����              
clc; clear; close all;
vnx=368;
nx=368;
nz=200;
dx=10.0;
dz=10.0;
ns_hcp=19;%19
ds_hcp=20;
fs_hcp=1;
zs_hcp=1;
dt=0.5;
dtout=0.5;
tmax=3.0;
npd=50;
npd1=50;
frequency=20.0;
ns_start=1;%?
mm=5;
ntrace=1;
stype=2; 	
wtype=1;
zrec=2;
pi=4.0*atan(1.0);  
pfac=1.0;


c=cal_c(mm);

fid=fopen('vpsmooth.dat','rb');  
      [vp0,count]=fread(fid,[nz,nx],'float');
    vp0=vp0';
      fclose(fid);
      
       for i=1:nx
         for j=1:nz
             rho0(i,j)=1.0;
         end
     end
      
fid=fopen('ref.dat','rb');  
     [ref_true,count1]=fread(fid,[nz,nx],'float');
  
     ref_true=ref_true';
     fclose(fid);

[dt,dtout,ndtt,nt,ntout]=get_constant(vp0,dx,dz,vnx,nz,frequency,dt,tmax,dtout);


tic
% add source 
ts=2.0/frequency; 
nts=floor(ts/dt); 
source=get_source(nt,dt,ts,nts,frequency,pi,pfac);

ref0=zeros(nx,nz);

cal1=zeros(ns_hcp,nx,nt);
parfor is=1:ns_hcp
 cal=zeros(nx,nt);
    res0=0;
    [xsn,vp,rho]=current_shot1(vp0,rho0,nx,nz,npd,npd1,vnx,fs_hcp,ds_hcp,is,stype);

    vp=pad_vv(nx,nz,npd,npd1,vp);

    rho=pad_vv(nx,nz,npd,npd1,rho);

    vp=rho.*vp.^2;
    rho=1.0./rho;
    
 [cal] = linear_modeling(vp,rho,mm,c,cal,ref_true,source,xsn,zs_hcp,zrec,dt,tmax,npd1,npd,nx,nz,dx,dz,nt,ndtt);
    
 cal1(is,:,:)=cal(:,:);


% 
end
cal2=zeros(nt,ns_hcp*nx);
for k=1:ns_hcp
   for j=1:nt
	    for i=1:nx
            cal2(j,i+(k-1)*nx)=cal1(k,i,j);
        end
   end
end


fid=fopen('shot_born.dat','wb');
fwrite(fid,cal2,'float');
fclose(fid);

toc









