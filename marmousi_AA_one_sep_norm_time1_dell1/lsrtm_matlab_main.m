         
clc; clear; close all;
vnx=368;
nx=368;
nz=200;
dx=10.0;
dz=10.0;
ns_hcp=19;
ds_hcp=20;
fs_hcp=1;
zs_hcp=1;
dt=0.5;
dtout=0.5;
tmax=3.0;
npd=50;
npd1=50;
frequency=20.0;
ns_start=1;
mm=5;
ntrace=1;
stype=2; 	
wtype=1;
zrec=2;
pi=4.0*atan(1.0);  
pfac=1.0;
nit=30;
mem_size=10;
beta = 1;



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


    fid=fopen('shot_born.dat','rb');  
     [shot,count1]=fread(fid,[nt,nx*ns_hcp],'float');
  
     shot=shot';
     fclose(fid);

shot1=zeros(ns_hcp,nx,nt);
for k=1:ns_hcp
    for i=1:nx
        for j=1:nt
         shot1(k,i,j)=shot(i+(k-1)*nx,j);   
        end
    end
end
% shot2(1:nx,1:nt)=shot1(1,1:nx,1:nt);
% imagesc(shot2');
tic;
% add source 
ts=2.0/frequency; 
nts=floor(ts/dt); 
source=get_source(nt,dt,ts,nts,frequency,pi,pfac);
ref0=zeros(nx,nz);



xc=zeros(nx*nz,1);
Fx1=zeros(nx*nz,1);
Fx0=zeros(nx*nz,1);
g1=zeros(nx*nz,1);
g0=zeros(nx*nz,1);
  Smem = [];
    Ymem = [];
    
for iter=1:nit
    res0=0;
    energy_all=zeros(nx,nz);
	grad0_all=zeros(nx,nz);
    energy=zeros(nx,nz);
	grad0=zeros(nx,nz);
    
       if iter>1
        if iter == 2
            x1 = Fx0;
        else
            x1 = x0 - g0 - (Smem - Ymem) * ((Smem'*Ymem) \ (Smem' * g0));
        end
        ref0=reshape(x1,nx,nz);
       end
    
     parfor is=1:ns_hcp
   
    [xsn,vp,rho,ref]=current_shot(vp0,rho0,ref0,nx,nz,npd,npd1,vnx,fs_hcp,ds_hcp,is,stype);
    vp=pad_vv(nx,nz,npd,npd1,vp);
    rho=pad_vv(nx,nz,npd,npd1,rho);
    vp=rho.*vp.^2;
    rho=1.0./rho;
    obs=zeros(nx,nt);
    cal=zeros(nx,nt);
    obs(1:nx,1:nt)=shot1(is,1:nx,1:nt);
    if(iter==1)
       deta=-obs; 
 [left_u,right_u,top_w,bottom_w,top_p,bottom_p,left_p,right_p,p_end,u_end,w_end] = modeling(vp,rho,mm,c,source,xsn,zs_hcp,zrec,dt,tmax,npd1,npd,nx,nz,dx,dz,nt,ndtt);
         
    else
 [cal,left_u,right_u,top_w,bottom_w,top_p,bottom_p,left_p,right_p,p_end,u_end,w_end] = linear_modeling(vp,rho,mm,c,cal,ref,source,xsn,zs_hcp,zrec,dt,tmax,npd1,npd,nx,nz,dx,dz,nt,ndtt);
	  deta=cal-obs;
    end
    for i=1:nx
	  for it=1:ntout
	    res0=res0+deta(i,it)^2;
      end
    end
    
     [energy,grad0] =rtm(vp,rho,mm,c,deta,source,xsn,zs_hcp,zrec,dt,tmax,npd1,npd,nx,nz,dx,dz,nt,ndtt,left_u,right_u,top_w,bottom_w,top_p,bottom_p,left_p,right_p,p_end,u_end,w_end);	
%     fprintf('grad=%f\n',grad0(100,100));
     energy_all=energy_all+energy;
    grad0_all= grad0_all+grad0;
    end
    
     disp(iter);
    disp(res0);
       
    fid=fopen('residual.txt','a');
a='iter=,res0=';
fprintf(fid,'%s %f %f\n',a,iter*1.0,res0);
fclose(fid);


%        if(iter==1)
       dk=grad0_all;
    for i=1:nx
	 for j=1:nz
	    if(energy_all(i,j)~=0) 
            dk(i,j)=grad0_all(i,j)/energy_all(i,j);
        end
     end
    end


  fid=fopen('dk_norm.txt','a');
    fprintf(fid,'iter=%d %.6e\n',iter,norm(dk));
 fclose(fid);


%        else
%            
%            dk=grad0_all;
%     for i=1:nx
% 	 for j=1:nz
% 	    if(energy_all(i,j)~=0) 
%             dk(i,j)=grad0_all(i,j)/energy_all(i,j);
%         end
%      end
%     end
% 
%         gg=0;
%         dgg=0;
%      for i=1:nx
% 	 for j=1:nz
%             gg=gg+grad0_all(i,j)*dk(i,j);
%             dgg=dgg+grad1(i,j)*dk1(i,j);
%         end
%      end
%     beta=gg/dgg;
%     dk=dk+beta.*dk1;



%        end
fid1=fopen(['gradient',num2str(iter),'.dat'],'wb');
fwrite(fid1, dk','float');
fclose(fid1);
    
       
       
dgg_all=0;
  if(iter==1)
 parfor is=1:ns_hcp

    cal=zeros(nx,nt);
    [xsn,vp,rho,dk2]=current_shot(vp0,rho0,dk,nx,nz,npd,npd1,vnx,fs_hcp,ds_hcp,is,stype);
    vp=pad_vv(nx,nz,npd,npd1,vp);
    rho=pad_vv(nx,nz,npd,npd1,rho);
    vp=rho.*vp.^2;
    rho=1.0./rho;
[cal,left_u,right_u,top_w,bottom_w,top_p,bottom_p,left_p,right_p,p_end,u_end,w_end] = linear_modeling(vp,rho,mm,c,cal,dk2,source,xsn,zs_hcp,zrec,dt,tmax,npd1,npd,nx,nz,dx,dz,nt,ndtt);
for i=1:nx
	  for it=1:ntout
	    dgg_all=dgg_all+cal(i,it)^2;
      end
   end 
   
   
end


    gg=0;
	for i=1:nx
	  for j=1:nz
	    gg=gg+grad0_all(i,j)*dk(i,j);
      end
    end
	alfa=gg/dgg_all;
    
end
     fid=fopen('alfa.txt','a'); 
a='iter=,alfa=';
fprintf(fid,'%s %f %.6e\n',a,iter*1.0,alfa);
fclose(fid);


   
    ref0(:,:)=ref0(:,:)-alfa*dk(:,:);
    
    
    if iter==1
       Fx0=reshape(ref0,nx*nz,1) ;
    else
        Fx1=reshape(ref0,nx*nz,1) ;
    end
    
% Apply Anderson acceleration.     
% [xc,Fx1,Fx0,g1,g0,Q,R,DF,mAA]=aaqr(nx,nz,dx,dz,iter,mem_size,dk_reshape,alfa,xc,Fx1,Fx0,g1,g0,Q,R,DF,mAA,droptol,beta);
  if iter>1
     if iter-1 <= mem_size - 1
            Smem(:, iter-1) = x1 - x0;
            Ymem(:, iter-1) = x1 - Fx1 - g0;
        else
            Smem = [Smem(:, 2:end), x1 - x0];
            Ymem = [Ymem(:, 2:end), x1 - Fx1 - g0];
     end
    end    
    
    if(iter==1)
        x0=zeros(nx*nz,1);
    g0 = x0 - Fx0;
    else
       x0 = x1;
        g0 = x0 - Fx1;    
    end
% ref0=reshape(xc,nx,nz);


  fid=fopen('ref0_norm.txt','a');
    fprintf(fid,'iter=%d %.6e\n',iter,norm(ref0));
 fclose(fid);


misfit=0;
    for i=1:nx
	  for j=1:nz
	    misfit=misfit+(ref0(i,j)-ref_true(i,j))^2;
      end
    end
    
     fid5=fopen('misfit.txt','a');
    fprintf(fid5,'iter=%d %.6e\n',iter,misfit);
 fclose(fid5);
fid2=fopen(['mig_lsrtm',num2str(iter),'.dat'],'wb');
fwrite(fid2, ref0','float');
fclose(fid2);
toc;
    fid4=fopen('time.txt','a');
    fprintf(fid4,'iter=%d %f\n',iter,toc);
 fclose(fid4); 
end






   
  

   
 
   


   




