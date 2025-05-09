% Clear workspace variables and close all figures
clc; clear; close all;

% 
vnx=368;              % Number of grid points in x-direction for velocity model
nx=368;               % Number of grid points in x-direction for computation
nz=200;               % Number of grid points in z-direction
dx=10.0;              % Grid spacing in x-direction (m)
dz=10.0;              % Grid spacing in z-direction (m)
ns_hcp=19;            % Number of shots
ds_hcp=20;            % Shot spacing (grid points)
fs_hcp=1;             % First shot location (grid points)
zs_hcp=1;           % Source depth (grid points)
dt=0.5;             % Time step (ms)
dtout=0.5;          % Output time interval (ms)
tmax=3.0;           % Total simulation time (s)
npd=50;            % Absorbing boundary thickness (grid points)
npd1=50;           % Another absorbing boundary thickness (grid points)
frequency=20.0;      % Source dominant frequency (Hz)
ns_start=1;          % Starting shot number
mm=5;              % Staggered grid parameter
ntrace=1;            % Number of traces
stype=2;             % Source type
wtype=1;            % Wavefield type
zrec=2;              % Receiver depth (grid points)
pi=4.0*atan(1.0);      % pi
pfac=1.0;            % Source factor
nit=30;              % Number of iterations
mem_size=10;        % Memory size for storing history information
beta = 1;             % Inversion parameter

% Calculate staggered grid coefficients
c=cal_c(mm);

% Read velocity model
fid=fopen('vpsmooth.dat','rb');  
[vp0,count]=fread(fid,[nz,nx],'float');
vp0=vp0';  % Transpose to match MATLAB's row-column representation
fclose(fid);

% Initialize density model
for i=1:nx
    for j=1:nz
        rho0(i,j)=1.0;
    end
end

% Read true reflectivity model
fid=fopen('ref.dat','rb');  
[ref_true,count1]=fread(fid,[nz,nx],'float');
ref_true=ref_true';
fclose(fid);

% Calculate time-related parameters
[dt,dtout,ndtt,nt,ntout]=get_constant(vp0,dx,dz,vnx,nz,frequency,dt,tmax,dtout);

% Read observed data
fid=fopen('shot_born.dat','rb');  
[shot,count1]=fread(fid,[nt,nx*ns_hcp],'float');
shot=shot';
fclose(fid);

% Reorganize shot gather data
shot1=zeros(ns_hcp,nx,nt);
for k=1:ns_hcp
    for i=1:nx
        for j=1:nt
            shot1(k,i,j)=shot(i+(k-1)*nx,j);   
        end
    end
end

% Start timing
tic;

% Generate source function
ts=2.0/frequency; 
nts=floor(ts/dt); 
source=get_source(nt,dt,ts,nts,frequency,pi,pfac);

% Initialize reflectivity model
ref0=zeros(nx,nz);

% Initialize variables
xc=zeros(nx*nz,1);
Fx1=zeros(nx*nz,1);
Fx0=zeros(nx*nz,1);
g1=zeros(nx*nz,1);
g0=zeros(nx*nz,1);
Smem = [];  % Store historical model changes
Ymem = [];  % Store historical gradient changes

% Iterative inversion starts
for iter=1:nit
    res0=0;  % Initialize residual
    energy_all=zeros(nx,nz);  % Store energy
    grad0_all=zeros(nx,nz);   % Store gradient
    energy=zeros(nx,nz);
    grad0=zeros(nx,nz);
    
    % Anderson acceleration algorithm (to optimize inversion process)
    if iter>1
        if iter == 2
            x1 = Fx0;
        else
            x1 = x0 - g0 - (Smem - Ymem) * ((Smem'*Ymem) \ (Smem' * g0));
        end
        ref0=reshape(x1,nx,nz);
    end
    
    % Process each shot in parallel
    parfor is=1:ns_hcp
        % Prepare model and parameters for current shot
        [xsn,vp,rho,ref]=current_shot(vp0,rho0,ref0,nx,nz,npd,npd1,vnx,fs_hcp,ds_hcp,is,stype);
        vp=pad_vv(nx,nz,npd,npd1,vp);
        rho=pad_vv(nx,nz,npd,npd1,rho);
        vp=rho.*vp.^2;
        rho=1.0./rho;
        
        % Initialize observed and calculated data
        obs=zeros(nx,nt);
        cal=zeros(nx,nt);
        obs(1:nx,1:nt)=shot1(is,1:nx,1:nt);
        
        % Forward modeling for the first iteration
        if(iter==1)
            deta=-obs; 
            [left_u,right_u,top_w,bottom_w,top_p,bottom_p,left_p,right_p,p_end,u_end,w_end] = modeling(vp,rho,mm,c,source,xsn,zs_hcp,zrec,dt,tmax,npd1,npd,nx,nz,dx,dz,nt,ndtt);
        % Linearized modeling for subsequent iterations
        else
            [cal,left_u,right_u,top_w,bottom_w,top_p,bottom_p,left_p,right_p,p_end,u_end,w_end] = linear_modeling(vp,rho,mm,c,cal,ref,source,xsn,zs_hcp,zrec,dt,tmax,npd1,npd,nx,nz,dx,dz,nt,ndtt);
            deta=cal-obs;
        end
        
        % Calculate residual energy
        for i=1:nx
            for it=1:ntout
                res0=res0+deta(i,it)^2;
            end
        end
        
        % Reverse time migration to compute gradient and energy
        [energy,grad0] = rtm(vp,rho,mm,c,deta,source,xsn,zs_hcp,zrec,dt,tmax,npd1,npd,nx,nz,dx,dz,nt,ndtt,left_u,right_u,top_w,bottom_w,top_p,bottom_p,left_p,right_p,p_end,u_end,w_end);	
        
        % Accumulate energy and gradient from all shots
        energy_all=energy_all+energy;
        grad0_all=grad0_all+grad0;
    end
    
    % Display iteration number and residual
    disp(iter);
    disp(res0);
    
    % Write residual to file
    fid=fopen('residual.txt','a');
    a='iter=,res0=';
    fprintf(fid,'%s %f %f\n',a,iter*1.0,res0);
    fclose(fid);
    
    % Normalize gradient by energy
    dk=grad0_all;
    for i=1:nx
        for j=1:nz
            if(energy_all(i,j)~=0) 
                dk(i,j)=grad0_all(i,j)/energy_all(i,j);
            end
        end
    end
    
    % Write gradient norm to file
    fid=fopen('dk_norm.txt','a');
    fprintf(fid,'iter=%d %.6e\n',iter,norm(dk));
    fclose(fid);
    
    % Write gradient to file
    fid1=fopen(['gradient',num2str(iter),'.dat'],'wb');
    fwrite(fid1, dk','float');
    fclose(fid1);
    
    % Compute step length (alpha) for the first iteration
    dgg_all=0;
    if(iter==1)
        parfor is=1:ns_hcp
            cal=zeros(nx,nt);
            [xsn,vp,rho,dk2]=current_shot(vp0,rho0,dk,nx,nz,npd,npd1,vnx,fs_hcp,ds_hcp,is,stype);
            vp=pad_vv(nx,nz,npd,npd1,vp);
            rho=pad_vv(nx,nz,npd,npd1,rho);
            vp=rho.*vp.^2;
            rho=1.0./rho;
            
            % Linearized modeling for gradient
            [cal,left_u,right_u,top_w,bottom_w,top_p,bottom_p,left_p,right_p,p_end,u_end,w_end] = linear_modeling(vp,rho,mm,c,cal,dk2,source,xsn,zs_hcp,zrec,dt,tmax,npd1,npd,nx,nz,dx,dz,nt,ndtt);
            
            % Calculate step length denominator
            for i=1:nx
                for it=1:ntout
                    dgg_all=dgg_all+cal(i,it)^2;
                end
            end 
        end
        
        % Calculate step length numerator
        gg=0;
        for i=1:nx
            for j=1:nz
                gg=gg+grad0_all(i,j)*dk(i,j);
            end
        end
        
        % Compute the step length
        alfa=gg/dgg_all;
    end
    
    % Write the step length to file
    fid=fopen('alfa.txt','a'); 
    a='iter=,alfa=';
    fprintf(fid,'%s %f %.6e\n',a,iter*1.0,alfa);
    fclose(fid);
    
    % Update reflectivity model
    ref0(:,:)=ref0(:,:)-alfa*dk(:,:);
    
    % Vectorize reflectivity model
    if iter==1
        Fx0=reshape(ref0,nx*nz,1);
    else
        Fx1=reshape(ref0,nx*nz,1);
    end
    
    % Update Anderson acceleration memory
    if iter>1
        if iter-1 <= mem_size - 1
            Smem(:, iter-1) = x1 - x0;
            Ymem(:, iter-1) = x1 - Fx1 - g0;
        else
            Smem = [Smem(:, 2:end), x1 - x0];
            Ymem = [Ymem(:, 2:end), x1 - Fx1 - g0];
        end
    end    
    
    % Update previous model and gradient
    if(iter==1)
        x0=zeros(nx*nz,1);
        g0 = x0 - Fx0;
    else
        x0 = x1;
        g0 = x0 - Fx1;    
    end
    
    % Write reflectivity model norm to file
    fid=fopen('ref0_norm.txt','a');
    fprintf(fid,'iter=%d %.6e\n',iter,norm(ref0));
    fclose(fid);
    
    % Calculate misfit with true model
    misfit=0;
    for i=1:nx
        for j=1:nz
            misfit=misfit+(ref0(i,j)-ref_true(i,j))^2;
        end
    end
    
    % Write misfit to file
    fid5=fopen('misfit.txt','a');
    fprintf(fid5,'iter=%d %.6e\n',iter,misfit);
    fclose(fid5);
    
    % Write migrated image to file
    fid2=fopen(['mig_lsrtm',num2str(iter),'.dat'],'wb');
    fwrite(fid2, ref0','float');
    fclose(fid2);
    
    % Record iteration time
    toc;
    fid4=fopen('time.txt','a');
    fprintf(fid4,'iter=%d %f\n',iter,toc);
    fclose(fid4); 
end
