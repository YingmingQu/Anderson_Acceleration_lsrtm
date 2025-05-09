% Born modeling
clc; clear; close all;

% Model parameters
vnx=368;              % Number of grid points in x-direction for velocity model
nx=368;               % Number of grid points in x-direction for computation
nz=200;               % Number of grid points in z-direction
dx=10.0;              % Grid spacing in x-direction (m)
dz=10.0;              % Grid spacing in z-direction (m)
ns_hcp=19;            % Number of shots (originally 19)
ds_hcp=20;            % Shot spacing (grid points)
fs_hcp=1;             % First shot location (grid points)
zs_hcp=1;             % Source depth (grid points)
dt=0.5;               % Time step (ms)
dtout=0.5;            % Output time interval (ms)
tmax=3.0;             % Total simulation time (s)
npd=50;               % Absorbing boundary thickness (grid points)
npd1=50;              % Another absorbing boundary thickness (grid points)
frequency=20.0;       % Source dominant frequency (Hz)
ns_start=1;           % Starting shot number (purpose uncertain)
mm=5;                 % Staggered grid parameter
ntrace=1;             % Number of traces
stype=2;              % Source type
wtype=1;              % Wavefield type
zrec=2;               % Receiver depth (grid points)
pi=4.0*atan(1.0);     % Compute pi
pfac=1.0;             % Source factor

% Calculate staggered grid coefficients
c=cal_c(mm);

% Read velocity model from file
fid=fopen('vpsmooth.dat','rb');  
[vp0,count]=fread(fid,[nz,nx],'float');
vp0=vp0';  % Transpose to match MATLAB's row-column representation
fclose(fid);

% Initialize density model (uniform density)
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

% Start timing
tic

% Generate source function
ts=2.0/frequency; 
nts=floor(ts/dt); 
source=get_source(nt,dt,ts,nts,frequency,pi,pfac);

% Initialize reflectivity model
ref0=zeros(nx,nz);

% Initialize array to store modeled data for all shots
cal1=zeros(ns_hcp,nx,nt);

% Parallel processing for each shot
parfor is=1:ns_hcp
    cal=zeros(nx,nt);  % Initialize modeled data for current shot
    res0=0;            % Initialize residual
    
    % Prepare model parameters for current shot
    [xsn,vp,rho]=current_shot1(vp0,rho0,nx,nz,npd,npd1,vnx,fs_hcp,ds_hcp,is,stype);
    
    % Apply padding for absorbing boundaries
    vp=pad_vv(nx,nz,npd,npd1,vp);
    rho=pad_vv(nx,nz,npd,npd1,rho);
    
    % Convert to Lame parameters
    vp=rho.*vp.^2;
    rho=1.0./rho;
    
    % Perform linearized modeling using true reflectivity
    [cal] = linear_modeling(vp,rho,mm,c,cal,ref_true,source,xsn,zs_hcp,zrec,dt,tmax,npd1,npd,nx,nz,dx,dz,nt,ndtt);
    
    % Store modeled data for current shot
    cal1(is,:,:)=cal(:,:);
end

% Reshape modeled data into a single matrix
cal2=zeros(nt,ns_hcp*nx);
for k=1:ns_hcp
    for j=1:nt
        for i=1:nx
            cal2(j,i+(k-1)*nx)=cal1(k,i,j);
        end
    end
end

% Write modeled data to file
fid=fopen('shot_born.dat','wb');
fwrite(fid,cal2,'float');
fclose(fid);

% End timing and display elapsed time
toc
