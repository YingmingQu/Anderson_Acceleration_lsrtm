clc; clear; close all;
nx=368;
nz=6001;
ns=19;
fid=fopen('E:\程序\AA-QR\marmousi_AA_QR_sep_norm_time\forward_modeling\shot_born.dat','rb');  
     [shot]=fread(fid,[nz,nx*ns],'float');
  fclose(fid);
     shot=shot';
     sho1=zeros(ns,nx,nz);
     for i=1:ns
         for j=1:nx
             for k=1:nz
             shot1(i,j,k)=shot((i-1)*nx+j,k);    
             end
         end
     end
     shot2=zeros(nx,nz);
     shot2(1:nx,1:nz)=shot1(10,1:nx,1:nz);
     single=shot2(nx/2,1:nz);
     shot3=shot2';
     figure(1);
     imagesc(shot3);
       figure(2);
       
       
        single_noise_total=zeros(ns,nx,nz);
        single_noise_total1=zeros(ns*nx,nz);
       
for i=1:ns    
 for  j=1:nx  
      single=zeros(1,nz);
   single(1,1:nz)= shot1(i,j,1:nz);
     noise=rand(1,nz);
     noise1= (noise-0.5)/4;
     single_noise1=single+noise1;
     
     if(i==10&&j==(nx/2))
      hold on;
        plot(single_noise1);
          plot(single);
     end
      single_noise_total(i,j,1:nz)=single_noise1;   
 end  

end



 for i=1:ns
         for j=1:nx
             for k=1:nz
             single_noise_total1((i-1)*nx+j,k)=single_noise_total(i,j,k);    
             end
         end
     end
  single_noise_total2= single_noise_total1';
  
 snr=0;
         I=shot; %initial signal
         In= single_noise_total1;           %signal with noise
Ps=sum(sum((I-mean(mean(I))).^2));%signal power
Pn=sum(sum((I-In).^2)); %noise power
snr=10*log10(Ps/Pn); 