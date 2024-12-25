function [source]=get_source(nt,dt,ts,nts,frequency,pi,pfac)
 pi2=pi*pi;
  t=0.0;
    source=zeros(nt,1);
  for i=1:nts 
    t=t+dt;
    x=frequency*(t-ts/2); 
    xx=x*x;	 
    source(i)=(1-2*pi2*xx)*exp(-(pi2*xx))*pfac;
  end 
end