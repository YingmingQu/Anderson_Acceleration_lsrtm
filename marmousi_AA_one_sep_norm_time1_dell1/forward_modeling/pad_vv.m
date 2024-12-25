function [ee]=pad_vv(nx,nz,npd,npd1,ee) 
%pad left side
for j=(1+npd1):(nz+npd1)
    for i=(0+npd1):(-1):(-npd+1+npd1)
        ee(i,j)=ee((1+npd1),j);
    end
end
%pad right side
for j=(1+npd1):(nz+npd1)
    for i=(nx+1+npd1):(nx+npd+npd1)
         ee(i,j)=ee(nx+npd1,j);
    end
end
%pad upper side
for j=(0+npd1):(-1):(-npd+1+npd1)
     for i=(-npd+1+npd1):(nx+npd+npd1)
          ee(i,j)=ee(i,1+npd1);
     end
end
%lower side
for j=(nz+1+npd1):(nz+npd+npd1)
     for i=(-npd+1+npd1):(nx+npd+npd1)
           ee(i,j)=ee(i,nz+npd1);
     end
end
end