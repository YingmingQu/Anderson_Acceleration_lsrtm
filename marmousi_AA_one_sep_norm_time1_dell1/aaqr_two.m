 function [xc,Fx1,Fx0,g1,g0,Q,R,DF,mAA]=aaqr(nx,nz,dx,dz,iter,mem_size,dk_reshape,alfa,xc,Fx1,Fx0,g1,g0,Q,R,DF,mAA,droptol,beta)

% Update the df vector and the DG array.
    if(iter>1)
dg=g1-g0;

  if mAA < mem_size
  DF = [DF Fx1-Fx0];
  else
  DF = [DF(:,2:mAA) Fx1-Fx0];
  end
  mAA=mAA+1;
    end

    Fx0=Fx1;
g0=g1;

if(mAA==0)
    xc=Fx1;
else
if(mAA==1)
% form the initial QR decomposition.
R(1,1) = norm(dg);
Q = R(1,1)\dg;

else

 if(mAA > mem_size)
[Q,R] = qrdelete(Q,R,1);
mAA=mAA-1;
  if size(R,1) ~= size(R,2)
  Q = Q(:,1:mAA-1); R = R(1:mAA-1,:);
  end
  
 end

for j = 1:mAA-1
    QQ=Q(:,j);
R(j,mAA) = QQ'*dg;
dg = dg - R(j,mAA)*Q(:,j);
end
R(mAA,mAA) = norm(dg);
Q = [Q,R(mAA,mAA)\dg];

end

if droptol > 0
% Drop residuals to improve conditioning if necessary.
condDF = cond(R);
fprintf(' condDF %.4e \n', condDF);

while condDF > droptol && mAA > 1
    
fprintf(' cond(D) = %e, reducing iter to %d \n', condDF, mAA-1);
[Q,R] = qrdelete(Q,R,1);
DF = DF(:,2:mAA);
mAA=mAA-1;
% The following treats the qrdelete quirk described above.
if size(R,1) ~= size(R,2)
Q = Q(:,1:mAA); R = R(1:mAA,:);
end
condDF = cond(R);
end
end

% Solve the least-squares problem.
gamma = R\(Q'*g1);
% Update the approximate solution.
xc = Fx1 - DF*gamma;
% Apply damping if beta is a function handle or if beta > 0
% (and beta ~= 1).
if isa(beta,'function_handle')
xc = xc - (1-beta(k))*(Fx1- Q*R*gamma);
else
if beta > 0 && beta ~= 1
xc = xc - (1-beta)*(Fx1 - Q*R*gamma);
end
end
% xc=(1-beta)*xc+beta*Fx1;
end
end