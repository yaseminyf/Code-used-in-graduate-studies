% Written: 04.02.1999
% Updated: 23.03.2003
%
% The residual received by the slaves is g lagged and the p 
% vector is zero. linesearch is used for computing the relaxation
% parameter. the relaxation parameter is not limited between zero
% and 2 to check if it is equivalent to planar search in this manner.

matdir='D:\Matlab\';         

format long e;

files(1,:)=  'ash219';
files(5,:)=  'abb313';
files(3,:)=  'ash331';
files(4,:)=  'ash608';
files(2,:)=  'ash958';
files(6,:)=  'wl1033';
files(7,:)=  'wl1252';
files(8,:)=  'wl1364';
files(9,:)=  'wl1641';
files(10,:)= 'wl1850';
files(11,:)= 'wl1991';
files(12,:)= 'wl2808';
 
% matfile is the input given. the test matrix to be used

string=['load ' matdir files(matfile,:)];
eval(string);

[m,n]=size(A);

% choose number of blocks

group=15;
nc=fix(n/group);
last=n-nc*(group-1);

% create the submatrices A1..Agroup

for i=1:group
   j=(i-1)*nc;
   if i~=group 
      k=nc;
   else 
      k=last;
   end
   eval(['A' num2str(i) '=A(:,j+1:j+k);']);
end
   
% read in the x-values

xfiles(1,:)=  'ash219x';
xfiles(5,:)=  'abb313x';
xfiles(3,:)=  'ash331x';
xfiles(4,:)=  'ash608x';
xfiles(2,:)=  'ash958x';
xfiles(6,:)=  'wl1033x';
xfiles(7,:)=  'wl1252x';
xfiles(8,:)=  'wl1364x';
xfiles(9,:)=  'wl1641x';
xfiles(10,:)=  'wl1850x';
xfiles(11,:)=  'wl1991x';
xfiles(12,:)=  'wl2808x';

string=['load ' xfiles(matfile,:)];
eval(string);

if matfile == 1
   x=ash219x;
elseif matfile == 2
   x=ash958x;
elseif matfile == 3
   x=ash331x;
elseif matfile == 4
   x=ash608x;
elseif matfile == 5
   x=abb313x;
elseif matfile == 6
   x=wl1033x;
elseif matfile == 7
   x=wl1252x;
elseif matfile == 8
   x=wl1364x;
elseif matfile == 9
   x=wl1641x;
elseif matfile == 10
   x=wl1850x;
elseif matfile == 11
   x=wl1991x;
elseif matfile == 12
   x=wl2808x;
end       

% compute the right-hand-side - zero residual case

rhs=A*x;

% compute the right-hand-side - nonzero residual case
%
% rhs=rand(m,1);

% assign the initial x vector to be all zeros

x_guess=zeros(n,1);

% assign the initial guess to be a random vector
%
% x_guess=rand(n,1);

% compute the initial residual

residual=A*x_guess-rhs;

% do the cholesky factorization of all blocks
% Bn = An'*An
% Rn = chol(Bn)

for i=1:group
   aa=eval(['A' num2str(i)]);
   aatrans=aa';
   eval(['B' num2str(i) '= aatrans*aa;']);
   eval(['R' num2str(i) '=chol(B' num2str(i) ');']);
end   

resnorms=[];
res_zero_norm=norm(residual);
res_norm=norm(residual);
resnorms=[resnorms res_norm];
err=res_norm/res_zero_norm;
if matfile == 6 
   epsilon=1.0e-3;
elseif  matfile == 10
   epsilon=1.0e-3;  
else
   epsilon=1.0e-5;
end   
iter=0;

% initialize the working residual vector

res_work=[];
for i=1:group
   res_work=[res_work residual];
end

% begin the main loop

converged=0;

while converged==0 & iter <30000
   
   for i=1:group
      iter=iter+1;

      aa=eval(['A' num2str(i)]);
      RR=eval(['R' num2str(i)]);
      atrans=aa';
      Rtrans=RR';
      eval(['y = Rtrans \ (-atrans*res_work(:,group));']); 
      eval(['s' num2str(i) '=RR \ y;']);   % solve for s
      ss=eval(['s' num2str(i)]);
      eval(['As' num2str(i) '=aa*ss;']);
      ass=eval(['As' num2str(i)]);
      check=ass'*residual;
      dotnorm=ass'*ass;
      omega=-check/dotnorm;
      if (omega > 0) & (omega < 2.0)
         residual=residual+omega*ass;               % update the residual
         if i == group
            x_guess((i-1)*nc+1:n)=x_guess((i-1)*nc+1:n)+omega*ss;
         else
            x_guess((i-1)*nc+1:i*nc)=x_guess((i-1)*nc+1:i*nc)+omega*ss;
         end
         res_norm=norm(residual);
         err=res_norm/res_zero_norm;
         resnorms=[resnorms;res_norm];

         if err < epsilon
            converged=1;
            break;
         else
            res_work=[residual res_work(:,1:group-1)];
         end
      end
   end
   
end
      


   
                                