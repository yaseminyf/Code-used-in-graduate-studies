% Written: 26.02.2002
% Updated: 06.04.2002
%
% This routine is for checking the Jacobi with planar search
% approach with p=0. It is to simulate the behaviour of the synchronous algorithm.
% The residual received by the slaves is the same and the p
% vector is zero. Planar search is done after "kcount" updates which is given when the function is called.
% The number of groups g is 8.
% The identity of the matrix, "matfile",is also given outside at call time.
% Residual norm is taken only after kcount*g updates on the residual.
%
% All the processors use the same residual in their computations.
% Iteration counter is updated after each loop through g blocks.
% The updates d_i computed by the processors in kcount loops are accumulated, so as the A_id_i.
% The residual value at the end of a planar search is used in the next planar search.
% Between two planar searches a working residual is updated and used in the computations.

matdir='D:\matlab\';

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

% choose the coefficient matrix

% matfile is the input given. the test matrix to be used

string=['load ' matdir files(matfile,:)];
eval(string);

[m,n]=size(A);

% choose number of blocks

group=8;
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

% kcount is given outside. it is to decide when to do a planar search

% begin the main loop

residual_work=residual;
converged=0;
j=0;
hatA=zeros(m,group);
d=zeros(n,1);
while converged==0 & iter <10000
   residual_update=[];

   for i=1:group
      aa=eval(['A' num2str(i)]);
      RR=eval(['R' num2str(i)]);
      atrans=aa';
      Rtrans=RR';
      eval(['y = Rtrans \ (-atrans*residual_work);']);
      eval(['d' num2str(i) '=RR \ y;']);   % solve for d
      dd=eval(['d' num2str(i)]);
      eval(['Ad' num2str(i) '=aa*dd;']);
      add=eval(['Ad' num2str(i)]);
      hatA(:,i)=hatA(:,i)+add;                     % add A*d to hatA
      residual_update=[residual_update add];
      if i == group                      % accumulate the d vectors to update x vector later
         d(((i-1)*nc)+1:n)=d(((i-1)*nc)+1:n)+dd;
      else
         d(((i-1)*nc)+1:i*nc)=d(((i-1)*nc)+1:i*nc)+dd;
      end
   end
   iter=iter+1;
   for i=1:group
      residual_work=residual_work+residual_update(:,i); % update the working residual
   end
   j=j+1;
   if j==kcount     % decide if planar search time
      j=0;
      Bighat=hatA'*hatA;
      Rhat=chol(Bighat);
      y=Rhat' \ (-hatA'*residual);
      s=Rhat \ y;                     % compute s with the last residual
      xupdate=[];
      for i=1:group          % update residual and x vector
         residual=residual+s(i)*hatA(:,i);
         if i == group
            x_guess(((i-1)*nc)+1:n)=x_guess(((i-1)*nc)+1:n)+s(i)*d(((i-1)*nc)+1:n);
         else
            x_guess(((i-1)*nc)+1:i*nc)=x_guess(((i-1)*nc)+1:i*nc)+s(i)*d(((i-1)*nc)+1:i*nc);
         end
      end

% check if converged

      res_norm=norm(residual);
      err=res_norm/res_zero_norm;
      resnorms=[resnorms res_norm];
      if err < epsilon
         converged=1;
         numplan=iter/kcount;
      else
         residual_work=residual;
         hatA=zeros(m,group);
         d=zeros(n,1);
      end
   end

end





