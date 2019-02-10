
% Written: 28.03.2002
% Updated: 05.04.2002
%
% This routine is for checking the Jacobi with planar search
% approach with p set to Ds. It is to simulate the behaviour of the asynchronous algorithm.
% The residual received by the slaves is the same and the p
% vector is set to Ds. Planar search is done after "kcount" updates which is given when the
% function is called.
% The number of groups g is 8.
% The identity of the matrix, "matfile", is also given outside at call time.
% Residual norm is taken only after kcount*g updates on the residual.
%
% This is the timelagged simulation of asynchronous operation.
% The timelag is equal to the number of groups meaning each group is assigned to one processor.
% The order of the updates is fixed, cyclic from group one up to group g.
%
% All the processors use a different residual in their computations.
% Iteration counter is updated after each block sends his update.
% The updates d_i computed by the processors in kcount loops are accumulated, so as the A_id_i.
% The residual value at the end of a planar search is used in the next planar search.
% Between two planar searches a working residual is updated and used in the computations.
%
% The parameter "deltasin" is used to decide whether the delta values are used in the updates
% or are thrown away. This parameter receives its value outside of the function call.
% Even if the deltasin is set to 1, the delta values are not used outside
% of the planar search because they are incomplete.

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

% choose the coefficient matrix

% matfile is the input given. the test matrix to be used

string=['load ' matdir files(matfile,:)];
eval(string);

[m,n]=size(A);

% choose number of blocks

group=30;
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

% choose the starting p vector - all ones

p=ones(n,1);

% choose the p vector - FM values

%k=1;
%p=zeros(n,1);
%for i=1:group
%   aa=eval(['A' num2str(i)]);
%   value=aa'*aa;
%   [mm,nn]=size(aa);
%   evec=ones(nn,1);
%   s=value*evec;
%   for j=1:nn
%      p(k)=1/s(j);
%      k=k+1;
%   end
%end

% compute the F matrix

F=[];

for i=1:group
   aa=eval(['A' num2str(i)]);
   if i==group
      ap=aa*p((i-1)*nc+1:n);
   else
      ap=aa*p((i-1)*nc+1:i*nc);
   end
   F=[F ap];
end

% compute tildeA matrices and factorize them
% Bn = An'*An
% Rn = chol(Bn)

for i=1:group
   aa=eval(['A' num2str(i)]);
   eval(['tildeA' num2str(i) '=[F(:,1:i-1) aa F(:,i+1:group)];']);
   tA=eval(['tildeA' num2str(i)]);
   tAtrans=tA';
   eval(['B' num2str(i) '= tAtrans*tA;']);
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

% initialize the working residual vectors

residual_work=[];
for i=1:group
   residual_work=[residual_work residual];
end
residual_current=residual;

% residual_current is used during no planar search phase
% the original residual right after the planar search is saved to be used
% in the next planar search
tmp_res=[];
tmp_res=[tmp_res norm(residual_current)];

converged=0;
j=0;
hatA=zeros(m,group);
d=zeros(n,1);  % to accumulate d vectors
if deltasin == 1
   delta=zeros(group,1);   % to accumulate deltas
end

while converged==0 & iter <30000
   
   for i=1:group
      iter=iter+1;

      aa=eval(['tildeA' num2str(i)]);
      RR=eval(['R' num2str(i)]);
      atrans=aa';
      Rtrans=RR';
      eval(['y = Rtrans \ (-atrans*residual_work(:,group));']);
      eval(['d' num2str(i) '=RR \ y;']);   % solve for d
      dtilde=eval(['d' num2str(i)]);

% separate d and deltas and store for future use and compute Ad's

      aa=eval(['A' num2str(i)]);
      if i==group
         eval(['Ad' num2str(i) '=aa*dtilde(i:i+last-1);']);
         d(((i-1)*nc)+1:n)=d(((i-1)*nc)+1:n)+dtilde(i:i+last-1);
      else
         eval(['Ad' num2str(i) '=aa*dtilde(i:i+nc-1);']);
         d(((i-1)*nc)+1:i*nc)=d(((i-1)*nc)+1:i*nc)+dtilde(i:i+nc-1);
      end
      add=eval(['Ad' num2str(i)]);
      hatA(:,i)=hatA(:,i)+add;                     % add A*d to hatA
      if deltasin == 1
         for jj=1:group
            if jj~=i
               if jj<i
                  delta(jj)=delta(jj)+dtilde(jj);
               elseif j>i
                  delta(jj)=delta(jj)+dtilde(jj+nc-1);
               end
            end
         end
      end
      residual_current=residual_current+add;
      residual_work=[residual_current residual_work(:,1:group-1)];
      tmp_res=[tmp_res; norm(residual_current)];

   end
   j=j+1;
   if j==kcount     % decide if planar search time
      j=0;
      if deltasin == 1 % add the delta part to hatA
         for i=1:group
	        hatA(:,i)=hatA(:,i)+delta(i)*F(:,i);
         end
      end
     
      Bighat=hatA'*hatA;
      Rhat=chol(Bighat);
      y=Rhat' \ (-hatA'*residual);
      s=Rhat \ y;                     % compute s with the last residual
     
      for i=1:group          % update residual and x vector
         residual=residual+s(i)*hatA(:,i);
         F(:,i)=s(i)*hatA(:,i);

% we do not need to form the new p vector because we only need the F matrix
% the rest is done in pieces according to the value of deltasin variable
% p vector is needed for updating the x vector if deltasin is set to one

         if i == group
            x_guess(((i-1)*nc)+1:n)=x_guess(((i-1)*nc)+1:n)+s(i)*d(((i-1)*nc)+1:n);
	         if deltasin == 1
               x_guess(((i-1)*nc)+1:n)=x_guess(((i-1)*nc)+1:n)+s(i)*delta(i)*p(((i-1)*nc)+1:n);
               p(((i-1)*nc)+1:n)=s(i)*(d(((i-1)*nc)+1:n)+delta(i)*p(((i-1)*nc)+1:n));
	         end
         else
            x_guess(((i-1)*nc)+1:i*nc)=x_guess(((i-1)*nc)+1:i*nc)+s(i)*d(((i-1)*nc)+1:i*nc);
            if deltasin == 1
               x_guess(((i-1)*nc)+1:i*nc)=x_guess(((i-1)*nc)+1:i*nc)+s(i)*delta(i)*p(((i-1)*nc)+1:i*nc);
               p(((i-1)*nc)+1:i*nc)=s(i)*(d(((i-1)*nc)+1:i*nc)+delta(i)*p(((i-1)*nc)+1:i*nc));
	         end
         end
      end

% check if converged

      res_norm=norm(residual);
      err=res_norm/res_zero_norm;
      resnorms=[resnorms; res_norm];
      if err < epsilon
         converged=1;
         numplan=iter/(group*kcount);
      else

% form and factorize the new tildeA matrices

         for i=1:group
            aa=eval(['A' num2str(i)]);
            eval(['tildeA' num2str(i) '=[F(:,1:i-1) aa F(:,i+1:group)];']);
            tA=eval(['tildeA' num2str(i)]);
            tAtrans=tA';
            eval(['B' num2str(i) '= tAtrans*tA;']);
            eval(['R' num2str(i) '=chol(B' num2str(i) ');']);
         end
         residual_work=[];
         for i=1:group
            residual_work=[residual_work residual];
         end
         residual_current=residual;
         hatA=zeros(m,group);
         d=zeros(n,1);
	      if deltasin == 1
	         delta=zeros(group,1);
	      end
      end
   end

end





