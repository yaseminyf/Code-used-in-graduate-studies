% Written: 17.03.1999
% Updated and changed: 03.06.1999
%
% This routine is for checking the time-lagged planar search
% approach with synchronization.
% The residual received by the slaves is g lagged and the p 
% vector is p=Ds. Planar search is done after 
% j=g,2g,3g,4g updates. The D values are accumulated between two consecutive
% planar searches.
% In the planar search the 'old' value of the residual is used.
% In between two synchronizations a temporary residual is used.
% This temporary residual is updated using the 'short' version of updates 
% received from each block.
% Residual norm is computed only after each planar search.
% p vector is also 'short'
%
% this is the above algorithm with p=0

matdir='/Home/jute/yasemin/matlab/';         

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

group=4;
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

xfiles(1,:)=  '/Home/jute/yasemin/programs/sequential/ash219x';
xfiles(5,:)=  '/Home/jute/yasemin/programs/sequential/abb313x';
xfiles(3,:)=  '/Home/jute/yasemin/programs/sequential/ash331x';
xfiles(4,:)=  '/Home/jute/yasemin/programs/sequential/ash608x';
xfiles(2,:)=  '/Home/jute/yasemin/programs/sequential/ash958x';
xfiles(6,:)=  '/Home/jute/yasemin/matlab/wl1033x             ';
xfiles(7,:)=  '/Home/jute/yasemin/matlab/wl1252x             ';
xfiles(8,:)=  '/Home/jute/yasemin/matlab/wl1364x             ';
xfiles(9,:)=  '/Home/jute/yasemin/matlab/wl1641x             ';
xfiles(10,:)=  '/Home/jute/yasemin/matlab/wl1850x             ';
xfiles(11,:)=  '/Home/jute/yasemin/matlab/wl1991x             ';
xfiles(12,:)=  '/Home/jute/yasemin/matlab/wl2808x             ';

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


% compute tildeA matrices and factorize them
% Bn = An'*An
% Rn = chol(Bn)

for i=1:group
   aa=eval(['A' num2str(i)]);
   aatrans=aa';
   eval(['B' num2str(i) '= aatrans*aa;']);
   eval(['R' num2str(i) '=chol(B' num2str(i) ');']);
end      

% initialize the working residual vector

res_work=[];
for i=1:group
   res_work=[res_work residual];
end

resnorms=[];
res_zero_norm=norm(residual);
res_norm=norm(residual);
err=res_norm/res_zero_norm;
resnorms=[resnorms;res_norm];
if matfile == 6 
   epsilon=1.0e-3;
elseif  matfile == 10
   epsilon=1.0e-3;  
else
   epsilon=1.0e-5;
end   
iter=0;

% begin the main loop

% kcount is given outside. it is to decide when to do a planar search

converged=0;
kount=0;
Ad=zeros(m,group);
D=zeros(n,1);
[mm,nn]=size(D);
v=residual;
numplan=0;
while converged==0 & iter <30000
   for i=1:group
      iter=iter+1;
      aa=eval(['A' num2str(i)]);
      RR=eval(['R' num2str(i)]);
      atrans=aa';
      Rtrans=RR';
      eval(['y = Rtrans \ (-atrans*res_work(:,group));']); 
      eval(['d' num2str(i) '=RR \ y;']);        % solve for d
      dtilde=eval(['d' num2str(i)]);
      aa=eval(['A' num2str(i)]);
      eval(['Ad' num2str(i) '=aa*dtilde;']);
      ad=eval(['Ad' num2str(i)]);
      Ad(:,i)=Ad(:,i)+ad;                         % add A*d to Ad
      v=v+ad;                   % update the temporary residual
      res_work=[v res_work(:,1:group-1)];   % update working residual
      if i == group
        D((i-1)*nc+1:n)=D((i-1)*nc+1:n)+dtilde;
      else
        D((i-1)*nc+1:i*nc)=D((i-1)*nc+1:i*nc)+dtilde;
      end;
   end
   kount=kount+1;
   if kount==kcount         % decide if planar search time
      kount=0;
      numplan=numplan+1;   
      hatA=[];                     % form hatA
      for k=1:group
         hatA=[hatA Ad(:,k)]; 
      end
          
      

% do the planar search
      
      Bighat=hatA'*hatA;
      Rhat=chol(Bighat);
      y=Rhat' \ (-hatA'*residual);
      s=Rhat \ y;                  % compute s with the last residual
      xupdate=[];
      for i=1:group                % update residual and x vector
         residual=residual+s(i)*hatA(:,i);
         if i==group
            sss=s(i)*D((i-1)*nc+1:n);
         else
            sss=s(i)*D((i-1)*nc+1:i*nc);
         end
         xupdate=[xupdate;sss];
      end
      x_guess=x_guess+xupdate;
     
      Ad=zeros(m,group);
      D=zeros(n,1);
      res_work=[];
      for i=1:group
         res_work=[res_work residual];
      end
      res_norm=norm(residual);
      err=res_norm/res_zero_norm;
      resnorms=[resnorms; res_norm];
      v=residual;
      if err < epsilon
         converged=1;
      end
   end
  
end
      


   
                                