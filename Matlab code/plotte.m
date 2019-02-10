 [m,n]=size(resnorms);
x=[1:group:(m)*group];
%drawchar is given from outside
semilogy(x,resnorms,'o-');