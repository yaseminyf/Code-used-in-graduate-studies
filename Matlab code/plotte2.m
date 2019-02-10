[m,n]=size(resnorms_async);
x=[0:1:m-1];
semilogy(x,resnorms_async)
[m,n]=size(tmp_res);
x=[0:1:m-1];
semilogy(x,tmp_res,':')
x=[0:group*kcount:m-1]
semilogy(x,resnorms_l1,'x')