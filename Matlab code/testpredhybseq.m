fid=fopen('testpredhybseq.results','w');

for jjjj=1:6
  kcount=jjjj;
  fprintf(fid,'kcount=%g\n',kcount);
  for iiii=1:12
     matfile=iiii;
     predict_hybrid_seq
     [hh,ll]=size(resnorms);
     ii=iter;
     fprintf(fid,'matfile=%g\n',matfile);
     fprintf(fid,'number of planar searches for r %g\n',numplanr);
     fprintf(fid,'number of planar searches for p %g\n',numplanp);
     fprintf(fid,'number of iterations %g\n',ii);
     fprintf(fid,'resnorm at zero %e\n',resnorms(1));
     fprintf(fid,'resnorm at finish %e\n\n',resnorms(hh));
  end  
end   
fclose(fid);
