fid=fopen('jacobi_sync_pds_nodelta.results','a');
deltasin=0;
for jjjj=1:6
  kcount=jjjj;
  fprintf(fid,'kcount=%g\n',kcount);
  for iiii=1:12
     matfile=iiii;
     jacobi_sync_pds_planar
     [hh,ll]=size(resnorms);
     ii=iter;
     numplan=iter/kcount;
     fprintf(fid,'matfile=%g\n',matfile);
     fprintf(fid,'number of iterations %g\n',ii);
     fprintf(fid,'number of planar searches %g\n',numplan);
     fprintf(fid,'resnorm at zero %e\n',resnorms(1));
     fprintf(fid,'resnorm at finish %e\n\n',resnorms(ll));
  end  
end   
fclose(fid);
