%
kcount=2;
%fid=fopen('jacobi_sync_p0_planar.results','a');
%fprintf(fid,'kcount=%g\n',kcount);
%fclose(fid);
for iiii=8:12
   fid=fopen('jacobi_sync_p0_planar.results','a');
   matfile=iiii;
   jacobi_sync_p0_planar
   [hh,ll]=size(resnorms);
   ii=iter;
   numplan=ii/kcount;
   fprintf(fid,'matfile=%g\n',matfile);
   fprintf(fid,'number of iterations %g\n',ii);
   fprintf(fid,'number of planar searches %g\n',numplan);
   fprintf(fid,'resnorm at zero %e\n',resnorms(1));
   fprintf(fid,'resnorm at finish %e\n\n',resnorms(ll));
   fclose(fid);
end
%
kcount=3;
fid=fopen('jacobi_sync_p0_planar.results','a');
fprintf(fid,'kcount=%g\n',kcount);
fclose(fid);
for iiii=1:6
   fid=fopen('jacobi_sync_p0_planar.results','a');
   matfile=iiii;
   jacobi_sync_p0_planar
   [hh,ll]=size(resnorms);
   ii=iter;
   numplan=ii/kcount;
   fprintf(fid,'matfile=%g\n',matfile);
   fprintf(fid,'number of iterations %g\n',ii);
   fprintf(fid,'number of planar searches %g\n',numplan);
   fprintf(fid,'resnorm at zero %e\n',resnorms(1));
   fprintf(fid,'resnorm at finish %e\n\n',resnorms(ll));
   fclose(fid);
end
%
kcount=3;
for iiii=7:12
   fid=fopen('jacobi_sync_p0_planar.results','a');
   matfile=iiii;
   jacobi_sync_p0_planar
   [hh,ll]=size(resnorms);
   ii=iter;
   numplan=ii/kcount;
   fprintf(fid,'matfile=%g\n',matfile);
   fprintf(fid,'number of iterations %g\n',ii);
   fprintf(fid,'number of planar searches %g\n',numplan);
   fprintf(fid,'resnorm at zero %e\n',resnorms(1));
   fprintf(fid,'resnorm at finish %e\n\n',resnorms(ll));
   fclose(fid);
end
%
kcount=4;
fid=fopen('jacobi_sync_p0_planar.results','a');
fprintf(fid,'kcount=%g\n',kcount);
fclose(fid);
for iiii=1:6
   fid=fopen('jacobi_sync_p0_planar.results','a');
   matfile=iiii;
   jacobi_sync_p0_planar
   [hh,ll]=size(resnorms);
   ii=iter;
   numplan=ii/kcount;
   fprintf(fid,'matfile=%g\n',matfile);
   fprintf(fid,'number of iterations %g\n',ii);
   fprintf(fid,'number of planar searches %g\n',numplan);
   fprintf(fid,'resnorm at zero %e\n',resnorms(1));
   fprintf(fid,'resnorm at finish %e\n\n',resnorms(ll));
   fclose(fid);
end
%
kcount=4;
for iiii=7:12
   fid=fopen('jacobi_sync_p0_planar.results','a');
   matfile=iiii;
   jacobi_sync_p0_planar
   [hh,ll]=size(resnorms);
   ii=iter;
   numplan=ii/kcount;
   fprintf(fid,'matfile=%g\n',matfile);
   fprintf(fid,'number of iterations %g\n',ii);
   fprintf(fid,'number of planar searches %g\n',numplan);
   fprintf(fid,'resnorm at zero %e\n',resnorms(1));
   fprintf(fid,'resnorm at finish %e\n\n',resnorms(ll));
   fclose(fid);
end
%
kcount=5;
fid=fopen('jacobi_sync_p0_planar.results','a');
fprintf(fid,'kcount=%g\n',kcount);
fclose(fid);
for iiii=1:6
   fid=fopen('jacobi_sync_p0_planar.results','a');
   matfile=iiii;
   jacobi_sync_p0_planar
   [hh,ll]=size(resnorms);
   ii=iter;
   numplan=ii/kcount;
   fprintf(fid,'matfile=%g\n',matfile);
   fprintf(fid,'number of iterations %g\n',ii);
   fprintf(fid,'number of planar searches %g\n',numplan);
   fprintf(fid,'resnorm at zero %e\n',resnorms(1));
   fprintf(fid,'resnorm at finish %e\n\n',resnorms(ll));
   fclose(fid);
end
%
kcount=5;
for iiii=7:12
   fid=fopen('jacobi_sync_p0_planar.results','a');
   matfile=iiii;
   jacobi_sync_p0_planar
   [hh,ll]=size(resnorms);
   ii=iter;
   numplan=ii/kcount;
   fprintf(fid,'matfile=%g\n',matfile);
   fprintf(fid,'number of iterations %g\n',ii);
   fprintf(fid,'number of planar searches %g\n',numplan);
   fprintf(fid,'resnorm at zero %e\n',resnorms(1));
   fprintf(fid,'resnorm at finish %e\n\n',resnorms(ll));
   fclose(fid);
end
%
kcount=6;
fid=fopen('jacobi_sync_p0_planar.results','a');
fprintf(fid,'kcount=%g\n',kcount);
fclose(fid);
for iiii=1:6
   fid=fopen('jacobi_sync_p0_planar.results','a');
   matfile=iiii;
   jacobi_sync_p0_planar
   [hh,ll]=size(resnorms);
   ii=iter;
   numplan=ii/kcount;
   fprintf(fid,'matfile=%g\n',matfile);
   fprintf(fid,'number of iterations %g\n',ii);
   fprintf(fid,'number of planar searches %g\n',numplan);
   fprintf(fid,'resnorm at zero %e\n',resnorms(1));
   fprintf(fid,'resnorm at finish %e\n\n',resnorms(ll));
   fclose(fid);
end
%
kcount=6;
for iiii=7:12
   fid=fopen('jacobi_sync_p0_planar.results','a');
   matfile=iiii;
   jacobi_sync_p0_planar
   [hh,ll]=size(resnorms);
   ii=iter;
   numplan=ii/kcount;
   fprintf(fid,'matfile=%g\n',matfile);
   fprintf(fid,'number of iterations %g\n',ii);
   fprintf(fid,'number of planar searches %g\n',numplan);
   fprintf(fid,'resnorm at zero %e\n',resnorms(1));
   fprintf(fid,'resnorm at finish %e\n\n',resnorms(ll));
   fclose(fid);
end
%
