deltasin=0;
kcount=1;
%fid=fopen('jacobi_sync_pconst_nodelta.results','a');
%fprintf(fid,'kcount=%g\n',kcount);
%fclose(fid);
%for iiii=1:6
%   fid=fopen('jacobi_sync_pconst_nodelta.results','a');
%   matfile=iiii;
%   jacobi_sync_pconst_planar
%   [hh,ll]=size(resnorms);
%   ii=iter;
%   numplan=iter/kcount;
%   fprintf(fid,'matfile=%g\n',matfile);
%   fprintf(fid,'number of iterations %g\n',ii);
%   fprintf(fid,'number of planar searches %g\n',numplan);
%   fprintf(fid,'resnorm at zero %e\n',resnorms(1));
%   fprintf(fid,'resnorm at finish %e\n\n',resnorms(ll));
%   fclose(fid);
%end
%
%deltasin=0;
%kcount=1;
%for iiii=6:12
%   fid=fopen('jacobi_sync_pconst_nodelta.results','a');
%   matfile=iiii;
%   jacobi_sync_pconst_planar
%   [hh,ll]=size(resnorms);
%   ii=iter;
%   numplan=iter/kcount;
%   fprintf(fid,'matfile=%g\n',matfile);
%   fprintf(fid,'number of iterations %g\n',ii);
%   fprintf(fid,'number of planar searches %g\n',numplan);
%   fprintf(fid,'resnorm at zero %e\n',resnorms(1));
%   fprintf(fid,'resnorm at finish %e\n\n',resnorms(ll));
%   fclose(fid);
%end
%
%deltasin=0;
%kcount=2;
%fid=fopen('jacobi_sync_pconst_nodelta.results','a');
%fprintf(fid,'kcount=%g\n',kcount);
%fclose(fid);
%for iiii=1:6
%   fid=fopen('jacobi_sync_pconst_nodelta.results','a');
%   matfile=iiii;
%   jacobi_sync_pconst_planar
%   [hh,ll]=size(resnorms);
%   ii=iter;
%   numplan=iter/kcount;
%   fprintf(fid,'matfile=%g\n',matfile);
%   fprintf(fid,'number of iterations %g\n',ii);
%   fprintf(fid,'number of planar searches %g\n',numplan);
%   fprintf(fid,'resnorm at zero %e\n',resnorms(1));
%   fprintf(fid,'resnorm at finish %e\n\n',resnorms(ll));
%   fclose(fid);
%end
%
deltasin=0;
kcount=2;
for iiii=8:12
   fid=fopen('jacobi_sync_pconst_nodelta.results','a');
   matfile=iiii;
   jacobi_sync_pconst_planar
   [hh,ll]=size(resnorms);
   ii=iter;
   numplan=iter/kcount;
   fprintf(fid,'matfile=%g\n',matfile);
   fprintf(fid,'number of iterations %g\n',ii);
   fprintf(fid,'number of planar searches %g\n',numplan);
   fprintf(fid,'resnorm at zero %e\n',resnorms(1));
   fprintf(fid,'resnorm at finish %e\n\n',resnorms(ll));
   fclose(fid);
end
%
deltasin=0;
kcount=3;
fid=fopen('jacobi_sync_pconst_nodelta.results','a');
fprintf(fid,'kcount=%g\n',kcount);
fclose(fid);
for iiii=1:6
   fid=fopen('jacobi_sync_pconst_nodelta.results','a');
   matfile=iiii;
   jacobi_sync_pconst_planar
   [hh,ll]=size(resnorms);
   ii=iter;
   numplan=iter/kcount;
   fprintf(fid,'matfile=%g\n',matfile);
   fprintf(fid,'number of iterations %g\n',ii);
   fprintf(fid,'number of planar searches %g\n',numplan);
   fprintf(fid,'resnorm at zero %e\n',resnorms(1));
   fprintf(fid,'resnorm at finish %e\n\n',resnorms(ll));
   fclose(fid);
end
%
deltasin=0;
kcount=3;
for iiii=7:12
   fid=fopen('jacobi_sync_pconst_nodelta.results','a');
   matfile=iiii;
   jacobi_sync_pconst_planar
   [hh,ll]=size(resnorms);
   ii=iter;
   numplan=iter/kcount;
   fprintf(fid,'matfile=%g\n',matfile);
   fprintf(fid,'number of iterations %g\n',ii);
   fprintf(fid,'number of planar searches %g\n',numplan);
   fprintf(fid,'resnorm at zero %e\n',resnorms(1));
   fprintf(fid,'resnorm at finish %e\n\n',resnorms(ll));
   fclose(fid);
end
%
deltasin=0;
kcount=4;
fid=fopen('jacobi_sync_pconst_nodelta.results','a');
fprintf(fid,'kcount=%g\n',kcount);
fclose(fid);
for iiii=1:6
   fid=fopen('jacobi_sync_pconst_nodelta.results','a');
   matfile=iiii;
   jacobi_sync_pconst_planar
   [hh,ll]=size(resnorms);
   ii=iter;
   numplan=iter/kcount;
   fprintf(fid,'matfile=%g\n',matfile);
   fprintf(fid,'number of iterations %g\n',ii);
   fprintf(fid,'number of planar searches %g\n',numplan);
   fprintf(fid,'resnorm at zero %e\n',resnorms(1));
   fprintf(fid,'resnorm at finish %e\n\n',resnorms(ll));
   fclose(fid);
end
%
deltasin=0;
kcount=4;
for iiii=7:12
   fid=fopen('jacobi_sync_pconst_nodelta.results','a');
   matfile=iiii;
   jacobi_sync_pconst_planar
   [hh,ll]=size(resnorms);
   ii=iter;
   numplan=iter/kcount;
   fprintf(fid,'matfile=%g\n',matfile);
   fprintf(fid,'number of iterations %g\n',ii);
   fprintf(fid,'number of planar searches %g\n',numplan);
   fprintf(fid,'resnorm at zero %e\n',resnorms(1));
   fprintf(fid,'resnorm at finish %e\n\n',resnorms(ll));
   fclose(fid);
end
%
deltasin=0;
kcount=5;
fid=fopen('jacobi_sync_pconst_nodelta.results','a');
fprintf(fid,'kcount=%g\n',kcount);
fclose(fid);
for iiii=1:6
   fid=fopen('jacobi_sync_pconst_nodelta.results','a');
   matfile=iiii;
   jacobi_sync_pconst_planar
   [hh,ll]=size(resnorms);
   ii=iter;
   numplan=iter/kcount;
   fprintf(fid,'matfile=%g\n',matfile);
   fprintf(fid,'number of iterations %g\n',ii);
   fprintf(fid,'number of planar searches %g\n',numplan);
   fprintf(fid,'resnorm at zero %e\n',resnorms(1));
   fprintf(fid,'resnorm at finish %e\n\n',resnorms(ll));
   fclose(fid);
end
%
deltasin=0;
kcount=5;
for iiii=7:12
   fid=fopen('jacobi_sync_pconst_nodelta.results','a');
   matfile=iiii;
   jacobi_sync_pconst_planar
   [hh,ll]=size(resnorms);
   ii=iter;
   numplan=iter/kcount;
   fprintf(fid,'matfile=%g\n',matfile);
   fprintf(fid,'number of iterations %g\n',ii);
   fprintf(fid,'number of planar searches %g\n',numplan);
   fprintf(fid,'resnorm at zero %e\n',resnorms(1));
   fprintf(fid,'resnorm at finish %e\n\n',resnorms(ll));
   fclose(fid);
end
%
deltasin=0;
kcount=6;
fid=fopen('jacobi_sync_pconst_nodelta.results','a');
fprintf(fid,'kcount=%g\n',kcount);
fclose(fid);
for iiii=1:6
   fid=fopen('jacobi_sync_pconst_nodelta.results','a');
   matfile=iiii;
   jacobi_sync_pconst_planar
   [hh,ll]=size(resnorms);
   ii=iter;
   numplan=iter/kcount;
   fprintf(fid,'matfile=%g\n',matfile);
   fprintf(fid,'number of iterations %g\n',ii);
   fprintf(fid,'number of planar searches %g\n',numplan);
   fprintf(fid,'resnorm at zero %e\n',resnorms(1));
   fprintf(fid,'resnorm at finish %e\n\n',resnorms(ll));
   fclose(fid);
end
%
deltasin=0;
kcount=6;
for iiii=7:12
   fid=fopen('jacobi_sync_pconst_nodelta.results','a');
   matfile=iiii;
   jacobi_sync_pconst_planar
   [hh,ll]=size(resnorms);
   ii=iter;
   numplan=iter/kcount;
   fprintf(fid,'matfile=%g\n',matfile);
   fprintf(fid,'number of iterations %g\n',ii);
   fprintf(fid,'number of planar searches %g\n',numplan);
   fprintf(fid,'resnorm at zero %e\n',resnorms(1));
   fprintf(fid,'resnorm at finish %e\n\n',resnorms(ll));
   fclose(fid);
end
%
%
%
%
deltasin=1;
kcount=1;
fid=fopen('jacobi_sync_pconst_delta.results','a');
fprintf(fid,'kcount=%g\n',kcount);
fclose(fid);
for iiii=1:6
   fid=fopen('jacobi_sync_pconst_delta.results','a');
   matfile=iiii;
   jacobi_sync_pconst_planar
   [hh,ll]=size(resnorms);
   ii=iter;
   numplan=iter/kcount;
   fprintf(fid,'matfile=%g\n',matfile);
   fprintf(fid,'number of iterations %g\n',ii);
   fprintf(fid,'number of planar searches %g\n',numplan);
   fprintf(fid,'resnorm at zero %e\n',resnorms(1));
   fprintf(fid,'resnorm at finish %e\n\n',resnorms(ll));
   fclose(fid);
end
%
deltasin=1;
kcount=1;
for iiii=7:12
   fid=fopen('jacobi_sync_pconst_delta.results','a');
   matfile=iiii;
   jacobi_sync_pconst_planar
   [hh,ll]=size(resnorms);
   ii=iter;
   numplan=iter/kcount;
   fprintf(fid,'matfile=%g\n',matfile);
   fprintf(fid,'number of iterations %g\n',ii);
   fprintf(fid,'number of planar searches %g\n',numplan);
   fprintf(fid,'resnorm at zero %e\n',resnorms(1));
   fprintf(fid,'resnorm at finish %e\n\n',resnorms(ll));
   fclose(fid);
end
%
deltasin=1;
kcount=2;
fid=fopen('jacobi_sync_pconst_delta.results','a');
fprintf(fid,'kcount=%g\n',kcount);
fclose(fid);
for iiii=1:6
   fid=fopen('jacobi_sync_pconst_delta.results','a');
   matfile=iiii;
   jacobi_sync_pconst_planar
   [hh,ll]=size(resnorms);
   ii=iter;
   numplan=iter/kcount;
   fprintf(fid,'matfile=%g\n',matfile);
   fprintf(fid,'number of iterations %g\n',ii);
   fprintf(fid,'number of planar searches %g\n',numplan);
   fprintf(fid,'resnorm at zero %e\n',resnorms(1));
   fprintf(fid,'resnorm at finish %e\n\n',resnorms(ll));
   fclose(fid);
end
%
deltasin=1;
kcount=2;
for iiii=7:12
   fid=fopen('jacobi_sync_pconst_delta.results','a');
   matfile=iiii;
   jacobi_sync_pconst_planar
   [hh,ll]=size(resnorms);
   ii=iter;
   numplan=iter/kcount;
   fprintf(fid,'matfile=%g\n',matfile);
   fprintf(fid,'number of iterations %g\n',ii);
   fprintf(fid,'number of planar searches %g\n',numplan);
   fprintf(fid,'resnorm at zero %e\n',resnorms(1));
   fprintf(fid,'resnorm at finish %e\n\n',resnorms(ll));
   fclose(fid);
end
%
deltasin=1;
kcount=3;
fid=fopen('jacobi_sync_pconst_delta.results','a');
fprintf(fid,'kcount=%g\n',kcount);
fclose(fid);
for iiii=1:6
   fid=fopen('jacobi_sync_pconst_delta.results','a');
   matfile=iiii;
   jacobi_sync_pconst_planar
   [hh,ll]=size(resnorms);
   ii=iter;
   numplan=iter/kcount;
   fprintf(fid,'matfile=%g\n',matfile);
   fprintf(fid,'number of iterations %g\n',ii);
   fprintf(fid,'number of planar searches %g\n',numplan);
   fprintf(fid,'resnorm at zero %e\n',resnorms(1));
   fprintf(fid,'resnorm at finish %e\n\n',resnorms(ll));
   fclose(fid);
end
%
deltasin=1;
kcount=3;
for iiii=7:12
   fid=fopen('jacobi_sync_pconst_delta.results','a');
   matfile=iiii;
   jacobi_sync_pconst_planar
   [hh,ll]=size(resnorms);
   ii=iter;
   numplan=iter/kcount;
   fprintf(fid,'matfile=%g\n',matfile);
   fprintf(fid,'number of iterations %g\n',ii);
   fprintf(fid,'number of planar searches %g\n',numplan);
   fprintf(fid,'resnorm at zero %e\n',resnorms(1));
   fprintf(fid,'resnorm at finish %e\n\n',resnorms(ll));
   fclose(fid);
end
%
deltasin=1;
kcount=4;
fid=fopen('jacobi_sync_pconst_delta.results','a');
fprintf(fid,'kcount=%g\n',kcount);
fclose(fid);
for iiii=1:6
   fid=fopen('jacobi_sync_pconst_delta.results','a');
   matfile=iiii;
   jacobi_sync_pconst_planar
   [hh,ll]=size(resnorms);
   ii=iter;
   numplan=iter/kcount;
   fprintf(fid,'matfile=%g\n',matfile);
   fprintf(fid,'number of iterations %g\n',ii);
   fprintf(fid,'number of planar searches %g\n',numplan);
   fprintf(fid,'resnorm at zero %e\n',resnorms(1));
   fprintf(fid,'resnorm at finish %e\n\n',resnorms(ll));
   fclose(fid);
end
%
deltasin=1;
kcount=4;
for iiii=7:12
   fid=fopen('jacobi_sync_pconst_delta.results','a');
   matfile=iiii;
   jacobi_sync_pconst_planar
   [hh,ll]=size(resnorms);
   ii=iter;
   numplan=iter/kcount;
   fprintf(fid,'matfile=%g\n',matfile);
   fprintf(fid,'number of iterations %g\n',ii);
   fprintf(fid,'number of planar searches %g\n',numplan);
   fprintf(fid,'resnorm at zero %e\n',resnorms(1));
   fprintf(fid,'resnorm at finish %e\n\n',resnorms(ll));
   fclose(fid);
end
%
deltasin=1;
kcount=5;
fid=fopen('jacobi_sync_pconst_delta.results','a');
fprintf(fid,'kcount=%g\n',kcount);
fclose(fid);
for iiii=1:6
   fid=fopen('jacobi_sync_pconst_delta.results','a');
   matfile=iiii;
   jacobi_sync_pconst_planar
   [hh,ll]=size(resnorms);
   ii=iter;
   numplan=iter/kcount;
   fprintf(fid,'matfile=%g\n',matfile);
   fprintf(fid,'number of iterations %g\n',ii);
   fprintf(fid,'number of planar searches %g\n',numplan);
   fprintf(fid,'resnorm at zero %e\n',resnorms(1));
   fprintf(fid,'resnorm at finish %e\n\n',resnorms(ll));
   fclose(fid);
end
%
deltasin=1;
kcount=5;
for iiii=7:12
   fid=fopen('jacobi_sync_pconst_delta.results','a');
   matfile=iiii;
   jacobi_sync_pconst_planar
   [hh,ll]=size(resnorms);
   ii=iter;
   numplan=iter/kcount;
   fprintf(fid,'matfile=%g\n',matfile);
   fprintf(fid,'number of iterations %g\n',ii);
   fprintf(fid,'number of planar searches %g\n',numplan);
   fprintf(fid,'resnorm at zero %e\n',resnorms(1));
   fprintf(fid,'resnorm at finish %e\n\n',resnorms(ll));
   fclose(fid);
end
%
deltasin=1;
kcount=6;
fid=fopen('jacobi_sync_pconst_delta.results','a');
fprintf(fid,'kcount=%g\n',kcount);
fclose(fid);
for iiii=1:6
   fid=fopen('jacobi_sync_pconst_delta.results','a');
   matfile=iiii;
   jacobi_sync_pconst_planar
   [hh,ll]=size(resnorms);
   ii=iter;
   numplan=iter/kcount;
   fprintf(fid,'matfile=%g\n',matfile);
   fprintf(fid,'number of iterations %g\n',ii);
   fprintf(fid,'number of planar searches %g\n',numplan);
   fprintf(fid,'resnorm at zero %e\n',resnorms(1));
   fprintf(fid,'resnorm at finish %e\n\n',resnorms(ll));
   fclose(fid);
end
%
deltasin=1;
kcount=6;
for iiii=7:12
   fid=fopen('jacobi_sync_pconst_delta.results','a');
   matfile=iiii;
   jacobi_sync_pconst_planar
   [hh,ll]=size(resnorms);
   ii=iter;
   numplan=iter/kcount;
   fprintf(fid,'matfile=%g\n',matfile);
   fprintf(fid,'number of iterations %g\n',ii);
   fprintf(fid,'number of planar searches %g\n',numplan);
   fprintf(fid,'resnorm at zero %e\n',resnorms(1));
   fprintf(fid,'resnorm at finish %e\n\n',resnorms(ll));
   fclose(fid);
end
%
%
%
%
deltasin=0;
kcount=1;
fid=fopen('jacobi_sync_pds_nodelta.results','a');
fprintf(fid,'kcount=%g\n',kcount);
fclose(fid);
for iiii=1:6
   fid=fopen('jacobi_sync_pds_nodelta.results','a');
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
   fclose(fid);
end
%
deltasin=0;
kcount=1;
for iiii=7:12
   fid=fopen('jacobi_sync_pds_nodelta.results','a');
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
   fclose(fid);
end
%
deltasin=0;
kcount=2;
fid=fopen('jacobi_sync_pds_nodelta.results','a');
fprintf(fid,'kcount=%g\n',kcount);
fclose(fid);
for iiii=1:6
   fid=fopen('jacobi_sync_pds_nodelta.results','a');
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
   fclose(fid);
end
%
deltasin=0;
kcount=2;
for iiii=7:12
   fid=fopen('jacobi_sync_pds_nodelta.results','a');
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
   fclose(fid);
end
%
deltasin=0;
kcount=3;
fid=fopen('jacobi_sync_pds_nodelta.results','a');
fprintf(fid,'kcount=%g\n',kcount);
fclose(fid);
for iiii=1:6
   fid=fopen('jacobi_sync_pds_nodelta.results','a');
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
   fclose(fid);
end
%
deltasin=0;
kcount=3;
for iiii=7:12
   fid=fopen('jacobi_sync_pds_nodelta.results','a');
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
   fclose(fid);
end
%
deltasin=0;
kcount=4;
fid=fopen('jacobi_sync_pds_nodelta.results','a');
fprintf(fid,'kcount=%g\n',kcount);
fclose(fid);
for iiii=1:6
   fid=fopen('jacobi_sync_pds_nodelta.results','a');
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
   fclose(fid);
end
%
deltasin=0;
kcount=4;
for iiii=7:12
   fid=fopen('jacobi_sync_pds_nodelta.results','a');
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
   fclose(fid);
end
%
deltasin=0;
kcount=5;
fid=fopen('jacobi_sync_pds_nodelta.results','a');
fprintf(fid,'kcount=%g\n',kcount);
fclose(fid);
for iiii=1:6
   fid=fopen('jacobi_sync_pds_nodelta.results','a');
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
   fclose(fid);
end
%
deltasin=0;
kcount=5;
for iiii=7:12
   fid=fopen('jacobi_sync_pds_nodelta.results','a');
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
   fclose(fid);
end
%
deltasin=0;
kcount=6;
fid=fopen('jacobi_sync_pds_nodelta.results','a');
fprintf(fid,'kcount=%g\n',kcount);
fclose(fid);
for iiii=1:6
   fid=fopen('jacobi_sync_pds_nodelta.results','a');
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
   fclose(fid);
end
%
deltasin=0;
kcount=6;
for iiii=7:12
   fid=fopen('jacobi_sync_pds_nodelta.results','a');
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
   fclose(fid);
end
%
%
%
%
deltasin=1;
kcount=1;
fid=fopen('jacobi_sync_pds_delta.results','a');
fprintf(fid,'kcount=%g\n',kcount);
fclose(fid);
for iiii=1:6
   fid=fopen('jacobi_sync_pds_delta.results','a');
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
   fclose(fid);
end
%
deltasin=1;
kcount=1;
for iiii=7:12
   fid=fopen('jacobi_sync_pds_delta.results','a');
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
   fclose(fid);
end
%
deltasin=1;
kcount=2;
fid=fopen('jacobi_sync_pds_delta.results','a');
fprintf(fid,'kcount=%g\n',kcount);
fclose(fid);
for iiii=1:6
   fid=fopen('jacobi_sync_pds_delta.results','a');
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
   fclose(fid);
end
%
deltasin=1;
kcount=2;
for iiii=7:12
   fid=fopen('jacobi_sync_pds_delta.results','a');
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
   fclose(fid);
end
%
deltasin=1;
kcount=3;
fid=fopen('jacobi_sync_pds_delta.results','a');
fprintf(fid,'kcount=%g\n',kcount);
fclose(fid);
for iiii=1:6
   fid=fopen('jacobi_sync_pds_delta.results','a');
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
   fclose(fid);
end
%
deltasin=1;
kcount=3;
for iiii=7:12
   fid=fopen('jacobi_sync_pds_delta.results','a');
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
   fclose(fid);
end
%
deltasin=1;
kcount=4;
fid=fopen('jacobi_sync_pds_delta.results','a');
fprintf(fid,'kcount=%g\n',kcount);
fclose(fid);
for iiii=1:6
   fid=fopen('jacobi_sync_pds_delta.results','a');
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
   fclose(fid);
end
%
deltasin=1;
kcount=4;
for iiii=7:12
   fid=fopen('jacobi_sync_pds_delta.results','a');
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
   fclose(fid);
end
%
deltasin=1;
kcount=5;
fid=fopen('jacobi_sync_pds_delta.results','a');
fprintf(fid,'kcount=%g\n',kcount);
fclose(fid);
for iiii=1:6
   fid=fopen('jacobi_sync_pds_delta.results','a');
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
   fclose(fid);
end
%
deltasin=1;
kcount=5;
for iiii=7:12
   fid=fopen('jacobi_sync_pds_delta.results','a');
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
   fclose(fid);
end
%
deltasin=1;
kcount=6;
fid=fopen('jacobi_sync_pds_delta.results','a');
fprintf(fid,'kcount=%g\n',kcount);
fclose(fid);
for iiii=1:6
   fid=fopen('jacobi_sync_pds_delta.results','a');
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
   fclose(fid);
end
%
deltasin=1;
kcount=6;
for iiii=7:12
   fid=fopen('jacobi_sync_pds_delta.results','a');
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
   fclose(fid);
end

