function crater_5eqns_eta_movie(N)           
km = 1e3;
n0 = 10000;
fname = ['fort.',num2str(n0)];
fname(6) = 'g';
fid = fopen(fname);
mx = fscanf(fid,'%d',1);  
my = fscanf(fid,'%d',1);     
bgrid = fscanf(fid,'%g %g',[2 inf]);
status = fclose(fid);
bgrid = bgrid';
%
x = reshape(bgrid(:,1),mx,my);
y = reshape(bgrid(:,2),mx,my);
%x = x/1.0e3;
%y = y/1.0e3;
%
for j=0:N
    clf

    n1 = j+10000;
    fname = ['fort.',num2str(n1)];
    fname(6) = 't';
    fid  = fopen(fname);
    t1   = fscanf(fid,'%g',1);      fscanf(fid,'%s',1);
    meqn = fscanf(fid,'%d',1);      fscanf(fid,'%s',1);
    ngrids = fscanf(fid,'%d',1);    fscanf(fid,'%s',1);
    fclose(fid);
%
    fname(6) = 'c';
    fid    = fopen(fname);
    data_5eqns  = fscanf(fid,'%g',[3 inf]);
    status = fclose(fid);
    data_5eqns = data_5eqns';
    size(data_5eqns)
%
    plot(data_5eqns(:,1)/km,data_5eqns(:,3),'b-',...
         'LineWidth',1)
%
    title(['time t=', num2str(t1),' seconds after impact ($DC=1000$m)'],...
       'fontsize',20,'interpreter','latex')
%
    legend('$2$-phase flow solution',...
       'fontsize',20,'interpreter','latex',...
       'Location','NorthWest',...
       'box','off')
%
    ylim([-50 50])
%
    xlabel('radial distance (km)','fontsize',20,'interpreter','latex')
    ylabel('water height (m)','fontsize',20,'interpreter','latex')
    set(gca,'TickLabelInterpreter','latex',...
        'fontsize',20)
%
    grid on
    F(j+1) = getframe(gcf);
end
%
v = VideoWriter('crater_5eqns_DC1000_eta.mp4','MPEG-4');
open(v);
writeVideo(v,F);
close(v);
end
