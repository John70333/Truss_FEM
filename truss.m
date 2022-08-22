clc; clear; nD = 2;
ui1 = [0,0]; ui2 = [2.5,2.5*3^.5]; ui3 = [5,0]; ui4 = [7.5,2.5*3^.5]; ui5 = [10,0];
ui = [ui1,ui2,ui3,ui4,ui5]'; conn = {[2,3],[1,3,4],[1,2,4,5],[2,3,5],[3,4]}';
Ra = [0,-10,0,-10,0,-10,0]'*1e3; ai = 3:9; ub = [0,0,0]'; bi = [1,2,10];
A = zeros(length(conn)); E = A; nu = .3; ns = 2;
for i = 1:length(conn)
    for j = conn{i}
        A(i,j) = 1e-4;
        E(i,j) = 1e10;
    end
end
[R,K,u,uf,lci,lcf,eps] = tfem(nD,ui,conn,A,E,nu,Ra,ai,ub,bi,ns);
draw(nD,ui,uf(:,end),conn);
%%
function [R,K,u,uf,lci,lcf,eps] = tfem(nD,ui,conn,A,E,nu,Ra,ai,ub,bi,ns)
n = length(conn); uc = ui; Ras = Ra./ns;
g = @(i)nD*(i-1)+1:nD*i; l = @(u)sum(u.^2)^.5;
a = {@(u)1,@(u)atan2(u(2),u(1)),@(u)[atan2(u(2),u(1)), atan2((u(1).^2+u(2).^2)^.5,u(3))]};
c = {@(a)1,@(a)[cos(a),sin(a)],@(a)[cos(a(1)).*sin(a(2)),sin(a(1)).*sin(a(2)),cos(a(2))]};
R = zeros(nD*n,1,ns); K = R; u = R; uf = R; lci = zeros(n,n,ns); lcf = lci; eps = lci;
for pg = 1:ns
    for i = 1:n
        s = zeros(nD);
        for j = conn{i}
            lc = l(uc(g(j))-uc(g(i)));
            ca = c{nD}(a{nD}(uc(g(j))-uc(g(i))));
            if pg == 1
                lci(i,j) = lc; pac = 1;
            else
                pac = (1-nu*eps(i,j,pg-1))^2;
            end
            k = E(i,j)*A(i,j)*pac/lc*(ca'*ca);
            K(g(i),g(j),pg) = -k; s = s+k;
        end
        K(g(i),g(i),pg) = s;
    end
    Kaa = K(ai,ai,pg); Kab = K(ai,bi,pg); Kba = K(bi,ai,pg); Kbb = K(bi,bi,pg);
    ua = Kaa\(Ras-Kab*ub); Rb = Kba*ua+Kbb*ub;
    u(ai,:,pg) = ua; u(bi,:,pg) = ub; R(ai,:,pg) = Ras; R(bi,:,pg) = Rb;
    uf(:,:,pg) = uc+u(:,:,pg); uc = uf(:,:,pg);
    for i = 1:n
        for j = conn{i}
            lcf(i,j,pg) = l(uf(g(i),:,pg)-uf(g(j),:,pg));
            eps(i,j,pg) = lcf(i,j,pg)./lci(i,j)-1;
        end
    end
end
end
%%
function draw(nD,ui,uf,conn)
figure('windowstate','maximized'); hold on; grid on;
xi = ui(1:nD:end); xf = uf(1:nD:end,:,end);
xlim([min([xi;xf])-1,max([xi;xf])+1]);
if nD >= 2
    yi = ui(2:nD:end); yf = uf(2:nD:end,:,end);
    ylim([min([yi;yf])-1,max([yi;yf])+1]);
    if nD == 3
        zi = ui(3:nD:end); zf = uf(3:nD:end,:,end);
        zlim([min([zi;zf])-1,max([zi;zf])+1]); view(3);
    end
end
for i = 1:length(conn)
    for j = conn{i}
        if j < i
            switch nD
                case 1
                    plot([xi(i),xi(j)],[0,0],'b-');
                    plot([xf(i),xf(j)],[0,0],'r-');
                case 2
                    plot([xi(i),xi(j)],[yi(i),yi(j)],'b-');
                    plot([xf(i),xf(j)],[yf(i),yf(j)],'r-');
                case 3
                    plot3([xi(i),xi(j)],[yi(i),yi(j)],[zi(i),zi(j)],'b-');
                    plot3([xf(i),xf(j)],[yf(i),yf(j)],[zf(i),zf(j)],'r-');
            end
        end
    end
end
saveas(gcf,'truss.png');
end
