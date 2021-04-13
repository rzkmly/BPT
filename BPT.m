clc;clear;
%% Forward Modeling

% gridding
x=linspace(0,100,6);y=linspace(0,100,6);
dx=abs(x(2)-x(1));dy=abs(y(2)-y(1));
xt=(min(x)+dx/2:dx:max(x)-dx/2);
yt=(min(y)+dy/2:dy:max(y)-dy/2);
nx=length(x)-1;
ny=length(y)-1;
kec=[1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 2000 2000 2000 2000 2000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000]';
V=reshape(kec,[nx ny])';
figure()
imagesc(xt,yt,V)
set(gca,'YDir','reverse')
title('Model Kecepatan & Ray path')
xlabel(colorbar(), 'V (m/s)')
hold on

% posisi source dan receiver
rs=[zeros(ny,1) yt'];
rr=[ones(length(yt),1)*max(x) yt'];
L=zeros(1,nx*ny);
k=1;
for i=1:length(rs(:,1))
    xs=rs(i,1);
    ys=rs(i,2);
    for j=1:length(rr(:,1))
        xr=rr(j,1);
        yr=rr(j,2);
        yx=(yr-ys)/(xr-xs);
        cond=1;
        x1=xs;
        y1=ys;
        plot([xs xr],[ys yr],'k--')
        p1=plot(xs,ys,'*r','linewidth',1,'markerfacecolor','r')
        p2=plot(xr,yr,'vb','linewidth',1,'markerfacecolor','b')
        while cond==1
            ii=floor(abs((x1-min(x))/dx))+1;
            jj=floor(abs((y1-min(y))/dy))+1;
            if x1==max(x);
                ii=ii-1;
            elseif y1==max(y);
                jj=jj-1;
            end
            blok=ii+(jj-1)*nx;
            batasx=x(abs(x-x1)>0 & abs(x-x1)<abs(xr-x1) & abs(xr-x)<abs(xr-x1));
            batasy=y(abs(y-y1)>0 & abs(y-y1)<abs(yr-y1) & abs(yr-y)<abs(yr-y1));
            x_next=batasx(abs(batasx-x1)==min(abs(batasx-x1)) & abs(batasx-x1)>0); 
            y_next=batasy(abs(batasy-y1)==min(abs(batasy-y1)) & abs(batasy-y1)>0);
            if length(x_next)==1 && length(y_next)==1
                x_cal=x1+(y_next-y1)/yx;
                y_cal=y1+yx*(x_next-x1);
                if (abs(x_cal-x1) <= abs(x_next-x1)) 
                    l=sqrt((x_cal-x1)^2+(y_next-y1)^2);
                    x1=x_cal;y1=y_next;
                else
                    l=sqrt((x_next-x1)^2+(y_cal-y1)^2);
                    x1=x_next;y1=y_cal;
                end
            elseif length(x_next)==1 && isempty(y_next)
                y_cal=y1+yx*(x_next-x1);
                l=sqrt((x_next-x1)^2+(y_cal-y1)^2);
                x1=x_next;y1=y_cal;
            elseif isempty(x_next) && length(y_next)==1
                x_cal=x1+(y_next-y1)/yx;
                l=sqrt((x_cal-x1)^2+(y_next-y1)^2);
                x1=x_cal;y1=y_next;
            else
                l=sqrt((xr-x1)^2+(yr-y1)^2);
                cond=0;
            end
            plot(x1,y1,'or','MarkerSize',3)
            L(k,blok)=l;
            x_next=[];y_next=[];
        end
        k=k+1;   
    end
end
legend([p1,p2],'source','receiver','location','southeast')
hold off
Tobs=L*(1./kec);
nray=k-1;

%% Back Projection 
S_avg=zeros(nray,1);wL=L*0;
for i=1:nray
    S_avg(i)=Tobs(i)/sum(L(i,:));
    wL(i,:)=L(i,:)/sum(L(i,:));
end
S_i=zeros(length(kec),1);
n=100;
for i=1:length(kec)
    S_i(i)=(sum(L(:,i).*S_avg)/sum(L(:,i)));
end
Tcal=L*S_i;
Vcal=reshape(1./S_i,[nx ny])';
figure()
imagesc(xt,yt,Vcal)
set(gca,'YDir','reverse')
title('Model Kecepatan BPT')
xlabel(colorbar(), 'V (m/s)')

%% Inversion
Sinv=(L'*L)\L'*Tobs;
Vinv=reshape(1./Sinv,[nx ny])';
figure()
imagesc(xt,yt,Vinv)
set(gca,'YDir','reverse')
title('Hasil Inversi Kecepatan')
xlabel(colorbar(), 'V (m/s)')

%% ERROR
error=sqrt((Vinv-V).^2)
percent=error./max(error)*100
figure()
imagesc(xt,yt,percent)
set(gca,'YDir','reverse')
title('Error RMS dalam Persen')
xlabel(colorbar(), 'V (m/s)')
