%%
datafile='asl_raw_data';
addpath /opt/fmrib/fsl/etc/matlab
[y,~,scales]=read_avw(datafile);
m=read_avw('asl_mask1');
y=reshape(y,numel(m),size(y,4));
y=y(~~m,:);
yd=detrend(y')';
ystd=std(yd,0,2);
x=zeros(size(m));
x(~~m)=ystd;
save_avw(x,'data_std1','f',scales);
unix(['fslcpgeom ' datafile ' data_std1 -d']);
