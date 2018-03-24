%%
datafile='asl_raw_data_corrected_mcf';
addpath /opt/fmrib/fsl/etc/matlab
[y,~,scales]=read_avw(datafile);
m=read_avw('asl_mask2');
y=reshape(y,numel(m),size(y,4));
y=y(~~m,:);
yd=detrend(y')';
ystd=std(yd,0,2);
x=zeros(size(m));
x(~~m)=ystd;
save_avw(x,'data_std','f',scales);
unix(['fslcpgeom ' datafile ' data_std -d']);
