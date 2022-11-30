%% test hsbound
E=[];
s=[];
for i=0:0.001:1
    [b,s]=hsbound(1e6,0.3,1-i);
    E(end+1)=b(1);
end
figure;plot(E);
%%
