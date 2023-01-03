function [bulk, shear] = hsbound(youngs,psratio,vol)
%     bulks=[youngs,youngs*1e-9];
%     shears=[youngs/(2*(1+psratio)),youngs*1e-9/(2*(1+psratio))];
%     [bulk,shear]=hsbound_imp(bulks,shears,vol);
      [c,s]=CS(youngs,psratio);
      K = sum(sum(c(1:3,1:3)))/9;
      G = c(4,4);
      bulk = 4*K*G*vol/(3*K*(1-vol)+4*G);
      shear=vol*G*(9*K+8*G)/(K*(15-6*vol)+G*(20-16*vol));
end

function [bulk, shear] = hsbound_imp(bulks, shears,ratio0)
    bulk=[];
    shear=[];
    phi=1-ratio0;
    k2=bulks(1);k1=bulks(2);
    mu2=shears(1);mu1=shears(2);
    bulk(1)=k2+phi/((k1-k2).^(-1)+(1-phi)*(k2+4/3*mu2).^(-1));
    shear(1)=mu2+phi/((mu1-mu2)^(-1)+2*(1-phi)*(k2+2*mu2)/(5*mu2*(k2+4/3*mu2)));
    k2=bulks(2);k1=bulks(1);
    mu2=shears(2);mu1=shears(1);
    bulk(2)=k2+phi/((k1-k2).^(-1)+(1-phi)*(k2+4/3*mu2).^(-1));
    shear(2)=mu2+phi/((mu1-mu2)^(-1)+2*(1-phi)*(k2+2*mu2)/(5*mu2*(k2+4/3*mu2)));
end