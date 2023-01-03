function [c,s]=CS(E,mu)
    s = [1/E, -mu/E, -mu/E, 0, 0, 0;
         -mu/E, 1/E, -mu/E, 0, 0, 0;
         -mu/E,-mu/E,  1/E, 0, 0, 0;
         0, 0, 0, 2*(1+mu)/E,0, 0;
         0, 0, 0, 0,2*(1+mu)/E, 0;
         0, 0, 0, 0, 0,2*(1+mu)/E;];
    c=inv(s);
end