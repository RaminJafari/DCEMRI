function [A,B] = cg_A_B_mats(aorta_enhancement,pv_enhancement,enhancement,time_sim)
    
t_s = time_sim(1,1);
[p,l] = size(enhancement);
m = sqrt(l);
n = sqrt(l);
p=p-1;

enhancement(enhancement <= 0 )= 1e-6;
aorta_enhancement(aorta_enhancement <=0) = 1e-6;
pv_enhancement(pv_enhancement <=0) = 1e-6;

b = zeros(m*n*p,1);
cols = zeros(m*n*3*p,1);
rows = zeros(m*n*3*p,1);
vals = zeros(m*n*3*p,1);       
    
start = -1;
r = 0;
ind = 1;     
for c = 1:3:m*n*3;
        start = start + 1; 
        r = 1 + start*p;

        
         cols(ind:ind+p-1) = c;
         rows(ind:ind+p-1) = [r:r+p-1];
         vals(ind:ind+p-1) = cumtrapz(aorta_enhancement(1:end-1,start+1));
         
         ind = ind+p;
         
         cols(ind:ind+p-1) = c+1;
         rows(ind:ind+p-1) = [r:r+p-1];
         vals(ind:ind+p-1) = cumtrapz(pv_enhancement(1:end-1,start+1));
      
         ind = ind+p;
         
         cols(ind:ind+p-1) = c+2;
         rows(ind:ind+p-1) = [r:r+p-1];
         vals(ind:ind+p-1) = cumtrapz(enhancement(1:end-1,start+1));
         
         ind = ind+p;

         b(r:r+p-1,1) = enhancement(1:end-1,start+1);
end

a = sparse(rows,cols,vals,m*n*p,3*m*n,m*n*3*p);
   
b(b==0)= 1e-6;
b(isnan(b))= 1e-6;
a(isinf(a))= 1e-6;
a(isnan(a))= 1e-6;
    
A=a'*a;
B=a'*b;   
end

