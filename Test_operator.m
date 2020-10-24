clear all
m=4;
n=2^m;
X=[0,1;1,0];
Y=[0,-1i;1i,0];
Z=[1,0;0,-1];
I=eye(2);
II=eye(n);
S={X Y Z I};
for ii=1:(m-1)      %generate R gate
    R{2*ii-1}=1;
    R{2*ii}=1;
    for jj=1:m
        if jj==ii
            R{2*ii-1}=kron(R{2*ii-1},S{2});
            R{2*ii}=kron(R{2*ii},S{2});
        else
            if jj==ii+1
                R{2*ii-1}=kron(R{2*ii-1},S{2});
                R{2*ii}=kron(R{2*ii},S{1});
            else
                R{2*ii-1}=kron(R{2*ii-1},S{4});
                R{2*ii}=kron(R{2*ii},S{4});
            end
        end
    end
    R{2*ii-1}=1/sqrt(2)*(II+1i*R{2*ii-1});
    R{2*ii}=1/sqrt(2)*(II-1i*R{2*ii});
end
     
U{1,1}=II;    %generate U gate
for ii=1:m-1
    for jj=1:3^(ii-1)
        U{ii+1,3*jj-2}=U{ii,jj};
        U{ii+1,3*jj-1}=U{ii,jj}*inv(R{2*ii-1});
        U{ii+1,3*jj}=U{ii,jj}*inv(R{2*ii});
    end
end

for ii=0:2^m-1   %initial operators
    Q{ii+1}=1;
    for jj=1:m
        n(jj)=mod(floor(ii/2^(m-jj)),2);
        Q{ii+1}=kron(Q{ii+1},(-Z)^n(jj));
    end
end
 
for ii=1:2^m      %generate operators
    for jj=1:3^(m-1)
        UQ{3^(m-1)*(ii-1)+jj}=U{m,jj}'*Q{ii}*U{m,jj};
    end
end

for ii=1:m       %generate MF operator
    C{2*ii-1}=1;
    C{2*ii}=1;
    for jj=1:m
        if jj==ii
             C{2*ii-1}=kron(C{2*ii-1},X);
             C{2*ii}=kron(C{2*ii},Y);
        else
            if jj<ii
                C{2*ii-1}=kron(C{2*ii-1},I);
                C{2*ii}=kron(C{2*ii},I);
            else
                C{2*ii-1}=kron(C{2*ii-1},Z);
                C{2*ii}=kron(C{2*ii},Z);
            end
        end
    end
end

for ii=0:4^m-1       %combine MF operator
    CC{ii+1}=1;
    for jj=1:2*m
        n(jj)=mod(floor(ii/2^(2*m-jj)),2);
        CC{ii+1}=CC{ii+1}*C{jj}^n(jj);
    end
end

for ii=1:4^m     %calculate coefficients
    for jj=1:2^m*3^(m-1)
        F(ii,jj)=trace(CC{ii}*UQ{jj}/2^m);
    end
end

for ii=1:2^m*3^(m-1)     %test how many coefficients is 0 or 1 in each column and where it is.
    Sum0(ii)=0;
    Sum1(ii)=0;
    L=1;
    for jj=1:4^m
        if abs(F(jj,ii))<1e-8
            Sum0(ii)=Sum0(ii)+1;
        end
        if abs(F(jj,ii))>1-1e-8&&abs(F(jj,ii))<1+1e-8
            Sum1(ii)=Sum1(ii)+1;
            L1(L,ii)=jj;
            L=L+1;
        end
    end
end      %Sum0 is all 4^m-1 and Sum1 is all 1

a=0;
b=0;
for ii=0:4^m-1     %test all the coefficients 1 is in even space and every operator in even space have at least one coefficient 1
    HW=0;
    for jj=1:2*m
        n(jj)=mod(floor(ii/2^(2*m-jj)),2);
        HW=HW+n(jj);
    end
    if mod(HW,2)==0
        a=a+1;
        TFeven(a,:)=F(ii+1,:);    %even space
    end
    if mod(HW,2)==1
        b=b+1;
        TFodd(b,:)=F(ii+1,:);     %odd space
    end
end
SUMeven=sum(abs(TFeven.'));
SUModd=sum(abs(TFodd.'));
t1=0;
t2=0;
for ii=1:2^(2*m-1)
    if abs(SUMeven(ii))<1e-8
        t1=t1+1;
    end
end
for ii=1:2^(2*m-1)
    if abs(SUModd(ii))>1e-8
        t2=t2+1;
    end
end
% t1 and t2 is zero is right