clear all
m=3;
n=2^m;
X=[0,1;1,0];
Y=[0,-1i;1i,0];
Z=[1,0;0,-1];
I=eye(2);
II=eye(n);
S={X Y Z I};
for ii=1:(m-1)       %generate R gate
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

G{1,1}=II;       %generate G gate
for ii=1:m-1
    for jj=1:4^(ii-1)
        G{ii+1,4*jj-3}=G{ii,jj};
        G{ii+1,4*jj-2}=R{2*ii-1}*G{ii,jj};
        G{ii+1,4*jj-1}=R{2*ii}*G{ii,jj};
        G{ii+1,4*jj}=R{2*ii-1}^2*G{ii,jj};
    end
end

irho=1;
for ii=1:m      %initial state
    irho=1/2*kron(irho,I-Z);
end
for ii=1:4^(m-1)    %generate state
    rho{ii}=G{m,ii}*irho*G{m,ii}';
end

for ii=1:m      %generate MF operator
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

for ii=0:4^m-1    %combine MF operator
    CC{ii+1}=1;
    for jj=1:2*m
        n(jj)=mod(floor(ii/2^(2*m-jj)),2);
        CC{ii+1}=CC{ii+1}*C{jj}^n(jj);
    end
end

for ii=1:4^m    %calculate coefficients
    for jj=1:4^(m-1)
        F(ii,jj)=trace(CC{ii}*rho{jj})/2^m;
    end
end

a=0;
for ii=0:4^m-1    %Check for linear independence
    HW=0;
    for jj=1:2*m
        n(jj)=mod(floor(ii/2^(2*m-jj)),2);
        HW=HW+n(jj);
    end
    if mod(HW,2)==0&&n(2*m)==0
        a=a+1;
        TF(a,:)=F(ii+1,:);
    end
end
rank(TF)    %4^(m-1) is right