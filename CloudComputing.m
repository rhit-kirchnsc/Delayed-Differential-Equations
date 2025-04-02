%Exercise 7 - Cloud Computing



function dQdt = SIR_DDE_eqns(t,y,Z)

    

    dQdt(1,1)=lambda-
    
    dQdt=zeros(n,1);
    mu=ones(n,1);
    r=ones(n,1);
    k=ones(n,1);
    tauI=ones(n,1);
    tauJ=ones(n,1);

    for i=2:n+1
        dQdt(i,1) = -mu(i,1)+k(i,1)*Qm(t-tauI(i,1))-r(i,1)()
    end

    dSdt=-ItHat;          %susceptible
    dIdt=ItHat-ItauHat;   %infected
    dRdt=ItauHat;         %recovered


    dQdt=[dSdt;dIdt;dRdt];
end
