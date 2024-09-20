using AbstractAlgebra
print("\e[2J");#iiordanis
#"-----------------------------------------------------
function CreatePolRing(xyz,order)
    #Δημιουργία πολυωνυμικού δακτυλίου
    ord=Symbol(order)
    le=length(xyz)
    s=string("R,(",xyz[1])
    for i=2:le; s=string(s,",",xyz[i]);end
    global Pr=PolynomialRing(QQ, xyz, ordering=ord)
    s=string(s,")=Pr")
    eval(Meta.parse(s))
end
#"-----------------------------------------------------
function CreatePoly(pol)
    #κατασκευή συνόλου πολυωνύμων F
    s=string("F=[", pol[1])
    for i=2:length(pol);s=string(s,",",pol[i]);  end
    s=string(s,"]")
    eval(Meta.parse(s))
    F
end
#"------------------------------------------------------------"
function SPoly(f,g)
    #Κατασκευή S-πολυωνύμου
    n1=leading_term(f);v1 = exponent_vector(n1, 1)
    n2=leading_term(g);v2 = exponent_vector(n2, 1)
    v=max.(v1,v2)
    h=1;for i=1:length(gens(R));h=h*gens(R)[i]^v[i];end
    sp=div(h,n1)*f-div(h,n2)*g
end
#"------------------------------------------------------------"
function Groebasis(F)
    #Κατασκευή βάσης Groebner από αλγόριθμο Buchberger
    i=1
    while(i<length(F))
        i+=1
        for j=1:i-1
            sp=SPoly(F[i],F[j])
            _,rem=divrem(sp,F);#println();println(i,",",j);println(length(F))
            if rem!=0 ; F=vcat(F,rem); end;# println(length(F)); println("--------------")
        end
    end
    F
end
function ProMinGroebasis(G)
    #Έλεγχος ύπαρξης σταθερού πολυωνύμου στη βάση Groebner G
    s=0
    for i=1:length(G)
        GG=div(G[i],leading_coefficient(G[i])*gens(R)[1]^0)
        if GG==1;s+=1;end
    end
s
end
#"-------------------------------------------------------------"
function MinGroebasis(G)
    #Κατασκευή ελαχιστικής βάσης Groebner
    MO=deepcopy(G);MG=deepcopy(G);i=1
    for i=1:length(G);MO[i]=leading_term(G[i]);end
    while i<=length(MG)
        b=deleteat!(deepcopy(MO),i);_,rem=divrem(MO[i],b)
        if rem==0
            MG=deleteat!(deepcopy(MG),[i]);
            MO=deleteat!(deepcopy(MO),[i]);
            i-=1
        else
            MG[i]=div(MG[i],leading_coefficient(MG[i])*gens(R)[1]^0)
        end
        i+=1
    end
    MG
end
#"-------------------------------------------------------------"
function RedGroebasis(MG)
    #Κατασκευή ανηγμένης βάσης Groebner
    RG=deepcopy(MG)
    for i=1:length(RG)
        b=deleteat!(deepcopy(RG),[i])
        _,RG[i]=divrem(RG[i],b)
    end
    RG
end
#"--------------------- ----------------------------------------"
function CheckGroebasis(F)
    #Έλεγχος βάσης Groebner
    s=0; n=length(F)
    for i=1:n-1
        for j=i+1:n
            sp=SPoly(F[i],F[j]);_,rem=divrem(sp,F)
            if rem!=0;s+=1;end
        end
    end
    if s==0;println("Έλεγχος: Είναι Βάση Groebner")
    else; println("'Έλεγχος: Δεν είναι Βάση Groebner")
    end
    println("----------------------------------------")
end
#"-------------------------------------------------------------"
function Createxi(n)
    #Δημιουργία μεταβλητών xi
    xy=["x1"]
    for i=2:n ; xy=vcat(xy,string("x",i)) ;end
    xy
end
#"---------------------------------------------------"
function PrintGroebner()
    #Εκτύπωση αποτελεσμάτων
    println("")
    println("----------------------------------------")
    println();println("Πολυώνυμα F")
    print();println();for i=1:length(F);println(F[i]);end
    println("----------------------------------------")
    println();println("Βάση Groebner από αλγόριθμο Buchberg")
    print();println();for i=1:length(G);println(G[i]);end
    println("----------------------------------------")
    println();println("Ελαχιστική Βάση Groebner")
    println();println();for i=1:length(MinG);println(MinG[i]);end
    println("----------------------------------------")
    println();println("Ανηγμένη Βάση Groebner")
    println();for i=1:length(RedG);println(RedG[i]);end
    println("----------------------------------------")
end
#"---------------------------------------------------"
function go(xyz, order, pol )
    #Κύριο πρόγραμμα
    CreatePolRing(xyz,order)
    F=CreatePoly(pol)
    global G=Groebasis(F)
    CG=ProMinGroebasis(G)
    global MinG
    global RedG
    if CG==1;MinG=CG;else;MinG=MinGroebasis(G);end
    if CG==1;RedG=CG;else;RedG=RedGroebasis(MinG);end
    PrintGroebner()
    CheckGroebasis(RedG)
end
#"---------------------------------------------------"
#Δεδομένα
order="lex"  #είδος διάταξης μονωνύμων ("lex" ,"deglex",¨degrevlex")

xyz=[
"x",
"y"
]

pol=[
"x+y-8","x^2+x*y^2-84",
#"x-2"
]
#'----------------------------------'
go(xyz,order,pol)
