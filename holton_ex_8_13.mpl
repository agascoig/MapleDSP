                            "maple in ~/.mapleinit"

Warning, on line 1, incomplete string; use " to end the string
# Holton p. 560 Example 8.13 Type-II discrete-time inverse Chebyshev lowpass filter
# Design a discrete-time inverse (Type-II) Chebyshev lowpass filter using the bilinear-transformation method using the same parameters as the previous example: \317\211
# p
# = 0.5\317\200, \316\264
# p
# = -3 dB,
# \317\211
# s
# = 0.7\317\200 and \316\264
# s
# = -40 dB.
restart: 
assume(z,real): 
omega__p := `*`(.5,Pi): 
omega__s := `*`(.7,Pi): 
d__p := -3: 
d__s := -40: 
WpT := 2*tan(`/`(omega__p,2)): 
WsT := 2*tan(`/`(omega__s,2)): 
e__s := evalf(sqrt(10^(-`/`(d__s,10))-1)): 
e__p := evalf(sqrt(10^(-`/`(d__p,10))-1)): 
x := arccosh(e__s/e__p): 
N := ceil(x/arccosh(WsT/WpT)): 
WcT := WsT: 
k := Vector(N,i -> i-1): 
psi := `/`((`+`~(2.*k,1))*Pi,2*N): 
phi := arcsinh(e__s)/N: 
pkT := `/`~(WcT,-sin~(psi)*sinh(phi)+I*cos~(psi)*cosh(phi)): 
psi2 := Concatenate(1,psi[1 .. floor(`/`(N,2))],psi[floor(`/`(N,2)) +~ 2 .. ()]): 
zkT := I . WcT /~ cos~(psi2): 
pls := (`+`~(2,pkT)) /~ (`-`~(2,pkT)): # fine;
zrs := Concatenate(1,(`+`~(2,zkT)) /~ (`-`~(2,zkT)),-1): # fine;
b := expand(Re~(mul((zkT[i]-2)/zkT[i],i = 1 .. NumElems(zkT))*mul(pkT[i]/(pkT[i]-2),i = 1 .. NumElems(pkT)))*Re~(mul(z-zr,zr in zrs))): 
a := Re~(expand(mul(z-p,p in pls))): 
b2 := Vector([seq(coeff(b,z,i),i = degree(b,z) .. 0,-1)]): 
a2 := Vector([seq(coeff(a,z,i),i = degree(a,z) .. 0,-1)]): 
HF := FilterFrequencyResponse([b2],[a2],size = 1024,output = record,fftnormalization = none): 
HF[magnitudeplot];
HF[phaseplot];
freq := HF[frequencies]: 
response := HF[response]: 
plot(freq,abs(response),title = "Absolute response magnitude",color = blue,gridlines = true,labels = ["Normalized Frequency (&times; &pi; rad/sample)", typeset("|H(",omega,")|")],labeldirections = [horizontal, vertical]);

