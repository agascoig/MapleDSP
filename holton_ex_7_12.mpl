                            "maple in ~/.mapleinit"

Warning, on line 1, incomplete string; use " to end the string
# Holton DSP, p.454 Example 7.12 : Linear system frequency-sampled filter
# Type-1 FIR lowpass filter of Ex 7.10 using simultanous equations
restart: 
N := 17: 
eps := .1e-3: 
wc := `/`(Pi,2): # center frequency;
Nh := iquo(N-1,2): 
k := Vector(Nh+1,i -> i-1)^%T: # row vector;
w := 2*Pi*k^%H/N: # frequency samples, column vector;
Ak := Vector(Nh+1,i -> `if`(is(w[i] < wc or 2*Pi-w[i] < wc),1,eps)): # amplitude samples, column vector;
C := cos~(evalf~(w . k)): 
a := LinearSolve(C,Ak): # outer product;
h := Concatenate(1,.5*Reverse(a)[1 .. NumElems(a)-1],a[1],.5*a[2 .. ()]): # impulse response;
H := FilterFrequencyResponse(h,size = 1024,output = record,fftnormalization = none): #noncausal response;
H[magnitudeplot];
H[phaseplot];
freq := H[frequencies]: 
response := H[response]: 
plot(freq,abs(response),title = "Absolute response magnitude",color = blue,gridlines = true,labels = ["Normalized Frequency (&times; &pi; rad/sample)", typeset("|H(",omega,")|")],labeldirections = [horizontal, vertical]);

