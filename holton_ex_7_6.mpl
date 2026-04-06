                            "maple in ~/.mapleinit"

Warning, on line 1, incomplete string; use " to end the string
# Holton Example 7.6: Causal highpass filter
# Design a
# causal
# highpass filter with \317\211
# c
# = 0.25 and N = 25 using a rectangular-window lowpass prototype, raised
cosine, hamming, hann, blackman, kaiser.
restart: 
# filter parameters
f__c := .25: 
omega__c := 2*Pi*f__c: 
eps := .1e-3: 
N := 25: 
N2 := iquo(N-1,2): 
# generate truncated sinc function from ideal filter
w := Vector(N,i -> 2*Pi*(i-1)/N): # frequency samples:
w__A := Vector(N,i -> `if`(is(w[i] < omega__c or 2*Pi-w[i] < omega__c),1,eps)): # lowpass amplitude samples;
h__A := FFTShift(Re(InverseFFT(w__A,normalization = full))): # this should be ideal filter, shifted to middle;
# causal windows
win__rect := Vector(N,1): 
win__hann := Vector(N,i -> .5+(-1)*.5*cos(2.0*Pi*(i-1)/(N-1))): 
win__hanning := Vector(N,n -> .50000+(-1)*.50000*cos(2.0*Pi*(n-1)/(N-1))): 
win__blackman := Vector(N,n -> .42+(-1)*.5*cos(2.0*Pi*(n-1)/(N-1))+.8e-1*cos(4.0*Pi*(n-1)/(N-1))): 
(h__rect, h__hann, h__hanning, h__blackman) := `*`~(h__A,win__rect), `*`~(h__A,win__hann), `*`~(h__A,win__hanning), `*`~(h__A,win__blackman): 
# causal frequency response
ffr_options := ptmp: 
H__rect := FilterFrequencyResponse(h__rect,ffr_options("Rectangular Window Magnitude","Rectangular Window Phase")); 
H__hann := FilterFrequencyResponse(h__hann,ffr_options("Hann Window Magnitude","Hann Window Phase")); 
H__hanning := FilterFrequencyResponse(h__hanning,ffr_options("Hanning Window Magnitude","Hanning Phase")); 
H__blackman := FilterFrequencyResponse(h__blackman,ffr_options("Blackman Window Magnitude","Blackman Window Phase")); 

