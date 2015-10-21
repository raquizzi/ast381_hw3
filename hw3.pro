readcol,'/Users/ram/Dropbox/ut/2015fa/ast381/hw3/as_1e-1.txt',w,qa_s,qs_s,g_s
readcol,'/Users/ram/Dropbox/ut/2015fa/ast381/hw3/as_1e+0.txt',w,qa_m,qs_m,g_m
readcol,'/Users/ram/Dropbox/ut/2015fa/ast381/hw3/as_1e+1.txt',w,qa_l,qs_l,g_l

w=reverse(w)                    ;wavelength array from Draine table is in descending order
qa_s=reverse(qa_s)
qa_m=reverse(qa_m)
qa_l=reverse(qa_l)

foma1=part1(8590.0,1.842,10.0,w)
foma2=part1(8590.0,1.842,130.0,w)

set_plot,'ps'
device,file='/Users/ram/Dropbox/ut/2015fa/ast381/hw3/figures/hw3_1.eps',$
/landscape,/inches,xsize=11,ysize=8.5,/color;,yoffset=11,xoffset=0.25

plot,foma1[0,*],foma1[1,*],/xlog,/ylog,xrange=[1d-1,1d3],xtitle='Wavelength (microns)',ytitle='Flux Density (Jy)',charthick=4,thick=6,$
ythick=4,xthick=4,charsize=2
oplot,foma2[0,*],foma2[1,*],thick=6,linestyle=1
al_legend,['10 AU','130 AU'],linestyle=[0,1],colors=[cgcolor('black'),cgcolor('black')],position=[30,1d13],thick=[6,6],charthick=4,$
charsize=1.7,linsize=[0.5,0.5]

!p.multi=0
device,/close
set_plot,'X'

qa_p=dblarr(n_elements(qa_s))+1.0

p_s1=part2(w,foma1[1,*],qa_s,1e-1)
p_m1=part2(w,foma1[1,*],qa_m,1e+0)
p_l1=part2(w,foma1[1,*],qa_l,1e+1)
p_p1=part2(w,foma1[1,*],qa_p,1e+3)

p_s2=part2(w,foma2[1,*],qa_s,1e-1)
p_m2=part2(w,foma2[1,*],qa_m,1e+0)
p_l2=part2(w,foma2[1,*],qa_l,1e+1)
p_p2=part2(w,foma2[1,*],qa_p,1e+3)

print,'Powers absorbed:'
print,p_s1,p_s2
print,p_m1,p_m2
print,p_l1,p_l2
print,p_p1,p_p2

t_arr=(findgen(1000)+1.0)

p3s1=part3(p_s1,qa_s,t_arr,w,1e-1)
p3m1=part3(p_m1,qa_m,t_arr,w,1e+0)
p3l1=part3(p_l1,qa_l,t_arr,w,1e+1)
p3p1=part3(p_p1,qa_p,t_arr,w,1e+3)

p3s2=part3(p_s2,qa_s,t_arr,w,1e-1)
p3m2=part3(p_m2,qa_m,t_arr,w,1e+0)
p3l2=part3(p_l2,qa_l,t_arr,w,1e+1)
p3p2=part3(p_p2,qa_p,t_arr,w,1e+3)

t_s1=p3s1[0]
t_m1=p3m1[0]
t_l1=p3l1[0]
t_p1=p3p1[0]

t_s2=p3s2[0]
t_m2=p3m2[0]
t_l2=p3l2[0]
t_p2=p3p2[0]

l_s1=p3s1[1]
l_m1=p3m1[1]
l_l1=p3l1[1]
l_p1=p3p1[1]

l_s2=p3s2[1]
l_m2=p3m2[1]
l_l2=p3l2[1]
l_p2=p3p2[1]

s_s1=p3s1[2:*]
s_m1=p3m1[2:*]
s_l1=p3l1[2:*]
s_p1=p3p1[2:*]

s_s2=p3s2[2:*]
s_m2=p3m2[2:*]
s_l2=p3l2[2:*]
s_p2=p3p2[2:*]

print,'Grain temperatures:'
print,t_s1,t_s2
print,t_m1,t_m2
print,t_l1,t_l2
print,t_p1,t_p2

print,'Powers emitted:'
print,l_s1,l_s2
print,l_m1,l_m2
print,l_l1,l_l2
print,l_p1,l_p2

dis=4*!pi*(7.7*3.08567758d16)^2.0
sr_s=(1d-1*1d-6)^2.0/((7.7*3.08567758d16)^2.0)
sr_m=(1d+0*1d-6)^2.0/((7.7*3.08567758d16)^2.0)
sr_l=(1d+1*1d-6)^2.0/((7.7*3.08567758d16)^2.0)
sr_p=(1d+3*1d-6)^2.0/((7.7*3.08567758d16)^2.0)

set_plot,'ps'
device,file='/Users/ram/Dropbox/ut/2015fa/ast381/hw3/figures/hw3_3.eps',$
/landscape,/inches,xsize=11,ysize=8.5,/color ;,yoffset=11,xoffset=0.25

;plot,w,s_s1*1d29/dis,/xlog,/ylog,xrange=[1d-3,1d3],yrange=[1d-30,1d-10],ysty=1
;oplot,w,s_m1*1d29/dis
;oplot,w,s_l1*1d29/dis
;oplot,w,s_p1*1d29/dis
;oplot,w,s_s2*1d29/dis,linestyle=2
;oplot,w,s_m2*1d29/dis,linestyle=2
;oplot,w,s_l2*1d29/dis,linestyle=2
;oplot,w,s_p2*1d29/dis,linestyle=2

plot,w,s_s1,/xlog,/ylog,xrange=[1d-1,1d3],yrange=[1d-20,1d-4],ysty=1,xtitle='Wavelength (microns)',$
ytitle='Spectral Luminosity (J s!E-1!N Hz!E-1!N)',charthick=4,ythick=4,xthick=4,thick=6,charsize=1.5
oplot,w,s_m1,thick=6,color=cgcolor('red')
oplot,w,s_l1,thick=6,color=cgcolor('blue')
oplot,w,s_p1,thick=6,color=cgcolor('dark green')
oplot,w,s_s2,linestyle=1,thick=6
oplot,w,s_m2,linestyle=1,thick=6,color=cgcolor('red')
oplot,w,s_l2,linestyle=1,thick=6,color=cgcolor('blue')
oplot,w,s_p2,linestyle=1,thick=6,color=cgcolor('dark green')
al_legend,['a!Ias!N=0.1 !4l!3m; 10 AU','a!Ias!N=0.1 !4l!3m; 130 AU','a!Ias!N=1.0 !4l!3m; 10 AU','a!Ias!N=1.0 !4l!3m; 130 AU','a!Ias!N=10 !4l!3m; 10 AU','a!Ias!N=10 !4l!3m; 130 AU','a!Ipa!N=1000 !4l!3m; 10 AU','a!Ipa!N=1000 !4l!3m; 130 AU'],linestyle=[0,1,0,1,0,1,0,1],colors=[cgcolor('black'),cgcolor('black'),cgcolor('red'),cgcolor('red'),cgcolor('blue'),cgcolor('blue'),cgcolor('dark green'),cgcolor('dark green')],position=[0.13,1d-5],thick=[6,6,6,6,6,6,6,6],charthick=3,charsize=1.5,linsize=[0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25],box=-1

!p.multi=0
device,/close
set_plot,'X'

print,'Number of dust grains:'
print,7d2/max(s_s1*1d29/dis),1d4/max(s_s2*1d29/dis)
print,7d2/max(s_m1*1d29/dis),1d4/max(s_m2*1d29/dis)
print,7d2/max(s_l1*1d29/dis),1d4/max(s_l2*1d29/dis)
print,7d2/max(s_p1*1d29/dis),1d4/max(s_p2*1d29/dis)

print,'Number of dust grains (new):'
print,7d2/max(s_s1*1d29*sr_s),1d4/max(s_s2*1d29*sr_s)
print,7d2/max(s_m1*1d29*sr_m),1d4/max(s_m2*1d29*sr_m)
print,7d2/max(s_l1*1d29*sr_l),1d4/max(s_l2*1d29*sr_l)
print,7d2/max(s_p1*1d29*sr_p),1d4/max(s_p2*1d29*sr_p)

vol_s=4*!pi*(1d-1*1d-6)^3.0/3.0
vol_m=4*!pi*(1d+0*1d-6)^3.0/3.0
vol_l=4*!pi*(1d+1*1d-6)^3.0/3.0
vol_p=4*!pi*(1d+3*1d-6)^3.0/3.0

;Density of grains in kg m^-3 is 2.0 * 1d-3 *1d6 = 2d3

print,'Dust masses:'
print,7d2*vol_s*2d3/max(s_s1*1d29/dis),1d4*vol_s*2d3/max(s_s2*1d29/dis)
print,7d2*vol_m*2d3/max(s_m1*1d29/dis),1d4*vol_m*2d3/max(s_m2*1d29/dis)
print,7d2*vol_l*2d3/max(s_l1*1d29/dis),1d4*vol_l*2d3/max(s_l2*1d29/dis)
print,7d2*vol_p*2d3/max(s_p1*1d29/dis),1d4*vol_p*2d3/max(s_p2*1d29/dis)

print,'Dust masses (new):'
print,7d2*vol_s*2d3/max(s_s1*1d29*sr_s),1d4*vol_s*2d3/max(s_s2*1d29*sr_s)
print,7d2*vol_m*2d3/max(s_m1*1d29*sr_m),1d4*vol_m*2d3/max(s_m2*1d29*sr_m)
print,7d2*vol_l*2d3/max(s_l1*1d29*sr_l),1d4*vol_l*2d3/max(s_l2*1d29*sr_l)
print,7d2*vol_p*2d3/max(s_p1*1d29*sr_p),1d4*vol_p*2d3/max(s_p2*1d29*sr_p)

frad_s1=p_s1/2.99792458d8
frad_m1=p_m1/2.99792458d8
frad_l1=p_l1/2.99792458d8
frad_p1=p_p1/2.99792458d8
frad_s2=p_s2/2.99792458d8
frad_m2=p_m2/2.99792458d8
frad_l2=p_l2/2.99792458d8
frad_p2=p_p2/2.99792458d8

v1=sqrt(6.67408d-11*1.92*1.989d30/(10.0*1.496d11))
v2=sqrt(6.67408d-11*1.92*1.989d30/(130.0*1.496d11))

fpr_s1=p_s1*v1/(2.99792458d8)^2
fpr_m1=p_m1*v1/(2.99792458d8)^2
fpr_l1=p_l1*v1/(2.99792458d8)^2
fpr_p1=p_p1*v1/(2.99792458d8)^2
fpr_s2=p_s2*v2/(2.99792458d8)^2
fpr_m2=p_m2*v2/(2.99792458d8)^2
fpr_l2=p_l2*v2/(2.99792458d8)^2
fpr_p2=p_p2*v2/(2.99792458d8)^2

print,'Radiation forces:'
print,frad_s1,frad_s2
print,frad_m1,frad_m2
print,frad_l1,frad_l2
print,frad_p1,frad_p2

print,'Poynting Drag:'
print,fpr_s1,fpr_s2
print,fpr_m1,fpr_m2
print,fpr_l1,fpr_l2
print,fpr_p1,fpr_p2

t_s1=sqrt(2d3*vol_s*(10.0*1.496d11)/fpr_s1)/31557600.0
t_m1=sqrt(2d3*vol_m*(10.0*1.496d11)/fpr_m1)/31557600.0
t_l1=sqrt(2d3*vol_l*(10.0*1.496d11)/fpr_l1)/31557600.0
t_p1=sqrt(2d3*vol_p*(10.0*1.496d11)/fpr_p1)/31557600.0
t_s2=sqrt(2d3*vol_s*(130.0*1.496d11)/fpr_s2)/31557600.0
t_m2=sqrt(2d3*vol_m*(130.0*1.496d11)/fpr_m2)/31557600.0
t_l2=sqrt(2d3*vol_l*(130.0*1.496d11)/fpr_l2)/31557600.0
t_p2=sqrt(2d3*vol_p*(130.0*1.496d11)/fpr_p2)/31557600.0

print,'Timescales:'
print,t_s1,t_s2
print,t_m1,t_m2
print,t_l1,t_l2
print,t_p1,t_p2

beta_s1=(p_s1/2.99792458d8)/(6.67408d-11*1.92*1.989d30*2d3*vol_s/(10.0*1.496d11)^2)
beta_m1=(p_m1/2.99792458d8)/(6.67408d-11*1.92*1.989d30*2d3*vol_m/(10.0*1.496d11)^2)
beta_l1=(p_l1/2.99792458d8)/(6.67408d-11*1.92*1.989d30*2d3*vol_l/(10.0*1.496d11)^2)
beta_p1=(p_p1/2.99792458d8)/(6.67408d-11*1.92*1.989d30*2d3*vol_p/(10.0*1.496d11)^2)
beta_s2=(p_s2/2.99792458d8)/(6.67408d-11*1.92*1.989d30*2d3*vol_s/(130.0*1.496d11)^2)
beta_m2=(p_m2/2.99792458d8)/(6.67408d-11*1.92*1.989d30*2d3*vol_m/(130.0*1.496d11)^2)
beta_l2=(p_l2/2.99792458d8)/(6.67408d-11*1.92*1.989d30*2d3*vol_l/(130.0*1.496d11)^2)
beta_p2=(p_p2/2.99792458d8)/(6.67408d-11*1.92*1.989d30*2d3*vol_p/(130.0*1.496d11)^2)

t_s1_pap=400.0*10.0^2/beta_s1
t_m1_pap=400.0*10.0^2/beta_m1
t_l1_pap=400.0*10.0^2/beta_l1
t_p1_pap=400.0*10.0^2/beta_p1
t_s2_pap=400.0*130.0^2/beta_s2
t_m2_pap=400.0*130.0^2/beta_m2
t_l2_pap=400.0*130.0^2/beta_l2
t_p2_pap=400.0*130.0^2/beta_p2

print,'Paper Timescales:'
print,t_s1_pap,t_s2_pap
print,t_m1_pap,t_m2_pap
print,t_l1_pap,t_l2_pap
print,t_p1_pap,t_p2_pap

end
