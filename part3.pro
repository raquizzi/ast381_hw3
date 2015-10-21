function part3,p,q,t_arr,w,a;power absorbed in J s^-1, na, temperature array in K, wavelength array in microns, dust grain radius in microns

  h=6.6260704d-34               ;planck constant in J s
  c=2.99792458d8                ;speed of light in m s^-1
  k=1.3806488d-23               ;boltzmann constant J K^-1
  
  l=w*1d-6                      ;convert wavelength array to meters
  nu=c/reverse(l)               ;frequency array
  q_r=reverse(q)

  resid=dblarr(n_elements(t_arr))
  spec=dblarr(n_elements(t_arr),n_elements(w))

  for i0=0,n_elements(t_arr)-1 do begin
     
     for i1=0,n_elements(w)-2 do begin
        
        flux_a=2*h*nu[i1]^3*!pi*q_r[i1]/(c^2*(exp(h*nu[i1]/(k*t_arr[i0]))-1))
        flux_b=2*h*nu[i1+1]^3*!pi*q_r[i1+1]/(c^2*(exp(h*nu[i1+1]/(k*t_arr[i0]))-1))
        spec[i0,i1]=(nu[i1+1]-nu[i1])*(flux_a+flux_b)/2.0
        
     endfor
     
     resid[i0]=abs(p-(4*!pi*(a*1d-6)^2)*total(spec[i0,*]))
     
  endfor

  min_resid=min(resid,tg)
  spec_fin=(4*!pi*(a*1d-6)^2)*reform(spec[tg,*]) ;output spectral luminosity J s^-1 Hz^-1
  lum=total(spec_fin) ;output power emitted in J s^-1

  output=dblarr(n_elements(w)+2)
  output[0]=t_arr[tg]
  output[1]=lum
  output[2:*]=reverse(spec_fin)
  
  return,output

end
