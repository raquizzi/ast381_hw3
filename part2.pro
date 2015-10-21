function part2,lam,f_nu,qa,a  ;wavelength array in microns, flux density from star in Jy, na, radius of dust grain in microns

  c=2.99792458d8                ;speed of light in m s^-1
  nu=c/(reverse(lam)*1d-6)      ;frequency array in Hz
  f_nu_r=reverse(f_nu)
  qa_r=reverse(qa)

  flux=dblarr(n_elements(lam)-1)
  for i=0,n_elements(lam)-2 do begin
     flux[i]=(nu[i+1]-nu[i])*(f_nu_r[i]*qa_r[i]+f_nu_r[i+1]*qa_r[i+1])*1d-26/2.0;use trapezoid rule to calculate flux density in J s^-1 m^-2 Hz^-1
  endfor
  
  output=!pi*(a*1d-6)^2*total(flux) ;sum flux density over frequency to get total flux and multiply by grain cross-section to get power in J s^-1

  return,output

end
