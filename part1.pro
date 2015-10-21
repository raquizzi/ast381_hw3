function part1,t_s,r_s,r_o,w      ;star temperature in K, star radius in solar radii, orbital radius in au, wavelength array in microns

  h=6.6260704d-34               ;planck constant in J s
  c=2.99792458d8                ;speed of light in m s^-1
  k=1.3806488d-23               ;boltzmann constant in J K^-1
  au2m=1.496d11                 ;au to meters conversion (m)
  rsol=6.96d8                   ;solar radius to meters conversion (m)
  
  l=w*1d-6                      ;using wavelength array from Draine tables, convert it to meters
  nu=c/reverse(l)                        ;frequency array in Hz  
  f_nu=2*h*nu^3*!pi*atan(r_s*rsol/(r_o*au2m))^2*1d26/(c^2*(exp(h*nu/(k*t_s))-1)) ;calculate flux density in Jy from planck function and scaling by stellar radius and distance from star

  output=dblarr(2,n_elements(l))
  output[0,*]=l*1d6             ;output wavelength array in microns
  output[1,*]=reverse(f_nu)              ;Jy
  
  return,output

end
