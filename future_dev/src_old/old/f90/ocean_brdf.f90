! This code is a Fortran 90 version of SHDOM/ocean_brdf.f (F77 by Evans 1996)
!
! 2022-08-12
! by Hong Chen (hong.chen.cu@gmail.com) at Laboratory for Atmospheric and Space Physics

subroutine ocean_brdf_sw(pws, xsal, pcl, pwl, xmuo, xmu, xphi, xpaw, brdfx)

  ! INPUT:  pws=wind speed (in m/s) [value not allowed below 0.25 m/s]
  !         xsal=salinity (in ppt) [if xsal<0 then 34.3ppt by default]
  !         pcl=pigment concentration (in mg/m^3)
  !         pwl=wavelength of the computation (in micrometer)
  !         xmuo=cosine of the zenith angle of the incident radiation
  !         xmu=cosine of the zenith angle of the reflected radiation
  !         xphi=relative azimuth angle (between 0 and 2*pi)
  !         xpaw=azimuth angle of sun - azimuth angle of wind (in radian)
  !              For the direct solar radiation, the wind direction
  !              is fixed at the azimuth angle of 0 degree so that
  !              the azimuth angle of the sun = xpaw.
  !
  ! OUTPUT: brdfx=the total reflectance of the sea water

  implicit none

  real pws, xsal, pcl, pwl, xmuo, xmu, xphi, xpaw
  real brdfx

  real phi
  integer iwl, iang, isz1, isz2, ivz1, ivz2, iws, iws1, iws2

  real paw, wl, wspd, c, azw

  real ref(39)
  real pi , fac , nr , ni , n12
  real tetas , w , wlp , ref_i , rwc , rw , tds
  real tetav , tdv , fi , rog , a , rwb
  real angbnd(5) , wsbnd(6) , tdsbnd(5,6) , tdvbnd(5,6)

  ! effective reflectance of the whitecaps (Koepke, 1984)
  data ref/0.220 , 0.220 , 0.220 , 0.220 , 0.220 , 0.220 , 0.215 , &
         & 0.210 , 0.200 , 0.190 , 0.175 , 0.155 , 0.130 , 0.080 , &
         & 0.100 , 0.105 , 0.100 , 0.080 , 0.045 , 0.055 , 0.065 , &
         & 0.060 , 0.055 , 0.040 , 0.000 , 0.000 , 0.000 , 0.000 , &
         & 0.000 , 0.000 , 0.000 , 0.000 , 0.000 , 0.000 , 0.000 , &
         & 0.000 , 0.000 , 0.000 , 0.000/

  data angbnd/0.0 , 45.0 , 60.0 , 75.0 , 85.0/

  data wsbnd/1.0 , 3.0 , 5.0 , 7.0 , 9.0 , 20.0/

  data ((tdsbnd(iang,iws),iws=1,6),iang=1,5)/0.9787803 , 0.9787738 ,&
      & 0.9787626 , 0.9787467 , 0.9787264 , 0.9785573 , 0.9706900 , &
      & 0.9698871 , 0.9691746 , 0.9685547 , 0.9680276 , 0.9666586 , &
      & 0.9479931 , 0.9404608 , 0.9385692 , 0.9381815 , 0.9384519 , &
      & 0.9430056 , 0.9690591 , 0.9275920 , 0.9058769 , 0.8951812 , &
      & 0.8899654 , 0.8892645 , 0.9980542 , 0.9602273 , 0.9114283 , &
      & 0.8713799 , 0.8417820 , 0.7800314/

  data ((tdvbnd(iang,iws),iws=1,6),iang=1,5)/0.9787764 , 0.9787535 ,&
      & 0.9787106 , 0.9786453 , 0.9785548 , 0.9775019 , 0.9692680 , &
      & 0.9637051 , 0.9564344 , 0.9495727 , 0.9438773 , 0.9288712 , &
      & 0.9225163 , 0.9069787 , 0.9044844 , 0.9052351 , 0.9068328 , &
      & 0.9153687 , 0.8048478 , 0.8479503 , 0.8678726 , 0.8797889 , &
      & 0.8878716 , 0.9091171 , 0.7294627 , 0.8137348 , 0.8453338 , &
      & 0.8629867 , 0.8745421 , 0.9036854/

  pi  = atan(1.0)*4.0
  fac = pi/180.0

  ! paw= azim. of sun - azim. of wind (in deg.)
  paw = xpaw/fac

  ! conversion angle parameter
  if (xphi .lt. 0.0) then
    phi = -xphi
  elseif (xphi .ge. 2.0*pi) then
    phi = xphi - 2.0*pi
  else
    phi = xphi
  endif

  if (xmuo .le. 0.028) then
    tetas = acos(0.028)/fac
  else
    tetas = acos(xmuo)/fac
  endif

  if (xmu .le. 0.028) then
    tetav = acos(0.028)/fac
  else
    tetav = acos(xmu)/fac
  endif

  if (pwl .lt. 0.4) then
    wl = 0.4
  elseif (pwl .gt. 4.0) then
    wl = 4.0
  else
    wl = pwl
  endif

  fi = 180.0 - phi/fac

  c = pcl
  wspd = max(0.25, pws)
  azw = paw

  ! compute index of water
  call indwat(wl, Xsal, nr, ni)
  n12 = sqrt(nr*nr + ni*ni)

  ! compute whitecaps reflectance (lambertian)
  w = 2.95e-06 * (wspd**3.52)
  iwl = 1 + int((wl-0.2)/0.1)
  wlp = 0.5 + (iwl-1)*0.1
  ref_i = ref(iwl+1) + (wl-wlp)/0.1*(ref(iwl)-ref(iwl+1))
  rwc = w*ref_i

  ! compute backscattered reflectance from the sea water (lambertian)
  ! water reflectance below the sea surface
  call morcasiwat(wl, c, rw)

  call getbound(wsbnd , 1, 6, wspd , iws1, iws2)
  call getbound(angbnd, 1, 5, tetas, isz1, isz2)
  call getbound(angbnd, 1, 5, tetav, ivz1, ivz2)

  tds = tdsbnd(isz1, iws1)
  tdv = tdvbnd(ivz1, iws1)

  ! sunglint reflectance
  call sunglint(wspd, nr, ni, azw, tetas, tetav, fi, rog)

  ! water reflectance above the sea surface
  ! for explanation on value of a see OCEAALBE.f
  a = 0.485
  ! add change in solid angle from under to above to surface
  ! that account for 1/(n12*n12) decrease in sea water directional
  ! reflectance
  rwb = (1/(n12*n12))*tds*tdv*rw/(1-a*rw)

  ! total reflectance of sea water
  brdfx = rwc + (1-w)*rog + (1-rwc)*rwb

end subroutine




subroutine morcasiwat(wl, c, r2)
  ! Spectral diffuse attenuation coefficient of Case I Waters as Predicted
  ! by MOREL within the spectral range 400-700nm (1988, Journal of Geophysical
  ! Research, Vol.93, No C9, pp 10749-10768)
  !
  ! input parameters:
  !     wl: wavelength (IN MICROMETERS)
  !      c: pigment concentration
  !
  ! output parameter:
  !     r2:  reflectance of water
  !
  ! According Morel,1988, we use:
  !
  ! Kd   : spectral value of the attenuation coefficient for
  !        downwelling irradiance
  !        with Kd=Kw+Xc*C**e
  ! Kw   : spectral value of the diffuse attenuation coefficient
  !        for pure oceanic water
  ! Xc, e: spectral coefficients to compute the diffuse attenuation
  !        coefficient for pigment
  ! bb   : total backscattering coefficient
  !        with bb=0.5*bw+bbt*b
  ! bw   : spectral value of the molecular scattering coefficient of water
  ! bbt,b: parameters to compute the scattering coefficients of pigments
  ! R2   : reflectance of water below the surface
  !        with R2=(0.33/u)*(bb/Kd) where u is depending of R2

  implicit none
  real kw , kd
  real tkw(61) , txc(61) , te(61) , tbw(61)
  real wl , c , r2 , xc , e , bw , bb , b , bbt , u1 , r1 , u2 , err
  integer iwl

  data tkw/0.0209 , 0.0200 , 0.0196 , 0.0189 , 0.0183 , 0.0182 , &
         & 0.0171 , 0.0170 , 0.0168 , 0.0166 , 0.0168 , 0.0170 , &
         & 0.0173 , 0.0174 , 0.0175 , 0.0184 , 0.0194 , 0.0203 , &
         & 0.0217 , 0.0240 , 0.0271 , 0.0320 , 0.0384 , 0.0445 , &
         & 0.0490 , 0.0505 , 0.0518 , 0.0543 , 0.0568 , 0.0615 , &
         & 0.0640 , 0.0640 , 0.0717 , 0.0762 , 0.0807 , 0.0940 , &
         & 0.1070 , 0.1280 , 0.1570 , 0.2000 , 0.2530 , 0.2790 , &
         & 0.2960 , 0.3030 , 0.3100 , 0.3150 , 0.3200 , 0.3250 , &
         & 0.3300 , 0.3400 , 0.3500 , 0.3700 , 0.4050 , 0.4180 , &
         & 0.4300 , 0.4400 , 0.4500 , 0.4700 , 0.5000 , 0.5500 , 0.6500/

  data txc/0.1100 , 0.1110 , 0.1125 , 0.1135 , 0.1126 , 0.1104 , &
         & 0.1078 , 0.1065 , 0.1041 , 0.0996 , 0.0971 , 0.0939 , &
         & 0.0896 , 0.0859 , 0.0823 , 0.0788 , 0.0746 , 0.0726 , &
         & 0.0690 , 0.0660 , 0.0636 , 0.0600 , 0.0578 , 0.0540 , &
         & 0.0498 , 0.0475 , 0.0467 , 0.0450 , 0.0440 , 0.0426 , &
         & 0.0410 , 0.0400 , 0.0390 , 0.0375 , 0.0360 , 0.0340 , &
         & 0.0330 , 0.0328 , 0.0325 , 0.0330 , 0.0340 , 0.0350 , &
         & 0.0360 , 0.0375 , 0.0385 , 0.0400 , 0.0420 , 0.0430 , &
         & 0.0440 , 0.0445 , 0.0450 , 0.0460 , 0.0475 , 0.0490 , &
         & 0.0515 , 0.0520 , 0.0505 , 0.0440 , 0.0390 , 0.0340 , 0.0300/

  data te/0.668 , 0.672 , 0.680 , 0.687 , 0.693 , 0.701 , 0.707 , &
        & 0.708 , 0.707 , 0.704 , 0.701 , 0.699 , 0.700 , 0.703 , &
        & 0.703 , 0.703 , 0.703 , 0.704 , 0.702 , 0.700 , 0.700 , &
        & 0.695 , 0.690 , 0.685 , 0.680 , 0.675 , 0.670 , 0.665 , &
        & 0.660 , 0.655 , 0.650 , 0.645 , 0.640 , 0.630 , 0.623 , &
        & 0.615 , 0.610 , 0.614 , 0.618 , 0.622 , 0.626 , 0.630 , &
        & 0.634 , 0.638 , 0.642 , 0.647 , 0.653 , 0.658 , 0.663 , &
        & 0.667 , 0.672 , 0.677 , 0.682 , 0.687 , 0.695 , 0.697 , &
        & 0.693 , 0.665 , 0.640 , 0.620 , 0.600/

  data tbw/0.0076 , 0.0072 , 0.0068 , 0.0064 , 0.0061 , 0.0058 , &
         & 0.0055 , 0.0052 , 0.0049 , 0.0047 , 0.0045 , 0.0043 , &
         & 0.0041 , 0.0039 , 0.0037 , 0.0036 , 0.0034 , 0.0033 , &
         & 0.0031 , 0.0030 , 0.0029 , 0.0027 , 0.0026 , 0.0025 , &
         & 0.0024 , 0.0023 , 0.0022 , 0.0022 , 0.0021 , 0.0020 , &
         & 0.0019 , 0.0018 , 0.0018 , 0.0017 , 0.0017 , 0.0016 , &
         & 0.0016 , 0.0015 , 0.0015 , 0.0014 , 0.0014 , 0.0013 , &
         & 0.0013 , 0.0012 , 0.0012 , 0.0011 , 0.0011 , 0.0010 , &
         & 0.0010 , 0.0010 , 0.0010 , 0.0009 , 0.0008 , 0.0008 , &
         & 0.0008 , 0.0007 , 0.0007 , 0.0007 , 0.0007 , 0.0007 , 0.0007/

  if (wl .lt. 0.400 .or. wl .gt. 0.700) then
    r2 = 0.000
    return
  endif

  iwl = 1 + nint((wl-0.400)/0.005)
  kw  = tkw(iwl)
  xc  = txc(iwl)
  e   = te(iwl)
  bw  = tbw(iwl)

  if (abs(c) .lt. 0.0001) then
     bb = 0.5*bw
     kd = kw
  else
     b = 0.30*c**0.62
     bbt = 0.002 + 0.02*(0.5-0.25*alog10(c))*0.550/wl
     bb = 0.5*bw + bbt*b
     kd = kw + xc*c**e
  endif

  u1 = 0.75
  r1 = 0.33*bb/u1/kd

  do
    u2 = 0.90*(1.-r1)/(1.+2.25*r1)
    r2 = 0.33*bb/u2/kd
    err = abs((r2-r1)/r2)

    if (err .ge. 0.0001) then
      r1 = r2
    else
      exit
    endif
  enddo

end subroutine






subroutine indwat(wl, xsal, nr, ni)

  ! input parameters:  wl=wavelength (in micrometers)
  !                    xsal=salinity (in ppt), if xsal<0 then 34.3ppt by default
  ! output parameters: nr=index of refraction of sea water
  !                    ni=extinction coefficient of sea water

  implicit none

  real twl(62) , tnr(62) , tni(62)
  real nr, ni, wl, xwl, yr, yi, nrc, nic, xsal
  integer i

  ! Indices of refraction for pure water from Hale and Querry,
  ! Applied Optique, March 1973, Vol. 12,  No. 3, pp. 555-563
  data twl/0.250 , 0.275 , 0.300 , 0.325 , 0.345 , 0.375 , 0.400 , &
         & 0.425 , 0.445 , 0.475 , 0.500 , 0.525 , 0.550 , 0.575 , &
         & 0.600 , 0.625 , 0.650 , 0.675 , 0.700 , 0.725 , 0.750 , &
         & 0.775 , 0.800 , 0.825 , 0.850 , 0.875 , 0.900 , 0.925 , &
         & 0.950 , 0.975 , 1.000 , 1.200 , 1.400 , 1.600 , 1.800 , &
         & 2.000 , 2.200 , 2.400 , 2.600 , 2.650 , 2.700 , 2.750 , &
         & 2.800 , 2.850 , 2.900 , 2.950 , 3.000 , 3.050 , 3.100 , &
         & 3.150 , 3.200 , 3.250 , 3.300 , 3.350 , 3.400 , 3.450 , &
         & 3.500 , 3.600 , 3.700 , 3.800 , 3.900 , 4.000/

  data tnr/1.362 , 1.354 , 1.349 , 1.346 , 1.343 , 1.341 , 1.339 , &
         & 1.338 , 1.337 , 1.336 , 1.335 , 1.334 , 1.333 , 1.333 , &
         & 1.332 , 1.332 , 1.331 , 1.331 , 1.331 , 1.330 , 1.330 , &
         & 1.330 , 1.329 , 1.329 , 1.329 , 1.328 , 1.328 , 1.328 , &
         & 1.327 , 1.327 , 1.327 , 1.324 , 1.321 , 1.317 , 1.312 , &
         & 1.306 , 1.296 , 1.279 , 1.242 , 1.219 , 1.188 , 1.157 , &
         & 1.142 , 1.149 , 1.201 , 1.292 , 1.371 , 1.426 , 1.467 , &
         & 1.483 , 1.478 , 1.467 , 1.450 , 1.432 , 1.420 , 1.410 , &
         & 1.400 , 1.385 , 1.374 , 1.364 , 1.357 , 1.351/

  data tni/3.35e-08 , 2.35e-08 , 1.60e-08 , 1.08e-08 , 6.50e-09 , &
         & 3.50e-09 , 1.86e-09 , 1.30e-09 , 1.02e-09 , 9.35e-10 , &
         & 1.00e-09 , 1.32e-09 , 1.96e-09 , 3.60e-09 , 1.09e-08 , &
         & 1.39e-08 , 1.64e-08 , 2.23e-08 , 3.35e-08 , 9.15e-08 , &
         & 1.56e-07 , 1.48e-07 , 1.25e-07 , 1.82e-07 , 2.93e-07 , &
         & 3.91e-07 , 4.86e-07 , 1.06e-06 , 2.93e-06 , 3.48e-06 , &
         & 2.89e-06 , 9.89e-06 , 1.38e-04 , 8.55e-05 , 1.15e-04 , &
         & 1.10e-03 , 2.89e-04 , 9.56e-04 , 3.17e-03 , 6.70e-03 , &
         & 1.90e-02 , 5.90e-02 , 1.15e-01 , 1.85e-01 , 2.68e-01 , &
         & 2.98e-01 , 2.72e-01 , 2.40e-01 , 1.92e-01 , 1.35e-01 , &
         & 9.24e-02 , 6.10e-02 , 3.68e-02 , 2.61e-02 , 1.95e-02 , &
         & 1.32e-02 , 9.40e-03 , 5.15e-03 , 3.60e-03 , 3.40e-03 , &
         & 3.80e-03 , 4.60e-03/

  i = 2
  do
    if (wl .ge. twl(i)) then
      exit
    else
      if (i .lt. 62 ) then
        i = i + 1
      else
        exit
      endif
    endif
  enddo

  xwl = twl(i) - twl(i-1)
  yr  = tnr(i) - tnr(i-1)
  yi  = tni(i) - tni(i-1)
  nr  = tnr(i-1) + (wl-twl(i-1))*yr/xwl
  ni  = tni(i-1) + (wl-twl(i-1))*yi/xwl

  ! Correction to be applied to the index of refraction and to the extinction
  ! coefficients of the pure water to obtain the ocean water one (see for
  ! example Friedman). By default, a typical sea water is assumed
  ! (Salinity=34.3ppt, Chlorinity=19ppt) as reported by Sverdrup.
  ! In that case there is no correction for the extinction coefficient between
  ! 0.25 and 4 microns. For the index of refraction, a correction of +0.006
  ! has to be applied (McLellan). For a chlorinity of 19.0ppt the correction
  ! is a linear function of the salt concentration. Then, in 6S users are able
  ! to enter the salt concentration (in ppt).
  ! REFERENCES:
  ! Friedman D., Applied Optics, 1969, Vol.8, No.10, pp.2073-2078.
  ! McLellan H.J., Elements of physical Oceanography, Pergamon Press, Inc.,
  !        New-York, 1965, p 129.
  ! Sverdrup H.V. et al., The Oceans (Prentice-Hall, Inc., Englewood Cliffs,
  !        N.J., 1942, p 173.

  nrc = 0.006
  nic = 0.000
  if (xsal .ge. 0.0) then
     nr = nr + nrc*(xsal/34.3)
     ni = ni + nic*(xsal/34.3)
  else
     nr = nr + nrc
     ni = ni + nic
   endif

end subroutine



subroutine sunglint(wspd, nr, ni, azw, ts, tv, fi, rog)

  ! input parameters:   wspd=speed of the wind (in m/s)
  !                     nr=index of refraction of the sea water
  !                     ni=extinction coefficient of the sea water
  !                     azw=azim. of the sun - azim. of the wind (in deg.)
  !                     ts=solar zenith angle (in deg.)
  !                     tv=view zenith angle (in deg.)
  !                     fi=relative azimuth (sun-satellite)
  ! output parameters:  rog=reflectance of the sun glint

  implicit none
  real pi, fac
  real wspd, nr, ni, ts, tv, fi, rog, azw, phw
  real cs, cv, ss, sv, phi, zx, zy, tantilt, tilt, proba, xe, xn, xe2, xn2
  real coef, cos2chi, coschi, sinchi
  real r1, sigmac, sigmau, c21, c03, c40, c04, c22

  pi  = atan(1.0)*4.0
  fac = pi/180.0
  phw = azw*fac
  cs  = cos(ts*fac)
  cv  = cos(tv*fac)
  ss  = cos(ts*fac)
  sv  = cos(tv*fac)
  phi = fi*fac
  zx  = -sv*sin(phi)/(cs+cv)
  zy  = (ss+sv*cos(phi))/(cs+cv)
  tantilt = sqrt(zx*zx+zy*zy)
  tilt    = atan(tantilt)

  ! Anisotropic Gaussian distribution
  ! phw=phi_sun-phi_wind
  sigmac = 0.003 + 0.00192*wspd
  sigmau = 0.00316*wspd
  c21 = 0.01 - 0.0086*wspd
  c03 = 0.04 - 0.033*wspd
  c40 = 0.40
  c22 = 0.12
  c04 = 0.23
  xe  = (cos(phw)*zx+sin(phw)*zy)/sqrt(sigmac)
  xn  = (-sin(phw)*zx+cos(phw)*zy)/sqrt(sigmau)
  xe2 = xe*xe
  xn2 = xn*xn
  coef  = 1.0 - c21/2.0*(xe2-1.0)*xn - c03/6.0*(xn2-3.0)*xn
  coef  = coef + c40/24.0*(xe2*xe2-6.0*xe2+3.0)
  coef  = coef + c04/24.0*(xn2*xn2-6.0*xn2+3.0)
  coef  = coef + c22/4.0*(xe2-1.0)*(xn2-1.0)
  proba = coef/2.0/pi/sqrt(sigmau)/sqrt(sigmac)*exp(-(xe2+xn2)/2.0)

  ! Compute Fresnel's coefficient R1
  cos2chi = cv*cs + sv*ss*COS(phi)
  if (cos2chi .gt.  1.0) cos2chi =  0.99999999999
  if (cos2chi .lt. -1.0) cos2chi = -0.99999999999

  coschi = sqrt(0.5*(1+cos2chi))
  sinchi = sqrt(0.5*(1-cos2chi))
  call fresnel(nr, ni, coschi, sinchi, r1)

  ! Compute Reflectance of the sun glint
  rog = pi*r1*proba/4.0/cs/cv/(cos(tilt)**4.0)

end subroutine




subroutine fresnel(nr, ni, coschi, sinchi, r1)

  ! to compute the Fresnel's coefficient of reflection (see for
  ! example M. Born and E. Wolf, Principles of Optics, Pergamon Press, fifth
  ! edition, 1975, pp 628
  ! input parameters: nr=index of refraction of the sea water
  !                   ni=extinction coefficient of the sea water
  !                   coschi & sinchi=cosine and sine of the incident radiation
  !                                   with respect of the wave facet normal.
  ! output parameter: R1=Fresnel's coefficient for reflection

  implicit none
  real nr, ni, a1, a2, u, v, rr2, rl2, b1, b2, r1, coschi, sinchi

  ! absolute value for a1 to get v=0 when ni=0
  a1  = abs(nr*nr-ni*ni-sinchi*sinchi)
  a2  = sqrt((nr*nr-ni*ni-sinchi*sinchi)**2.0+4.0*nr*nr*ni*ni)
  u   = sqrt(0.5*(a1+a2))
  v   = sqrt(max(0.0, 0.5*(-a1+a2)))
  rr2 = ((coschi-u)**2+v*v)/((coschi+u)**2+v*v)
  b1  = (nr*nr-ni*ni)*coschi
  b2  = 2.0*nr*ni*coschi
  rl2 = ((b1-u)**2.0+(b2-v)**2.0)/((b1+u)**2.0+(b2+v)**2.0)
  r1  = (rr2+rl2)/2.0

end subroutine



subroutine glitable(wspd, nr, ni, azw, rge)

  ! To compute the spherical albedo of the sea water. See for example
  ! Masuda et al., Remote Sens. Environ., 24, 313-329, 1988.
  !
  ! input parameters: wsp=wind of speed
  !                   nr=index of refraction of the sea water
  !                   ni=extinction coefficient of the sea water
  !                   azw=azim. of sun - azim. of wind (in deg.)
  ! output parameter: rge=spherical albedo of the sun glint

  implicit none
  integer nfa, nta, ntb

  real Nr, Ni, Azw, phw, Rge, q, Wspd, prefl, proba, pr, pp, pi, fac
  real sigma, sigmac, sigmau, c21, c03, c40, c04, c22
  real costt, hta, htb, hfa, cotb, cota, cofa, diff, coef
  real phin, cosphin, sinphin, costet, tet, sintet
  real costetn, sintetn, tantetn, coschi, sinchi
  real zx, zy, xe, xn, xe2, xn2, fonc0, pond, r1
  integer km, i, j

  pi     = atan(1.0)*4.0
  fac    = pi/180.0
  sigma  = 0.003 + 0.00512*wspd
  sigmac = 0.003 + 0.00192*wspd
  sigmau = 0.00316*wspd
  c21 = 0.01 - 0.0086*wspd
  c03 = 0.04 - 0.033*wspd
  c40 = 0.40
  c22 = 0.12
  c04 = 0.23

  ! costt to minimize the time of the computation
  ! integration between 1 and costt instead of 1 and 0
  q = 50
  costt = 1.0/sqrt(1.0+q*sigma/4.0)
  phw = azw*fac

  prefl = 0.0
  proba = 0.0

  ntb = 31
  htb = 1.0/float(ntb-1.0)

  ! loops on the zenith angle of the emitted radiation
  do km = 1 , ntb
    costet = (km-1)*htb
    tet = acos(costet)
    sintet = sin(tet)
    tet = tet/fac

    ! Simpson's rules for the angle of the emitted radiation teta
    cotb = 2.0
    diff = abs(km/2 - km/2.)
    if (diff .lt. 0.00001) cotb = 4.0
    if (km .eq. 1 .OR. km.EQ.ntb ) cotb = 1.0
!  loops step for phiN and tetaN (N is the facet unit normal vector)
         IF ( tet.LT.91 ) nta = 801
         IF ( tet.LT.81 ) nta = 301
         IF ( tet.LT.75 ) nta = 101
         IF ( tet.LT.65 ) nta = 31
         nfa = nta
         hta = (1.-costt)/FLOAT(nta-1)
         hfa = pi/FLOAT(nfa-1)
! loops on phiN (azimuth angle of the facet normal vector)
         pr = 0.
         pp = 0.
         DO i = 1 , nfa
            phin = (i-1)*hfa
            cosphin = COS(phin)
            sinphin = SIN(phin)
!  Simpson's rules for phin
            cofa = 2.
            diff = ABS(i/2-i/2.)
            IF ( diff.LT.0.00001 ) cofa = 4.
            IF ( i.EQ.1 .OR. i.EQ.nfa ) cofa = 1.0
! loops on tetaN (zenith angle of the facet normal vector)
            DO j = 1 , nta
               costetn = costt + (j-1)*hta
               sintetn = SQRT(1-costetn*costetn)
               tantetn = sintetn/costetn
!  Simpson's rules for tetaN
               cota = 2.
               diff = ABS(j/2-j/2.)
               IF ( diff.LT.0.00001 ) cota = 4.
               IF ( j.EQ.1 .OR. j.EQ.nta ) cota = 1.0
! Fresnel's reflection coefficient R1
               coschi = costet*costetn + sintet*sintetn*cosphin
               IF ( coschi*coschi.GT.1.0 ) coschi = 0.99999999999
               sinchi = SQRT(1-coschi*coschi)
               IF ( coschi.LT.0.0 ) THEN
                  r1 = 0.
                  cota = 0.
               ELSE
                  CALL FRESNEL(Nr,Ni,coschi,sinchi,r1)
               ENDIF
!  Anisotropic Gaussian distribution for wave facets slopes
               zx = -tantetn*cosphin
               zy = -tantetn*sinphin
               xe = (COS(phw)*zx+SIN(phw)*zy)/SQRT(sigmac)
               xn = (-SIN(phw)*zx+COS(phw)*zy)/SQRT(sigmau)
               xe2 = xe*xe
               xn2 = xn*xn
               coef = 1 - c21/2.*(xe2-1)*xn - c03/6.*(xn2-3)*xn
               coef = coef + c40/24.*(xe2*xe2-6*xe2+3)
               coef = coef + c04/24.*(xn2*xn2-6*xn2+3)
               coef = coef + c22/4.*(xe2-1)*(xn2-1)
               fonc0 = 0.5*coschi*coef*EXP(-(xe2+xn2)/2.)/(costetn**4)
               pr = pr + r1*fonc0*cofa*cota*cotb
               pp = pp + fonc0*cofa*cota*cotb
            ENDDO
         ENDDO
!
         pond = 2.*hta*hfa*htb/pi/SQRT(sigmac)/SQRT(sigmau)/3./3./3.
         prefl = prefl + pr*pond
         proba = proba + pp*pond
   enddo

  rge = prefl/proba

end subroutine



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     PROGRAM: GETBOUND.FOR
!
!     DESCRIPTION: FINDS THE BOUNDS XVALS(IND1) AND XVALS(IND2)
!                  FOR WHICH X IS BOUNDED.
!
      SUBROUTINE GETBOUND(Xvals,Ifirst,Ilast,X,Ind1,Ind2)
      IMPLICIT NONE
!*--GETBOUND555
      REAL Xvals(*) , X
      INTEGER Ind1 , Ind2 , Ifirst , Ilast , i , imid

      imid = Ilast/2 + 1

      IF ( Xvals(Ilast).GT.Xvals(Ifirst) ) THEN
         IF ( X.GT.Xvals(imid) ) THEN
            DO i = imid , Ilast - 1
               IF ( Xvals(i).LE.X .AND. Xvals(i+1).GE.X ) THEN
                  Ind1 = i
                  Ind2 = i + 1
                  GOTO 99999
               ENDIF
            ENDDO
         ELSE
            DO i = Ifirst , imid
               IF ( Xvals(i).LE.X .AND. Xvals(i+1).GE.X ) THEN
                  Ind1 = i
                  Ind2 = i + 1
                  GOTO 99999
               ENDIF
            ENDDO
         ENDIF
         IF ( X.LT.Xvals(Ifirst) ) THEN
            Ind1 = Ifirst
            Ind2 = Ifirst + 1
         ELSE
            Ind1 = Ilast - 1
            Ind2 = Ilast
         ENDIF
!
      ELSE
!
         IF ( X.LT.Xvals(imid) ) THEN
            DO i = imid , Ilast - 1
               IF ( Xvals(i).GE.X .AND. Xvals(i+1).LE.X ) THEN
                  Ind1 = i
                  Ind2 = i + 1
                  GOTO 99999
               ENDIF
            ENDDO
         ELSE
            DO i = Ifirst , imid
               IF ( Xvals(i).GE.X .AND. Xvals(i+1).LE.X ) THEN
                  Ind1 = i
                  Ind2 = i + 1
                  GOTO 99999
               ENDIF
            ENDDO
         ENDIF
         IF ( X.GT.Xvals(Ifirst) ) THEN
            Ind1 = Ifirst
            Ind2 = Ifirst + 1
         ELSE
            Ind1 = Ilast - 1
            Ind2 = Ilast
         ENDIF
      ENDIF
99999 END
