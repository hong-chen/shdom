; IDL script for plotting "run_polarized_example" files.

;   The following plot types are available:
; plottype='xz_lwcre'   X-Z  liquid water content and effective radius
; plottype='xz_extalb'  X-Z  extinction and single scattering albedo
; plottype='radiance'  radiance for several directions vs X
plottype='xz_lwcre'


; This script outputs Postscript formats: 
;   Postscript ('ps'), Encapsulated Postscript ('eps')
out='ps'
psfile='pol_example_'+plottype+'.'+out

if (out eq 'ps') then begin
  ;  For Postscript use device fonts (pick Times Roman)
  set_plot, 'ps'
  device, filename=psfile, /color, /portrait, /times, $
   	 xsize=18, ysize=24, xoffset=0.0, yoffset=0.0
endif
if (out eq 'eps') then begin
  set_plot, 'ps'
  device, filename=psfile, /color, /portrait, /encapsulated, /times, $
   	xsize=18.0, ysize=24.0, xoffset=1.0, yoffset=1.0
endif
nlev=40  ; number of color levels
cola=(0.95*!d.n_colors/(1.0*nlev))*(indgen(nlev)+0.5)
!p.color = 0
!p.font=0
!p.charsize=1.2

;  Select the color table
;loadct, 0    ; gray scale
loadct, 39    ; blue to red color

; Set other global plotting variables
!p.thick=1.0  &  !x.thick=1.0  &  !y.thick=1.0
!x.ticklen = -0.013  &  !y.ticklen = -0.010




;  Plot X-Z slice of liquid water content and effective radius for one
;   particle type (component).    Reads a particle properties file.
if (plottype eq 'xz_lwcre') then begin
  lwcfile='aero0822nh15c_t13_iy249.lwc'  ; LWC file
  partype=-1                  ; flag to read LWC file instead of particle file
  parfile='les_cumulus_dust.part'   ; particle file 
  partype=1                         ; particle type to extract from file

  zmin=0.0  &  zmax=5.0          ; height range to plot
  if (abs(partype) eq 1) then begin
    parname='Water Droplets'     ; particle type name
    lwcmax=2.0                   ; max liquid water content
    lwcform='(F4.2)'
    reffmin=4.0 & reffmax=24.0   ; min and max effective radius
  endif
  if (partype eq 2) then begin
    parname='Dust Aerosol'       ; particle type name
    lwcmax=1.0E-4                ; max liquid water content
    lwcform='(E8.1)'
    reffmin=0.0 & reffmax=1.0    ; min and max effective radius
  endif
  iy=0                         ; Y index slice (0 for 2D field)
  title='Example LWC/r!Deff!N: '+parname   ; title for plot

  if (partype eq -1) then begin
    read_lwc, lwcfile, nx,ny,nz,xgrid,ygrid,zgrid, temp,lwc,reff
  endif else begin
    read_particle, parfile, partype, nx,ny,nz,xgrid,ygrid,zgrid, temp,lwc,reff
  endelse

  plotdata1 = reform(lwc[*,iy,*])
  plotdata2 = reform(reff[*,iy,*])

  xyouts, 0.50,0.95, title, /normal, size=1.2*!p.charsize, alignment=0.5

  ;  Calculate size of plot
  x1plot=0.15
  delxplot=0.70
  y1plot=0.47+0.18
  delyplot = 0.25

  ;  Make color contour plot with color bar of LWC field
  contlev=0.0+((lwcmax-0.0)/nlev)*findgen(nlev+1)
  plotdata1 = (0.999*lwcmax < plotdata1)
  lablev=[contlev(0),contlev(nlev/4),contlev(nlev/2),contlev(3*nlev/4),contlev(nlev)]
  contour, plotdata1, xgrid, zgrid, $
	level=contlev,c_color=cola,/fill, /closed,/follow, /noerase, $
        /xstyle, yrange=[zmin,zmax],/ystyle, $
	position=[x1plot,y1plot,x1plot+delxplot,y1plot+delyplot], $
	title='', xtitle='X (km)', ytitle='Z (km)'
  makekey, x1plot,y1plot-0.10, delxplot,0.02, 0.0,-0.025, $
        labels=string(lablev,format=lwcform), colors=cola, charsize=!p.charsize
  xyouts, 0.50,y1plot-0.15, 'Liquid Water Content (g/m!U3!N)', size=1.2*!p.charsize, $
	/normal, alignment=0.5


  ;  Make color contour plot with color bar of effective radius field
  y1plot=0.18
  contlev=reffmin+((reffmax-reffmin)/nlev)*findgen(nlev+1)
  plotdata2 = 0.999*reffmin > (0.999*reffmax < plotdata2)
  lablev=[contlev(0),contlev(nlev/4),contlev(nlev/2),contlev(3*nlev/4),contlev(nlev)]
  contour, plotdata2, xgrid, zgrid, $
	level=contlev,c_color=cola,/fill, /closed,/follow, /noerase, $
        /xstyle, yrange=[zmin,zmax],/ystyle, $
	position=[x1plot,y1plot,x1plot+delxplot,y1plot+delyplot], $
	title='', xtitle='X (km)', ytitle='Z (km)'
  makekey, x1plot,y1plot-0.10, delxplot,0.02, 0.0,-0.025, $
        labels=string(lablev,format='(F5.2)'), colors=cola, charsize=!p.charsize
  xyouts, 0.50,y1plot-0.15, 'Effective Radius (!9m!7m)', size=1.2*!p.charsize, $
	/normal, alignment=0.5

endif





;  Plot X-Z slice of extinction and single scattering albedo 
;      Reads SHDOM input property file
if (plottype eq 'xz_extalb') then begin
  prpfile='les_cumulus_dust_w0.646.prp'  ; SHDOM property file to plot
  zmin=0.0  &  zmax=5.0             ; height range to plot
  extinctmin=-1.8                   ; minimum log10 extinction to plot
  extinctmax=2.2                    ; maximum log10 extinction to plot
  ssalbmin=0.84                     ; minimum single scattering albedo
  ssalbmax=1.00                     ; maximum single scattering albedo
  iy=0                              ; Y index slice (0 for 2D field)
  title='LES cloud/dust example: 0.646 !9m!7m'   ; title for plot

  ;  Read in property file
  read_prp, prpfile, nx,ny,nz, xgrid,ygrid,zgrid, temps,extinct,ssalb
  plotdata1 = alog10(1.0E-4 > reform(extinct(*,iy,*)))
  plotdata2 = reform(ssalb(*,iy,*))

  ;  Plot title
  xyouts, 0.50,0.96, title, /normal, size=1.4*!p.charsize, alignment=0.5

  ;  Calculate size of plot: X and Z are to scale if 4:3 aspect ratio output
  x1plot=0.15
  delxplot=0.70
  y1plot=0.47+0.18
  delyplot = 0.25

  ;  Make color contour plot with color bar of extinction field
  contlev=extinctmin+((extinctmax-extinctmin)/nlev)*findgen(nlev+1)
;  plotdata1 = 1.001*extinctmax > (0.999*extinctmax < plotdata1)
  lablev=[contlev(0),contlev(nlev/4),contlev(nlev/2),contlev(3*nlev/4),contlev(nlev)]
  contour, plotdata1, xgrid, zgrid, $
	level=contlev,c_color=cola,/fill, /closed,/follow, /noerase, $
        /xstyle, yrange=[zmin,zmax],/ystyle, $
	position=[x1plot,y1plot,x1plot+delxplot,y1plot+delyplot], $
	title='', xtitle='X (km)', ytitle='Z (km)'
  makekey, x1plot,y1plot-0.10, delxplot,0.02, 0.0,-0.025, $
        labels=string(lablev,format='(F4.1)'), colors=cola, charsize=!p.charsize
  xyouts, 0.50,y1plot-0.15, 'Log!D10!N Extinction (km!U-1!N)', size=1.2*!p.charsize, $
	/normal, alignment=0.5


  ;  Make color contour plot with color bar of single scattering albedo
  y1plot=0.18
  contlev=ssalbmin+((ssalbmax-ssalbmin)/nlev)*findgen(nlev+1)
  plotdata2 = 1.001*ssalbmin > (0.999*ssalbmax < plotdata2)
  lablev=[contlev(0),contlev(nlev/4),contlev(nlev/2),contlev(3*nlev/4),contlev(nlev)]
  contour, plotdata2, xgrid, zgrid, $
	level=contlev,c_color=cola,/fill, /closed,/follow, /noerase, $
        /xstyle, yrange=[zmin,zmax],/ystyle, $
	position=[x1plot,y1plot,x1plot+delxplot,y1plot+delyplot], $
	title='', xtitle='X (km)', ytitle='Z (km)'
  makekey, x1plot,y1plot-0.10, delxplot,0.02, 0.0,-0.025, $
        labels=string(lablev,format='(F4.2)'), colors=cola, charsize=!p.charsize
  xyouts, 0.50,y1plot-0.15, 'Single scattering albedo', size=1.2*!p.charsize, $
	/normal, alignment=0.5

endif






;  Plot I, Q, U radiances for several directions vs X
;       Reads SHDOM radiance file
if (plottype eq 'radiance') then begin
  radfile='les_cumulus_dust_w0.646_sfcW_ns3.rad' ; SHDOM radiance output file
  nstokes=3                    ; number of Stokes parameters in radiance file
  iy=0                         ; Y index slice (0 for 2D field)
  idir=[1,4,7]                 ; indices for directions
  radmin=0.0                   ; minimum I radiance
  radmax=0.35                  ; maximum I radiance
  qumin=-0.02                  ; minimum Q, U radiance
  qumax=0.02                   ; maximum Q, U radiance
  title1='LES Cloud/Dust Example 0.646 !9m!7m Radiances at 5 km'  ; title for plot

  read_radpol, radfile, nstokes, nx,ny,ndir, xgrid2,ygrid,zlev, mu, phi, radiance
  nrad=n_elements(idir)
  Irad=fltarr(nx,nrad)
  Qrad=fltarr(nx,nrad)
  Urad=fltarr(nx,nrad)
  for i = 0, nrad-1 do begin
    Irad(*,i) = radiance[0,*,iy,idir[i]]
    Qrad(*,i) = radiance[1,*,iy,idir[i]]
    Urad(*,i) = radiance[2,*,iy,idir[i]]
  endfor

  ;  Plot the I radiances vs X
  cols=FIX([0.25,0.60,0.95]*!d.n_colors)
  plot, xgrid2, Irad[*,0], /xstyle, /ystyle, yrange=[radmin,radmax], $
	xtitle='X (km)', ytitle='I Radiance', /noerase, /nodata, $
	position=[0.14,0.55,0.97,0.93], title=title2
  for i = 0, nrad-1 do $
    oplot, xgrid2, Irad[*,i], linestyle=0, color=cols[i]

  ;  Plot legend for the radiance directions
  smu='!9m!7'   &  sphi='!9f!7'
  dx=xgrid2[nx-1]-xgrid2[0]
  x0=xgrid2[0]+0.0*dx
  x1=x0+0.02*dx  &  x2=x0+0.08*dx  &  x3=x0+0.10*dx
  y0=radmax      &  dy=radmax-radmin  & ey=0.02*dy
  for i = 0, nrad-1 do begin
    y1=y0-(0.06*i+0.05)*dy
    oplot, [x1,x2], [y1,y1], linestyle=0, color=cols[i]
    xyouts, x3, y1-ey, 'Direction: ' $
	+smu+'='+string(mu(idir(i)),format='(F5.3)')+'  ' $
	+sphi+'='+string(fix(phi(idir(i))),format='(I3)'), $
	charsize=1.0*!p.charsize
  endfor

  ;  Plot the Q and U radiances vs X
  plot, xgrid2, Qrad[*,0], /xstyle, /ystyle, yrange=[qumin,qumax], $
	xtitle='X (km)', ytitle='Q, U Radiance', /noerase, /nodata, $
	position=[0.14,0.09,0.97,0.47], title=''
  for i = 0, nrad-1 do $
    oplot, xgrid2, Qrad[*,i], linestyle=0, color=cols[i]
  for i = 0, nrad-1 do $
    oplot, xgrid2, Urad[*,i], linestyle=2, color=cols[i]

  ;  Plot legend indicating Q and U
  dx=xgrid2[nx-1]-xgrid2[0]
  x0=xgrid2[0]+0.0*dx
  x1=x0+0.02*dx  &  x2=x0+0.08*dx  &  x3=x0+0.10*dx
  y0=qumax      &  dy=qumax-qumin  & ey=0.02*dy
  y1=y0-0.05*dy
  oplot, [x1,x2], [y1,y1], linestyle=0, color=0
  xyouts, x3, y1-ey, 'Q', charsize=1.0*!p.charsize
  y1=y0-0.11*dy
  oplot, [x1,x2], [y1,y1], linestyle=2, color=0
  xyouts, x3, y1-ey, 'U', charsize=1.0*!p.charsize

endif


device, /close

end
