; IDL script for plotting SHDOM radiance output files from run_brdf.

BRDF = 'L'

radmin=0.0                   ; minimum I radiance
if BRDF eq 'L' then begin
  radfile='brdf_L1r.out'     ; SHDOM radiance output file
  radmax=0.10                ; maximum I radiance
  title='SHDOM BRDF Example: Lambertian Surface'
  psfile='brdf_exampleL.ps'
endif
if BRDF eq 'R' then begin
  radfile='brdf_R1r.out'
  radmax=0.20
  title='SHDOM BRDF Example: RPV Surface'
  psfile='brdf_exampleR.ps'
endif
if BRDF eq 'O' then begin
  radfile='brdf_O1r.out'
  radmax=0.24
  title='SHDOM BRDF Example: Ocean Surface'
  psfile='brdf_exampleO.ps'
endif
if BRDF eq 'W' then begin
  radfile='brdf_W1r.out'
  radmax=0.20                  ; maximum I radiance
  qmin=-0.2                    ; minimum Q radiance
  qmax=0.0                    ; maximum Q radiance
  title='SHDOM BRDF Example: Wave-Fresnel Surface'  ; title for plot
  psfile='brdf_exampleW.ps'
endif
if BRDF eq 'D' then begin
  radfile='brdf_D1r.out'
  radmax=0.16
  qmin=-0.07
  qmax=0.01
  title='SHDOM BRDF Example: Diner Surface'
  psfile='brdf_exampleD.ps'
endif


; This script outputs Postscript formats: 
;   Postscript ('ps'), Encapsulated Postscript ('eps')
out='ps'

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


read_radpol, radfile, nstokes, nx,ny,ndir, xgrid,ygrid,zlev, mu, phi, radiance

zamin=-90 &  zamax=90  ; zenith angle range to plot (<0 for phi=180, >0 for phi=0)
iy=0                   ; Y index slice (0 for 2D field)
 ; Convert mu to zenith angle and negate for phi=180
zenang = 180*acos(mu)/acos(-1.0)
for i = 0, ndir-1 do begin
  if phi[i] gt 90. then zenang[i] = -zenang[i]
endfor

plotdata1 = reform(radiance[0,*,iy,*])

xyouts, 0.50,0.96, title, /normal, size=1.2*!p.charsize, alignment=0.5

 ; Calculate size of plot
x1plot=0.15
delxplot=0.70
y1plot=0.64
delyplot = 0.30

  ; Make color contour plot with color bar of I radiance field
  contlev=radmin+((radmax-radmin)/nlev)*findgen(nlev+1)
  plotdata1 = 0.999*radmin > (0.999*radmax < plotdata1)
  lablev=[contlev(0),contlev(nlev/4),contlev(nlev/2),contlev(3*nlev/4),contlev(nlev)]
  contour, plotdata1, xgrid, zenang, $
	level=contlev,c_color=cola,/fill, /closed,/follow, /noerase, $
        /xstyle, yrange=[zamin,zamax],/ystyle, $
	position=[x1plot,y1plot,x1plot+delxplot,y1plot+delyplot], $
	title='', xtitle='BRDF Parameter Interpolation Fraction', ytitle='Zenith Angle'
  makekey, x1plot,y1plot-0.09, delxplot,0.02, 0.0,-0.025, $
        labels=string(lablev,format='(F5.3)'), colors=cola, charsize=!p.charsize
  xyouts, 0.50,y1plot-0.14, 'I Radiance', size=!p.charsize, $
	/normal, alignment=0.5

if nstokes gt 1 then begin
  plotdata2 = reform(radiance[1,*,iy,*])
   ;  Make color contour plot with color bar of Q radiance field
  y1plot=0.16
  contlev=qmin+((qmax-qmin)/nlev)*findgen(nlev+1)
  plotdata2 = 0.999*qmin > (0.999*qmax < plotdata2)
  lablev=[contlev(0),contlev(nlev/4),contlev(nlev/2),contlev(3*nlev/4),contlev(nlev)]
  contour, plotdata2, xgrid, zenang, $
	level=contlev,c_color=cola,/fill, /closed,/follow, /noerase, $
        /xstyle, yrange=[zamin,zamax],/ystyle, $
	position=[x1plot,y1plot,x1plot+delxplot,y1plot+delyplot], $
	title='', xtitle='BRDF Parameter Interpolation Fraction', ytitle='Zenith Angle'
  makekey, x1plot,y1plot-0.09, delxplot,0.02, 0.0,-0.025, $
        labels=string(lablev,format='(F6.3)'), colors=cola, charsize=!p.charsize
  xyouts, 0.50,y1plot-0.14, 'Q Radiance', size=!p.charsize, $
	/normal, alignment=0.5
endif

device, /close

end
