; IDL procedure  -  read_radpol.pro
;
;****************************************************************************
;
; NAME:
;      read_radpol
;
; PURPOSE:
;      Read an SHDOM-produced polarized radiance file.
;
; CATEGORY:
;      SHDOM utility
;
; CALLING SEQUENCE:
;      read_radpol, filename, nstokes, nx, ny, ndir, xgrid, ygrid, zlev, 
;                             mu, phi, radiance
;
; INPUTS:
;      filename:  a string with the directory path and file name
;
; OUTPUTS:
;      nstokes:    number of Stokes components to radiance
;      nx, ny:     output grid size
;      ndir:       number of output radiance directions
;      xgrid:      X grid locations  [fltarr(nx)]
;      ygrid:      Y grid locations  [fltarr(ny)]
;      zlev:       height of Z level
;      mu:         array of mu directions (cosine zenith angle)
;      phi:        array of phi directions (azimuth in degrees)
;      radiance:   radiance array [fltarr(nstokes,nx,ny,ndir)]
;
;
; MODIFICATION HISTORY:
;      Written by:   Frank Evans, July 2013
;
; NOTES:
;      For use with files created by the Spherical Harmonic Discrete
;        Ordinate Method (SHDOM) (Evans, 1998) radiative transfer model.
;
;
;****************************************************************************

PRO read_radpol, filename, nstokes, nx, ny, ndir, xgrid, ygrid, zlev, $
                           mu, phi, radiance

;****************************************************************************

;  Open file and get the base grid size from the header
ss = ' '
openr,lun, filename, /get_lun
readf,lun, ss
readf,lun, ss
readf,lun, ss
nstokes = fix(strmid(ss,strpos(ss,'NSTOKES=')+8,1))
;nx = fix(strmid(ss,strpos(ss,'NX=')+3,4))
;ny = fix(strmid(ss,strpos(ss,'NY=')+3,4))
;nz = fix(strmid(ss,strpos(ss,'NZ=')+3,4))

;  Find the Z level and number of radiance outputs
repeat begin
  readf,lun, ss
endrep until (strpos(ss,'RADIANCE AT') ne -1)
zlev = float(strmid(ss,strpos(ss,'Z=')+2,7))
nx = fix(strmid(ss,strpos(ss,'NXO=')+4,4))
ny = fix(strmid(ss,strpos(ss,'NYO=')+4,4))
ndir = fix(strmid(ss,strpos(ss,'NDIR=')+5,4))
readf,lun, ss

;  Read each block (direction) of radiances 
mu = fltarr(ndir)
phi = fltarr(ndir)
radiance = fltarr(nstokes,nx,ny,ndir)
n = long(nx)*long(ny)
buffer = fltarr(2+nstokes,n)
for i = 0, ndir-1 do begin
    readf,lun, ss
    mu(i) = float(strmid(ss,3,8))
    phi(i) = float(strmid(ss,12,6))
    readf,lun, buffer
    temp=reform(buffer(0,*),nx,ny)
    xgrid=reform(temp(*,0),nx)
    temp=reform(buffer(1,*),nx,ny)
    ygrid=reform(temp(0,*),ny)
    radiance[*,*,*,i] = reform(buffer(2:1+nstokes,*),nstokes,nx,ny)
endfor

free_lun, lun

END

