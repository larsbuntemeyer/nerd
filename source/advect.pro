pro advect,x,xi,q,ui,dt,fluxlim,nghost
nx     = n_elements(xi)-1
if n_elements(x) ne nx then stop
if n_elements(xi) ne nx+1 then stop
if n_elements(q) ne nx then stop
if n_elements(ui) ne nx+1 then stop
if nghost lt 1 then stop
;
; Determine the r_{i-1/2} for the flux limiter
;
r = dblarr(nx+1)
for i=2,nx-2 do begin
dq = (q[i]-q[i-1])
if abs(dq) gt 0.d0 then begin
if(ui[i] ge 0.d0) then begin
r[i] = (q[i-1]-q[i-2])/dq
endif else begin
r[i] = (q[i+1]-q[i])/dq
endelse
endif
endfor
;
; Determine the flux limiter
; (many other flux limiters can be implemented here!)
;
case fluxlim of
   'donor-cell': begin
     phi = dblarr(nx+1)
   end
   'superbee': begin
     phi = dblarr(nx+1)
     for i=1,nx-1 do begin
       a = min([1.d0,2.d0*r[i]])
       b = min([2.d0,r[i]])
       phi[i] = max([0.d0,a,b])
     endfor
   end
   else: stop
endcase
;
; Now construct the flux
;
flux = dblarr(nx+1)
for i=1,nx-1 do begin
   if ui[i] ge 0.d0 then begin
      flux[i] = ui[i] * q[i-1]
   endif else begin
      flux[i] = ui[i] * q[i]
   endelse
   flux[i] = flux[i] + 0.5 * abs(ui[i]) * $
                      (1-abs(ui[i]*dt/(x[i]-x[i-1]))) * $ 
                       phi[i] * (q[i]-q[i-1])
endfor
;
; Update the cells, except the ghost cells
;
for i=nghost,nx-1-nghost do begin
q[i] = q[i] - dt*( flux[i+1]-flux[i] ) / ( xi[i+1] - xi[i] )
endfor
;
end

pro boundary,rho,rhou,periodic=periodic,mirror=mirror
;
; Get the number of grid points including the ghost cells ;
nx = n_elements(rho)
;
; If periodic, then install periodic BC, using two ghost cells
; on each side (two are required for the non-lin flux limiter) ;
if keyword_set(periodic) then begin
   rho[0]     = rho[nx-4]
   rho[1]     = rho[nx-3]
   rho[nx-2]  = rho[2]
   rho[nx-1]  = rho[3]
   rhou[0]    = rhou[nx-4]
   rhou[1]    = rhou[nx-3]
   rhou[nx-2] = rhou[2]
   rhou[nx-1] = rhou[3]
endif ;
; If mirror symmetry, then install mirror BC, using two ghost cells
; on each side (two are required for the non-lin flux limiter) ;
if keyword_set(mirror) then begin
   rho[0]     = rho[3]
   rho[1]     = rho[2]
   rho[nx-2]  = rho[nx-3]
   rho[nx-1]  = rho[nx-4]
   rhou[0]    = -rhou[3]
   rhou[1]    = -rhou[2]
   rhou[nx-2] = -rhou[nx-3]
   rhou[nx-1] = -rhou[nx-4]
endif
end

pro hydrostep,x,xi,rho,rhou,e,gamma,dt,periodic=periodic, $
    mirror=mirror,fluxlim=fluxlim,nrvisc=nrvisc
;
; Check for conflicting settings
;
if keyword_set(mirror) and keyword_set(periodic) then stop ;
; Use 2 ghost cells on each side
;
nghost = 2
;
; If not defined, install default flux limiter
;
if not keyword_set(fluxlim) then fluxlim='donor-cell'
;
; Get the number of grid points including the ghost cells ;
nx = n_elements(x)
;
; Impose boundary conditions
;
boundary,rho,rhou,periodic=periodic,mirror=mirror
;
; Compute the velocity at the cell interfaces
;
ui = dblarr(nx+1)
for ix=1,nx-1 do begin
  ui[ix] = 0.5 * ( rhou[ix]/rho[ix] + rhou[ix-1]/rho[ix-1] ) 
endfor
;
; Advect rho
;
advect,x,xi,rho,ui,dt,fluxlim,nghost
;
; Advect rho u
;
advect,x,xi,rhou,ui,dt,fluxlim,nghost
;
; Re-impose boundary conditions
;
boundary,rho,rhou,periodic=periodic,mirror=mirror
;
; Compute the pressure
;
p = (gamma-1.d0)*rho*e
;
; Now add the pressure force, for all cells except the ghost cells 
;
for ix=2,nx-3 do begin
  rhou[ix] = rhou[ix] - dt*(p[ix+1]-p[ix-1])/(x[ix+1]-x[ix-1])
endfor
;
; Re-impose boundary conditions a last time (not ; strictly necessary)
; boundary,rho,rhou,periodic=periodic,mirror=mirror ;
; Done
;
end



pro hydroiso_cen,x,xi,rho,rhou,e,gamma,dt
nx = n_elements(x)
;
; Compute the velocity at the cell interfaces 
;
ui = dblarr(nx+1)
for ix=1,nx-1 do begin
  ui[ix] = 0.5 * ( rhou[ix]/rho[ix] + rhou[ix-1]/rho[ix-1] ) 
endfor
;
; Compute the flux for rho
;
fluxrho = dblarr(nx+1)
for ix=1,nx-1 do begin
  if ui[ix] gt 0. then begin 
    fluxrho[ix] = rho[ix-1] * ui[ix]
  endif else begin
    fluxrho[ix] = rho[ix] * ui[ix]
  endelse 
end
;
; Update the density
;
for ix=0,nx-1 do begin
  rho[ix] = rho[ix] - (dt/(xi[ix+1]-xi[ix])) * $ 
            (fluxrho[ix+1]-fluxrho[ix])
endfor
;
; Compute the flux for rho u ;
fluxrhou = dblarr(nx+1)
for ix=1,nx-1 do begin
  if ui[ix] gt 0. then begin
    fluxrhou[ix] = rhou[ix-1]^2 / rho[ix-1]
  endif else begin
    fluxrhou[ix] = rhou[ix]^2 / rho[ix]
  endelse
end
;
; Update the momentum
;
for ix=0,nx-1 do begin
  rhou[ix] = rhou[ix] - (dt/(xi[ix+1]-xi[ix])) * $ 
              (fluxrhou[ix+1]-fluxrhou[ix])
endfor
;
; Compute the pressure
;
p    = (gamma-1.d0)*rho*e
;
; Now add the pressure force, for all cells
; except the ones near the boundary ;
for ix=1,nx-2 do begin
  rhou[ix] = rhou[ix] - dt*(p[ix+1]-p[ix-1])/(x[ix+1]-x[ix-1]) 
endfor
;
; Now do the boundary cells, assuming mirror
; symmetry in the boundaries
;
rhou[0] = rhou[0] - 0.5*dt*(p[1]-p[0])/(x[1]-x[0])
rhou[nx-1] = rhou[nx-1] - 0.5*dt*(p[nx-1]-p[nx-2])/(x[nx-1]-x[nx-2]) 
;
; Done
;
end




nx = 100
nt = 1000
x0 = 0.d0
x1 = 100.d0
xmid = 0.5 * (x0+x1)
dt = 0.25
cfl = 0.5
x = x0 + (x1-x0)*(dindgen(nx)/(nx-1.d0)) 
gamma = 7./5.
rho = dblarr(nx,nt+1)
rhou = dblarr(nx,nt+1)
e = dblarr(nx)+1.d0
time = dblarr(nt+1)
dg = 0.1*(x1-x0)
rho[*,0] = 1.d0+0.3*exp(-(x-xmid)^2/dg^2) 
;rho[*,0] = 1.d0
;rho[50:nx-1,0] = 0.d0


;
; Now some additional arrays are set up
;
xi = dblarr(nx+1)
xi[1:nx-1] = 0.5 * ( x[1:nx-1] + x[0:nx-2] ) 
xi[0] = 2*xi[1] - xi[2]
xi[nx] = 2*xi[nx-1] - xi[nx-2]
dx = ( xi[1:nx] - xi[0:nx-1] )

;
; Now the hydro is done
;
for it=1,nt do begin
  qrho = rho[*,it-1]
  qrhou = rhou[*,it-1]
  cs = sqrt(gamma*(gamma-1)*e)
  dum = dx/(cs+abs(qrhou/qrho))
  dt = cfl*min(dum)
  time[it] = time[it-1]+dt
  ;;
  print,'Time step ',it,', Time = ',time[it],', Dt = ',dt ;;
  hydrostep,x,xi,qrho,qrhou,e,gamma,dt,periodic=1, $
    fluxlim='superbee',nrvisc=1
  ;hydroiso_cen,x,xi,qrho,qrhou,e,gamma,dt
  plot,  x, qrho, psym=-6;,yrange=[0.0,dens_left*1.1]
  ;;
  rho[*,it] = qrho
  rhou[*,it] = qrhou
endfor


END
