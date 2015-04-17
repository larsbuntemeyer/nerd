;
; a simple real time roe solver
; =====================================
; I simply wrote this script to
; understand the basics of
; a characteristics Riemann solver.
;
; simply run it with .r roe_fl.pro should
; do the job. have a look at basic
; parameters further down!
; most important to tweak are probably
; the CFL number and the flux limiter...
; Enjoy!
;
;   /\_/\
;  ( o.o )
;   > ^ <
;
; 04/13 - Lars Buntemeyer 
;

function SIGN, X1,X2

on_error, 1				; On error return to main level
if n_elements(X1) eq 0 or n_elements(X2) eq 0 then message, 'Both arguments must exist

A = n_elements(X1)  
B = n_elements(X2)

if A ne B and B ne 1 then message, 'Dimension error
OUT = abs(X1)
A = where(X2 lt 0)
if A(0) ne -1 then if B ne 1 then OUT(A) = -OUT(A) else OUT = -OUT
return, OUT  
end


;simple gamma eos
pro eos,dens,pres,ener,eint,enth,velx,vely,velz,gamma
   eint = ener - 0.5*(velx^2+vely^2+velz^2)
   pres = dens*eint*(gamma-1.0)
   ener = eint + 0.5*(velx^2+vely^2+velz^2)
   enth = ener+pres/dens
end

;Roe's average states
pro average,q,qa,w,n
   for i=1,n-1 do begin
      qa(i) = (q(i-1)*w(i-1)+q(i)*w(i))/(w(i-1)+w(i)) 
   endfor
   qa(0) = qa(1)
   qa(n) = qa(n-1)
end

;minmod
function minmod,a,b
   c = 0.d0
   if(a*b gt 0.d0) then begin
      if(abs(a) lt abs(b)) then begin
         c = a
      endif else begin
         c = b
      endelse
   endif else begin
      c = 0.d0
   endelse
   return, c 
end

;hyperbee
function hyperbee,r
   fl_hyperbee = 0.d0
   a = 0.1d0
   eps = 0.d0
   if(r le eps) then begin
      fl_hyperbee = 0.d0
   endif else if(abs(r-1.d0) le eps) then begin
      fl_hyperbee = 1.d0
   endif else begin
      fl_hyperbee = 2.d0*r/(a*(1.d0-a)) * $
          (a*(r-1.d0)+(1.d0-r^a))/(r-1.d0)^2
   endelse
;;   case r of
;;      0.d0 : fl_hyperbee = 0.d0
;;      1.d0 : fl_hyperbee = 1.d0
;;      else : fl_hyperbee = 2.d0*r/(a*(1.d0-a)) * $
;;              (a*(r-1.d0)+(1.d0-r^a))/(r-1.d0)^2
;;   endcase
   return, fl_hyperbee
end

;flux limiter
function phi,fl,r
   limiter = 0.d0
   case fl of
   'donor-cell'  : limiter = 0.d0
   'Lax-Wendroff': limiter = 1.d0
   'Beam-Warming': limiter = r
   'Fromm'       : limiter = 0.5d0*(1.d0+r)
   'minmod'      : limiter = minmod(1.d0,r)
   'superbee'    : limiter = max([0.d0,min([1.d0,2.d0*r]),min([2.d0,r])])
   'hyperbee'    : limiter = hyperbee(r) 
   'MC'          : limiter = max([0.d0,min([0.5d0*(1.d0+r),2.d0,2.d0*r])])
   'van Leer'    : limiter = (r+abs(r))/(1.d0+abs(r)) 
   'van Albada 1': limiter = (r^2+r)/(r^2+1.d0) 
   'van Albada 2': limiter = (2.d0*r)/(r^2+1.d0) 
   endcase
   return, limiter 
end

;==================
;main program
;==================
; some parameters
dt = 0.1
dt_max = 200.0
dt_min = 0.01
t0 = 0.0
tmax = 0.2
timesteps = 10000
x0 = 0.0
x1 = 1.0
xc = 0.3*(x1-x0)
n  = 100
k  = 1.0
gamma = 7.0/5.0
; have a look above at the available flux limiter 
fl = 'superbee'
cfl = 0.8
; set ps = 1 to write a ps file
; set wf = 1 to write a data file of the end results
; if both are 0, we plot in real time 
ps = 0
wf = 0

;Toro Test 1
dens_left = 1.0
dens_right = 0.125
pres_left = 1.0
pres_right = 0.1
velx_left = 0.75
velx_right = 0.0

;Toro Test 2
;dens_left = 1.0
;dens_right = 1.0
;pres_left = 0.4
;pres_right = 0.4
;velx_left = -2.0
;velx_right = 2.0

;Toro Test 3
;dens_left = 1.0
;dens_right = 1.0
;pres_left = 1000.0
;pres_right = 0.01
;velx_left = 0.0
;velx_right = 0.0

;Toro Test 4
;dens_left = 5.99924
;dens_right = 5.99242
;pres_left = 460.894
;pres_right = 46.0950
;velx_left = 19.5975
;velx_right = -6.19633

;Toro Test 5
;dens_left = 1.0 
;dens_right = 1.0 
;pres_left = 1000.0 
;pres_right = 0.01
;velx_left = -19.59745
;velx_right = -19.59745
;==================

;dx
dx = (x1-x0)/n

;hydro variables at cell centers
dens = dindgen(n)
pres = dindgen(n)
velx = dindgen(n)
vely = dindgen(n)
velz = dindgen(n)
ener = dindgen(n)
eint = dindgen(n)
enth = dindgen(n)

;averaged hydro variables at cell faces
dens_a = dindgen(n+1)
pres_a = dindgen(n+1)
velx_a = dindgen(n+1)
vely_a = dindgen(n+1)
velz_a = dindgen(n+1)
ener_a = dindgen(n+1)
eint_a = dindgen(n+1)
enth_a = dindgen(n+1)

;conserved quantities at cell centers
q1 = dindgen(n)
q2 = dindgen(n)
q3 = dindgen(n)
q4 = dindgen(n)
q5 = dindgen(n)

;jumps in conserved quantities at cell faces
dq1 = dindgen(n+1)
dq2 = dindgen(n+1)
dq3 = dindgen(n+1)
dq4 = dindgen(n+1)
dq5 = dindgen(n+1)

;averaged states at cell interfaces
q1a = dindgen(n+1)
q2a = dindgen(n+1)
q3a = dindgen(n+1)
q4a = dindgen(n+1)
q5a = dindgen(n+1)

;flux at cell interfaces
f1 = dindgen(n+1)
f2 = dindgen(n+1)
f3 = dindgen(n+1)
f4 = dindgen(n+1)
f5 = dindgen(n+1)
f1l = dindgen(n+1)
f2l = dindgen(n+1)
f3l = dindgen(n+1)
f4l = dindgen(n+1)
f5l = dindgen(n+1)
f1r = dindgen(n+1)
f2r = dindgen(n+1)
f3r = dindgen(n+1)
f4r = dindgen(n+1)
f5r = dindgen(n+1)

;eigenvelocities
l1 = dindgen(n+1)
l2 = dindgen(n+1)
l3 = dindgen(n+1)
l4 = dindgen(n+1)
l5 = dindgen(n+1)

;eigenvelocities
a1 = dindgen(n+1)
a2 = dindgen(n+1)
a3 = dindgen(n+1)
a4 = dindgen(n+1)
a5 = dindgen(n+1)

;eigenvectors
K1 = dindgen(n+1,5)
K2 = dindgen(n+1,5)
K3 = dindgen(n+1,5)
K4 = dindgen(n+1,5)
K5 = dindgen(n+1,5)

;flux limiter
phi = dindgen(n+1,5)
r1 = dindgen(n+1)
r2 = dindgen(n+1)
r3 = dindgen(n+1)
r4 = dindgen(n+1)
r5 = dindgen(n+1)
t1 = dindgen(n+1)
t2 = dindgen(n+1)
t3 = dindgen(n+1)
t4 = dindgen(n+1)
t5 = dindgen(n+1)
e1 = dindgen(n+1)
e2 = dindgen(n+1)
e3 = dindgen(n+1)
e4 = dindgen(n+1)
e5 = dindgen(n+1)
p1 = dindgen(n+1)
p2 = dindgen(n+1)
p3 = dindgen(n+1)
p4 = dindgen(n+1)
p5 = dindgen(n+1)
t1v = dindgen(n+1)
t2v = dindgen(n+1)
t3v = dindgen(n+1)
t4v = dindgen(n+1)
t5v = dindgen(n+1)

dti = dindgen(n+1)

;dissipative flux
diss = dindgen(n+1,5)

;velocity at cell interaces
uf = dindgen(n+1)
;x-grid
x = dindgen(n)
;average weighting
w = dindgen(n)

;prepare x-grid
x(0) = x0
for i=1,n-1 do begin
   x(i) = x(i-1) + dx
endfor

;initial conditions
for i=0,n-1 do begin
   if(x(i) lt xc) then begin
      dens(i) = dens_left
      pres(i) = pres_left
      velx(i) = velx_left
   endif else begin 
      dens(i) = dens_right
      pres(i) = pres_right
      velx(i) = velx_right
   endelse
   vely(i) = 0.0
   velz(i) = 0.0
   eint(i) = pres(i)/(dens(i)*(gamma-1.0))
   ener(i) = eint(i) + 0.5*(velx(i)^2+vely(i)^2+velz(i)^2) 
endfor

;eos
eos,dens,pres,ener,eint,enth,velx,vely,velz,gamma

time = 0.0

for t=1,timesteps do begin

   print, 'step',t,'          time',time,'              timestep',dt

   ;conserved quantities
   for i=0,n-1 do begin
      q1(i) = dens(i) 
      q2(i) = dens(i)*velx(i) 
      q3(i) = dens(i)*vely(i) 
      q4(i) = dens(i)*velz(i) 
      q5(i) = ener(i)*dens(i)
   endfor
   ;jumps in conserved quantities
   dq1(0) = 0.0 
   dq2(0) = 0.0 
   dq3(0) = 0.0 
   dq4(0) = 0.0 
   dq5(0) = 0.0 
   dq1(n) = 0.0 
   dq2(n) = 0.0 
   dq3(n) = 0.0 
   dq4(n) = 0.0 
   dq5(n) = 0.0 
   for i=1,n-1 do begin
      dq1(i) = q1(i)-q1(i-1) 
      dq2(i) = q2(i)-q2(i-1) 
      dq3(i) = q3(i)-q3(i-1) 
      dq4(i) = q4(i)-q3(i-1) 
      dq5(i) = q5(i)-q5(i-1) 
   endfor

   ;average 
   average,velx,velx_a,sqrt(dens),n
   average,vely,vely_a,sqrt(dens),n
   average,velz,velz_a,sqrt(dens),n
   average,enth,enth_a,sqrt(dens),n
   vel_a = sqrt(velx_a^2+vely_a^2+velz_a^2)
   cs_a = sqrt((gamma-1.0)*((enth_a-0.5*vel_a^2)))

   ;average eigenvalues
   l1 = velx_a - cs_a
   l2 = velx_a
   l3 = velx_a
   l4 = velx_a
   l5 = velx_a + cs_a

   ;determine timestep from eigenvalues
   for i=0,n-1 do begin
      lmax = max([0.d0,l1(i),l2(i),l3(i),l4(i),l5(i)]) 
      lmin = min([0.d0,l1(i),l2(i),l3(i),l4(i),l5(i)])
      if((lmax-lmin) gt 0.d0) then begin
          dti(i) = dx/(lmax-lmin)
      endif else begin
          dti(i) = dt_max
      endelse
   endfor

   dt = cfl * min(dti)
   if(dt gt dt_max) then dt = dt_max


   ;average eigenvectors
   K1(*,0) = 1.0  
   K1(*,1) = velx_a - cs_a  
   K1(*,2) = vely_a
   K1(*,3) = velz_a
   K1(*,4) = enth_a - velx_a*cs_a

   K2(*,0) = 1.0  
   K2(*,1) = velx_a
   K2(*,2) = vely_a
   K2(*,3) = velz_a
   K2(*,4) = 0.5*vel_a^2 

   K3(*,0) = 0.0  
   K3(*,1) = 0.0 
   K3(*,2) = 1.0 
   K3(*,3) = 0.0 
   K3(*,4) = vely_a 

   K4(*,0) = 0.0  
   K4(*,1) = 0.0 
   K4(*,2) = 0.0 
   K4(*,3) = 1.0 
   K4(*,4) = velz_a 

   K5(*,0) = 1.0  
   K5(*,1) = velx_a + cs_a 
   K5(*,2) = vely_a 
   K5(*,3) = velz_a
   K5(*,4) = enth_a + velx_a*cs_a 

   ;wavestrengths
   a3 = dq3-vely_a*dq1
   a4 = dq4-velz_a*dq1
   dq5a = dq5-(dq3-vely_a*dq1)*vely_a-(dq4-velz_a*dq1)*velz_a
   a2 = (gamma-1.0)/(cs_a^2) * (dq1*(enth_a-velx_a^2)+velx_a*dq2-dq5a)
   a1 = 1.0/(2.0*cs_a) * (dq1*(velx_a+cs_a)-dq2-cs_a*a2)
   a5 = dq1-(a1+a2)

   ;
   ; flip flop function
   for i=0,n do begin
      t1v(i) = sign(1.0,velx_a(i))
      t2v(i) = sign(1.0,velx_a(i))
      t3v(i) = sign(1.0,velx_a(i))
      t4v(i) = sign(1.0,velx_a(i))
      t5v(i) = sign(1.0,velx_a(i))
      t1v(i) = 0.0 
      t2v(i) = 0.0 
      t3v(i) = 0.0 
      t4v(i) = 0.0 
      t5v(i) = 0.0 
   endfor 
   for i=1,n-1 do begin
      f1l(i) = (1.d0+t1v(i))*q1(i-1)*velx(i-1)
      f2l(i) = (1.d0+t1v(i))*(q2(i-1)*velx(i-1)+pres(i-1)) 
      f3l(i) = (1.d0+t1v(i))*q3(i-1)*velx(i-1)
      f4l(i) = (1.d0+t1v(i))*q4(i-1)*velx(i-1)
      f5l(i) = (1.d0+t1v(i))*(q5(i-1)+pres(i-1))*velx(i-1)
      f1r(i) = (1.d0-t1v(i))*q1(i)*velx(i)
      f2r(i) = (1.d0-t1v(i))*(q2(i)*velx(i)+pres(i)) 
      f3r(i) = (1.d0-t1v(i))*q3(i)*velx(i)
      f4r(i) = (1.d0-t1v(i))*q4(i)*velx(i)
      f5r(i) = (1.d0-t1v(i))*(q5(i)+pres(i))*velx(i)
   endfor
   f1l(0) = (1.d0+t1v(0))*q1(0)*velx(0)
   f2l(0) = (1.d0+t1v(0))*(q2(0)*velx(0)+pres(0)) 
   f3l(0) = (1.d0+t1v(0))*q3(0)*velx(0)
   f4l(0) = (1.d0+t1v(0))*q4(0)*velx(0)
   f5l(0) = (1.d0+t1v(0))*(q5(0)+pres(0))*velx(0)
   f1r(0) = (1.d0-t1v(0))*q1(0)*velx(0)
   f2r(0) = (1.d0-t1v(0))*(q2(0)*velx(0)+pres(0)) 
   f3r(0) = (1.d0-t1v(0))*q3(0)*velx(0)
   f4r(0) = (1.d0-t1v(0))*q4(0)*velx(0)
   f5r(0) = (1.d0-t1v(0))*(q5(0)+pres(0))*velx_a(0)
   f1l(n) = (1.d0+t1v(n))*q1(n-1)*velx(n-1)
   f2l(n) = (1.d0+t1v(n))*(q2(n-1)*velx(n-1)+pres(n-1)) 
   f3l(n) = (1.d0+t1v(n))*q3(n-1)*velx(n-1)
   f4l(n) = (1.d0+t1v(n))*q4(n-1)*velx(n-1)
   f5l(n) = (1.d0+t1v(n))*(q5(n-1)+pres(n-1))*velx(n-1)
   f1r(n) = (1.d0-t1v(n))*q1(n-1)*velx(n-1)
   f2r(n) = (1.d0-t1v(n))*(q2(n-1)*velx(n-1)+pres(n-1)) 
   f3r(n) = (1.d0-t1v(n))*q3(n-1)*velx(n-1)
   f4r(n) = (1.d0-t1v(n))*q4(n-1)*velx_a(n-1)
   f5r(n) = (1.d0-t1v(n))*(q5(n-1)+pres(n-1))*velx_a(n-1)

   ;flux limiter
   for i=2,n-2 do begin
      if(abs(a1(i)) gt 0.d0) then begin
         if(l1(i) gt 0.d0) then begin
            r1(i) = a1(i-1)/a1(i) 
         endif else begin
            r1(i) = a1(i+1)/a1(i) 
         endelse
      endif else begin
         r1(i) = 0.d0
      endelse
      if(abs(a2(i)) gt 0.d0) then begin
         if(l2(i) gt 0.d0) then begin
            r2(i) = a2(i-1)/a2(i) 
         endif else begin
            r2(i) = a2(i+1)/a2(i) 
         endelse
      endif else begin
         r2(i) = 0.d0
      endelse
      if(abs(dq3(i)) gt 0.d0) then begin
         if(l3(i) gt 0.d0) then begin
            r3(i) = a3(i-1)/a3(i) 
         endif else begin
            r3(i) = a3(i+1)/a3(i) 
         endelse
      endif else begin
         r3(i) = 0.d0
      endelse
      if(abs(dq4(i)) gt 0.d0) then begin
         if(l4(i) gt 0.d0) then begin
            r4(i) = a4(i-1)/a4(i) 
         endif else begin
            r4(i) = a4(i-1)/a4(i) 
         endelse
      endif else begin
         r4(i) = 0.d0
      endelse
      if(abs(dq5(i)) gt 0.d0) then begin
         if(l5(i) gt 0.d0) then begin
            r5(i) = a5(i-1)/a5(i) 
         endif else begin
            r5(i) = a5(i+1)/a5(i) 
         endelse
      endif else begin
         r5(i) = 0.d0
      endelse
   endfor
      r1(0) = r1(2)
      r1(1) = r1(2)
      r1(n) = r1(n-2)
      r1(n-1) = r1(n-2)
      r2(0) = r2(2)
      r2(1) = r2(2)
      r2(n) = r2(n-2)
      r2(n-1) = r2(n-2)
      r3(0) = r3(2)
      r3(1) = r3(2)
      r3(n) = r3(n-2)
      r3(n-1) = r3(n-2)
      r4(0) = r4(2)
      r4(1) = r4(2)
      r4(n) = r4(n-2)
      r4(n-1) = r4(n-2)
      r5(0) = r5(2)
      r5(1) = r5(2)
      r5(n) = r5(n-2)
      r5(n-1) = r5(n-2)

   ; flip flop function
   for i=0,n do begin
      t1(i) = sign(1.0,l1(i))
      t2(i) = sign(1.0,l2(i))
      t3(i) = sign(1.0,l3(i))
      t4(i) = sign(1.0,l4(i))
      t5(i) = sign(1.0,l5(i))
      p1(i) = phi(fl,r1(i))
      p2(i) = phi(fl,r2(i))
      p3(i) = phi(fl,r3(i))
      p4(i) = phi(fl,r4(i))
      p5(i) = phi(fl,r5(i))
   endfor 

   e1 = (l1) * dt/dx
   e2 = (l2) * dt/dx
   e3 = (l3) * dt/dx
   e4 = (l4) * dt/dx
   e5 = (l5) * dt/dx

   ;compute fluxes, FINALLY!
   f1(0) = 0.0
   f2(0) = 0.0 
   f3(0) = 0.0 
   f4(0) = 0.0 
   f5(0) = 0.0 
   f1(n) = 0.0
   f2(n) = 0.0 
   f3(n) = 0.0 
   f4(n) = 0.0 
   f5(n) = 0.0 


   for k=0,4 do begin
      diss(*,k) = (a1*(l1)*K1(*,k)*(t1+p1*(e1-t1)) + a2*(l2)*K2(*,k)*(t2+p2*(e2-t2))+ $
                   a3*(l3)*K3(*,k)*(t3+p3*(e3-t3)) + a4*(l4)*K4(*,k)*(t4+p4*(e4-t4))+ $
                   a5*(l5)*K5(*,k)*(t5+p5*(e5-t5)))
   endfor
   ;
   f1 = 0.5*(f1l+f1r) - 0.5 * diss(*,0) 
   f2 = 0.5*(f2l+f2r) - 0.5 * diss(*,1) 
   f3 = 0.5*(f3l+f3r) - 0.5 * diss(*,2) 
   f4 = 0.5*(f4l+f4r) - 0.5 * diss(*,3) 
   f5 = 0.5*(f5l+f5r) - 0.5 * diss(*,4)

   ;save old timestep quantities
   q1_old = q1
   q2_old = q2
   q3_old = q3
   q4_old = q4
   q5_old = q5

   ;advect
   for i=0,n-1 do begin
      
       q1(i) = q1_old(i) + dt/dx * (f1(i)-f1(i+1))
       q2(i) = q2_old(i) + dt/dx * (f2(i)-f2(i+1))
       q3(i) = q3_old(i) + dt/dx * (f3(i)-f3(i+1))
       q4(i) = q4_old(i) + dt/dx * (f4(i)-f4(i+1))
       q5(i) = q5_old(i) + dt/dx * (f5(i)-f5(i+1))
      
   endfor
    
   ;hydro quantities
   dens = q1
   velx = q2/q1
   vely = q3/q1
   velz = q4/q1
   ener = q5/q1

   ;eos
   eos,dens,pres,ener,eint,enth,velx,vely,velz,gamma

   time = time+dt
   if(time gt tmax) then break
 
   if (NOT ps and NOT wf) then plot,  x, dens, color=0, psym=-6;,yrange=[0.0,dens_left*1.1]

endfor


; plot to psfile

if(ps) then begin

   xsize = 7
   ysize = 5
   
   vars = [[dens],[velx],[pres],[enth]]
   var_names = ['dens','velx','pres','enth']
   fnum  = string(tmax,format='(i4.4)')
   
   for i=0,3 do begin
     psfile = var_names(i) + '_' + fl + '_' + fnum + '.ps' 
     PS_Start, psfile, xsize=xsize, ysize=ysize, /encapsulated 
     xtitle = 'x'
     ytitle = var_names(i)
     title = ytitle
     print, 'plotting ',ytitle
     cgplot, x, vars(*,i), xst=1, xtit=xtitle, title=title  $
          , yst=1, ytit=ytitle, yrange=[0.0,1.1*max(vars(*,i))], psym=-6, symsize=0.6, $
          charsize=1.0
     PS_End
   endfor

endif

if(wf) then begin
  data_file = 'roe_fl.out'
  openw,1,data_file
  for i=0,n-1 do begin
    printf,1,x(i),dens(i),velx(i),pres(i),eint(i)
  endfor
  close,1
  print,'written data to ',data_file
endif

print, 'SUCCESS!'

END

