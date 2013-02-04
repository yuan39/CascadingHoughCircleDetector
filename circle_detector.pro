;+
; :Author: Yuan Yuan
; :History: 20100713
; :Description: Find out the radius and center location of a solar disk in a halpha image
;-
pro circle_detector,filename,minr,maxr,cenx,ceny,radius,timespan
starttime = systime(1)
img = readfits(filename); readfits('c:\bbso_halph_fr_20021228_192204.fts');
WRITE_GIF, filename+'_original.gif', bytscl(img)
sz = SIZE(img, /DIMENSIONS)
if maxr EQ 0 then begin
  if( sz[0] LT sz[1]) then begin
    minr = sz[0]/4
    maxr = sz[0]/2
  endif else begin
    minr = sz[1]/4
    maxr = sz[1]/2
  endelse
endif
;img = BytScl(img, MIN=mean(img)-3*stdev(img), MAX=mean(img)+3*stdev(img))
mean_img = mean(img);
stdev_img = stdev(img);
tra_ind = where(img GT mean_img+3*stdev_img,count)
if(count GT 0) then begin
  img(tra_ind) = mean_img+3*stdev_img;
endif
inf_ind = where(img LT mean_img-3*stdev_img,count)
if(count GT 0) then begin
  img(inf_ind) = mean_img-3*stdev_img;
endif
;img = congrid(img,800,800)
shrink_num = 5;
S = SIZE(img)
NX = S[1] & NY = S[2]
shrinked = CONGRID(img, fix(NX/shrink_num), fix(NY/shrink_num));
S = SIZE(shrinked)
NX = S[1] & NY = S[2]
WINDOW, 1, XSIZE = NX, YSIZE = NY
TVSCL, shrinked
smoothed_shrinked = MEDIAN(shrinked, 20)

WINDOW, 2, XSIZE = NX, YSIZE = NY
TVSCL, smoothed_shrinked

shrinked_edge = SOBEL(smoothed_shrinked)

;WINDOW, 3, XSIZE = NX, YSIZE = NY
;TVSCL, shrinked_edge

filter_x = [[1,2,1], [0,0,0],[-1,-2,-1]]/8.0
;print,filter_x
filter_y = [[1,0,-1],[2,0,-2],[1,0,-1]]/8.0
;print,filter_y
filtered_x =  convol(smoothed_shrinked,filter_x,/edge_zero);
filtered_y =  convol(smoothed_shrinked,filter_y,/edge_zero);
edge_map = sqrt(filtered_x*filtered_x + filtered_y*filtered_y);
WINDOW, 3, XSIZE = NX, YSIZE = NY
TVSCL, edge_map
;thresh = sqrt(1000*mean(shrinked_edge));
;thresh = sqrt(1000*mean(edge_map));
thresh = 3*stdev(edge_map);max(edge_map)/2;
bin_edge = edge_map GE thresh
WINDOW, 4, XSIZE = NX, YSIZE = NY
TVSCL, bin_edge

thin_edge = bytarr(NX,NY)
;;;;;;;;;;;;;;;;;thinning the edge map
for i = 0,NX-1 do begin
  for j = 0,NY-1 do begin
    if( (i-1) LT 0 ) then b1 = 1 else b1 = edge_map(i-1,j) LE edge_map(i,j)
    if( (i+1) GT NX-1 ) then b2 = 1 else b2 = edge_map(i+1,j) LT edge_map(i,j)
    if( (j-1) LT 0 ) then b3 = 1 else b3 = edge_map(i,j-1) LE edge_map(i,j)
    if( (j+1) GT  NY-1 ) then b4 = 1 else b4 = edge_map(i,j+1) LT edge_map(i,j)
    thin_edge(i,j) =  bin_edge(i,j) and ( ( abs(filtered_x(i,j)) GE abs(filtered_y(i,j)) and  b1 and b2 ) or  ( abs(filtered_x(i,j)) LE abs(filtered_y(i,j)) and b3 and b4 ))

    ;thin_edge(i,j) = ( (edge_map(i,j) GT thresh) and abs(filtered_x(i,j)) GE abs(filtered_y(i,j)) and  b1 and b2 )
    ; or ( (edge_map(i,j)GT thresh) and abs(filtered_x(i,j)) LE abs(filtered_y(i,j)) and b3 and b4 )
  endfor
endfor
thin_edge(0:10,*) = 0;
thin_edge(*,0:10) = 0;
thin_edge(NX-11:NX-1,*) = 0;
thin_edge(*,NY-11:NY-1) = 0;
WINDOW, 5, XSIZE = NX, YSIZE = NY
TVSCL, thin_edge

tmpmaxval = 0L;
tmpHM = 0L;
tmpx_cen =0
tmpy_cen =0
tmpr = 0;
for r = round(minr/shrink_num),round(maxr/shrink_num) do begin
  yyHough,thin_edge,r,maxval,HM,x_cen,y_cen
  if maxval GT tmpmaxval then begin
    tmpmaxval = maxval
    tmpHM = HM
    tmpr = r
    tmpx_cen = x_cen
    tmpy_cen = y_cen
  endif
endfor
print,'the radius,  x-center location,  y-center location  is: '
print,tmpr,tmpx_cen,tmpy_cen
;FOR k=0,n_elements(tp)-1 DO Begin
;  x =
;  OPLOT,circles[k].cx+circles[k].R*COS(tp),circles[k].cy+circles[k].R*SIN(tp)
;endfor





;r2 = r^2;
;for i = 0,NX-1 do begin
;  for j = 0,fix(NY/shrink_num)-1 do begin
;    if thin_edge(i,j) GT 0 then begin
;      for m = 0,NX-1 do begin
;        for n = 0,fix(NY/shrink_num)-1 do begin
;          if abs((m-i)^2 + (n-j)^2 - r2) LT 4 then begin
;            HM(m,n) = HM(m,n) + 1;
;          endif
;        endfor
;      endfor
;    endif
;  endfor
;endfor

WINDOW, 6, XSIZE = NX, YSIZE = NY
TVSCL, tmpHM


;Plot the results
WINDOW,7, XSIZE=NX,YSIZE=NY,TITLE='Circles Found';,XPOS=305,YPOS=0
;tp=FINDGEN(201)*!PI/50
;PLOT,X,Y,psym=4,XRANGE=[0,5],YRANGE=[0,5],XTICKLEN=1,YTICKLEN=1,$
;  xmargin=[2,2],ymargin=[2,2]
tvscl,shrinked
PLOTS, CIRCLE(tmpx_cen, tmpy_cen, tmpr), /Device
;FOR k=0,N_ELEMENTS(circles)-1 DO $
;xcen=tp
;xcen[*]=tmpx_cen
;ycen=tp
;ycen[*]=tmpy_cen
;  OPLOT,xcen+tmpr*COS(tp),ycen+tmpr*SIN(tp)



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;search in original domain
S = SIZE(img)
NX = S[1] & NY = S[2]
smoothed = MEDIAN(img, 20)


window,8,xs=400,ys=400
tvscl,congrid(smoothed,400,400)

;ori_edge = SOBEL(smoothed)

;WINDOW, 3, XSIZE = NX, YSIZE = NY
;TVSCL, ori_edge

filter_x = [[1,2,1], [0,0,0],[-1,-2,-1]]/8.0
;print,filter_x   
filter_y = [[1,0,-1],[2,0,-2],[1,0,-1]]/8.0
;print,filter_y
filtered_x =  convol(smoothed,filter_x,/edge_zero);
filtered_y =  convol(smoothed,filter_y,/edge_zero);
ori_edge_map = sqrt(filtered_x*filtered_x + filtered_y*filtered_y);

window,9, xs=400,ys=400
tvscl,congrid(ori_edge_map,400,400)
;thresh = sqrt(1000*mean(shrinked_edge));
thresh = 3*stdev(ori_edge_map);sqrt(10000*mean(ori_edge_map));
ori_bin_edge = ori_edge_map GE thresh

window,10,xs=400,ys=400
tvscl,congrid(ori_bin_edge,400,400)

ori_thin_edge = bytarr(NX,NY)
;;;;;;;;;;;;;;;;;thinning the edge map
for i = 0,NX-1 do begin
  for j = 0,NY-1 do begin
    if( (i-1) LT 0 ) then b1 = 1 else b1 = ori_edge_map(i-1,j) LE ori_edge_map(i,j)
    if( (i+1) GT NX-1 ) then b2 = 1 else b2 = ori_edge_map(i+1,j) LT ori_edge_map(i,j)
    if( (j-1) LT 0 ) then b3 = 1 else b3 = ori_edge_map(i,j-1) LE ori_edge_map(i,j)
    if( (j+1) GT  NY-1 ) then b4 = 1 else b4 = ori_edge_map(i,j+1) LT ori_edge_map(i,j)
    ori_thin_edge(i,j) =  ori_bin_edge(i,j) and ( ( abs(filtered_x(i,j)) GE abs(filtered_y(i,j)) and  b1 and b2 ) or  ( abs(filtered_x(i,j)) LE abs(filtered_y(i,j)) and b3 and b4 ))
    
    ;thin_edge(i,j) = ( (edge_map(i,j) GT thresh) and abs(filtered_x(i,j)) GE abs(filtered_y(i,j)) and  b1 and b2 )
    ; or ( (edge_map(i,j)GT thresh) and abs(filtered_x(i,j)) LE abs(filtered_y(i,j)) and b3 and b4 )
  endfor
endfor
ori_thin_edge(0:40,*) = 0;
ori_thin_edge(*,0:40) = 0;
ori_thin_edge(NX-41:NX-1,*) = 0;
ori_thin_edge(*,NY-41:NY-1) = 0;
WRITE_PNG, filename+'_edge.png', bytscl(ori_thin_edge)
window,11,xs=400,ys=400
tvscl,congrid(ori_thin_edge,400,400)


o_tmpmaxval = 0;
o_tmpHM = 0;
o_tmpx_cen =0
o_tmpy_cen =0
o_tmpr = 0;

for ori_r = tmpr*shrink_num-1*shrink_num,tmpr*shrink_num+1*shrink_num do begin
  yyHough,ori_thin_edge,ori_r,maxval,HM,x_cen,y_cen
  if maxval GT o_tmpmaxval then begin
    o_tmpmaxval = maxval
    o_tmpHM = HM
    o_tmpr = ori_r
    o_tmpx_cen = x_cen
    o_tmpy_cen = y_cen
  endif
endfor
print,'the radius, x-center location, y-center location  is: '
print,o_tmpr,o_tmpx_cen,o_tmpy_cen
cenx = o_tmpx_cen
ceny = o_tmpy_cen
radius = o_tmpr
newimg = img(o_tmpx_cen-o_tmpr:o_tmpx_cen+o_tmpr,o_tmpy_cen-o_tmpr:o_tmpy_cen+o_tmpr)


window,12,xs=400,ys=400
tvscl,congrid(newimg,400,400)
WRITE_JPEG, filename+'_crop.jpg', bytscl(newimg)
writefits, filename+'_crop.fits', newimg

window,13,xs=400,ys=400
tvscl,congrid(o_tmpHM,400,400)
timespan = systime(1)- starttime
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro yyHough,thin_edge,r,maxval,HM,x_cen,y_cen
S = SIZE(thin_edge)
NX = S[1] & NY = S[2]
;r = 185L;
HM = lonarr(NX,NY)
;tp=FINDGEN(101)*!PI/50
ind = where(thin_edge EQ 1,cnt)


for c = 0L,cnt-1L do begin
  cx = ind(c) mod NX
  cy = ind(c) / NX
  x = 0L & y = r & F = 3L-2*r
  while( x LT y) do begin
    ;cirpot,cx,cy,x,y
    ;if (cx+x GE 0 and cx+x LT NX and $
    ;    cy+y GE 0 and cy+y LT fix(NY/shrink_num) )
    ;    HM(cx+x,cy+y) = HM(cx+x,cy+y)+1;
    ;endif
    inc,HM,cx,cy,x,y,NX,NY
    inc,HM,cx,cy,y,x,NX,NY
    inc,HM,cx,cy,y,-x,NX,NY
    inc,HM,cx,cy,x,-y,NX,NY
    inc,HM,cx,cy,-x,-y,NX,NY
    inc,HM,cx,cy,-y,-x,NX,NY
    inc,HM,cx,cy,-y,x,NX,NY
    inc,HM,cx,cy,-x,y,NX,NY
    if(F LT 0) then  begin
        F = F + 4*x + 6
    endif else begin
        F = F + 4 * (x-y) + 10;
        y = y - 1
    endelse
    x = x + 1
  endwhile
  if x EQ y then begin
  ;cirpot(cx,cy,x,y);
    inc,HM,cx,cy,x,y,NX,NY
    inc,HM,cx,cy,y,x,NX,NY
    inc,HM,cx,cy,y,-x,NX,NY
    inc,HM,cx,cy,x,-y,NX,NY
    inc,HM,cx,cy,-x,-y,NX,NY
    inc,HM,cx,cy,-y,-x,NX,NY
    inc,HM,cx,cy,-y,x,NX,NY
    inc,HM,cx,cy,-x,y,NX,NY
   endif
endfor
maxval = max(HM)
print,'the number of limb pixel, current investigating radius, and the maximum response value is: '
print,cnt,r,maxval
maxind = round(mean(where(HM EQ maxval)))
x_cen = maxind mod NX
y_cen = maxind / NX

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro inc,HM,cx,cy,x,y,xsize,ysize
if (cx+x GE 0 and cx+x LT xsize and  cy+y GE 0 and cy+y LT ysize )then begin
        HM(cx+x,cy+y) = HM(cx+x,cy+y)+1;
endif
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
FUNCTION CIRCLE, xcenter, ycenter, radius
   points = (2 * !PI / 99.0) * FINDGEN(100)
   x = xcenter + radius * COS(points )
   y = ycenter + radius * SIN(points )
   RETURN, TRANSPOSE([[x],[y]])
END
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;