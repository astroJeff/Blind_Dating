PRO CONTOUR_PLUS, x,y,XBINSIZE=xbinsize,YBINSIZE=ybinsize,XMIN=xmin,XMAX=xmax,YMIN=ymin,YMAX=ymax,PMINX=PMINX,PMAXX=pmaxx,PMINY=PMINY, PMAXY=PMAXY,LEVELS=levels,THRESHOLD=threshold,REVERSE=reverse,NLEVELS=nlevels,_EXTRA=extra,NOCONT=nocont,RELATIVE=RELATIVE, C_COLOR=c_color,INTERP=interp,SYMBOL=symbol,NODATA=nodata,NOCOLOR=nocolor

; PROGRAM: CONTOUR_PLUS
;
; WRITTEN: Anil Seth, 2005; adjusted by MA, 2012
;
; PURPOSE: Created for use in color-magnitude diagrams, contour_plus is
; a program which makes a contour plot above a certain data density
; threshold and plot individual data points below that threshold.
;
; ARGUMENTS:
;           X --> data values for the xaxis (e.g. V-I)
;           Y --> data values for the yaxis (e.g. V)
;
; OPTIONAL KEYWORDS:
;           XBINSIZE --> size of bins in X direction
;           YBINSIZE --> size of bins in Y direction
;           XMIN,XMAX,YMIN,YMAX --> boundaries of the plot/binning
;           PMINX,PMINY,PMAXX,PMAXY --> boundaries of plot, if they're different from binning
;           LEVELS --> density levels (points/bin) at which to draw contours
;             e.g. levels=[75,100,150,200,250,350,500,750,1000,1500,2000]
;           THRESHOLD --> density level above which to not draw points,
;              should be higher than lowest contour
;           REVERSE --> reverse the yaxis when plotting
;           NOCONT --> don't plot contours
;           _EXTRA --> applied to both PLOT and CONTOUR commands, can
;             include keywords such as C_COLOR, XMARGIN, etc.
;           INTERP --> interpolate the data to a grid 5 times smaller
;           SYMBOl --> plot objects outside of contours as blue stars
;           NOCOLOR --> plot objects outside of contours as black points
;           NODATA --> don't plot objects outside of contours 
;   MODIFICATIONS: Uses USERSYM rather than PSYM to avoid PDFing problems

IF NOT (KEYWORD_SET(xbinsize)) THEN xbinsize=.1
IF NOT (KEYWORD_SET(ybinsize)) THEN ybinsize=.1
;IF NOT (KEYWORD_SET(xmin)) THEN xmin=MIN(x)
;IF NOT (KEYWORD_SET(ymin)) THEN ymin=MIN(y)
;IF NOT (KEYWORD_SET(xmax)) THEN xmax=MAX(x)
;IF NOT (KEYWORD_SET(ymax)) THEN ymax=MAX(y)
IF NOT (KEYWORD_SET(PMINX)) THEN PMINX=xmin
IF NOT (KEYWORD_SET(PMINY)) THEN PMINY=ymin
IF NOT (KEYWORD_SET(PMAXX)) THEN PMAXX=xmax
IF NOT (KEYWORD_SET(PMAXY)) THEN PMAXY=ymax
IF NOT (KEYWORD_SET(threshold)) THEN threshold=75.
;IF NOT (KEYWORD_SET(levels)) THEN levels=[50,75,100,150,200,250,350,500,750]
IF NOT (KEYWORD_SET(c_color)) THEN c_color=[180,160,140,120,100,80,60,40,20]

IF (KEYWORD_SET(REVERSE)) THEN yrange=[PMAXY,PMINY] ELSE yrange=[PMINY,PMAXY]

nxbins=FLOOR((xmax-xmin) / xbinsize) + 1L
nybins=FLOOR((ymax-ymin) / ybinsize) + 1L
xbins=(FINDGEN(nxbins)*xbinsize)+xmin+xbinsize/2.
ybins=(FINDGEN(nybins)*ybinsize)+ymin+ybinsize/2.


hist=HIST_2D(x,y,bin1=xbinsize,bin2=ybinsize,min1=xmin,min2=ymin,max1=xmax,max2=ymax)

; To fix some interpolation issues
;FOR i=0,nxbins-1 DO BEGIN
;   FOR j=0,nybins-1 DO BEGIN
;      IF hist[i,j] LT 6 THEN hist[i,j]=0
;   ENDFOR
;ENDFOR

; To set limits to sigma contours
IF NOT (KEYWORD_SET(levels)) THEN BEGIN
   levels = FINDGEN(3)
   limits = FINDGEN(3)
   sigma = [0.003,0.0455,0.317]
   hist1d = reform(hist,N_ELEMENTS(hist))
   hist1d = hist1d[sort(hist1d)]
   total_val = 0
   FOR i=0,N_ELEMENTS(hist1d)-1,1 DO BEGIN
      total_val = total_val + hist1d[i]
   ENDFOR
   limits = sigma*total_val
   cumu_val = 0
   FOR i=0,N_ELEMENTS(hist1d)-1,1 DO BEGIN
      FOR j=0,2,1 DO BEGIN
         IF ((cumu_val LT limits[j]) AND (cumu_val+hist1d[i] GT limits[j])) THEN BEGIN
            levels[j] = hist1d[i]+1
         ENDIF
      ENDFOR
      
      cumu_val = cumu_val + hist1d[i]
   ENDFOR
   print,levels
ENDIF

IF (KEYWORD_SET(INTERP)) THEN BEGIN
   dims=size(hist, /dim)        ; get dimensions of 2d histogram
                                ;interpolate image to a grid 5 times smaller
   interpimage=min_curve_surf(hist, /regular, nx=(dims[0]-1)*5+1, ny=(dims[1]-1)*5+1)

   dims2=size(interpimage, /dim) ; get dimensions of the 2d histogram
 
;   FOR i=0,dims2[0]-1 DO BEGIN
;      FOR j=0,dims2[1]-1 DO BEGIN
;         IF interpimage[i,j] LT 6 THEN interpimage[i,j] = 0.
;      ENDFOR
;   ENDFOR

ENDIF

PLOTSYM, 0, .1, /FILL

IF (KEYWORD_SET(levels) AND KEYWORD_SET(relative)) THEN levels=levels*MAX(hist)
;IF NOT (KEYWORD_SET(levels)) THEN levels=[0.001,0.003,0.007,0.01,0.03,0.07,0.1,0.3,0.7]*MAX(hist)
IF NOT (KEYWORD_SET(threshold)) THEN threshold=levels[2]
IF (threshold LT MIN(hist)) THEN threshold=MIN(hist)*1.5

IF (MAX(hist) LT threshold*1.5) THEN BEGIN
   oplot,x,y,psym=8
;   plot,x,y,xrange=[PMINX,PMAXX],yrange=yrange,psym=8,_EXTRA=extra,xstyle=4,ystyle=4
;   COMMENTED TO GET USERSYM RATHER THAN PSYM=3
;   plot,x,y,xrange=[PMINX,PMAXX],yrange=yrange,psym=3,_EXTRA=extra
ENDIF ELSE BEGIN

   IF (KEYWORD_SET(SYMBOL)) THEN BEGIN
      PLOTSYM, 3, /FILL
      LOADCT, 13, /SILENT
   ENDIF

   IF (KEYWORD_SET(NODATA)) THEN BEGIN
;      plot,x[ind],y[ind],xrange=[PMINX,PMAXX],yrange=yrange,psym=8,_EXTRA=extra,xstyle=5,ystyle=5,symsize=1.1, /NODATA
   ENDIF ELSE BEGIN

;now create array for indices of stars that fall outside of contours
      ind=[0l]
      FOR i=0l,N_ELEMENTS(x)-1 DO BEGIN
         mindistx=MIN(ABS(xbins-x[i]),PMINXos)
         mindisty=MIN(ABS(ybins-y[i]),PMINYos)
         IF (hist[PMINXos,PMINYos] LT threshold) THEN ind=[ind,i]
      ENDFOR
      
      ind=ind[1:*]
      

      IF (KEYWORD_SET(NOCOLOR)) THEN BEGIN
         oplot,x[ind],y[ind],psym=8,symsize=1.1
;      plot,x[ind],y[ind],xrange=[PMINX,PMAXX],yrange=yrange,psym=8,_EXTRA=extra,xstyle=5,ystyle=5,symsize=1.1 
      ENDIF ELSE BEGIN
         LOADCT, 13, /SILENT
         oplot,x[ind],y[ind],psym=8,symsize=1.1,color=50
;;      plot,x[ind],y[ind],xrange=[PMINX,PMAXX],yrange=yrange,psym=8,_EXTRA=extra,xstyle=5,ystyle=5,symsize=1.1, color = 50
      ENDELSE
;   COMMENTED TO GET USERSYM RATHER THAN PSYM=3
;   plot,x,y,xrange=[PMINX,PMAXX],yrange=yrange,psym=3,_EXTRA=extra
      
   ENDELSE

   IF NOT (KEYWORD_SET(nlevels)) THEN nlevels=MAX(hist)/threshold
   IF NOT (KEYWORD_SET(levels)) THEN levels=threshold*FINDGEN(nlevels)+threshold
   
   LOADCT, 0, /SILENT
   IF (KEYWORD_SET(INTERP)) THEN BEGIN
      dims2=size(interpimage, /dim) ; get dimensions of the 2d histogram
      
                                ; find the corresponding x and y coordinates for each bin in 2d histogram... 
                                ; subtract off half a bin to set coordinates to the *center of the bin*
      
      xs=findgen(dims2[0])*xbinsize/5.+(xmin+xbinsize/2.)
      ys=findgen(dims2[1])*ybinsize/5.+(ymin+ybinsize/2.)
      contour,interpimage,xs,ys,yrange=yrange,levels=levels,/overplot,/FILL,_EXTRA=extra,C_COLOR=c_color 
;,C_COLOR=REPLICATE(!P.BACKGROUND,nlevels),xstyle=4,ystyle=4
 
      IF (NOT KEYWORD_SET(nocont)) THEN BEGIN
         contour,interpimage,xs,ys,yrange=yrange,levels=levels,/overplot,color=!P.BACKGROUND
      ENDIF ELSE BEGIN
         contour,hist,xbins,ybins,yrange=yrange,levels=levels,/overplot,/FILL,_EXTRA=extra,C_COLOR=c_color ;    
         contour,hist,xbins,ybins,yrange=yrange,levels=levels,/overplot,color=!P.BACKGROUND
         IF (NOT KEYWORD_SET(nocont)) THEN BEGIN
            print,"HERE"
;            contour,hist,xbins,ybins,yrange=yrange,levels=levels,/overplot,color=!P.BACKGROUND
         ENDIF
            
      ENDELSE
         
   ENDIF ELSE BEGIN

      contour,hist,xbins,ybins,yrange=yrange,levels=levels,/overplot,/FILL,_EXTRA=extra,C_COLOR=c_color ;    
      contour,hist,xbins,ybins,yrange=yrange,levels=levels,/overplot,color=!P.BACKGROUND

   ENDELSE


ENDELSE

END
