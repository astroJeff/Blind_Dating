 ; Output to postscript file or to screen

   output = 'screen'
   output = 'PS'

; Input and output files
;   infile = COMMAND_LINE_ARGS()

   dir = './'

   fin1 = dir + 'dist_07.dat'
   fin2 = dir + 'dist_07_NS.dat'
   fin3 = dir + 'dist_flat_NS.dat'

   fin4 = dir + 'posterior_07.dat'
   fin5 = dir + 'posterior_07_NS.dat'
   fin6 = dir + 'posterior_flat_NS.dat'

   fin7 = dir + 'dist_07_ana.dat'
   fin8 = dir + 'dist_07_NS_ana.dat'
   fin9 = dir + 'dist_flat_NS_ana.dat'

   fin10 = dir + 'hist_P_NS_ind_07.dat'
   fin11 = dir + 'hist_P_NS_ind_07_NS.dat'
   fin12 = dir + 'hist_P_NS_ind_flat_NS.dat'

   fin13 = dir + 'hist_trueNS_07_NS.dat'
   fin14 = dir + 'hist_trueNS_flat_NS.dat'

   fout = './model_post.eps'

; Default line and character thickness

 ; !P.THICK = 6.0
 ; !P.CHARTHICK = 6.0
 ; !X.THICK = 6.0
 ; !Y.THICK = 6.0

; Default character size

  ;!P.CHARSIZE = 2.0
  !P.FONT = 0
  !P.CHARTHICK = 2.0
  !X.THICK = 2.0
  !Y.THICK = 2.0
  !P.CHARSIZE = 1.1

; Define filled circle plotting symbol

  A = FINDGEN(17)*(!PI*2/16.)
  USERSYM, 1.0*COS(A), 1.0*SIN(A), /FILL

; Plot positions

  pos1 = [0.05,0.70,0.24,0.93]
  pos2 = [0.30,0.70,0.49,0.93]
  pos3 = [0.55,0.70,0.74,0.93]
  pos4 = [0.79,0.70,0.98,0.93]

  pos5 = [0.05,0.38,0.24,0.65]
  pos6 = [0.30,0.38,0.49,0.65]
  pos7 = [0.55,0.38,0.74,0.65]
  pos8 = [0.79,0.38,0.98,0.65]

  pos9 = [0.05,0.05,0.24,0.33]
  pos10 = [0.30,0.05,0.49,0.33]
  pos11 = [0.55,0.05,0.74,0.33]
  pos12 = [0.79,0.05,0.98,0.33]

; Setup color table

  red   = [0,255,200,162,108,54,255,255]
  green = [0,0,54,108,162,200,255,255]
  blue  = [0,0,54,108,162,200,255,255]

; Load color table

  TVLCT, [red], [green], [blue]

; Plot axes and ranges and labels


;   xlab = 'M!DHe!N (M' + SUNSYMBOL() +')'
;   ylab = 'Log (P!Dpre-SN!N/days)'


; Open ps-file for output

  IF (output EQ 'PS') THEN BEGIN
     SET_PLOT, 'PS'
     DEVICE, /LANDSCAPE,ENCAPSULATE=0, XSIZE=16,YSIZE=12,FILE=fout, SCALE=1.0, $
             /INCHES,/SCHOOLBOOK,/COLOR
  ENDIF ELSE BEGIN
     DEVICE, DECOMPOSED=0
     WINDOW, 1, XSIZE=600, YSIZE=600, RETAIN=2
  ENDELSE

; Set background color (polyfill is necessary for postscript output)

  POLYFILL, [1,1,0,0,1], [1,0,0,1,1], /NORMAL, COLOR=7

; Make plot



  READCOL, fin1, a1,a2,a3,a4,a5, FORMAT='F,F,F,F,F'
  READCOL, fin2, b1,b2,b3,b4,b5, FORMAT='F,F,F,F,F'
  READCOL, fin3, c1,c2,c3,c4,c5, FORMAT='F,F,F,F,F'


  READCOL, fin4, d1,d2,d3,d4,d5,d6,d7, FORMAT='F,F,F,F,F,F,F'
  READCOL, fin5, e1,e2,e3,e4,e5,e6,e7, FORMAT='F,F,F,F,F,F,F'
  READCOL, fin6, f1,f2,f3,f4,f5,f6,f7, FORMAT='F,F,F,F,F,F,F'

  
  READCOL, fin7, g1,g2, FORMAT='F,F'
  READCOL, fin8, h1,h2, FORMAT='F,F'
  READCOL, fin9, i1,i2, FORMAT='F,F'

  
  READCOL, fin10, j1,j2, FORMAT='F,F'
  READCOL, fin11, k1,k2, FORMAT='F,F'
  READCOL, fin12, l1,l2, FORMAT='F,F'
  

  READCOL, fin13, m1,m2, FORMAT='F,F'
  READCOL, fin14, n1,n2, FORMAT='F,F'


  xmin = 0.0
  xmax = 1.50
  ymin = 0.00
  ymax = 20.0

  xlab = 'Mass (M' + SUNSYMBOL() + ')'
  ylab = 'N'

  PLOT, [0], [0], XRANGE=[xmin,xmax], YRANGE=[ymin,ymax],     $
     XSTYLE=1, YSTYLE=1, XTITLE=xlab, YTITLE=ylab, $
     COLOR=0, POSITION=pos1, /NODATA, /NOERASE

  pdf_M2 = HISTOGRAM(a4,MIN=xmin,MAX=xmax,BINSIZE=0.1,LOCATIONS=pdf_M2_x)
  pdf_M2_min = HISTOGRAM(a5,MIN=xmin,MAX=xmax,BINSIZE=0.1,LOCATIONS=pdf_M2_min_x)

;  OPLOT, PDF_M2_x, PDF_M2, PSYM=10, LINESTYLE=2, THICK=4
  OPLOT, PDF_M2_min_x, PDF_M2_min, PSYM=10, THICK=4

;  pdf_M2_min = lindgen( N_elements()) * 0.1 + min(arr) 


  OPLOT, g1,2*g2, LINESTYLE=0, COLOR=1






  xmin = 0.05
  xmax = 1.00
  ymin = 0.05
  ymax = 0.50

  xlab = 'Mu'
  ylab = 'Sigma'

  PLOT, [0], [0], XRANGE=[xmin,xmax], YRANGE=[ymin,ymax],     $
     XSTYLE=1, YSTYLE=1, XTITLE=xlab, YTITLE=ylab,  $
     COLOR=0, POSITION=pos2, /NODATA, /NOERASE
  
  CONTOUR_PLUS,d3,d4,XBINSIZE=0.03,YBINSIZE=0.03,XMIN=xmin,XMAX=xmax,YMIN=ymin,YMAX=ymax,PMINY=ymin,PMAXY=ymax,MINX=xmin,PMAXX=xmax,THRESHOLD=1,/INTERP,/NODATA,/NOERASE,/NOCONT

  mu_true = FINDGEN(1)
  sd_true = FINDGEN(1)
  mu_true[0] = 0.70
  sd_true[0] = 0.20
  OPLOT,mu_true,sd_true,PSYM=7,SYMSIZE=2.0
;  CONTOUR_PLUS,d1,d2,XBINSIZE=0.03,YBINSIZE=0.02,XMIN=xmin,XMAX=xmax,YMIN=ymin,YMAX=ymax,PMINY=ymin,PMAXY=ymax,MINX=xmin,PMAXX=xmax,LEVELS=levels,THRESHOLD=1,/INTERP,/NODATA,/NOERASE,/NOCONT


;  levels=[2.0,14.0,245.0,2645.0]
;  CONTOUR_PLUS,f1,f_Porb,XBINSIZE=0.1,YBINSIZE=0.1,XMIN=x1min,XMAX=x1max,YMIN=y1min,YMAX=y1max,PMINX=x1min,PMAXX=x1max,PMINY=y1min,PMAXY=y1max,LEVELS=levels,THRESHOLD=1,/NODATA,/NOERASE,/NOCONT




  PLOT, [0], [0], XRANGE=[xmin,xmax], YRANGE=[ymin,ymax],     $
     XSTYLE=1, YSTYLE=1,  $
     COLOR=0, POSITION=pos2, /NODATA, /NOERASE



; NS Probability Posterior

  xmin = 0.0
  xmax = 0.3
  ymin = 0.00
  ymax = 0.10

  xlab = 'NS Fraction'
  ylab = 'N'

  PLOT, [0], [0], XRANGE=[xmin,xmax], YRANGE=[ymin,ymax],     $
     XSTYLE=1, YSTYLE=1, XTITLE=xlab, YTITLE=ylab,  $
     COLOR=0, POSITION=pos3, /NODATA, /NOERASE

  pdf_NS = float(HISTOGRAM(d6,MIN=xmin,MAX=xmax,BINSIZE=0.005,LOCATIONS=pdf_NS_x))
  pdf_NS = pdf_NS/N_ELEMENTS(d1)

  OPLOT, PDF_NS_x, PDF_NS, PSYM=10



; Individual System NS Probability

  xmin = 0.0
  xmax = 1.0
  ymin = 0.00
  ymax = 80.0

  xlab = 'NS Probability'
  ylab = 'N'

  PLOT, [0], [0], XRANGE=[xmin,xmax], YRANGE=[ymin,ymax],     $
     XSTYLE=1, YSTYLE=1, XTITLE=xlab, YTITLE=ylab,  $
     COLOR=0, POSITION=pos4, /NODATA, /NOERASE


;  pdf_NS = float(HISTOGRAM(j1,MIN=xmin,MAX=xmax,BINSIZE=0.02,LOCATIONS=pdf_NS_x))

  OPLOT, j1, j2, PSYM=10





  TVLCT, [red], [green], [blue]

; Random sample M2, M2_min 

  xmin = 0.0
  xmax = 1.50
  ymin = 0.00
  ymax = 20.0

  xlab = 'Mass (M' + SUNSYMBOL() + ')'
  ylab = 'N'


  PLOT, [0], [0], XRANGE=[xmin,xmax], YRANGE=[ymin,ymax],     $
     XSTYLE=1, YSTYLE=1, XTITLE=xlab, YTITLE=ylab,  $
     COLOR=0, POSITION=pos5, /NODATA, /NOERASE

  pdf_M2 = HISTOGRAM(b4,MIN=xmin,MAX=xmax,BINSIZE=0.1,LOCATIONS=pdf_M2_x)
  pdf_M2_min = HISTOGRAM(b5,MIN=xmin,MAX=xmax,BINSIZE=0.1,LOCATIONS=pdf_M2_min_x)

;  OPLOT, PDF_M2_x, PDF_M2, PSYM=10, LINESTYLE=2, THICK=4
  OPLOT, PDF_M2_min_x, PDF_M2_min, PSYM=10, THICK=4

  OPLOT, h1,2*h2, LINESTYLE=0, COLOR=1





; Gaussian parameters posterior

  xmin = 0.05
  xmax = 1.00
  ymin = 0.05
  ymax = 0.50

  xlab = 'Mu'
  ylab = 'Sigma'


  PLOT, [0], [0], XRANGE=[xmin,xmax], YRANGE=[ymin,ymax],     $
     XSTYLE=1, YSTYLE=1, XTITLE=xlab, YTITLE=ylab,  $
     COLOR=0, POSITION=pos6, /NODATA, /NOERASE

  CONTOUR_PLUS,e3,e4,XBINSIZE=0.03,YBINSIZE=0.03,XMIN=xmin,XMAX=xmax,YMIN=ymin,YMAX=ymax,PMINY=ymin,PMAXY=ymax,MINX=xmin,PMAXX=xmax,THRESHOLD=1,/INTERP,/NODATA,/NOERASE,/NOCONT

  OPLOT,mu_true,sd_true,PSYM=7,SYMSIZE=2.0

  PLOT, [0], [0], XRANGE=[xmin,xmax], YRANGE=[ymin,ymax],     $
     XSTYLE=1, YSTYLE=1, XTITLE=xlab, YTITLE=ylab,  $
     COLOR=0, POSITION=pos6, /NODATA, /NOERASE


; NS Probability Posterior

  xmin = 0.0
  xmax = 0.3
  ymin = 0.00
  ymax = 0.10

  xlab = 'NS Fraction'
  ylab = 'N'

  PLOT, [0], [0], XRANGE=[xmin,xmax], YRANGE=[ymin,ymax],     $
     XSTYLE=1, YSTYLE=1, XTITLE=xlab, YTITLE=ylab,  $
     COLOR=0, POSITION=pos7, /NODATA, /NOERASE

  pdf_NS = float(HISTOGRAM(e6,MIN=xmin,MAX=xmax,BINSIZE=0.005,LOCATIONS=pdf_NS_x))
  OPLOT, PDF_NS_x, PDF_NS/N_ELEMENTS(e1), PSYM=10



; Individual System NS Probability

  xmin = 0.0
  xmax = 1.0
  ymin = 0.00
  ymax = 60.0

  xlab = 'NS Probability'
  ylab = 'N'

  PLOT, [0], [0], XRANGE=[xmin,xmax], YRANGE=[ymin,ymax],     $
     XSTYLE=1, YSTYLE=1, XTITLE=xlab, YTITLE=ylab,  $
     COLOR=0, POSITION=pos8, /NODATA, /NOERASE

;  pdf_NS = float(HISTOGRAM(k1,MIN=xmin,MAX=xmax,BINSIZE=0.02,LOCATIONS=pdf_NS_x))
;  pdf_true_NS = float(HISTOGRAM(m1,MIN=xmin,MAX=xmax,BINSIZE=0.02,LOCATIONS=pdf_true_NS_x))

  

  OPLOT, k1, k2, PSYM=10
  dx = (m1[2]-m1[1])/2.0
  FOR i=0, N_ELEMENTS(m2)-1, 1 DO BEGIN
     POLYFILL, [m1[i]-dx,m1[i]-dx,m1[i]+dx,m1[i]+dx,m1[i]-dx], [0,m2[i],m2[i],0,0], COLOR=2
  ENDFOR








  TVLCT, [red], [green], [blue]


; Random sample M2, M2_min 

  xmin = 0.0
  xmax = 1.50
  ymin = 0.00
  ymax = 20.0

  xlab = 'Mass (M' + SUNSYMBOL() + ')'
  ylab = 'N'

  PLOT, [0], [0], XRANGE=[xmin,xmax], YRANGE=[ymin,ymax],     $
     XSTYLE=1, YSTYLE=1, XTITLE=xlab, YTITLE=ylab,  $
     COLOR=0, POSITION=pos9, /NODATA, /NOERASE

  pdf_M2 = HISTOGRAM(c4,MIN=xmin,MAX=xmax,BINSIZE=0.1,LOCATIONS=pdf_M2_x)
  pdf_M2_min = HISTOGRAM(c5,MIN=xmin,MAX=xmax,BINSIZE=0.1,LOCATIONS=pdf_M2_min_x)

;  OPLOT, PDF_M2_x, PDF_M2, PSYM=10, LINESTYLE=2, THICK=4
  OPLOT, PDF_M2_min_x, PDF_M2_min, PSYM=10, THICK=4

  OPLOT, i1,2*i2, LINESTYLE=0, COLOR=1



; Gaussian parameters posterior

  xmin = 0.05
  xmax = 1.00
  ymin = 0.05
  ymax = 0.80

  xlab = 'Mu'
  ylab = 'Sigma'

  PLOT, [0], [0], XRANGE=[xmin,xmax], YRANGE=[ymin,ymax],     $
     XSTYLE=1, YSTYLE=1, XTITLE=xlab, YTITLE=ylab,  $
     COLOR=0, POSITION=pos10, /NODATA, /NOERASE

  CONTOUR_PLUS,f3,f4,XBINSIZE=0.03,YBINSIZE=0.03,XMIN=xmin,XMAX=xmax,YMIN=ymin,YMAX=ymax,PMINY=ymin,PMAXY=ymax,MINX=xmin,PMAXX=xmax,THRESHOLD=1,/INTERP,/NODATA,/NOERASE,/NOCONT

  PLOT, [0], [0], XRANGE=[xmin,xmax], YRANGE=[ymin,ymax],     $
     XSTYLE=1, YSTYLE=1, XTITLE=xlab, YTITLE=ylab,  $
     COLOR=0, POSITION=pos10, /NODATA, /NOERASE



; NS Probability Posterior

  xmin = 0.0
  xmax = 0.3
  ymin = 0.00
  ymax = 0.10

  xlab = 'NS Fraction'
  ylab = 'N'

  PLOT, [0], [0], XRANGE=[xmin,xmax], YRANGE=[ymin,ymax],     $
     XSTYLE=1, YSTYLE=1, XTITLE=xlab, YTITLE=ylab,  $
     COLOR=0, POSITION=pos11, /NODATA, /NOERASE


  pdf_NS = float(HISTOGRAM(f6,MIN=xmin,MAX=xmax,BINSIZE=0.005,LOCATIONS=pdf_NS_x))
  OPLOT, PDF_NS_x, PDF_NS/N_ELEMENTS(f1), PSYM=10



; Individual System NS Probability

  xmin = 0.0
  xmax = 1.0
  ymin = 0.00
  ymax = 40.0

  xlab = 'NS Probability'
  ylab = 'N'

  PLOT, [0], [0], XRANGE=[xmin,xmax], YRANGE=[ymin,ymax],     $
     XSTYLE=1, YSTYLE=1, XTITLE=xlab, YTITLE=ylab,  $
     COLOR=0, POSITION=pos12, /NODATA, /NOERASE

;  pdf_NS = float(HISTOGRAM(l1,MIN=xmin,MAX=xmax,BINSIZE=0.02,LOCATIONS=pdf_NS_x))
;  pdf_true_NS = float(HISTOGRAM(m1,MIN=xmin,MAX=xmax,BINSIZE=0.02,LOCATIONS=pdf_true_NS_x))


  OPLOT, l1, l2, PSYM=10

  dx = (n1[2]-n1[1])/2.0
  FOR i=0, N_ELEMENTS(n2)-1, 1 DO BEGIN
     POLYFILL, [n1[i]-dx,n1[i]-dx,n1[i]+dx,n1[i]+dx,n1[i]-dx], [0,n2[i],n2[i],0,0], COLOR=2
  ENDFOR




; Close ps-file

   DEVICE, /CLOSE

END
