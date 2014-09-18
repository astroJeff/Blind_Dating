 ; Output to postscript file or to screen

   output = 'screen'
   output = 'PS'

; Input and output files
;   infile = COMMAND_LINE_ARGS()

   dir = './'

   fin1 = dir + 'dist_07.dat'
   fin2 = dir + 'dist_07_NS.dat'
   fin3 = dir + 'dist_flat_NS.dat'


   fin4 = dir + 'post_07.dat'
   fin5 = dir + 'post_07_NS.dat'
   fin6 = dir + 'post_flat_NS.dat'

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

  pos1 = [0.06,0.70,0.31,0.97]
  pos2 = [0.39,0.70,0.64,0.97]
  pos3 = [0.72,0.70,0.97,0.97]
  pos4 = [0.06,0.38,0.31,0.65]
  pos5 = [0.39,0.38,0.64,0.65]
  pos6 = [0.72,0.38,0.97,0.65]
  pos7 = [0.06,0.05,0.31,0.33]
  pos8 = [0.39,0.05,0.64,0.33]
  pos9 = [0.72,0.05,0.97,0.33]

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
     DEVICE, /LANDSCAPE,ENCAPSULATE=0, XSIZE=12,YSIZE=12,FILE=fout, SCALE=1.0, $
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


  READCOL, fin4, d1,d2,d3,d4,d5, FORMAT='F,F,F,F,F'
  READCOL, fin5, e1,e2,e3,e4,e5, FORMAT='F,F,F,F,F'
  READCOL, fin6, f1,f2,f3,f4,f5, FORMAT='F,F,F,F,F'



  xmin = 0.0
  xmax = 1.50
  ymin = 0.00
  ymax = 24.0

  xlab = 'Mass (M' + SUNSYMBOL() + ')'
  ylab = 'N'

  PLOT, [0], [0], XRANGE=[xmin,xmax], YRANGE=[ymin,ymax],     $
     XSTYLE=1, YSTYLE=1, XTITLE=xlab, YTITLE=ylab, $
     COLOR=0, POSITION=pos1, /NODATA, /NOERASE

  pdf_M2 = HISTOGRAM(a4,BINSIZE=0.1,LOCATIONS=pdf_M2_x)
  pdf_M2_min = HISTOGRAM(a5,BINSIZE=0.1,LOCATIONS=pdf_M2_min_x)

  OPLOT, PDF_M2_x, PDF_M2, PSYM=10, LINESTYLE=2
  OPLOT, PDF_M2_min_x, PDF_M2_min, PSYM=10


  xmin = 0.45
  xmax = 0.85
  ymin = 0.05
  ymax = 0.35

  xlab = 'Mu'
  ylab = 'Sigma'

  PLOT, [0], [0], XRANGE=[xmin,xmax], YRANGE=[ymin,ymax],     $
     XSTYLE=1, YSTYLE=1, XTITLE=xlab, YTITLE=ylab,  $
     COLOR=0, POSITION=pos2, /NODATA, /NOERASE
  
  CONTOUR_PLUS,d1,d2,XBINSIZE=0.02,YBINSIZE=0.01,XMIN=xmin,XMAX=xmax,YMIN=ymin,YMAX=ymax,PMINY=ymin,PMAXY=ymax,MINX=xmin,PMAXX=xmax,THRESHOLD=1,/INTERP,/NODATA,/NOERASE,/NOCONT

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
  xmax = 0.2
  ymin = 0.00
  ymax = 600

  xlab = 'NS Probability'
  ylab = 'N'

  PLOT, [0], [0], XRANGE=[xmin,xmax], YRANGE=[ymin,ymax],     $
     XSTYLE=1, YSTYLE=1, XTITLE=xlab, YTITLE=ylab,  $
     COLOR=0, POSITION=pos3, /NODATA, /NOERASE

  pdf_NS = HISTOGRAM(d4,BINSIZE=0.005,LOCATIONS=pdf_NS_x)
  OPLOT, PDF_NS_x, PDF_NS, PSYM=10





; Random sample M2, M2_min 

  xmin = 0.0
  xmax = 1.50
  ymin = 0.00
  ymax = 24.0

  xlab = 'Mass (M' + SUNSYMBOL() + ')'
  ylab = 'N'


  PLOT, [0], [0], XRANGE=[xmin,xmax], YRANGE=[ymin,ymax],     $
     XSTYLE=1, YSTYLE=1, XTITLE=xlab, YTITLE=ylab,  $
     COLOR=0, POSITION=pos4, /NODATA, /NOERASE

  pdf_M2 = HISTOGRAM(b4,BINSIZE=0.1,LOCATIONS=pdf_M2_x)
  pdf_M2_min = HISTOGRAM(b5,BINSIZE=0.1,LOCATIONS=pdf_M2_min_x)

  OPLOT, PDF_M2_x, PDF_M2, PSYM=10, LINESTYLE=2
  OPLOT, PDF_M2_min_x, PDF_M2_min, PSYM=10




; Gaussian parameters posterior

  xmin = 0.45
  xmax = 0.85
  ymin = 0.05
  ymax = 0.35

  xlab = 'Mu'
  ylab = 'Sigma'


  PLOT, [0], [0], XRANGE=[xmin,xmax], YRANGE=[ymin,ymax],     $
     XSTYLE=1, YSTYLE=1, XTITLE=xlab, YTITLE=ylab,  $
     COLOR=0, POSITION=pos5, /NODATA, /NOERASE

  CONTOUR_PLUS,e1,e2,XBINSIZE=0.02,YBINSIZE=0.01,XMIN=xmin,XMAX=xmax,YMIN=ymin,YMAX=ymax,PMINY=ymin,PMAXY=ymax,MINX=xmin,PMAXX=xmax,THRESHOLD=1,/INTERP,/NODATA,/NOERASE,/NOCONT

  OPLOT,mu_true,sd_true,PSYM=7,SYMSIZE=2.0



; NS Probability Posterior

  xmin = 0.0
  xmax = 0.2
  ymin = 0.00
  ymax = 600

  xlab = 'NS Probability'
  ylab = 'N'

  PLOT, [0], [0], XRANGE=[xmin,xmax], YRANGE=[ymin,ymax],     $
     XSTYLE=1, YSTYLE=1, XTITLE=xlab, YTITLE=ylab,  $
     COLOR=0, POSITION=pos6, /NODATA, /NOERASE

  pdf_NS = HISTOGRAM(e4,BINSIZE=0.005,LOCATIONS=pdf_NS_x)
  OPLOT, PDF_NS_x, PDF_NS, PSYM=10


; Random sample M2, M2_min 

  xmin = 0.0
  xmax = 1.50
  ymin = 0.00
  ymax = 24.0

  xlab = 'Mass (M' + SUNSYMBOL() + ')'
  ylab = 'N'

  PLOT, [0], [0], XRANGE=[xmin,xmax], YRANGE=[ymin,ymax],     $
     XSTYLE=1, YSTYLE=1, XTITLE=xlab, YTITLE=ylab,  $
     COLOR=0, POSITION=pos7, /NODATA, /NOERASE

  pdf_M2 = HISTOGRAM(c4,BINSIZE=0.1,LOCATIONS=pdf_M2_x)
  pdf_M2_min = HISTOGRAM(c5,BINSIZE=0.1,LOCATIONS=pdf_M2_min_x)

  OPLOT, PDF_M2_x, PDF_M2, PSYM=10, LINESTYLE=2
  OPLOT, PDF_M2_min_x, PDF_M2_min, PSYM=10




; Gaussian parameters posterior

  xmin = 0.05
  xmax = 0.85
  ymin = 0.15
  ymax = 0.60

  xlab = 'Mu'
  ylab = 'Sigma'

  PLOT, [0], [0], XRANGE=[xmin,xmax], YRANGE=[ymin,ymax],     $
     XSTYLE=1, YSTYLE=1, XTITLE=xlab, YTITLE=ylab,  $
     COLOR=0, POSITION=pos8, /NODATA, /NOERASE

  CONTOUR_PLUS,f1,f2,XBINSIZE=0.04,YBINSIZE=0.04,XMIN=xmin,XMAX=xmax,YMIN=ymin,YMAX=ymax,PMINY=ymin,PMAXY=ymax,MINX=xmin,PMAXX=xmax,THRESHOLD=1,/INTERP,/NODATA,/NOERASE,/NOCONT




; NS Probability Posterior

  xmin = 0.0
  xmax = 0.2
  ymin = 0.00
  ymax = 600

  xlab = 'NS Probability'
  ylab = 'N'

  PLOT, [0], [0], XRANGE=[xmin,xmax], YRANGE=[ymin,ymax],     $
     XSTYLE=1, YSTYLE=1, XTITLE=xlab, YTITLE=ylab,  $
     COLOR=0, POSITION=pos9, /NODATA, /NOERASE


  pdf_NS = HISTOGRAM(f4,BINSIZE=0.005,LOCATIONS=pdf_NS_x)
  OPLOT, PDF_NS_x, PDF_NS, PSYM=10





; Close ps-file

   DEVICE, /CLOSE

END
