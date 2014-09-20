 ; Output to postscript file or to screen

   output = 'screen'
   output = 'PS'

; Input and output files
;   infile = COMMAND_LINE_ARGS()

   dir = './'

   fin1 = dir + 'ELM_WD.dat'

   fin2 = dir + 'ELM_post.dat'

   fout = './ELM_post.eps'



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

  pos1 = [0.06,0.15,0.31,0.93]
  pos2 = [0.39,0.15,0.64,0.93]
  pos3 = [0.72,0.15,0.97,0.93]

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
     DEVICE, /LANDSCAPE,ENCAPSULATE=0, XSIZE=12,YSIZE=4,FILE=fout, SCALE=1.0, $
             /INCHES,/SCHOOLBOOK,/COLOR
  ENDIF ELSE BEGIN
     DEVICE, DECOMPOSED=0
     WINDOW, 1, XSIZE=600, YSIZE=600, RETAIN=2
  ENDELSE

; Set background color (polyfill is necessary for postscript output)

  POLYFILL, [1,1,0,0,1], [1,0,0,1,1], /NORMAL, COLOR=7

; Make plot



  READCOL, fin1, a1,a2,a3,a4,a5,a6,a7,a8,a9,a10, $
           a11,a12, FORMAT='A,F,F,F,F,F,F,F,F,F,F,F'

  READCOL, fin2, b1,b2,b3,b4,b5,b6, FORMAT='F,F,F,F,F,F'





  xmin = 0.001
  xmax = 1.50
  ymin = 0.00
  ymax = 12.0

  xlab = 'Mass (M' + SUNSYMBOL() + ')'
  ylab = 'N'

  PLOT, [0], [0], XRANGE=[xmin,xmax], YRANGE=[ymin,ymax],     $
     XSTYLE=1, YSTYLE=1, XTITLE=xlab, YTITLE=ylab, $
     COLOR=0, POSITION=pos1, /NODATA, /NOERASE

  pdf_M2_min = HISTOGRAM(a9,MIN=xmin,MAX=xmax,BINSIZE=0.1,LOCATIONS=pdf_M2_min_x)

  OPLOT, PDF_M2_min_x, PDF_M2_min, PSYM=10, THICK=4





  xmin = 0.05
  xmax = 1.00
  ymin = 0.05
  ymax = 0.70

  xlab = 'Mu'
  ylab = 'Sigma'

  PLOT, [0], [0], XRANGE=[xmin,xmax], YRANGE=[ymin,ymax],     $
     XSTYLE=1, YSTYLE=1, XTITLE=xlab, YTITLE=ylab,  $
     COLOR=0, POSITION=pos2, /NODATA, /NOERASE
  
  CONTOUR_PLUS,b2,b3,XBINSIZE=0.03,YBINSIZE=0.03,XMIN=xmin,XMAX=xmax,YMIN=ymin,YMAX=ymax,PMINY=ymin,PMAXY=ymax,MINX=xmin,PMAXX=xmax,THRESHOLD=1,/INTERP,/NODATA,/NOERASE,/NOCONT





  PLOT, [0], [0], XRANGE=[xmin,xmax], YRANGE=[ymin,ymax],     $
     XSTYLE=1, YSTYLE=1,  $
     COLOR=0, POSITION=pos2, /NODATA, /NOERASE



; NS Probability Posterior

  xmin = 0.0
  xmax = 0.3
  ymin = 0.00
  ymax = 0.07

  xlab = 'NS Probability'
  ylab = 'N'

  PLOT, [0], [0], XRANGE=[xmin,xmax], YRANGE=[ymin,ymax],     $
     XSTYLE=1, YSTYLE=1, XTITLE=xlab, YTITLE=ylab,  $
     COLOR=0, POSITION=pos3, /NODATA, /NOERASE

  pdf_NS = float(HISTOGRAM(b5,MIN=xmin,MAX=xmax,BINSIZE=0.005,LOCATIONS=pdf_NS_x))
  pdf_NS = pdf_NS/N_ELEMENTS(b1)

  OPLOT, PDF_NS_x, PDF_NS, PSYM=10







; Close ps-file

   DEVICE, /CLOSE

END
