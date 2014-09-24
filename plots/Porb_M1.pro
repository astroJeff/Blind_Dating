 ; Output to postscript file or to screen

   output = 'screen'
   output = 'PS'

; Input and output files
;   infile = COMMAND_LINE_ARGS()

   dir = './'

   fin1 = dir + 'ELM_WD.dat'
   fin2 = dir + 'NS_WD.dat'

   fout = './Porb_M1.eps'



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

  pos = [0.10,0.10,0.95,0.93]

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
     DEVICE, /LANDSCAPE,ENCAPSULATE=0, XSIZE=6,YSIZE=6,FILE=fout, SCALE=1.0, $
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

  READCOL, fin2, b1,b2,b3,b4,b5, FORMAT='A,F,F,F,A'


  xmin = 0.001
  xmax = 0.50
  ymin = 0.00
  ymax = 30.0

  xlab = 'Mass (M' + SUNSYMBOL() + ')'
  ylab = 'Orbital Period (hrs)'

  PLOT, [0], [0], XRANGE=[xmin,xmax], YRANGE=[ymin,ymax],     $
     XSTYLE=1, YSTYLE=1, XTITLE=xlab, YTITLE=ylab, $
     COLOR=0, POSITION=pos, /NODATA, /NOERASE

  xtemp = FINDGEN(1)
  ytemp = FINDGEN(1)


  M_1 = FINDGEN(55)
  Porb = FINDGEN(55)
  non_det_M1 = FINDGEN(6)
  non_det_Porb = FINDGEN(6)


  M_1 = a8[0:54]
  Porb = a2[0:54]*24.0

  non_det_M1 = a8[55:60]
  non_det_Porb = a2[55:60]*24.0

  PLOTSYM, 0
  OPLOT, M_1, Porb, PSYM=8, THICK=2

  PLOTSYM, 0, /FILL
  FOR i=0,N_ELEMENTS(M_1)-1 DO BEGIN
     IF (a1[i] EQ "J0651+2844") THEN BEGIN
        xtemp[0] = M_1[i]
        ytemp[0] = Porb[i]
        OPLOT, xtemp, ytemp, PSYM=8, THICK=2
     ENDIF
     IF (a1[i] EQ "J0106-1000") THEN BEGIN
        xtemp[0] = M_1[i]
        ytemp[0] = Porb[i]
        OPLOT, xtemp, ytemp, PSYM=8, THICK=2
     ENDIF
     IF (a1[i] EQ "J0345+1748") THEN BEGIN
        xtemp[0] = M_1[i]
        ytemp[0] = Porb[i]
        OPLOT, xtemp, ytemp, PSYM=8, THICK=2
     ENDIF
  ENDFOR
  

  FOR i=0,N_ELEMENTS(non_det_M1)-1 DO BEGIN 
    ARROW, non_det_M1[i],3.0,non_det_M1[i],0.0, THICK=1, /DATA 
  ENDFOR



  PLOTSYM, 4, /FILL

  OPLOT, b2, b4*24.0, PSYM=8, THICK=4


; Close ps-file

   DEVICE, /CLOSE

END
