# -*- coding: utf-8 -*-
if 64 - 64: i11iIiiIii
import shapely
import shapely . wkt
from osgeo import ogr
import numpy
from math import sqrt
import fileIO
from subprocess import call
if 65 - 65: O0 / iIii1I11I1II1 % OoooooooOO - i1IIi
if 73 - 73: II111iiii
def fillEdgeCells ( cellDict , wl ) :
 for IiII1IiiIiI1 , iIiiiI1IiI1I1 in cellDict . keys ( ) :
  wl [ IiII1IiiIiI1 , iIiiiI1IiI1I1 ] = wl [ cellDict [ ( IiII1IiiIiI1 , iIiiiI1IiI1I1 ) ] [ 0 ] , cellDict [ ( IiII1IiiIiI1 , iIiiiI1IiI1I1 ) ] [ 1 ] ]
  if 87 - 87: OoOoOO00
 return
 if 27 - 27: OOOo0 / Oo - Ooo00oOo00o . I1IiI
 if 73 - 73: OOooOOo / ii11ii1ii
def maskZeroPoly ( wktFileName , rasterFileName ) :
 if 94 - 94: OoOO + OoOO0ooOOoo0O + o0000oOoOoO0o * o00O0oo
 print wktFileName , rasterFileName
 if 97 - 97: oO0o0ooO0 - IIII / O0oO - OoOO
 iiI11iii111 = [ 'gdal_rasterize' ]
 iiI11iii111 += [ '-burn' , '0' ]
 if 37 - 37: o00O0oo
 iiI11iii111 += [ wktFileName ]
 iiI11iii111 += [ rasterFileName ]
 if 71 - 71: OoOO * o00O0oo . OoOO / IIII
 call ( iiI11iii111 )
 if 14 - 14: iIii1I11I1II1
 return
 if 56 - 56: OoOO0ooOOoo0O - i1IIi
 if 64 - 64: IIII + o00O0oo
def weirFlowCalc ( wl1 , wl2 , cl , len , wc , ml ) :
 if 10 - 10: i11iIiiIii / ii11ii1ii % II111iiii
 Ooo00O0 = wl1 - cl
 oo0 = wl2 - cl
 if 80 - 80: i1IIi * iIii1I11I1II1 / oO0o0ooO0 % ii11ii1ii + IIII . oO0o0ooO0
 if Ooo00O0 < 0 and oo0 < 0 :
  return 0
  if 98 - 98: i11iIiiIii * OoOoOO00 % o00O0oo * o00O0oo * II111iiii
 if wl1 > wl2 and ( oo0 / Ooo00O0 ) <= ml :
  return wc * len * ( Ooo00O0 ** 1.5 )
  if 79 - 79: oO0o0ooO0
 if wl1 > wl2 and ( oo0 / Ooo00O0 ) > ml :
  return wc * len * Ooo00O0 * ( ( ( Ooo00O0 - oo0 ) / ( 1 - ml ) ) ** 0.5 )
  if 86 - 86: Ooo00oOo00o % OoOoOO00
 if wl1 < wl2 and ( Ooo00O0 / oo0 ) <= ml :
  return - wc * len * ( oo0 ** 1.5 )
  if 80 - 80: OoooooooOO . OoOoOO00
 if wl1 < wl2 and ( Ooo00O0 / oo0 ) > ml :
  return - wc * len * oo0 * ( ( ( oo0 - Ooo00O0 ) / ( 1 - ml ) ) ** 0.5 )
  if 87 - 87: ii11ii1ii / O0oO + IIII - O0oO . O0oO / II111iiii
 return 0.0
 if 11 - 11: OoOoOO00 % I1IiI - OOOo0
def calcSpills ( wlGrid , spills , qX , qY ) :
 oo0O000OoO = 0.9
 if 34 - 34: OoOO0ooOOoo0O * OoOoOO00
 for ( IiII1IiiIiI1 , iIiiiI1IiI1I1 , dir ) in spills . keys ( ) :
  if 31 - 31: II111iiii + Oo . IIII
  OoOooOOOO = spills [ ( IiII1IiiIiI1 , iIiiiI1IiI1I1 , dir ) ] [ 0 ]
  i11iiII = spills [ ( IiII1IiiIiI1 , iIiiiI1IiI1I1 , dir ) ] [ 1 ]
  I1iiiiI1iII = spills [ ( IiII1IiiIiI1 , iIiiiI1IiI1I1 , dir ) ] [ 2 ]
  if 20 - 20: i1IIi + i11iIiiIii - o0000oOoOoO0o % Oo . OoooooooOO
  if dir == 'x' :
   Ooo00O00O0O0O = wlGrid [ IiII1IiiIiI1 - 1 , iIiiiI1IiI1I1 ]
   OooO0OO = wlGrid [ IiII1IiiIiI1 , iIiiiI1IiI1I1 ]
   if 28 - 28: II111iiii
   qX [ IiII1IiiIiI1 , iIiiiI1IiI1I1 ] += weirFlowCalc ( Ooo00O00O0O0O , OooO0OO , OoOooOOOO , i11iiII , I1iiiiI1iII , oo0O000OoO )
   if 28 - 28: iIii1I11I1II1 - i1IIi
   if 70 - 70: Oo . Oo - Oo / OOooOOo * OoOO
   if len ( spills [ ( IiII1IiiIiI1 , iIiiiI1IiI1I1 , dir ) ] ) > 3 :
    if 86 - 86: i11iIiiIii + o0000oOoOoO0o + O0oO * OoOO0ooOOoo0O + I1IiI
    OoOooOOOO = spills [ ( IiII1IiiIiI1 , iIiiiI1IiI1I1 , dir ) ] [ 3 ]
    i11iiII = spills [ ( IiII1IiiIiI1 , iIiiiI1IiI1I1 , dir ) ] [ 4 ]
    I1iiiiI1iII = spills [ ( IiII1IiiIiI1 , iIiiiI1IiI1I1 , dir ) ] [ 5 ]
    if 61 - 61: Oo / i11iIiiIii
    qX [ IiII1IiiIiI1 , iIiiiI1IiI1I1 ] += weirFlowCalc ( Ooo00O00O0O0O , OooO0OO , OoOooOOOO , i11iiII , I1iiiiI1iII , oo0O000OoO )
    if 34 - 34: OoooooooOO + iIii1I11I1II1 + i11iIiiIii - OOooOOo + i11iIiiIii
  if dir == 'y' :
   Ooo00O00O0O0O = wlGrid [ IiII1IiiIiI1 , iIiiiI1IiI1I1 ]
   OooO0OO = wlGrid [ IiII1IiiIiI1 , iIiiiI1IiI1I1 - 1 ]
   if 65 - 65: Ooo00oOo00o
   qY [ IiII1IiiIiI1 , iIiiiI1IiI1I1 ] += weirFlowCalc ( Ooo00O00O0O0O , OooO0OO , OoOooOOOO , i11iiII , I1iiiiI1iII , oo0O000OoO )
   if 6 - 6: OoOoOO00 / OOOo0 % o0000oOoOoO0o
   if 84 - 84: i11iIiiIii . I1IiI
   if len ( spills [ ( IiII1IiiIiI1 , iIiiiI1IiI1I1 , dir ) ] ) > 3 :
    OoOooOOOO = spills [ ( IiII1IiiIiI1 , iIiiiI1IiI1I1 , dir ) ] [ 3 ]
    i11iiII = spills [ ( IiII1IiiIiI1 , iIiiiI1IiI1I1 , dir ) ] [ 4 ]
    I1iiiiI1iII = spills [ ( IiII1IiiIiI1 , iIiiiI1IiI1I1 , dir ) ] [ 5 ]
    if 100 - 100: o0000oOoOoO0o - o0000oOoOoO0o - IIII
    qY [ IiII1IiiIiI1 , iIiiiI1IiI1I1 ] += weirFlowCalc ( Ooo00O00O0O0O , OooO0OO , OoOooOOOO , i11iiII , I1iiiiI1iII , oo0O000OoO )
    if 20 - 20: OoooooooOO
 return
 if 13 - 13: i1IIi - o0000oOoOoO0o % ii11ii1ii / iIii1I11I1II1 % o00O0oo
 if 97 - 97: i11iIiiIii
 if 32 - 32: OOOo0 * O0 % ii11ii1ii % o0000oOoOoO0o . oO0o0ooO0
 if 61 - 61: O0oO
def spillways ( xll , yll , cellSize , xsz , ysz , spills , shpFileName ) :
 oOOO00o = fileIO . readPointShapefile ( shpFileName )
 if 97 - 97: OoOO0ooOOoo0O % OoOO0ooOOoo0O + II111iiii * o00O0oo
 for o0o00o0 , iIi1ii1I1 , o0 in oOOO00o :
  IiII1IiiIiI1 = int ( ( o0o00o0 - xll ) / cellSize )
  iIiiiI1IiI1I1 = int ( ( iIi1ii1I1 - yll ) / cellSize )
  if 9 - 9: o0000oOoOoO0o + ii11ii1ii % o0000oOoOoO0o + i1IIi . OoOO
  III1i1i = IiII1IiiIiI1 * cellSize + xll
  iiI1 = iIiiiI1IiI1I1 * cellSize + yll
  if 19 - 19: OoOO0ooOOoo0O + O0oO
  if 53 - 53: OoooooooOO . i1IIi
  ii1I1i1I = [ ]
  if ( IiII1IiiIiI1 , iIiiiI1IiI1I1 , 'x' ) in spills . keys ( ) :
   ii1I1i1I . append ( 'w' )
  if ( IiII1IiiIiI1 + 1 , iIiiiI1IiI1I1 , 'x' ) in spills . keys ( ) :
   ii1I1i1I . append ( 'e' )
  if ( IiII1IiiIiI1 , iIiiiI1IiI1I1 , 'y' ) in spills . keys ( ) :
   ii1I1i1I . append ( 's' )
  if ( IiII1IiiIiI1 , iIiiiI1IiI1I1 + 1 , 'y' ) in spills . keys ( ) :
   ii1I1i1I . append ( 'n' )
   if 88 - 88: Oo + O0 / Ooo00oOo00o * o00O0oo
   if 41 - 41: ii11ii1ii
  if len ( ii1I1i1I ) == 0 :
   print "***** Warning - spillway at %f %f cannot be assigned to edge - ignoring ***" % ( o0o00o0 , iIi1ii1I1 )
   continue
   if 6 - 6: OOooOOo
   if 31 - 31: o0000oOoOoO0o . o0000oOoOoO0o - I1IiI / Oo + O0oO * OoOoOO00
  if len ( ii1I1i1I ) > 1 :
   O0ooOooooO = 1e20
   if 'n' in ii1I1i1I :
    o00O = iiI1 + cellSize - iIi1ii1I1
    if o00O < O0ooOooooO :
     O0ooOooooO = o00O
     OOO0OOO00oo = 'n'
   if 's' in ii1I1i1I :
    o00O = iIi1ii1I1 - iiI1
    if o00O < O0ooOooooO :
     O0ooOooooO = o00O
     OOO0OOO00oo = 's'
   if 'e' in ii1I1i1I :
    o00O = III1i1i + cellSize - o0o00o0
    if o00O < O0ooOooooO :
     O0ooOooooO = o00O
     OOO0OOO00oo = 'e'
   if 'w' in ii1I1i1I :
    o00O = o0o00o0 - III1i1i
    if o00O < O0ooOooooO :
     O0ooOooooO = o00O
     OOO0OOO00oo = 'w'
  else :
   OOO0OOO00oo = ii1I1i1I [ 0 ]
   if 31 - 31: II111iiii - OoOO . IIII % Ooo00oOo00o - O0
  if OOO0OOO00oo == 'n' :
   spills [ ( IiII1IiiIiI1 , iIiiiI1IiI1I1 + 1 , 'y' ) ] += [ o0 [ 'cl' ] , o0 [ 'width' ] , o0 [ 'wc' ] ]
   spills [ ( IiII1IiiIiI1 , iIiiiI1IiI1I1 + 1 , 'y' ) ] [ 1 ] -= o0 [ 'width' ]
  if OOO0OOO00oo == 's' :
   spills [ ( IiII1IiiIiI1 , iIiiiI1IiI1I1 , 'y' ) ] += [ o0 [ 'cl' ] , o0 [ 'width' ] , o0 [ 'wc' ] ]
   spills [ ( IiII1IiiIiI1 , iIiiiI1IiI1I1 , 'y' ) ] [ 1 ] -= o0 [ 'width' ]
  if OOO0OOO00oo == 'e' :
   spills [ ( IiII1IiiIiI1 + 1 , iIiiiI1IiI1I1 , 'x' ) ] += [ o0 [ 'cl' ] , o0 [ 'width' ] , o0 [ 'wc' ] ]
   spills [ ( IiII1IiiIiI1 + 1 , iIiiiI1IiI1I1 , 'x' ) ] [ 1 ] -= o0 [ 'width' ]
  if OOO0OOO00oo == 'w' :
   spills [ ( IiII1IiiIiI1 , iIiiiI1IiI1I1 , 'x' ) ] += [ o0 [ 'cl' ] , o0 [ 'width' ] , o0 [ 'wc' ] ]
   spills [ ( IiII1IiiIiI1 , iIiiiI1IiI1I1 , 'x' ) ] [ 1 ] -= o0 [ 'width' ]
   if 4 - 4: II111iiii / O0oO . o00O0oo
 return
 if 58 - 58: OoOO * i11iIiiIii / Ooo00oOo00o % IIII - OOooOOo / ii11ii1ii
 if 50 - 50: OoOoOO00
 if 34 - 34: OoOoOO00 * II111iiii % o00O0oo * Ooo00oOo00o - OoOoOO00
def overlaps ( xmin1 , xmax1 , ymin1 , ymax1 , xmin2 , xmax2 , ymin2 , ymax2 ) :
 return ( xmin1 < xmax2 and ymin1 < ymax2 and xmax1 > xmin2 and ymax1 > ymin2 )
 if 33 - 33: I1IiI + OoOO * Oo - OOOo0 / ii11ii1ii % o0000oOoOoO0o
def almostSame ( x1 , x2 ) :
 return abs ( ( x1 - x2 ) ) / ( 0.5 * ( x1 + x2 ) ) < 1e-6
 if 21 - 21: Oo * iIii1I11I1II1 % ii11ii1ii * i1IIi
def edgeSame ( e1 , e2 ) :
 return e1 == e2 or e1 == [ e2 [ 2 ] , e2 [ 3 ] , e2 [ 0 ] , e2 [ 1 ] ]
 if 16 - 16: O0 - IIII * iIii1I11I1II1 + o00O0oo
def addEdgeArray ( arrX , arrY , clX , clY , edge , xll , yll , dx ) :
 Ii11iII1 = edge [ 0 ]
 Oo0O0O0ooO0O = edge [ 1 ]
 IIIIii = edge [ 2 ]
 O0o0 = edge [ 3 ]
 if 71 - 71: OoOO + O0oO % i11iIiiIii + OOooOOo - oO0o0ooO0
 if Ii11iII1 == IIIIii :
  oO0OOoO0 = int ( ( Ii11iII1 - xll ) / dx )
  I111Ii111 = int ( ( 0.5 * ( Oo0O0O0ooO0O + O0o0 ) - yll ) / dx )
  if 4 - 4: ii11ii1ii
  arrX [ oO0OOoO0 , I111Ii111 ] = True
  if 93 - 93: Oo % ii11ii1ii . Oo * IIII % o0000oOoOoO0o . II111iiii
  clX [ oO0OOoO0 , I111Ii111 ] = dx
  if 38 - 38: I1IiI
 if Oo0O0O0ooO0O == O0o0 :
  oO0OOoO0 = int ( ( 0.5 * ( Ii11iII1 + IIIIii ) - xll ) / dx )
  I111Ii111 = int ( ( Oo0O0O0ooO0O - yll ) / dx )
  if 57 - 57: O0 / ii11ii1ii * IIII / Ooo00oOo00o . II111iiii
  arrY [ oO0OOoO0 , I111Ii111 ] = True
  clY [ oO0OOoO0 , I111Ii111 ] = dx
  if 26 - 26: o00O0oo
 return
 if 91 - 91: Oo . OOooOOo + Oo - o00O0oo / OoooooooOO
def rmEdgeArray ( arrX , arrY , edge , xll , yll , dx ) :
 Ii11iII1 = edge [ 0 ]
 Oo0O0O0ooO0O = edge [ 1 ]
 IIIIii = edge [ 2 ]
 O0o0 = edge [ 3 ]
 if 39 - 39: OOooOOo / O0oO - II111iiii
 if Ii11iII1 == IIIIii :
  oO0OOoO0 = int ( ( Ii11iII1 - xll ) / dx )
  I111Ii111 = int ( ( 0.5 * ( Oo0O0O0ooO0O + O0o0 ) - yll ) / dx )
  if 98 - 98: OOooOOo / OoOO0ooOOoo0O % ii11ii1ii . Ooo00oOo00o
  arrX [ oO0OOoO0 , I111Ii111 ] = False
  if 91 - 91: ii11ii1ii % OOOo0
 if Oo0O0O0ooO0O == O0o0 :
  oO0OOoO0 = int ( ( 0.5 * ( Ii11iII1 + IIIIii ) - xll ) / dx )
  I111Ii111 = int ( ( Oo0O0O0ooO0O - yll ) / dx )
  if 64 - 64: OoOO0ooOOoo0O % o00O0oo - IIII - ii11ii1ii
  arrY [ oO0OOoO0 , I111Ii111 ] = False
  if 31 - 31: OoOO0ooOOoo0O - II111iiii . OoOO0ooOOoo0O
 return
 if 18 - 18: I1IiI
def vectorProduct ( x1 , y1 , x2 , y2 ) :
 return x1 * y2 - y1 * x2
 if 98 - 98: o00O0oo * o00O0oo / o00O0oo + OoOO0ooOOoo0O
def modifyConvPar ( xll , yll , cellSize , xsz , ysz , convParX , convParY , shpFileName , csvEdgeFileName = None ) :
 if 34 - 34: O0oO
 I1111I1iII11 = ogr . Open ( shpFileName , 0 )
 if 59 - 59: iIii1I11I1II1 * i11iIiiIii / OOooOOo * i1IIi * O0
 if 83 - 83: Oo / IIII . Ooo00oOo00o / oO0o0ooO0 . Ooo00oOo00o . OoOO
 O0oOoOO = I1111I1iII11 . GetLayer ( )
 oO00o0 = [ shapely . wkt . loads ( OOoo0O . GetGeometryRef ( ) . ExportToWkt ( ) ) for OOoo0O in O0oOoOO ]
 if 67 - 67: i11iIiiIii - i1IIi % OOooOOo . O0
 O0oOoOO . ResetReading ( )
 o0oo = [ OOoo0O . GetField ( "cl" ) for OOoo0O in O0oOoOO ]
 if 91 - 91: oO0o0ooO0
 O0oOoOO . ResetReading ( )
 iiIii = [ OOoo0O . GetField ( "wc" ) for OOoo0O in O0oOoOO ]
 if 79 - 79: OoooooooOO / O0
 if 75 - 75: Ooo00oOo00o % I1IiI % I1IiI . IIII
 III1iII1I1ii = numpy . zeros ( ( xsz + 1 , ysz ) , dtype = bool )
 oOOo0 = numpy . zeros ( ( xsz , ysz + 1 ) , dtype = bool )
 if 54 - 54: O0 - oO0o0ooO0 % OoOO
 OOoO = { }
 iII = { }
 ii1ii11IIIiiI = { }
 if 67 - 67: OoOO0ooOOoo0O * ii11ii1ii * OOooOOo + OoOO / i1IIi
 I1I111 = 0
 if 82 - 82: i11iIiiIii - o00O0oo * OoooooooOO / OoOO0ooOOoo0O
 if 31 - 31: oO0o0ooO0 . Oo - iIii1I11I1II1
 ooOOO00Ooo = { }
 IiIIIi1iIi = [ ]
 if 68 - 68: i11iIiiIii % OOooOOo + i11iIiiIii
 for iii , II1I in enumerate ( oO00o0 ) :
  if 84 - 84: oO0o0ooO0 . i11iIiiIii . oO0o0ooO0 * OOooOOo - OoOO0ooOOoo0O
  I1I111 += II1I . length
  if 42 - 42: i11iIiiIii
  I11i1iIII = [ ]
  iiIiI = [ ]
  if 91 - 91: o00O0oo % i1IIi % iIii1I11I1II1
  for IIi1I11I1II in II1I . coords :
   I11i1iIII . append ( IIi1I11I1II [ 0 ] )
   iiIiI . append ( IIi1I11I1II [ 1 ] )
   if 63 - 63: OoooooooOO - Oo . II111iiii / I1IiI . Ooo00oOo00o / O0
  o0OOOO00O0Oo = int ( ( I11i1iIII [ 0 ] - xll ) / cellSize )
  ii = int ( ( iiIiI [ 0 ] - yll ) / cellSize )
  if 90 - 90: I1IiI % i1IIi / Oo
  IIi = True
  i1Iii1i1I = False
  OOoO00 = [ ( I11i1iIII [ 0 ] , iiIiI [ 0 ] ) ]
  IiI111111IIII = [ ]
  if 37 - 37: IIII / Ooo00oOo00o
  if 23 - 23: O0
  o00oO0oOo00 = 0
  if 81 - 81: Oo
  if 42 - 42: Oo / OoOO0ooOOoo0O / I1IiI + o00O0oo / Ooo00oOo00o
  while True :
   o00oO0oOo00 += 1
   if o00oO0oOo00 > 100 :
    break
    if 84 - 84: O0oO * II111iiii + OOOo0
    if 53 - 53: o00O0oo % II111iiii . oO0o0ooO0 - iIii1I11I1II1 - oO0o0ooO0 * II111iiii
    if 77 - 77: iIii1I11I1II1 * Oo
   III1i1i = xll + o0OOOO00O0Oo * cellSize
   Ii11iII1 = III1i1i + cellSize
   iiI1 = yll + ii * cellSize
   Oo0O0O0ooO0O = iiI1 + cellSize
   if 95 - 95: OoOoOO00 + i11iIiiIii
   I1Ii = 0.5 * ( III1i1i + Ii11iII1 )
   O0oo00o0O = 0.5 * ( iiI1 + Oo0O0O0ooO0O )
   if 1 - 1: II111iiii
   OOooooO0Oo = "LINESTRING(%f %f,%f %f,%f %f,%f %f,%f %f)" % ( III1i1i , iiI1 , Ii11iII1 , iiI1 , Ii11iII1 , Oo0O0O0ooO0O , III1i1i , Oo0O0O0ooO0O , III1i1i , iiI1 )
   OO = shapely . wkt . loads ( OOooooO0Oo )
   if 25 - 25: Oo
   if 62 - 62: OoOO + O0
   oO0OOOO0 = OO . intersection ( II1I )
   if oO0OOOO0 . type == 'MultiPoint' :
    iI1I11iiI1i = [ ]
    for IiII1IiiIiI1 in range ( len ( oO0OOOO0 ) ) :
     iI1I11iiI1i . append ( oO0OOOO0 [ IiII1IiiIiI1 ] . coords [ 0 ] )
     if 78 - 78: ii11ii1ii % O0 % o0000oOoOoO0o
     if 46 - 46: OoooooooOO . i11iIiiIii
    for OOo0oO00ooO00 , oOO0O00oO0Ooo in iI1I11iiI1i :
     if ( OOo0oO00ooO00 , oOO0O00oO0Ooo ) not in OOoO00 :
      OOoO00 . append ( ( OOo0oO00ooO00 , oOO0O00oO0Ooo ) )
      if 67 - 67: Oo - OoOO
   else :
    if 36 - 36: oO0o0ooO0
    iI1I11iiI1i = oO0OOOO0 . coords [ 0 ]
    if IIi :
     OOoO00 . append ( iI1I11iiI1i )
    else :
     OOoO00 . append ( ( I11i1iIII [ - 1 ] , iiIiI [ - 1 ] ) )
     i1Iii1i1I = True
     if 36 - 36: O0oO / O0 * OOOo0 - OoOO % iIii1I11I1II1 * ii11ii1ii
     if 79 - 79: O0
   iI1I11iiI1i = iI1I11iiI1i [ : 2 ]
   if 78 - 78: OOooOOo + OoOO - IIII
   if 38 - 38: I1IiI - ii11ii1ii + iIii1I11I1II1 / Ooo00oOo00o % OOOo0
   if 57 - 57: Oo / O0oO
   OOooooO0Oo = "POLYGON((%f %f,%f %f,%f %f,%f %f,%f %f))" % ( III1i1i , iiI1 , Ii11iII1 , iiI1 , Ii11iII1 , Oo0O0O0ooO0O , III1i1i , Oo0O0O0ooO0O , III1i1i , iiI1 )
   OO = shapely . wkt . loads ( OOooooO0Oo )
   Ii1I1Ii = OO . intersection ( II1I )
   i11iiII = Ii1I1Ii . length
   if 69 - 69: OoOoOO00 / I1IiI . oO0o0ooO0 * IIII % o0000oOoOoO0o - I1IiI
   if 13 - 13: o0000oOoOoO0o . i11iIiiIii
   oOOoo00O00o = Ii1I1Ii . buffer ( 1.0 )
   O0O00Oo = OO . difference ( oOOoo00O00o )
   if 97 - 97: O0 * OoooooooOO . OoooooooOO
   if O0O00Oo . type == 'MultiPolygon' and len ( O0O00Oo ) == 2 :
    I111iI , oOOo0II1I1iiIII = O0O00Oo [ 0 ] . centroid . coords [ 0 ]
    oOOo0O00o , iIiIi11 = O0O00Oo [ 1 ] . centroid . coords [ 0 ]
    if 87 - 87: OOOo0 . OoOoOO00 - II111iiii + O0 / OOOo0 / ii11ii1ii
    IiI , IIIii1I = Ii1I1Ii . coords [ 0 ]
    ooO0OO , O000OOO = Ii1I1Ii . coords [ - 1 ]
    if 20 - 20: I1IiI . II111iiii % OoOO * iIii1I11I1II1
    oO00oOOoooO = vectorProduct ( ooO0OO - IiI , O000OOO - IIIii1I , oOOo0O00o - I111iI , iIiIi11 - oOOo0II1I1iiIII )
    if 46 - 46: OoOoOO00 - OoooooooOO - OoOO0ooOOoo0O * II111iiii
    if oO00oOOoooO > 0 :
     if 34 - 34: OoOO0ooOOoo0O - o00O0oo / OoOO + OOooOOo * o0000oOoOoO0o
     IiIIIi1iIi . append ( O0O00Oo [ 1 ] )
    else :
     if 73 - 73: Ooo00oOo00o . o0000oOoOoO0o * OOooOOo % OOooOOo % OoooooooOO
     IiIIIi1iIi . append ( O0O00Oo [ 0 ] )
     if 63 - 63: iIii1I11I1II1 * i11iIiiIii % iIii1I11I1II1 * i11iIiiIii
     if 32 - 32: OoOO
     if 42 - 42: oO0o0ooO0 * O0 % i1IIi . OoOO / I1IiI
     if 32 - 32: OoOoOO00 * OOOo0
     if 78 - 78: OoOO - OoooooooOO - OOooOOo / O0oO / II111iiii
   iiI11ii1I1 = 0.5 * ( OOoO00 [ - 1 ] [ 0 ] + OOoO00 [ - 2 ] [ 0 ] )
   Ooo0OOoOoO0 = 0.5 * ( OOoO00 [ - 1 ] [ 1 ] + OOoO00 [ - 2 ] [ 1 ] )
   if 70 - 70: ii11ii1ii
   if 59 - 59: I1IiI % ii11ii1ii
   if almostSame ( OOoO00 [ - 1 ] [ 1 ] , iiI1 ) :
    if 6 - 6: iIii1I11I1II1 % i11iIiiIii % OOooOOo
    if 93 - 93: oO0o0ooO0 * OoooooooOO + O0oO
    if not IIi and iiI11ii1I1 > I1Ii :
     IiI111111IIII . append ( [ Ii11iII1 , Oo0O0O0ooO0O , Ii11iII1 , iiI1 ] )
     III1iII1I1ii [ o0OOOO00O0Oo + 1 , ii ] = True
     OOoO [ ( o0OOOO00O0Oo + 1 , ii ) ] = i11iiII
     if 33 - 33: O0 * I1IiI - IIII % IIII
     if 18 - 18: IIII / OOOo0 * IIII + IIII * i11iIiiIii * OOooOOo
    if not IIi and iiI11ii1I1 <= I1Ii :
     IiI111111IIII . append ( [ III1i1i , Oo0O0O0ooO0O , III1i1i , iiI1 ] )
     III1iII1I1ii [ o0OOOO00O0Oo , ii ] = True
     OOoO [ ( o0OOOO00O0Oo , ii ) ] = i11iiII
     if 11 - 11: O0oO / Ooo00oOo00o - oO0o0ooO0 * OoooooooOO + OoooooooOO . Ooo00oOo00o
     ooOOO00Ooo [ ( o0OOOO00O0Oo , ii ) ] = ( o0OOOO00O0Oo - 1 , ii )
     if 26 - 26: o0000oOoOoO0o % OOooOOo
     if 76 - 76: oO0o0ooO0 * o00O0oo
    ii -= 1
    if 52 - 52: OoOO
    if 19 - 19: OoOoOO00
   if almostSame ( OOoO00 [ - 1 ] [ 1 ] , Oo0O0O0ooO0O ) :
    if 25 - 25: o0000oOoOoO0o / O0oO
    if 31 - 31: OoOO . O0 % OoOoOO00 . I1IiI + oO0o0ooO0
    if not IIi and iiI11ii1I1 > I1Ii :
     IiI111111IIII . append ( [ Ii11iII1 , iiI1 , Ii11iII1 , Oo0O0O0ooO0O ] )
     III1iII1I1ii [ o0OOOO00O0Oo + 1 , ii ] = True
     OOoO [ ( o0OOOO00O0Oo + 1 , ii ) ] = i11iiII
     if 71 - 71: IIII . II111iiii
     ooOOO00Ooo [ ( o0OOOO00O0Oo , ii ) ] = ( o0OOOO00O0Oo + 1 , ii )
     if 62 - 62: OoooooooOO . OoOO0ooOOoo0O
     if 61 - 61: Ooo00oOo00o - OoOO - i1IIi
    if not IIi and iiI11ii1I1 <= I1Ii :
     IiI111111IIII . append ( [ III1i1i , iiI1 , III1i1i , Oo0O0O0ooO0O ] )
     III1iII1I1ii [ o0OOOO00O0Oo , ii ] = True
     OOoO [ ( o0OOOO00O0Oo , ii ) ] = i11iiII
     if 25 - 25: O0 * OoOO0ooOOoo0O + OOooOOo . I1IiI . I1IiI
    ii += 1
    if 58 - 58: OoOoOO00
    if 53 - 53: i1IIi
   if almostSame ( OOoO00 [ - 1 ] [ 0 ] , III1i1i ) :
    if 59 - 59: I1IiI
    if 81 - 81: Ooo00oOo00o - Ooo00oOo00o . o00O0oo
    if not IIi and Ooo0OOoOoO0 > O0oo00o0O :
     IiI111111IIII . append ( [ Ii11iII1 , Oo0O0O0ooO0O , III1i1i , Oo0O0O0ooO0O ] )
     oOOo0 [ o0OOOO00O0Oo , ii + 1 ] = True
     iII [ ( o0OOOO00O0Oo , ii + 1 ) ] = i11iiII
     if 73 - 73: OoOO0ooOOoo0O % i11iIiiIii - OoOoOO00
     ooOOO00Ooo [ ( o0OOOO00O0Oo , ii ) ] = ( o0OOOO00O0Oo , ii + 1 )
     if 7 - 7: O0 * i11iIiiIii * o0000oOoOoO0o + O0oO % Oo - O0oO
     if 39 - 39: OOOo0 * OoOO % OoOO - OoooooooOO + I1IiI - OoOO0ooOOoo0O
    if not IIi and Ooo0OOoOoO0 <= O0oo00o0O :
     IiI111111IIII . append ( [ Ii11iII1 , iiI1 , III1i1i , iiI1 ] )
     oOOo0 [ o0OOOO00O0Oo , ii ] = True
     iII [ ( o0OOOO00O0Oo , ii ) ] = i11iiII
     if 23 - 23: i11iIiiIii
    o0OOOO00O0Oo -= 1
    if 30 - 30: I1IiI - i1IIi % II111iiii + OoOO0ooOOoo0O * iIii1I11I1II1
   if almostSame ( OOoO00 [ - 1 ] [ 0 ] , Ii11iII1 ) :
    if 81 - 81: oO0o0ooO0 % i1IIi . iIii1I11I1II1
    if 4 - 4: i11iIiiIii % Oo % i1IIi / oO0o0ooO0
    if not IIi and Ooo0OOoOoO0 > O0oo00o0O :
     IiI111111IIII . append ( [ III1i1i , Oo0O0O0ooO0O , Ii11iII1 , Oo0O0O0ooO0O ] )
     oOOo0 [ o0OOOO00O0Oo , ii + 1 ] = True
     iII [ ( o0OOOO00O0Oo , ii + 1 ) ] = i11iiII
     if 6 - 6: o00O0oo / OoOoOO00 % OoOO - OoOoOO00
    if not IIi and Ooo0OOoOoO0 <= O0oo00o0O :
     IiI111111IIII . append ( [ III1i1i , iiI1 , Ii11iII1 , iiI1 ] )
     oOOo0 [ o0OOOO00O0Oo , ii ] = True
     iII [ ( o0OOOO00O0Oo , ii ) ] = i11iiII
     if 31 - 31: OoOO
     ooOOO00Ooo [ ( o0OOOO00O0Oo , ii ) ] = ( o0OOOO00O0Oo , ii - 1 )
     if 23 - 23: IIII . oO0o0ooO0
    o0OOOO00O0Oo += 1
    if 92 - 92: Ooo00oOo00o + IIII * o0000oOoOoO0o % OoOoOO00
    if 42 - 42: OOOo0
   IIi = False
   if i1Iii1i1I :
    break
    if 76 - 76: OoOoOO00 * o00O0oo % IIII
    if 57 - 57: iIii1I11I1II1 - i1IIi / IIII - O0 * OoooooooOO % II111iiii
  for IiII1IiiIiI1 , OOO0OOO00oo in enumerate ( IiI111111IIII [ : - 1 ] ) :
   if OOO0OOO00oo [ 2 ] != IiI111111IIII [ IiII1IiiIiI1 + 1 ] [ 0 ] or OOO0OOO00oo [ 3 ] != IiI111111IIII [ IiII1IiiIiI1 + 1 ] [ 1 ] :
    if 68 - 68: OoooooooOO * OoOO0ooOOoo0O % Ooo00oOo00o - oO0o0ooO0
    IiI111111IIII . append ( [ OOO0OOO00oo [ 2 ] , OOO0OOO00oo [ 3 ] , IiI111111IIII [ IiII1IiiIiI1 + 1 ] [ 0 ] , IiI111111IIII [ IiII1IiiIiI1 + 1 ] [ 1 ] ] )
    addEdgeArray ( III1iII1I1ii , oOOo0 , OOoO , iII , IiI111111IIII [ - 1 ] , xll , yll , cellSize )
    if 34 - 34: IIII . iIii1I11I1II1 * Ooo00oOo00o * ii11ii1ii / IIII / OOooOOo
    if 78 - 78: OOOo0 - I1IiI / Ooo00oOo00o
  for IiII1IiiIiI1 , I11IIIi in enumerate ( IiI111111IIII [ : ] ) :
   for iIIiiI1II1i11 in IiI111111IIII [ IiII1IiiIiI1 + 1 : ] :
    if edgeSame ( I11IIIi , iIIiiI1II1i11 ) :
     rmEdgeArray ( III1iII1I1ii , oOOo0 , I11IIIi , xll , yll , cellSize )
     IiI111111IIII . remove ( I11IIIi )
     IiI111111IIII . remove ( iIIiiI1II1i11 )
     if 65 - 65: o0000oOoOoO0o / OoOO0ooOOoo0O / Ooo00oOo00o
     if 92 - 92: O0 - o00O0oo . OoOO * o0000oOoOoO0o
 I1iI = 0
 for IiII1IiiIiI1 in range ( xsz ) :
  for iIiiiI1IiI1I1 in range ( ysz ) :
   if III1iII1I1ii [ IiII1IiiIiI1 , iIiiiI1IiI1I1 ] :
    convParX [ IiII1IiiIiI1 , iIiiiI1IiI1I1 , : ] = - 9999
    ii1ii11IIIiiI [ ( IiII1IiiIiI1 , iIiiiI1IiI1I1 , 'x' ) ] = [ o0oo [ iii ] , OOoO [ ( IiII1IiiIiI1 , iIiiiI1IiI1I1 ) ] , iiIii [ iii ] ]
    I1iI += OOoO [ ( IiII1IiiIiI1 , iIiiiI1IiI1I1 ) ]
   if oOOo0 [ IiII1IiiIiI1 , iIiiiI1IiI1I1 ] :
    convParY [ IiII1IiiIiI1 , iIiiiI1IiI1I1 , : ] = - 9999
    ii1ii11IIIiiI [ ( IiII1IiiIiI1 , iIiiiI1IiI1I1 , 'y' ) ] = [ o0oo [ iii ] , iII [ ( IiII1IiiIiI1 , iIiiiI1IiI1I1 ) ] , iiIii [ iii ] ]
    I1iI += iII [ ( IiII1IiiIiI1 , iIiiiI1IiI1I1 ) ]
    if 38 - 38: ii11ii1ii % Ooo00oOo00o + OOooOOo . i11iIiiIii
    if 53 - 53: i11iIiiIii * o00O0oo
 if csvEdgeFileName is not None :
  OooooO0oOOOO = open ( csvEdgeFileName , "w" )
  OooooO0oOOOO . write ( "id;wkt;level;length;wc\n" )
  o00oO0oOo00 = 0
  for IiII1IiiIiI1 in range ( xsz ) :
   for iIiiiI1IiI1I1 in range ( ysz ) :
    if III1iII1I1ii [ IiII1IiiIiI1 , iIiiiI1IiI1I1 ] :
     Ii11iII1 = IiII1IiiIiI1 * cellSize + xll
     Oo0O0O0ooO0O = iIiiiI1IiI1I1 * cellSize + yll
     IIIIii = Ii11iII1
     O0o0 = Oo0O0O0ooO0O + cellSize
     if 100 - 100: o00O0oo % OoOO
     OooooO0oOOOO . write ( "%i;LINESTRING(%f %f, %f %f);%f;%f;%f\n" % ( o00oO0oOo00 , Ii11iII1 , Oo0O0O0ooO0O , IIIIii , O0o0 , ii1ii11IIIiiI [ ( IiII1IiiIiI1 , iIiiiI1IiI1I1 , 'x' ) ] [ 0 ] , ii1ii11IIIiiI [ ( IiII1IiiIiI1 , iIiiiI1IiI1I1 , 'x' ) ] [ 1 ] , ii1ii11IIIiiI [ ( IiII1IiiIiI1 , iIiiiI1IiI1I1 , 'x' ) ] [ 2 ] ) )
     if 86 - 86: OOOo0 . O0 - OoooooooOO . Oo + o0000oOoOoO0o
     if 57 - 57: I1IiI . i1IIi . oO0o0ooO0 * i11iIiiIii + IIII . oO0o0ooO0
     o00oO0oOo00 += 1
     if 57 - 57: IIII
    if oOOo0 [ IiII1IiiIiI1 , iIiiiI1IiI1I1 ] :
     Ii11iII1 = IiII1IiiIiI1 * cellSize + xll
     Oo0O0O0ooO0O = iIiiiI1IiI1I1 * cellSize + yll
     IIIIii = Ii11iII1 + cellSize
     O0o0 = Oo0O0O0ooO0O
     if 32 - 32: o0000oOoOoO0o - OOOo0 % OoooooooOO . o00O0oo / oO0o0ooO0 + OoOoOO00
     OooooO0oOOOO . write ( "%i;LINESTRING(%f %f, %f %f);%f;%f;%f\n" % ( o00oO0oOo00 , Ii11iII1 , Oo0O0O0ooO0O , IIIIii , O0o0 , ii1ii11IIIiiI [ ( IiII1IiiIiI1 , iIiiiI1IiI1I1 , 'y' ) ] [ 0 ] , ii1ii11IIIiiI [ ( IiII1IiiIiI1 , iIiiiI1IiI1I1 , 'y' ) ] [ 1 ] , ii1ii11IIIiiI [ ( IiII1IiiIiI1 , iIiiiI1IiI1I1 , 'y' ) ] [ 2 ] ) )
     if 76 - 76: O0oO
     if 73 - 73: O0 * o00O0oo + o0000oOoOoO0o + O0oO
     o00oO0oOo00 += 1
  OooooO0oOOOO . close ( )
  if 40 - 40: II111iiii . Ooo00oOo00o * IIII + OoOO + OoOO
  if 9 - 9: OoOO0ooOOoo0O % OoooooooOO . ii11ii1ii % OoOO0ooOOoo0O
  if 32 - 32: i11iIiiIii
  if 31 - 31: iIii1I11I1II1 / Oo / OOooOOo
  if 41 - 41: OOOo0
  if 10 - 10: OOOo0 / OOOo0 / IIII . IIII
  if 98 - 98: OOOo0 / OoOoOO00 . O0 + Oo
  if 43 - 43: II111iiii . ii11ii1ii / OOooOOo
  if 20 - 20: OoOoOO00
  if 95 - 95: o00O0oo - OoOoOO00
  if 34 - 34: O0oO * OoOoOO00 . i1IIi * O0oO / O0oO
  if 30 - 30: OOooOOo + OOOo0 / OOOo0 % OOooOOo . OOooOOo
  if 55 - 55: O0oO - OoOO0ooOOoo0O + II111iiii + o00O0oo % o0000oOoOoO0o
  if 41 - 41: i1IIi - OoOO0ooOOoo0O - o0000oOoOoO0o
  if 8 - 8: Oo + IIII - I1IiI % OOOo0 % I1IiI * ii11ii1ii
  if 9 - 9: OOOo0 - i11iIiiIii - OoOO * o0000oOoOoO0o + O0oO
  if 44 - 44: II111iiii
  if 52 - 52: OOooOOo - OOOo0 + OOooOOo % I1IiI
  if 35 - 35: iIii1I11I1II1
  if 42 - 42: IIII . OoOoOO00 . i1IIi + Ooo00oOo00o + OoOO + OoOoOO00
  if 31 - 31: o00O0oo . OoOO - O0oO . OoooooooOO / OoooooooOO
  if 56 - 56: Oo / ii11ii1ii / i11iIiiIii + OoooooooOO - OOOo0 - OoOO0ooOOoo0O
  if 21 - 21: O0 % oO0o0ooO0 . OoOoOO00 / II111iiii + oO0o0ooO0
  if 53 - 53: ii11ii1ii - OoOoOO00 - ii11ii1ii * o00O0oo
  if 71 - 71: O0 - iIii1I11I1II1
  if 12 - 12: OoOO / I1IiI
  if 42 - 42: OOOo0
  if 19 - 19: ii11ii1ii % OOooOOo * iIii1I11I1II1 + OoOoOO00
  if 46 - 46: OOOo0
  if 1 - 1: o00O0oo
  if 97 - 97: OoOO + o00O0oo + O0 + i11iIiiIii
  if 77 - 77: I1IiI / OoooooooOO
  if 46 - 46: I1IiI % iIii1I11I1II1 . o00O0oo % o00O0oo + i11iIiiIii
 print "Total polyline length=" , I1I111
 print "Total crest length added=" , I1iI
 if 72 - 72: iIii1I11I1II1 * o0000oOoOoO0o % O0oO / Oo
 return ii1ii11IIIiiI , ooOOO00Ooo , IiIIIi1iIi # dd678faae9ac167bc83abf78e5cb2f3f0688d3a3
