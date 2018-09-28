import fileIO
import numpy
from bresenhamalgorithm import bresenham
import scipy . stats
import scipy . optimize
import random
import os
import sys
import time
import ctypes
from numpy . ctypeslib import ndpointer
import tempfile
from subprocess import call
import reservoirs
if 64 - 64: i11iIiiIii
if 65 - 65: O0 / iIii1I11I1II1 % OoooooooOO - i1IIi
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
arrayType = numpy . float64
if 73 - 73: II111iiii
def setPrecision32 ( ) :
 global arrayType
 arrayType = numpy . float32
 if 22 - 22: I1IiiI * Oo0Ooo / OoO0O00 . OoOoOO00 . o0oOOo0O0Ooo / I1ii11iIi11i
def setPrecision64 ( ) :
 global arrayType
 arrayType = numpy . float64
 if 48 - 48: oO0o / OOooOOo / I11i / Ii1I
def getPrecision ( ) :
 return arrayType
 if 48 - 48: iII111i % IiII + I1Ii111 / ooOoO0o * Ii1I
 if 46 - 46: ooOoO0o * I11i - OoooooooOO
 if 30 - 30: o0oOOo0O0Ooo - O0 % o0oOOo0O0Ooo - OoooooooOO * O0 * OoooooooOO
 if 60 - 60: iIii1I11I1II1 / i1IIi * oO0o - I1ii11iIi11i + o0oOOo0O0Ooo
 if 94 - 94: i1IIi % Oo0Ooo
def replaceArrayVals ( arr , val1 , val2 ) :
 o0oO0 , oo00 = arr . shape
 if 88 - 88: iII111i . oO0o % ooOoO0o
 for ooO0oooOoO0 in range ( oo00 ) :
  II11i = numpy . where ( arr [ : , ooO0oooOoO0 ] == val1 )
  arr [ : , ooO0oooOoO0 ] [ II11i ] = val2
  if 43 - 43: Ii1I . oO0o
 return
 if 27 - 27: OoO0O00 - O0 . I1Ii111 * iII111i - I1ii11iIi11i
 if 15 - 15: I1IiiI
 if 90 - 90: IiII * i1IIi / Ii1I . OoO0O00 * oO0o
Ct = 0
previousRainfall = 0.
rainfallProfileSummer = [ 0.024 , 0.036 , 0.054 , 0.087 , 0.154 , 0.291 , 0.154 , 0.087 , 0.054 , 0.036 , 0.024 ]
if 16 - 16: ooOoO0o * IiII % I11i . I1Ii111 / IiII % iII111i
rainfallProfileWinter = [ 0.022 , 0.037 , 0.061 , 0.101 , 0.165 , 0.229 , 0.165 , 0.101 , 0.061 , 0.037 , 0.022 ]
if 27 - 27: IiII . i1IIi * OoOoOO00 % Ii1I / i1IIi
if 3 - 3: IiII / ooOoO0o
def refhRunoff ( t , stormDuration , rainfallDepth , bfihost , propwet , summerProfile = False ) :
 global Ct , previousRainfall , rainfallProfile , Cmax
 if 28 - 28: ooOoO0o + I1Ii111 - ooOoO0o . OoooooooOO
 if 97 - 97: OoO0O00 . I11i
 if t <= 0. :
  Cmax = 596.7 * ( bfihost ** 0.95 ) * ( propwet ** - 0.24 )
  Ct = ( Cmax / 2 ) * ( 0.90 - 0.82 * bfihost - 0.43 * propwet )
  if summerProfile is not None and summerProfile :
   rainfallProfile = rainfallProfileSummer
  else :
   rainfallProfile = rainfallProfileWinter
  previousRainfall = 0.
  if 32 - 32: Oo0Ooo - II111iiii - i11iIiiIii % I1Ii111
 if t < 0 or t > stormDuration :
  return 0
  if 54 - 54: OOooOOo % O0 + I1IiiI - iII111i / I11i
 iIiiI1 = int ( 11. * t / stormDuration )
 if 68 - 68: I1IiiI - i11iIiiIii - OoO0O00 / OOooOOo - OoO0O00 + i1IIi
 if iIiiI1 >= 1 and t < stormDuration :
  IiiIII111ii = sum ( rainfallProfile [ : iIiiI1 ] )
 elif t >= stormDuration :
  IiiIII111ii = 1.0
 else :
  IiiIII111ii = 0.
  if 3 - 3: iII111i + O0
 if t < stormDuration :
  IiiIII111ii += ( t - iIiiI1 * stormDuration / 11. ) / ( stormDuration / 11. ) * rainfallProfile [ iIiiI1 ]
  if 42 - 42: OOooOOo / i1IIi + i11iIiiIii - Ii1I
  if 78 - 78: OoO0O00
 IiiIII111ii *= rainfallDepth
 if 18 - 18: O0 - iII111i / iII111i + ooOoO0o % ooOoO0o - IiII
 IiiIII111ii = min ( IiiIII111ii , rainfallDepth )
 if 62 - 62: iII111i - IiII - OoOoOO00 % i1IIi / oO0o
 OoooooOoo = IiiIII111ii - previousRainfall
 if 70 - 70: OoO0O00 . OoO0O00 - OoO0O00 / I1ii11iIi11i * OOooOOo
 if 86 - 86: i11iIiiIii + Ii1I + ooOoO0o * I11i + o0oOOo0O0Ooo
 if 61 - 61: OoO0O00 / i11iIiiIii
 IiIiIi = OoooooOoo * ( Ct / Cmax + 0.5 * OoooooOoo / Cmax )
 if 40 - 40: oO0o . OoOoOO00 . Oo0Ooo . i1IIi
 Ct += OoooooOoo
 previousRainfall = IiiIII111ii
 if 33 - 33: Ii1I + II111iiii % i11iIiiIii . ooOoO0o - I1IiiI
 if 66 - 66: Ii1I - OoooooooOO * OoooooooOO . OOooOOo . I1ii11iIi11i
 return IiIiIi
 if 22 - 22: OoooooooOO % I11i - iII111i . iIii1I11I1II1 * i11iIiiIii
 if 32 - 32: Oo0Ooo * O0 % oO0o % Ii1I . IiII
 if 61 - 61: ooOoO0o
 if 79 - 79: Oo0Ooo + I1IiiI - iII111i
rainfallProfileSummer = [ 0.024 , 0.036 , 0.054 , 0.087 , 0.154 , 0.291 , 0.154 , 0.087 , 0.054 , 0.036 , 0.024 ]
if 83 - 83: ooOoO0o
rainfallProfileWinter = [ 0.022 , 0.037 , 0.061 , 0.101 , 0.165 , 0.229 , 0.165 , 0.101 , 0.061 , 0.037 , 0.022 ]
if 64 - 64: OoO0O00 % ooOoO0o % iII111i / OoOoOO00 - OoO0O00
if 74 - 74: iII111i * O0
def pRunoff ( t , stormDuration , rainfallDepth , pr , summerProfile = False ) :
 if 89 - 89: oO0o + Oo0Ooo
 global previousRainfall
 if 3 - 3: i1IIi / I1IiiI % I11i * i11iIiiIii / O0 * I11i
 if t <= 0 :
  previousRainfall = 0.
  if 49 - 49: oO0o % Ii1I + i1IIi . I1IiiI % I1ii11iIi11i
 if summerProfile is not None and summerProfile :
  I1i1iii = rainfallProfileSummer
 else :
  I1i1iii = rainfallProfileWinter
  if 20 - 20: o0oOOo0O0Ooo
 if t < 0 or t > stormDuration :
  return 0
  if 77 - 77: OoOoOO00 / I11i
 iIiiI1 = int ( 11. * t / stormDuration )
 if 98 - 98: iIii1I11I1II1 / i1IIi / i11iIiiIii / o0oOOo0O0Ooo
 if iIiiI1 >= 1 and t < stormDuration :
  IiiIII111ii = sum ( I1i1iii [ : iIiiI1 ] )
 elif t >= stormDuration :
  IiiIII111ii = 1.0
 else :
  IiiIII111ii = 0.
  if 28 - 28: OOooOOo - IiII . IiII + OoOoOO00 - OoooooooOO + O0
 if t < stormDuration :
  IiiIII111ii += ( t - iIiiI1 * stormDuration / 11. ) / ( stormDuration / 11. ) * I1i1iii [ iIiiI1 ]
  if 95 - 95: OoO0O00 % oO0o . O0
  if 15 - 15: ooOoO0o / Ii1I . Ii1I - i1IIi
 IiiIII111ii *= rainfallDepth
 if 53 - 53: IiII + I1IiiI * oO0o
 IiiIII111ii = min ( IiiIII111ii , rainfallDepth )
 if 61 - 61: i1IIi * OOooOOo / OoooooooOO . i11iIiiIii . OoOoOO00
 OoooooOoo = IiiIII111ii - previousRainfall
 if 60 - 60: I11i / I11i
 IiIiIi = OoooooOoo * pr / 100.
 if 46 - 46: Ii1I * OOooOOo - OoO0O00 * oO0o - I1Ii111
 previousRainfall = IiiIII111ii
 if 83 - 83: OoooooooOO
 return IiIiIi
 if 31 - 31: II111iiii - OOooOOo . I1Ii111 % OoOoOO00 - O0
 if 4 - 4: II111iiii / ooOoO0o . iII111i
 if 58 - 58: OOooOOo * i11iIiiIii / OoOoOO00 % I1Ii111 - I1ii11iIi11i / oO0o
def writeGridCSV ( xll , yll , cellSize , xsz , ysz , sp , fileName ) :
 if 50 - 50: I1IiiI
 if 34 - 34: I1IiiI * II111iiii % iII111i * OoOoOO00 - I1IiiI
 II1III = open ( fileName + "t" , "w" )
 II1III . write ( "\"Integer\",\"String\",\"Real\"\n" )
 II1III . close ( )
 if 19 - 19: oO0o % i1IIi % o0oOOo0O0Ooo
 oo0OooOOo0 = open ( fileName , "w" )
 if 92 - 92: iII111i . I11i + o0oOOo0O0Ooo
 id = 0
 if 28 - 28: i1IIi * Oo0Ooo - o0oOOo0O0Ooo * IiII * Ii1I / OoO0O00
 if 94 - 94: II111iiii % I1ii11iIi11i / OoOoOO00 * iIii1I11I1II1
 oo0OooOOo0 . write ( "id;wkt;elevation\n" )
 if 54 - 54: o0oOOo0O0Ooo - I1IiiI + OoooooooOO
 for O0o0 in range ( xsz ) :
  for OO00Oo in range ( ysz ) :
   O0OOO0OOoO0O = xll + O0o0 * cellSize
   O00Oo000ooO0 = O0OOO0OOoO0O + cellSize
   OoO0O00IIiII = yll + OO00Oo * cellSize
   o0 = OoO0O00IIiII + cellSize
   if 62 - 62: iIii1I11I1II1 * OoOoOO00
   i1 = "POLYGON (("
   i1 += "%0.3f %0.3f," % ( O0OOO0OOoO0O , OoO0O00IIiII )
   i1 += "%0.3f %0.3f," % ( O00Oo000ooO0 , OoO0O00IIiII )
   i1 += "%0.3f %0.3f," % ( O00Oo000ooO0 , o0 )
   i1 += "%0.3f %0.3f," % ( O0OOO0OOoO0O , o0 )
   i1 += "%0.3f %0.3f" % ( O0OOO0OOoO0O , OoO0O00IIiII )
   i1 += "))"
   if 91 - 91: OoO0O00 . I1ii11iIi11i + OoO0O00 - iII111i / OoooooooOO
   oo0OooOOo0 . write ( "%i;%s;%f\n" % ( id , i1 , sp [ O0o0 , OO00Oo , 0 ] ) )
   if 39 - 39: I1ii11iIi11i / ooOoO0o - II111iiii
   id += 1
 oo0OooOOo0 . close ( )
 if 98 - 98: I1ii11iIi11i / I11i % oO0o . OoOoOO00
 return
 if 91 - 91: oO0o % Oo0Ooo
 if 64 - 64: I11i % iII111i - I1Ii111 - oO0o
def i1ii1iiI ( xyList , raster , cellSize , strict = False ) :
 O0o0O00Oo0o0 = [ ]
 if 87 - 87: ooOoO0o * Oo0Ooo % i11iIiiIii % OoOoOO00 - OOooOOo
 if 68 - 68: I1Ii111 % i1IIi . IiII . I1ii11iIi11i
 if 92 - 92: iII111i . I1Ii111
 if 31 - 31: I1Ii111 . OoOoOO00 / O0
 o000O0o = len ( xyList )
 iI1iII1 = 0.
 for O0o0 in range ( o000O0o - 1 ) :
  O0OOO0OOoO0O = xyList [ O0o0 ] [ 0 ]
  OoO0O00IIiII = xyList [ O0o0 ] [ 1 ]
  O00Oo000ooO0 = xyList [ O0o0 + 1 ] [ 0 ]
  o0 = xyList [ O0o0 + 1 ] [ 1 ]
  if 86 - 86: OOooOOo
  if O0OOO0OOoO0O == O00Oo000ooO0 :
   iI1iII1 += ( o0 - OoO0O00IIiII )
   O0o0O00Oo0o0 = list ( raster [ O0OOO0OOoO0O , OoO0O00IIiII : o0 ] )
   O0o0O00Oo0o0 . append ( raster [ O00Oo000ooO0 , o0 ] )
   if 55 - 55: Oo0Ooo + iIii1I11I1II1 / OoOoOO00 * oO0o - i11iIiiIii - Ii1I
  elif OoO0O00IIiII == o0 :
   iI1iII1 += ( O00Oo000ooO0 - O0OOO0OOoO0O )
   O0o0O00Oo0o0 = list ( raster [ O0OOO0OOoO0O : O00Oo000ooO0 , OoO0O00IIiII ] )
   O0o0O00Oo0o0 . append ( raster [ O00Oo000ooO0 , o0 ] )
   if 25 - 25: I1ii11iIi11i
  else :
   iI1iII1 += numpy . sqrt ( ( O00Oo000ooO0 - O0OOO0OOoO0O ) ** 2 + ( o0 - OoO0O00IIiII ) ** 2 )
   if 7 - 7: i1IIi / I1IiiI * I1Ii111 . IiII . iIii1I11I1II1
   if O0OOO0OOoO0O == O00Oo000ooO0 and OoO0O00IIiII == o0 :
    O0o0O00Oo0o0 . append ( raster [ O0OOO0OOoO0O , OoO0O00IIiII ] )
   else :
    if O0o0 == ( o000O0o - 2 ) :
     for iIii , ooo0O in bresenham ( [ O0OOO0OOoO0O , OoO0O00IIiII ] , [ O00Oo000ooO0 , o0 ] ) . path [ : ] :
      O0o0O00Oo0o0 . append ( raster [ iIii , ooo0O ] )
    else :
     for iIii , ooo0O in bresenham ( [ O0OOO0OOoO0O , OoO0O00IIiII ] , [ O00Oo000ooO0 , o0 ] ) . path [ : - 1 ] :
      O0o0O00Oo0o0 . append ( raster [ iIii , ooo0O ] )
      if 75 - 75: o0oOOo0O0Ooo % o0oOOo0O0Ooo . I1Ii111
      if 5 - 5: o0oOOo0O0Ooo * ooOoO0o + OoOoOO00 . OOooOOo + OoOoOO00
 oO = cellSize * iI1iII1 / len ( O0o0O00Oo0o0 )
 if 7 - 7: o0oOOo0O0Ooo - I1IiiI
 if 100 - 100: oO0o + I11i . OOooOOo * Ii1I
 if strict :
  O0o0O00Oo0o0 = O0o0O00Oo0o0 [ 0 : - 1 ]
  if 73 - 73: i1IIi + I1IiiI
 return oO , O0o0O00Oo0o0
 if 46 - 46: OoO0O00 . Oo0Ooo - OoooooooOO
 if 93 - 93: iII111i
def i1IIIiiII1 ( x , a , b ) :
 return a * x * x + b * x
 if 87 - 87: oO0o * I1ii11iIi11i + OOooOOo / iIii1I11I1II1 / iII111i
def quadFit ( xList , yList ) :
 I1111IIi = scipy . optimize . curve_fit ( i1IIIiiII1 , numpy . array ( xList ) , numpy . array ( yList ) )
 if 93 - 93: OoooooooOO / I1IiiI % i11iIiiIii + I1ii11iIi11i * OoO0O00
 if 15 - 15: I11i . OoO0O00 / Oo0Ooo + I11i
 return I1111IIi [ 0 ]
 if 78 - 78: O0 . oO0o . II111iiii % OOooOOo
 if 49 - 49: Ii1I / OoO0O00 . II111iiii
 if 68 - 68: i11iIiiIii % I1ii11iIi11i + i11iIiiIii
def calcStorageParametersChannel ( rect , z , xll , yll , dx , cellSize , cag , cagXll , cagYll , cagDx , chanExp , chanMult , chanAR , chanMaxD ) :
 if 31 - 31: II111iiii . I1IiiI
 if 1 - 1: Oo0Ooo / o0oOOo0O0Ooo % iII111i * IiII . i11iIiiIii
 III1Iiii1I11 = int ( ( rect [ 0 ] [ 0 ] - xll ) / dx )
 IIII = int ( ( rect [ 1 ] [ 0 ] - xll ) / dx )
 if 32 - 32: OoooooooOO / iIii1I11I1II1 - o0oOOo0O0Ooo
 o00oooO0Oo = int ( ( rect [ 0 ] [ 1 ] - yll ) / dx )
 o0O0OOO0Ooo = int ( ( rect [ 1 ] [ 1 ] - yll ) / dx )
 if 45 - 45: O0 / o0oOOo0O0Ooo
 i1IIIII11I1IiI = z [ III1Iiii1I11 : ( IIII + 1 ) , o00oooO0Oo : ( o0O0OOO0Ooo + 1 ) ]
 if 16 - 16: iIii1I11I1II1
 i1IIIII11I1IiI = i1IIIII11I1IiI . reshape ( ( - 1 ) )
 if 90 - 90: o0oOOo0O0Ooo % i1IIi / OoO0O00
 if 44 - 44: Oo0Ooo . OoO0O00 / I1ii11iIi11i + Ii1I
 i1IIIII11I1IiI . sort ( )
 if 65 - 65: O0
 if 68 - 68: OOooOOo % I1Ii111
 if 88 - 88: iIii1I11I1II1 - ooOoO0o + OOooOOo
 III1Iiii1I11 = int ( ( rect [ 0 ] [ 0 ] - cagXll ) / cagDx )
 IIII = int ( ( rect [ 1 ] [ 0 ] - cagXll ) / cagDx )
 if 40 - 40: I1IiiI * Ii1I + OOooOOo % iII111i
 o00oooO0Oo = int ( ( rect [ 0 ] [ 1 ] - cagYll ) / cagDx )
 o0O0OOO0Ooo = int ( ( rect [ 1 ] [ 1 ] - cagYll ) / cagDx )
 if 74 - 74: oO0o - Oo0Ooo + OoooooooOO + I1Ii111 / OoOoOO00
 if 23 - 23: O0
 o00oO0oOo00 , oO0oOo0 = cag . shape
 if 45 - 45: iII111i / iII111i + I1Ii111 + ooOoO0o
 if III1Iiii1I11 > 0 and III1Iiii1I11 < o00oO0oOo00 and IIII > 0 and IIII < o00oO0oOo00 and o00oooO0Oo > 0 and o00oooO0Oo < oO0oOo0 and o0O0OOO0Ooo > 0 and o0O0OOO0Ooo < oO0oOo0 :
  if 47 - 47: o0oOOo0O0Ooo + ooOoO0o
  if 82 - 82: II111iiii . IiII - iIii1I11I1II1 - IiII * II111iiii
  if 77 - 77: iIii1I11I1II1 * OoO0O00
  if 95 - 95: I1IiiI + i11iIiiIii
  I1Ii = cag [ III1Iiii1I11 : IIII , o00oooO0Oo : o0O0OOO0Ooo ]
  if 94 - 94: Ii1I - II111iiii . OOooOOo % I11i . i11iIiiIii + O0
  I1IiiiiI = I1Ii . max ( )
 else :
  I1IiiiiI = 0.
  if 80 - 80: I1Ii111 . i11iIiiIii - o0oOOo0O0Ooo
  if 25 - 25: OoO0O00
  if 62 - 62: OOooOOo + O0
  if 98 - 98: o0oOOo0O0Ooo
 OOOO0oo0 = chanMult * ( I1IiiiiI ** chanExp )
 I11iiI1i1 = min ( OOOO0oo0 / chanAR , chanMaxD )
 if 47 - 47: iII111i - Ii1I . II111iiii + OoooooooOO . i11iIiiIii
 if 94 - 94: o0oOOo0O0Ooo * Ii1I / Oo0Ooo / Ii1I
 oO0 = i1IIIII11I1IiI . min ( )
 O0OO0O = OOOO0oo0 * I11iiI1i1 * cellSize
 OO = OOOO0oo0 * cellSize
 OoOoO = int ( OO / ( dx * dx ) )
 if 43 - 43: i11iIiiIii + Oo0Ooo * II111iiii * I1Ii111 * O0
 if 64 - 64: OOooOOo % iIii1I11I1II1 * oO0o
 if OoOoO > 0 :
  i1IIIII11I1IiI = i1IIIII11I1IiI [ OoOoO : ]
  i1IIIII11I1IiI = list ( i1IIIII11I1IiI ) + [ oO0 - I11iiI1i1 ] * OoOoO
  if 79 - 79: O0
 oOO00O = min ( i1IIIII11I1IiI )
 OOOoo0OO = max ( i1IIIII11I1IiI )
 if 57 - 57: OoO0O00 / ooOoO0o
 if 29 - 29: iIii1I11I1II1 + OoOoOO00 * OoO0O00 * OOooOOo . I1IiiI * I1IiiI
 if 7 - 7: IiII * I1Ii111 % Ii1I - o0oOOo0O0Ooo
 if 13 - 13: Ii1I . i11iIiiIii
 if 56 - 56: I1ii11iIi11i % O0 - I1IiiI
 O00o0OO0 = [ ]
 if 35 - 35: oO0o % ooOoO0o / I1Ii111 + iIii1I11I1II1 . OoooooooOO . I1IiiI
 if 71 - 71: IiII * II111iiii * oO0o
 if OOOoo0OO > ( oO0 + 5 ) :
  oOOo0 = [ oOO00O , oO0 , oO0 + 1 , oO0 + 5 , OOOoo0OO ]
 else :
  oOOo0 = [ oOO00O , oO0 , oO0 + 1 , oO0 + 5 , oO0 + 5 ]
  OOOoo0OO = oO0 + 5
  if 16 - 16: oO0o % I1ii11iIi11i * i11iIiiIii % i11iIiiIii
 O0OOOOo0O = [ ]
 if 81 - 81: O0 / OoO0O00 . i1IIi + I1IiiI - I11i
 for OoOOoOooooOOo in oOOo0 :
  oOo0O = sum ( [ max ( 0 , OoOOoOooooOOo - oo0O0 ) for oo0O0 in i1IIIII11I1IiI ] )
  O0OOOOo0O . append ( oOo0O / len ( i1IIIII11I1IiI ) )
  if 22 - 22: OoOoOO00 . OOooOOo * OoOoOO00
 return [ oOO00O , oO0 , OOOoo0OO ] + O0OOOOo0O [ 1 : ]
 if 54 - 54: IiII + Ii1I % OoO0O00 + OoooooooOO - O0 - o0oOOo0O0Ooo
 if 77 - 77: OOooOOo * iIii1I11I1II1
 if 98 - 98: I1IiiI % Ii1I * OoooooooOO
 if 51 - 51: iIii1I11I1II1 . OoOoOO00 / oO0o + o0oOOo0O0Ooo
def calcStorageParameters ( rect , z , xll , yll , dx , plotName = None , csvOutput = None ) :
 if 33 - 33: ooOoO0o . II111iiii % iII111i + o0oOOo0O0Ooo
 III1Iiii1I11 = int ( ( rect [ 0 ] [ 0 ] - xll ) / dx )
 IIII = int ( ( rect [ 1 ] [ 0 ] - xll ) / dx )
 if 71 - 71: Oo0Ooo % OOooOOo
 o00oooO0Oo = int ( ( rect [ 0 ] [ 1 ] - yll ) / dx )
 o0O0OOO0Ooo = int ( ( rect [ 1 ] [ 1 ] - yll ) / dx )
 if 98 - 98: I11i % i11iIiiIii % ooOoO0o + Ii1I
 if 78 - 78: I1ii11iIi11i % oO0o / iII111i - iIii1I11I1II1
 if 69 - 69: I1Ii111
 i1IIIII11I1IiI = z [ III1Iiii1I11 : ( IIII + 1 ) , o00oooO0Oo : ( o0O0OOO0Ooo + 1 ) ]
 if 11 - 11: I1IiiI
 i1IIIII11I1IiI = i1IIIII11I1IiI . reshape ( ( - 1 ) )
 if 16 - 16: Ii1I + IiII * O0 % i1IIi . I1IiiI
 oOO00O = i1IIIII11I1IiI . min ( )
 if 67 - 67: OoooooooOO / I1IiiI * Ii1I + I11i
 OOOoo0OO = i1IIIII11I1IiI . max ( )
 if 65 - 65: OoooooooOO - I1ii11iIi11i / ooOoO0o / II111iiii / i1IIi
 if numpy . isnan ( oOO00O ) or numpy . isnan ( OOOoo0OO ) :
  print i1IIIII11I1IiI
  if 71 - 71: I1Ii111 + Ii1I
 O00o0OO0 = [ ]
 if 28 - 28: OOooOOo
 if 38 - 38: ooOoO0o % II111iiii % I11i / OoO0O00 + OoOoOO00 / i1IIi
 if OOOoo0OO > ( oOO00O + 5 ) and OOOoo0OO > ( oOO00O + 1 ) :
  oOOo0 = [ oOO00O , oOO00O + 1 , oOO00O + 5 , OOOoo0OO ]
 else :
  oOOo0 = [ oOO00O , oOO00O + 1 , oOO00O + 5 , oOO00O + 5 ]
  OOOoo0OO = oOO00O + 5
  if 54 - 54: iIii1I11I1II1 % I1ii11iIi11i - OOooOOo / oO0o - OoO0O00 . I11i
 O0OOOOo0O = [ ]
 if 11 - 11: I1ii11iIi11i . OoO0O00 * IiII * OoooooooOO + ooOoO0o
 for OoOOoOooooOOo in oOOo0 :
  oOo0O = sum ( [ max ( 0 , OoOOoOooooOOo - oo0O0 ) for oo0O0 in i1IIIII11I1IiI ] )
  O0OOOOo0O . append ( oOo0O / len ( i1IIIII11I1IiI ) )
  if 33 - 33: O0 * o0oOOo0O0Ooo - I1Ii111 % I1Ii111
  if 18 - 18: I1Ii111 / Oo0Ooo * I1Ii111 + I1Ii111 * i11iIiiIii * I1ii11iIi11i
  if 11 - 11: ooOoO0o / OoOoOO00 - IiII * OoooooooOO + OoooooooOO . OoOoOO00
 return oOO00O , OOOoo0OO , O0OOOOo0O [ 1 ] , O0OOOOo0O [ 2 ] , O0OOOOo0O [ 3 ]
 if 26 - 26: Ii1I % I1ii11iIi11i
 if 76 - 76: IiII * iII111i
 if 52 - 52: OOooOOo
 if 19 - 19: I1IiiI
 if 25 - 25: Ii1I / ooOoO0o
def conveyanceParameters ( profile , dx , n , plotName = None , csvOutput = None , nList = None ) :
 II = 10.
 if 70 - 70: iIii1I11I1II1
 i11ii1iI = min ( profile )
 i1I = max ( profile )
 if 42 - 42: o0oOOo0O0Ooo + i1IIi - Ii1I / IiII
 iiIiIIIiiI = 5
 iiI1IIIi = [ ]
 for O0o0 in range ( iiIiIIIiiI ) :
  II11IiIi11 = ( 10 ** ( - 2. + O0o0 * 2. / ( iiIiIIIiiI ) ) )
  iiI1IIIi . append ( i11ii1iI + II11IiIi11 * II )
  if 7 - 7: OoO0O00 . Ii1I % oO0o * ooOoO0o + IiII + I1Ii111
 IIIIiII1i = [ ]
 i1II1 = [ ]
 if 25 - 25: I1Ii111 / iIii1I11I1II1 % iII111i
 if nList is None :
  nList = [ n ] * len ( iiI1IIIi )
  if 42 - 42: i11iIiiIii * iIii1I11I1II1 / I1ii11iIi11i . i11iIiiIii % I11i
 for OoOOoOooooOOo in iiI1IIIi :
  i1iI = 0.
  if 29 - 29: I1IiiI % OOooOOo - I1IiiI / OOooOOo . i1IIi
  i1II1 . append ( dx * len ( profile ) * ( ( OoOOoOooooOOo - i11ii1iI ) ** 1.666667 ) / n )
  if 31 - 31: I1Ii111
  if 88 - 88: OoO0O00 - ooOoO0o + OOooOOo * I1IiiI % iIii1I11I1II1 + Oo0Ooo
  if 76 - 76: I1IiiI * iII111i % I1Ii111
  if 57 - 57: iIii1I11I1II1 - i1IIi / I1Ii111 - O0 * OoooooooOO % II111iiii
  if 68 - 68: OoooooooOO * I11i % OoOoOO00 - IiII
  if 34 - 34: I1Ii111 . iIii1I11I1II1 * OoOoOO00 * oO0o / I1Ii111 / I1ii11iIi11i
  for O0o0 , oOoOOo0O in enumerate ( profile ) :
   i1iI += dx * ( max ( 0 , OoOOoOooooOOo - oOoOOo0O ) ** 1.6667 ) / nList [ O0o0 ]
   if 84 - 84: OoO0O00 + i1IIi - II111iiii . I1ii11iIi11i * OoooooooOO + I1IiiI
  IIIIiII1i . append ( i1iI )
  if 38 - 38: OOooOOo + II111iiii % ooOoO0o % OoOoOO00 - Ii1I / OoooooooOO
  if 73 - 73: o0oOOo0O0Ooo * O0 - i11iIiiIii
 i1II1 = [ numpy . log ( O0O0o0oOOO ) for O0O0o0oOOO in i1II1 ]
 IIIIiII1i = [ numpy . log ( O0O0o0oOOO ) for O0O0o0oOOO in IIIIiII1i ]
 if 67 - 67: OoOoOO00 + I1ii11iIi11i . o0oOOo0O0Ooo . II111iiii
 I1111IIi = scipy . stats . linregress ( i1II1 , IIIIiII1i )
 o000ooooO0o = I1111IIi [ 0 ]
 iI1i11 = I1111IIi [ 1 ]
 if 66 - 66: O0 % I1ii11iIi11i + i11iIiiIii . OoOoOO00 / Ii1I + I1ii11iIi11i
 if 86 - 86: o0oOOo0O0Ooo
 if plotName is not None :
  plt . figure ( 1 , facecolor = 'w' , edgecolor = 'w' )
  if 5 - 5: IiII * OoOoOO00
  IIIIiII1i = [ numpy . exp ( O0O0o0oOOO ) for O0O0o0oOOO in IIIIiII1i ]
  try :
   plt . semilogx ( IIIIiII1i , iiI1IIIi , '+' , color = 'k' )
  except :
   print IIIIiII1i , iiI1IIIi
   assert False
   if 5 - 5: I1Ii111
  O0I11Iiii1I = [ numpy . exp ( o000ooooO0o * oo00O0oO0O0 + iI1i11 ) for oo00O0oO0O0 in i1II1 ]
  if 96 - 96: i11iIiiIii % ooOoO0o / OoOoOO00
  plt . semilogx ( O0I11Iiii1I , iiI1IIIi , '-' , color = 'k' )
  if 36 - 36: OOooOOo + O0 - Ii1I - O0 % I11i . oO0o
  if csvOutput :
   ooo = open ( plotName + '.csv' , "w" )
   ooo . write ( "WL,K,Kapprox\n" )
   for O0o0 , OoOOoOooooOOo in enumerate ( iiI1IIIi ) :
    ooo . write ( "%f,%f,%f\n" % ( OoOOoOooooOOo , IIIIiII1i [ O0o0 ] , O0I11Iiii1I [ O0o0 ] ) )
   ooo . close ( )
   if 36 - 36: OoooooooOO . OoO0O00
  oOIIiIi = 0.
  OOoOooOoOOOoo = 0
  for O0o0 , i1iI in enumerate ( IIIIiII1i ) :
   OoOOoOooooOOo = iiI1IIIi [ O0o0 ]
   if 25 - 25: OoooooooOO - I1IiiI . I1IiiI * oO0o
   if 81 - 81: iII111i + IiII
   if 98 - 98: I1IiiI
   if 95 - 95: ooOoO0o / ooOoO0o
   IIiI1Ii = numpy . exp ( ( numpy . log ( i1iI ) - iI1i11 ) / o000ooooO0o )
   if 57 - 57: OOooOOo - ooOoO0o - I11i + OoO0O00
   I1IIIiI11i1 = i11ii1iI + ( IIiI1Ii * n / ( dx * len ( profile ) ) ) ** 0.6
   if 48 - 48: I1Ii111 - o0oOOo0O0Ooo % Ii1I
   oOIIiIi += ( I1IIIiI11i1 - OoOOoOooooOOo ) ** 2
   if 36 - 36: oO0o - Ii1I . Oo0Ooo - i11iIiiIii - OOooOOo * Oo0Ooo
   OOoOooOoOOOoo += 1
   if 76 - 76: i11iIiiIii + o0oOOo0O0Ooo / I1ii11iIi11i - OoO0O00 - Ii1I + I1ii11iIi11i
  oOIIiIi = numpy . sqrt ( oOIIiIi / OOoOooOoOOOoo )
  if 51 - 51: iIii1I11I1II1 . ooOoO0o + iIii1I11I1II1
  plt . title ( plotName + " RMS Error=%fm" % oOIIiIi )
  plt . xlabel ( 'Conveyance (m3s-1)' )
  plt . ylabel ( 'WL (m)' )
  if 95 - 95: I1IiiI
  pylab . savefig ( plotName + '.png' , bbox_inches = 0 )
  plt . close ( 'all' )
  if 46 - 46: OoOoOO00 + OoO0O00
  if 70 - 70: iII111i / iIii1I11I1II1
 return i11ii1iI , o000ooooO0o , iI1i11
 if 85 - 85: OoooooooOO % i1IIi * OoooooooOO / I1ii11iIi11i
 if 96 - 96: OoooooooOO + oO0o
def gridFlowSetupTiled ( dtmFileName , xll , yll , cellSize , xsz , ysz , nChan , nFP ,
 nFileName = None ,
 plotNamePrefix = None , outputPrefix = None ,
 ndv = None , ndr = None , conveyanceFunc = None , storageFunc = None ) :
 if 44 - 44: oO0o
 if plotNamePrefix is None :
  plotNamePrefix = ""
  if 20 - 20: I11i + Ii1I / O0 % iIii1I11I1II1
 if outputPrefix is None :
  outputPrefix = ""
  if 88 - 88: OoOoOO00 / II111iiii
 OOOOO0O00 , Iii , iIIiIiI1I1 , ooO , ii , OO0O0Ooo = fileIO . readScalarGridObj ( dtmFileName )
 if 77 - 77: o0oOOo0O0Ooo / OoooooooOO
 if nFileName is not None :
  IIii11I1i1I , o0o0OO0O00o , O0Oooo , iiIi1i , I1i11111i1i11 , OOoOOO0 = fileIO . readScalarGridObj ( nFileName )
  if 10 - 10: I1Ii111 / ooOoO0o + i11iIiiIii / Ii1I
  if 74 - 74: OOooOOo + O0 + i1IIi - i1IIi + II111iiii
 ndv = arrayType ( ndv )
 if 83 - 83: I1ii11iIi11i - I1IiiI + OOooOOo
 iIi1Ii1i1iI = numpy . zeros ( ( xsz + 1 , ysz , 7 ) , dtype = arrayType ) - 9999.
 IIiI1 = numpy . zeros ( ( xsz , ysz + 1 , 7 ) , dtype = arrayType ) - 9999.
 if 17 - 17: OOooOOo / OOooOOo / I11i
 ii1 = numpy . zeros ( ( xsz , ysz , 5 ) , dtype = arrayType ) - 9999.
 if 1 - 1: ooOoO0o % iIii1I11I1II1 + Oo0Ooo . iIii1I11I1II1 % I1IiiI
 o0o0oOoOO0O = None
 if 16 - 16: IiII % iIii1I11I1II1 . Ii1I
 oooooOOO000Oo = numpy . zeros ( 7 , dtype = arrayType )
 Ooo00OoOOO = numpy . zeros ( 5 , dtype = arrayType )
 if 98 - 98: iIii1I11I1II1 * I1ii11iIi11i * OOooOOo + ooOoO0o % i11iIiiIii % O0
 i1OO0oOOoo = 0
 if 52 - 52: o0oOOo0O0Ooo % Oo0Ooo
 if xsz > 100 :
  Oo000ooOOO = int ( xsz / 100 )
 else :
  Oo000ooOOO = 1
  if 31 - 31: iIii1I11I1II1 % I11i % ooOoO0o . Ii1I - I11i
 for O0o0 in range ( xsz ) :
  if ( i1OO0oOOoo % Oo000ooOOO ) == 0 : print "%i%% ..." % ( 100. * i1OO0oOOoo / xsz ) ,
  sys . stdout . flush ( )
  i1OO0oOOoo += 1
  for OO00Oo in range ( ysz ) :
   if 17 - 17: Ii1I
   O0OOO0OOoO0O = xll + O0o0 * cellSize
   O00Oo000ooO0 = O0OOO0OOoO0O + cellSize
   OoO0O00IIiII = yll + OO00Oo * cellSize
   o0 = OoO0O00IIiII + cellSize
   if 27 - 27: i11iIiiIii % II111iiii % I11i . O0 - Oo0Ooo + OoOoOO00
   III1Iiii1I11 = int ( ( O0OOO0OOoO0O - ii ) / Iii )
   IIII = int ( ( O00Oo000ooO0 - ii ) / Iii )
   o00oooO0Oo = int ( ( OoO0O00IIiII - OO0O0Ooo ) / Iii )
   o0O0OOO0Ooo = int ( ( o0 - OO0O0Ooo ) / Iii )
   if 57 - 57: iIii1I11I1II1 / I11i - i1IIi
   ooOOo00O00Oo = IIII - III1Iiii1I11 + 1
   IiII1 = o0O0OOO0Ooo - o00oooO0Oo + 1
   if 18 - 18: ooOoO0o * OoOoOO00 . iII111i / I1ii11iIi11i / i11iIiiIii
   if III1Iiii1I11 >= 0 and IIII < iIIiIiI1I1 and o00oooO0Oo >= 1 and o0O0OOO0Ooo < ooO :
    IIIII = OOOOO0O00 . ReadAsArray ( xoff = III1Iiii1I11 , yoff = ooO - o0O0OOO0Ooo , xsize = ooOOo00O00Oo , ysize = IiII1 ) . transpose ( ) . copy ( )
    IIIII = numpy . array ( IIIII [ : , : : - 1 ] , arrayType )
    if 78 - 78: Ii1I * i1IIi
    IIIII [ numpy . where ( IIIII == ndv ) ] = ndr
    IIIII [ numpy . where ( IIIII < - 1e6 ) ] = ndr
    IIIII [ numpy . where ( numpy . isnan ( IIIII ) ) ] = ndr
    if 1 - 1: I1IiiI / IiII * ooOoO0o
    if IIIII . max ( ) == ndr :
     continue
     if 1 - 1: I11i * o0oOOo0O0Ooo . OoOoOO00 / O0
     if 100 - 100: I1Ii111 . o0oOOo0O0Ooo * Oo0Ooo % O0 * O0
    if nFileName is not None :
     IIIii1 = IIii11I1i1I . ReadAsArray ( xoff = III1Iiii1I11 , yoff = ooO - o0O0OOO0Ooo , xsize = ooOOo00O00Oo , ysize = IiII1 ) . transpose ( ) . copy ( )
     IIIii1 = IIIii1 [ : , : : - 1 ]
    else :
     IIIii1 = IIIII . copy ( )
     IIIii1 [ : , : ] = nFP
     if 71 - 71: II111iiii / i1IIi . I1ii11iIi11i % OoooooooOO . OoOoOO00
     if 41 - 41: i1IIi * II111iiii / OoooooooOO . OOooOOo
     if 83 - 83: iII111i . O0 / Oo0Ooo / OOooOOo - II111iiii
    conveyanceFunc ( 0 , 0 , 0 , o0O0OOO0Ooo - o00oooO0Oo , IIIII , Iii , ooOOo00O00Oo , IiII1 , nFP , oooooOOO000Oo , False , 0. , IIIII , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , IIIii1 )
    if 100 - 100: OoO0O00
    if 46 - 46: OoOoOO00 / iIii1I11I1II1 % iII111i . iIii1I11I1II1 * iII111i
    if 38 - 38: I1ii11iIi11i - iII111i / O0 . I1Ii111
    if 45 - 45: I1Ii111
    if 83 - 83: OoOoOO00 . OoooooooOO
    iIi1Ii1i1iI [ O0o0 , OO00Oo , : ] = oooooOOO000Oo [ : ]
    if 58 - 58: i11iIiiIii + OoooooooOO % OoooooooOO / IiII / i11iIiiIii
    if 62 - 62: OoO0O00 / I1ii11iIi11i
    conveyanceFunc ( 0 , 0 , IIII - III1Iiii1I11 , 0 , IIIII , Iii , IIII - III1Iiii1I11 + 1 , o0O0OOO0Ooo - o00oooO0Oo + 1 , nFP , oooooOOO000Oo , False , 0. , IIIII , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , IIIii1 )
    if 7 - 7: OoooooooOO . IiII
    if 53 - 53: Ii1I % Ii1I * o0oOOo0O0Ooo + OoOoOO00
    if 92 - 92: OoooooooOO + i1IIi / Ii1I * O0
    if 100 - 100: ooOoO0o % iIii1I11I1II1 * II111iiii - iII111i
    if 92 - 92: ooOoO0o
    IIiI1 [ O0o0 , OO00Oo , : ] = oooooOOO000Oo [ : ]
    if 22 - 22: Oo0Ooo % iII111i * I1ii11iIi11i / OOooOOo % i11iIiiIii * I11i
    if 95 - 95: OoooooooOO - IiII * I1IiiI + OoOoOO00
    if storageFunc is not None :
     storageFunc ( arrayType ( O0OOO0OOoO0O ) , arrayType ( OoO0O00IIiII ) , arrayType ( O00Oo000ooO0 ) , arrayType ( o0 ) , IIIII , ooOOo00O00Oo , IiII1 , arrayType ( O0OOO0OOoO0O ) , arrayType ( OoO0O00IIiII ) , arrayType ( Iii ) , Ooo00OoOOO , False , IIIII , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 )
     if 10 - 10: o0oOOo0O0Ooo / i11iIiiIii
     if 92 - 92: I11i . I1Ii111
     if 85 - 85: I1ii11iIi11i . I1Ii111
     if 78 - 78: ooOoO0o * I1Ii111 + iIii1I11I1II1 + iIii1I11I1II1 / I1Ii111 . Ii1I
     if 97 - 97: ooOoO0o / I1Ii111 % i1IIi % I1ii11iIi11i
     ii1 [ O0o0 , OO00Oo , : ] = Ooo00OoOOO
    else :
     ii1 [ O0o0 , OO00Oo , : ] = calcStorageParameters ( [ ( O0OOO0OOoO0O , OoO0O00IIiII ) , ( O00Oo000ooO0 , o0 ) ] , IIIII , O0OOO0OOoO0O , OoO0O00IIiII , Iii , plotName = o0o0oOoOO0O , csvOutput = False )
     if 18 - 18: iIii1I11I1II1 % I11i
     if 95 - 95: ooOoO0o + i11iIiiIii * I1Ii111 - i1IIi * I1Ii111 - iIii1I11I1II1
 print "Done."
 if 75 - 75: OoooooooOO * IiII
 fileIO . saveConveyanceParametersCSV ( iIi1Ii1i1iI , IIiI1 , xll , yll , cellSize , outputPrefix + "conveyanceParams.csv" )
 if 9 - 9: IiII - II111iiii + O0 / iIii1I11I1II1 / i11iIiiIii
 fileIO . saveStorageParametersCSV ( ii1 , xll , yll , cellSize , outputPrefix + "storageParams.csv" )
 if 39 - 39: IiII * Oo0Ooo + iIii1I11I1II1 - IiII + OOooOOo
 if 69 - 69: O0
 return (iIi1Ii1i1iI , IIiI1 , ii1)
 if 85 - 85: ooOoO0o / O0
 if 18 - 18: o0oOOo0O0Ooo % O0 * I1ii11iIi11i
 if 62 - 62: I1Ii111 . IiII . OoooooooOO
def gridFlowSetupTiled2 ( dtmFileName , xll , yll , cellSize , xsz , ysz , nChan , nFP ,
 plotNamePrefix = None , outputPrefix = None ,
 ndv = None , ndr = None , conveyanceFunc = None , storageFunc = None ) :
 if 11 - 11: OOooOOo / I11i
 if plotNamePrefix is None :
  plotNamePrefix = ""
  if 73 - 73: i1IIi / i11iIiiIii
 if outputPrefix is None :
  outputPrefix = ""
  if 58 - 58: Oo0Ooo . II111iiii + oO0o - i11iIiiIii / II111iiii / O0
 OOOOO0O00 , Iii , iIIiIiI1I1 , ooO , ii , OO0O0Ooo = fileIO . readScalarGridObj ( dtmFileName )
 if 85 - 85: OoOoOO00 + OOooOOo
 if 10 - 10: IiII / OoO0O00 + OoOoOO00 / i1IIi
 ndv = arrayType ( ndv )
 if 27 - 27: Ii1I
 if 67 - 67: I1IiiI
 iIi1Ii1i1iI = numpy . zeros ( ( xsz + 1 , ysz , 6 ) , dtype = arrayType ) - 9999.
 IIiI1 = numpy . zeros ( ( xsz , ysz + 1 , 6 ) , dtype = arrayType ) - 9999.
 if 55 - 55: I1ii11iIi11i - iII111i * o0oOOo0O0Ooo + OoOoOO00 * OoOoOO00 * O0
 ii1 = numpy . zeros ( ( xsz , ysz , 5 ) , dtype = arrayType ) - 9999.
 if 91 - 91: I1Ii111 - OOooOOo % iIii1I11I1II1 - OoooooooOO % ooOoO0o
 OO0 = 0
 o0o0oOoOO0O = None
 if 44 - 44: iII111i - I1Ii111 / O0 * Oo0Ooo + II111iiii / OoOoOO00
 oooooOOO000Oo = numpy . zeros ( 6 , dtype = arrayType )
 Ooo00OoOOO = numpy . zeros ( 5 , dtype = arrayType )
 OOOOoO000 = numpy . zeros ( 7 , dtype = arrayType )
 if 57 - 57: II111iiii
 i1OO0oOOoo = 0
 if 54 - 54: Oo0Ooo + oO0o + i11iIiiIii
 i1i1ii111 = 10
 if 3 - 3: II111iiii / OOooOOo + IiII . ooOoO0o . OoO0O00
 if xsz > 100 :
  Oo000ooOOO = int ( xsz / 100 )
 else :
  Oo000ooOOO = 1
  if 83 - 83: oO0o + OoooooooOO
 for O0o0 in range ( 0 , xsz , i1i1ii111 ) :
  print O0o0 , xsz
  if ( i1OO0oOOoo % Oo000ooOOO ) == 0 : print "%i%% ..." % ( 100. * i1OO0oOOoo / xsz ) ,
  sys . stdout . flush ( )
  i1OO0oOOoo += 1
  for OO00Oo in range ( 0 , ysz , i1i1ii111 ) :
   if 22 - 22: Ii1I % iII111i * OoooooooOO - o0oOOo0O0Ooo / iIii1I11I1II1
   if 86 - 86: OoooooooOO . iII111i % OoOoOO00 / I11i * iII111i / o0oOOo0O0Ooo
   if 64 - 64: i11iIiiIii
   if 38 - 38: IiII / I1IiiI - IiII . I11i
   O0OOO0OOoO0O = xll + O0o0 * cellSize
   O00Oo000ooO0 = O0OOO0OOoO0O + cellSize * i1i1ii111
   OoO0O00IIiII = yll + OO00Oo * cellSize
   o0 = OoO0O00IIiII + cellSize * i1i1ii111
   if 69 - 69: OoooooooOO + I1ii11iIi11i
   III1Iiii1I11 = int ( ( O0OOO0OOoO0O - ii ) / Iii )
   IIII = int ( ( O00Oo000ooO0 - ii ) / Iii )
   o00oooO0Oo = int ( ( OoO0O00IIiII - OO0O0Ooo ) / Iii )
   o0O0OOO0Ooo = int ( ( o0 - OO0O0Ooo ) / Iii )
   if 97 - 97: OOooOOo - OoO0O00 / Ii1I . i11iIiiIii % oO0o * oO0o
   if 1 - 1: I1IiiI % ooOoO0o
   IIII = min ( IIII , iIIiIiI1I1 - 1 )
   o0O0OOO0Ooo = min ( o0O0OOO0Ooo , ooO - 1 )
   if 65 - 65: I1IiiI + OoOoOO00 / OOooOOo
   ooOOo00O00Oo = IIII - III1Iiii1I11 + 1
   IiII1 = o0O0OOO0Ooo - o00oooO0Oo + 1
   if 83 - 83: o0oOOo0O0Ooo . iII111i - Oo0Ooo
   if 65 - 65: iIii1I11I1II1 / ooOoO0o . IiII - II111iiii
   IIIII = OOOOO0O00 . ReadAsArray ( xoff = III1Iiii1I11 , yoff = ooO - o0O0OOO0Ooo , xsize = ooOOo00O00Oo , ysize = IiII1 ) . transpose ( ) . copy ( )
   if 72 - 72: iIii1I11I1II1 / IiII % iII111i % OOooOOo - I11i % OOooOOo
   IIIII = numpy . array ( IIIII [ : , : : - 1 ] , arrayType )
   IIIII [ numpy . where ( IIIII == ndv ) ] = ndr
   IIIII [ numpy . where ( numpy . isnan ( IIIII ) ) ] = ndr
   if 100 - 100: Oo0Ooo + i11iIiiIii
   if 71 - 71: I11i / o0oOOo0O0Ooo / I1Ii111 % OOooOOo
   for O0oooo00o0Oo in range ( O0o0 , O0o0 + i1i1ii111 ) :
    for I1iii in range ( OO00Oo , OO00Oo + i1i1ii111 ) :
     if 86 - 86: I1ii11iIi11i * O0 * IiII
     if O0oooo00o0Oo >= xsz or I1iii >= ysz :
      continue
      if 51 - 51: II111iiii + IiII . i1IIi . I1ii11iIi11i + OoOoOO00 * I1IiiI
      if 72 - 72: oO0o + oO0o / II111iiii . OoooooooOO % Ii1I
      if 49 - 49: oO0o . OoO0O00 - Oo0Ooo * OoooooooOO . Oo0Ooo
     ii1Ii1IiIIi = O0OOO0OOoO0O + ( O0oooo00o0Oo - O0o0 ) * cellSize
     o0OO0 = ii1Ii1IiIIi + cellSize
     oOo00Oo0o0Oo = OoO0O00IIiII + ( I1iii - OO00Oo ) * cellSize
     I1 = oOo00Oo0o0Oo + cellSize
     if 26 - 26: ooOoO0o . OOooOOo - OOooOOo . OoO0O00
     Ii1IiII = int ( ( ii1Ii1IiIIi - O0OOO0OOoO0O ) / Iii )
     I1i = int ( ( o0OO0 - O0OOO0OOoO0O ) / Iii )
     oooii1iiIi1 = int ( ( oOo00Oo0o0Oo - OoO0O00IIiII ) / Iii )
     i111iiI1ii = int ( ( I1 - OoO0O00IIiII ) / Iii )
     if 24 - 24: OoOoOO00 / OoooooooOO . II111iiii . I1IiiI % O0 % Ii1I
     if 5 - 5: OoooooooOO - OoO0O00 + IiII - iII111i . OoO0O00 / ooOoO0o
     i1I1i1i1iII1 = IIIII [ Ii1IiII : I1i + 1 , oooii1iiIi1 : i111iiI1ii + 1 ]
     if 68 - 68: Oo0Ooo + i11iIiiIii
     if 69 - 69: iIii1I11I1II1 * iIii1I11I1II1 * i11iIiiIii + I1IiiI / OOooOOo % Ii1I
     if 58 - 58: OOooOOo * o0oOOo0O0Ooo + O0 % OOooOOo
     if 25 - 25: Oo0Ooo % I1ii11iIi11i * ooOoO0o
     if 6 - 6: iII111i . IiII * OoOoOO00 . i1IIi
     if 98 - 98: i1IIi
     if 65 - 65: OoOoOO00 / OoO0O00 % IiII
     if 45 - 45: OoOoOO00
     if 66 - 66: OoO0O00
     if 56 - 56: O0
     if 61 - 61: o0oOOo0O0Ooo / OOooOOo / Oo0Ooo * O0
     if 23 - 23: oO0o - OOooOOo + I11i
     if 12 - 12: I1IiiI / ooOoO0o % o0oOOo0O0Ooo / i11iIiiIii % OoooooooOO
     if 15 - 15: iIii1I11I1II1 % OoooooooOO - Oo0Ooo * Ii1I + I11i
     if 11 - 11: iII111i * Ii1I - OoOoOO00
     OOO , III1iI1iII1I = i1I1i1i1iII1 . shape
     if 39 - 39: Ii1I * ooOoO0o / OoOoOO00 * OoO0O00 . I11i % II111iiii
     if 71 - 71: I1Ii111 % i1IIi - II111iiii - OOooOOo + OOooOOo * ooOoO0o
     conveyanceFunc ( 0 , 0 , 0 , III1iI1iII1I - 1 , i1I1i1i1iII1 , Iii , OOO , III1iI1iII1I , nFP , oooooOOO000Oo , False , 0. , i1I1i1i1iII1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 )
     if 51 - 51: iIii1I11I1II1 / OoOoOO00 + OOooOOo - I11i + iII111i
     if 29 - 29: o0oOOo0O0Ooo % iIii1I11I1II1 . OoooooooOO % OoooooooOO % II111iiii / iII111i
     if 70 - 70: i11iIiiIii % iII111i
     if 11 - 11: IiII % I1ii11iIi11i % Ii1I / II111iiii % I1Ii111 - Oo0Ooo
     if 96 - 96: I1ii11iIi11i / II111iiii . Ii1I - iII111i * I11i * oO0o
     iIi1Ii1i1iI [ O0oooo00o0Oo , I1iii , : ] = oooooOOO000Oo [ : ]
     if 76 - 76: Ii1I - II111iiii * OOooOOo / OoooooooOO
     if 18 - 18: OoO0O00 + iIii1I11I1II1 - II111iiii - I1IiiI
     conveyanceFunc ( 0 , 0 , OOO - 1 , 0 , i1I1i1i1iII1 , Iii , OOO , III1iI1iII1I , nFP , oooooOOO000Oo , False , 0. , i1I1i1i1iII1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 )
     if 71 - 71: OoooooooOO
     if 33 - 33: I1Ii111
     if 62 - 62: I1ii11iIi11i + Ii1I + i1IIi / OoooooooOO
     if 7 - 7: o0oOOo0O0Ooo + i1IIi . I1IiiI / Oo0Ooo
     if 22 - 22: ooOoO0o - ooOoO0o % OOooOOo . I1Ii111 + oO0o
     IIiI1 [ O0oooo00o0Oo , I1iii , : ] = oooooOOO000Oo [ : ]
     if 63 - 63: I1IiiI % I1Ii111 * o0oOOo0O0Ooo + I1Ii111 / Oo0Ooo % iII111i
     if 45 - 45: IiII
     if 20 - 20: OoooooooOO * o0oOOo0O0Ooo * O0 . OOooOOo
     storageFunc ( arrayType ( ii1Ii1IiIIi ) , arrayType ( oOo00Oo0o0Oo ) , arrayType ( o0OO0 ) , arrayType ( I1 ) , i1I1i1i1iII1 , OOO , III1iI1iII1I , arrayType ( ii1Ii1IiIIi ) , arrayType ( oOo00Oo0o0Oo ) , arrayType ( Iii ) , Ooo00OoOOO , False , i1I1i1i1iII1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 )
     if 78 - 78: iIii1I11I1II1 + I11i - Ii1I * I1Ii111 - OoooooooOO % OoOoOO00
     if 34 - 34: O0
     if 80 - 80: i1IIi - Oo0Ooo / OoO0O00 - i11iIiiIii
     if 68 - 68: oO0o - I1ii11iIi11i % O0 % I1Ii111
     if 11 - 11: O0 / OoO0O00 % OOooOOo + o0oOOo0O0Ooo + iIii1I11I1II1
     ii1 [ O0oooo00o0Oo , I1iii , : ] = Ooo00OoOOO
     if 40 - 40: ooOoO0o - OOooOOo . Ii1I * Oo0Ooo % I1Ii111
     if 56 - 56: i11iIiiIii . o0oOOo0O0Ooo - I1IiiI * I11i
     if 91 - 91: oO0o + OoooooooOO - i1IIi
     if 84 - 84: Ii1I / IiII
     if 86 - 86: OoOoOO00 * II111iiii - O0 . OoOoOO00 % iIii1I11I1II1 / OOooOOo
     if 11 - 11: I1IiiI * oO0o + I1ii11iIi11i / I1ii11iIi11i
 print "Done."
 if 37 - 37: i11iIiiIii + i1IIi
 fileIO . saveConveyanceParametersCSV ( iIi1Ii1i1iI , IIiI1 , xll , yll , cellSize , outputPrefix + "conveyanceParams.csv" )
 if 23 - 23: iII111i + I11i . OoOoOO00 * I1IiiI + I1ii11iIi11i
 fileIO . saveStorageParametersCSV ( ii1 , xll , yll , cellSize , outputPrefix + "storageParams.csv" )
 if 18 - 18: IiII * o0oOOo0O0Ooo . IiII / O0
 if 8 - 8: o0oOOo0O0Ooo
 return iIi1Ii1i1iI , IIiI1 , ii1
 if 4 - 4: I1ii11iIi11i + I1ii11iIi11i * ooOoO0o - OoOoOO00
 if 78 - 78: Ii1I / II111iiii % OoOoOO00
 if 52 - 52: OOooOOo - iII111i * oO0o
def gridFlowSetup ( dtmFileName , xll , yll , cellSize , xsz , ysz , nChan , nFP ,
 plotNamePrefix = None , outputPrefix = None ,
 ndv = None , ndr = None , conveyanceFunc = None , storageFunc = None , catchmentAreaGridFile = None , chanExp = None , chanMult = None , chanAR = None , chanMaxD = None ) :
 if 17 - 17: OoooooooOO + OOooOOo * I11i * OoOoOO00
 if 36 - 36: O0 + Oo0Ooo
 if plotNamePrefix is None :
  plotNamePrefix = ""
  if 5 - 5: Oo0Ooo * OoOoOO00
 if outputPrefix is None :
  outputPrefix = ""
  if 46 - 46: ooOoO0o
 I11iIiII , Iii , iIIiIiI1I1 , ooO , ii , OO0O0Ooo = fileIO . readScalarGrid ( dtmFileName , dataType = arrayType )
 if 66 - 66: Oo0Ooo - o0oOOo0O0Ooo * IiII + OoOoOO00 + o0oOOo0O0Ooo - iIii1I11I1II1
 if 17 - 17: oO0o
 if catchmentAreaGridFile is not None :
  i1ii11 , ii1i , o00oO0oOo00 , oO0oOo0 , IIi , oo0OO = fileIO . readScalarGrid ( catchmentAreaGridFile , dataType = arrayType )
  if 2 - 2: II111iiii - OoO0O00 . IiII * iII111i / oO0o
  if 80 - 80: OOooOOo / I11i / OoOoOO00 + i1IIi - Oo0Ooo
  i1ii11 [ numpy . where ( numpy . isnan ( i1ii11 ) ) ] = 0.
  if 11 - 11: o0oOOo0O0Ooo * OoO0O00
  if 15 - 15: OoOoOO00
  if 62 - 62: Ii1I
 ndv = arrayType ( ndv )
 if 51 - 51: OoOoOO00
 if ndv is not None and ndr is not None :
  replaceArrayVals ( I11iIiII , ndv , ndr )
  if 14 - 14: IiII % oO0o % Oo0Ooo - i11iIiiIii
  if 53 - 53: Ii1I % Oo0Ooo
 I11iIiII [ numpy . where ( numpy . isnan ( I11iIiII ) ) ] = ndr
 if 59 - 59: OOooOOo % iIii1I11I1II1 . i1IIi + II111iiii * IiII
 iIi1Ii1i1iI = numpy . zeros ( ( xsz + 1 , ysz , 6 ) , dtype = arrayType )
 IIiI1 = numpy . zeros ( ( xsz , ysz + 1 , 6 ) , dtype = arrayType )
 if 41 - 41: Ii1I % I1ii11iIi11i
 if catchmentAreaGridFile is None :
  ii1 = numpy . zeros ( ( xsz , ysz , 5 ) , dtype = arrayType )
 else :
  ii1 = numpy . zeros ( ( xsz , ysz , 7 ) , dtype = arrayType )
  if 12 - 12: OOooOOo
 ooOo0O = numpy . zeros ( ( xsz , ysz ) , dtype = arrayType )
 if 37 - 37: Ii1I % OoO0O00
 OO0 = 0
 o0o0oOoOO0O = None
 if 79 - 79: I1ii11iIi11i + I1IiiI / I1IiiI
 oooooOOO000Oo = numpy . zeros ( 6 , dtype = arrayType )
 if 71 - 71: OOooOOo * OoO0O00 % OoooooooOO % OoO0O00 / I1IiiI
 if 56 - 56: OoooooooOO % i11iIiiIii * iIii1I11I1II1 . OoO0O00 * O0
 i1OO0oOOoo = 0
 print "Calculating X-conveyance" ,
 for O0o0 in range ( xsz + 1 ) :
  if ( i1OO0oOOoo % int ( xsz / 10 ) ) == 0 : print "." ,
  i1OO0oOOoo += 1
  for OO00Oo in range ( ysz ) :
   if O0o0 == 0 or O0o0 == xsz :
    iIi1Ii1i1iI [ O0o0 , OO00Oo , : ] = - 9999.
   else :
    iIii = xll + O0o0 * cellSize
    OoO0O00IIiII = yll + OO00Oo * cellSize
    o0 = OoO0O00IIiII + cellSize
    if 23 - 23: i11iIiiIii
    II1I11IIi = int ( ( iIii - ii ) / Iii )
    o00oooO0Oo = int ( ( OoO0O00IIiII - OO0O0Ooo ) / Iii )
    o0O0OOO0Ooo = int ( ( o0 - OO0O0Ooo ) / Iii )
    if 66 - 66: i11iIiiIii / o0oOOo0O0Ooo - OoooooooOO / i1IIi . i11iIiiIii
    if II1I11IIi < 0 or II1I11IIi >= iIIiIiI1I1 or o00oooO0Oo < 0 or o0O0OOO0Ooo > ooO :
     iIi1Ii1i1iI [ O0o0 , OO00Oo , : ] = - 9999.
     continue
     if 16 - 16: Oo0Ooo % I1ii11iIi11i + I11i - O0 . iII111i / I1Ii111
    if conveyanceFunc is None :
     oO , O0o0O00Oo0o0 = i1ii1iiI ( [ ( II1I11IIi , o00oooO0Oo ) , ( II1I11IIi , o0O0OOO0Ooo ) ] , I11iIiII , ii , OO0O0Ooo , Iii )
     if 35 - 35: oO0o / I1Ii111 / II111iiii - iIii1I11I1II1 + II111iiii . I1Ii111
     if catchmentAreaGridFile is not None :
      O0O00O000OOO = int ( ( iIii - IIi ) / ii1i )
      iI = int ( ( OoO0O00IIiII - oo0OO ) / ii1i )
      Oo0O = int ( ( o0 - oo0OO ) / ii1i )
      if 1 - 1: O0 / iII111i % I1Ii111 . Oo0Ooo + IiII
      if O0O00O000OOO > 0 and O0O00O000OOO < o00oO0oOo00 and iI > 0 and Oo0O < oO0oOo0 and Oo0O > 0 and Oo0O < oO0oOo0 :
       if 27 - 27: I1Ii111 % OoooooooOO + IiII % i1IIi / oO0o / OoooooooOO
       if 11 - 11: OOooOOo % Ii1I - i11iIiiIii - oO0o + ooOoO0o + IiII
       if 87 - 87: I1Ii111 * i1IIi / I1ii11iIi11i
       IIII1i1 , oooOoOoOoO = i1ii1iiI ( [ ( O0O00O000OOO , iI ) , ( O0O00O000OOO , Oo0O ) ] , i1ii11 , IIi , oo0OO , ii1i , strict = True )
       if 62 - 62: i1IIi + Oo0Ooo % IiII
       if 28 - 28: I1ii11iIi11i . i1IIi
       OOOO0oo0 = chanMult * ( max ( oooOoOoOoO ) ** chanExp )
      else :
       OOOO0oo0 = 0.
       if 10 - 10: OoO0O00 / Oo0Ooo
      I11iiI1i1 = min ( OOOO0oo0 / chanAR , chanMaxD )
      if 15 - 15: iII111i . OoOoOO00 / iII111i * I11i - I1IiiI % I1ii11iIi11i
      oo0OOOOOO0 = O0o0O00Oo0o0 [ : ]
      i11 = Ii1I1I11I ( O0o0O00Oo0o0 , oO , OOOO0oo0 , I11iiI1i1 , nChan , nFP )
      if 29 - 29: OoooooooOO . I1IiiI % I1ii11iIi11i - iII111i
      if 8 - 8: i1IIi
     else :
      i11 = None
      if 32 - 32: oO0o / II111iiii
     iIi1Ii1i1iI [ O0o0 , OO00Oo , : ] = conveyanceParameters ( O0o0O00Oo0o0 , oO , nFP , plotName = o0o0oOoOO0O , csvOutput = False , nList = i11 )
     if 45 - 45: I1ii11iIi11i + OoO0O00 * i11iIiiIii / OOooOOo % I11i * O0
     if 17 - 17: O0
    else :
     if catchmentAreaGridFile is not None :
      if 88 - 88: Oo0Ooo . O0 % OoooooooOO / OOooOOo
      O0O00O000OOO = int ( ( iIii - IIi ) / ii1i )
      iI = int ( ( OoO0O00IIiII - oo0OO ) / ii1i )
      Oo0O = int ( ( o0 - oo0OO ) / ii1i )
      if 89 - 89: II111iiii / oO0o
      conveyanceFunc ( II1I11IIi , o00oooO0Oo , II1I11IIi , o0O0OOO0Ooo , I11iIiII , Iii , iIIiIiI1I1 , ooO , nFP , oooooOOO000Oo ,
      # iII111i
      # OOooOOo
 True , nChan , i1ii11 , o00oO0oOo00 , oO0oOo0 ,
 O0O00O000OOO , iI , O0O00O000OOO , Oo0O ,
 chanMult , chanExp , chanAR , chanMaxD )
      if 87 - 87: ooOoO0o + o0oOOo0O0Ooo
     else :
      if 28 - 28: OOooOOo * I1ii11iIi11i / oO0o
      conveyanceFunc ( II1I11IIi , o00oooO0Oo , II1I11IIi , o0O0OOO0Ooo , I11iIiII , Iii , iIIiIiI1I1 , ooO , nFP , oooooOOO000Oo , False , 0. , I11iIiII , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 )
      if 64 - 64: oO0o - I1IiiI / iII111i - OoO0O00
      if 37 - 37: i11iIiiIii / iII111i
      if 85 - 85: i11iIiiIii + I1Ii111 * OoOoOO00
      if 1 - 1: i1IIi / Oo0Ooo . OoO0O00
      if 57 - 57: I11i . Oo0Ooo + II111iiii
     iIi1Ii1i1iI [ O0o0 , OO00Oo , : ] = oooooOOO000Oo [ : ]
     if 43 - 43: I1Ii111 % iII111i
 print " Done"
 if 69 - 69: iII111i % OoO0O00
 if 86 - 86: oO0o / oO0o
 if 28 - 28: i11iIiiIii / o0oOOo0O0Ooo . iIii1I11I1II1 / II111iiii
 if 72 - 72: OoooooooOO / I1IiiI + Ii1I / OoOoOO00 * Ii1I
 OO0 = 0
 i1OO0oOOoo = 0
 print "Calculating Y-conveyance" ,
 for O0o0 in range ( xsz ) :
  if ( i1OO0oOOoo % int ( xsz / 10 ) ) == 0 : print "." ,
  i1OO0oOOoo += 1
  for OO00Oo in range ( ysz + 1 ) :
   if OO00Oo == 0 or OO00Oo == ysz :
    IIiI1 [ O0o0 , OO00Oo , : ] = - 9999.
   else :
    O0OOO0OOoO0O = xll + O0o0 * cellSize
    O00Oo000ooO0 = O0OOO0OOoO0O + cellSize
    ooo0O = yll + OO00Oo * cellSize
    if 34 - 34: O0 * O0 % OoooooooOO + iII111i * iIii1I11I1II1 % Ii1I
    III1Iiii1I11 = int ( ( O0OOO0OOoO0O - ii ) / Iii )
    IIII = int ( ( O00Oo000ooO0 - ii ) / Iii )
    I1iI1I1 = int ( ( ooo0O - OO0O0Ooo ) / Iii )
    if 48 - 48: I1IiiI / i11iIiiIii - o0oOOo0O0Ooo * oO0o / OoooooooOO
    if I1iI1I1 < 0 or I1iI1I1 >= ooO or III1Iiii1I11 < 0 or IIII > iIIiIiI1I1 :
     IIiI1 [ O0o0 , OO00Oo , : ] = - 9999.
     continue
     if 89 - 89: iIii1I11I1II1 / I1IiiI - II111iiii / Ii1I . i11iIiiIii . Ii1I
    if conveyanceFunc is None :
     oO , O0o0O00Oo0o0 = i1ii1iiI ( [ ( III1Iiii1I11 , I1iI1I1 ) , ( IIII , I1iI1I1 ) ] , I11iIiII , ii , OO0O0Ooo , Iii )
     if 48 - 48: O0 + O0 . I1Ii111 - ooOoO0o
     if catchmentAreaGridFile is not None :
      o00oo0000 = int ( ( O0OOO0OOoO0O - IIi ) / ii1i )
      iIi1IIi1ii = int ( ( O00Oo000ooO0 - IIi ) / ii1i )
      I11Ii = int ( ( ooo0O - oo0OO ) / ii1i )
      if 16 - 16: Oo0Ooo / i11iIiiIii
      if o00oo0000 > 0 and o00oo0000 < o00oO0oOo00 and iIi1IIi1ii > 0 and iIi1IIi1ii < o00oO0oOo00 and I11Ii > 0 and I11Ii < oO0oOo0 :
       if 64 - 64: i11iIiiIii / Ii1I * i1IIi
       if 73 - 73: Oo0Ooo - OoOoOO00 - oO0o - I1IiiI
       if 65 - 65: o0oOOo0O0Ooo
       IIII1i1 , oooOoOoOoO = i1ii1iiI ( [ ( o00oo0000 , I11Ii ) , ( iIi1IIi1ii , I11Ii ) ] , i1ii11 , IIi , oo0OO , ii1i , strict = True )
       if 7 - 7: IiII . OoOoOO00 / I1ii11iIi11i . OOooOOo * I11i - II111iiii
       if 37 - 37: I1Ii111 . OoOoOO00 / O0 * iII111i
       OOOO0oo0 = chanMult * ( max ( oooOoOoOoO ) ** chanExp )
      else :
       OOOO0oo0 = 0.
       if 7 - 7: OoO0O00 * I11i + II111iiii % i11iIiiIii
      I11iiI1i1 = min ( OOOO0oo0 / chanAR , chanMaxD )
      i11 = Ii1I1I11I ( O0o0O00Oo0o0 , oO , OOOO0oo0 , I11iiI1i1 , nChan , nFP )
     else :
      i11 = None
      if 8 - 8: ooOoO0o * O0
      if 73 - 73: o0oOOo0O0Ooo / oO0o / I11i / OoO0O00
      if 11 - 11: OoOoOO00 + IiII - OoooooooOO / OoO0O00
      if 34 - 34: ooOoO0o
      if 45 - 45: ooOoO0o / Oo0Ooo / Ii1I
      if 44 - 44: I1ii11iIi11i - Ii1I / II111iiii * OoO0O00 * Oo0Ooo
      if 73 - 73: o0oOOo0O0Ooo - I1IiiI * i1IIi / i11iIiiIii * OOooOOo % II111iiii
     IIiI1 [ O0o0 , OO00Oo , : ] = conveyanceParameters ( O0o0O00Oo0o0 , oO , nFP , plotName = o0o0oOoOO0O , csvOutput = False , nList = i11 )
     if 56 - 56: OoooooooOO * Oo0Ooo . Oo0Ooo . I1ii11iIi11i
     if 24 - 24: Oo0Ooo . I11i * Ii1I % iII111i / OOooOOo
    else :
     if 58 - 58: I1IiiI - I1ii11iIi11i % O0 . I1IiiI % OoO0O00 % IiII
     if catchmentAreaGridFile is not None :
      if 87 - 87: oO0o - i11iIiiIii
      o00oo0000 = int ( ( O0OOO0OOoO0O - IIi ) / ii1i )
      iIi1IIi1ii = int ( ( O00Oo000ooO0 - IIi ) / ii1i )
      I11Ii = int ( ( ooo0O - oo0OO ) / ii1i )
      if 78 - 78: i11iIiiIii / iIii1I11I1II1 - o0oOOo0O0Ooo
      conveyanceFunc ( III1Iiii1I11 , I1iI1I1 , IIII , I1iI1I1 , I11iIiII , arrayType ( Iii ) , iIIiIiI1I1 , ooO , nFP , oooooOOO000Oo , True , nChan , i1ii11 , o00oO0oOo00 , oO0oOo0 , o00oo0000 , I11Ii , iIi1IIi1ii , I11Ii , chanMult , chanExp , chanAR , chanMaxD )
      if 23 - 23: I11i
      if 40 - 40: o0oOOo0O0Ooo - II111iiii / Oo0Ooo
      if 14 - 14: I1ii11iIi11i
      if 5 - 5: o0oOOo0O0Ooo . iIii1I11I1II1 % iIii1I11I1II1
      if 56 - 56: OoooooooOO - I11i - i1IIi
      if 8 - 8: I1Ii111 / OOooOOo . I1IiiI + I1ii11iIi11i / i11iIiiIii
      if 31 - 31: ooOoO0o - iIii1I11I1II1 + iII111i . Oo0Ooo / IiII % iIii1I11I1II1
      if 6 - 6: IiII * i11iIiiIii % iIii1I11I1II1 % i11iIiiIii + o0oOOo0O0Ooo / i1IIi
     else :
      if 53 - 53: I11i + iIii1I11I1II1
      conveyanceFunc ( III1Iiii1I11 , I1iI1I1 , IIII , I1iI1I1 , I11iIiII , arrayType ( Iii ) , iIIiIiI1I1 , ooO , nChan , oooooOOO000Oo , False , 0 , I11iIiII , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 )
      if 70 - 70: I1ii11iIi11i
      if 67 - 67: OoooooooOO
      if 29 - 29: O0 - i11iIiiIii - II111iiii + OOooOOo * IiII
      if 2 - 2: i1IIi - ooOoO0o + I1IiiI . o0oOOo0O0Ooo * o0oOOo0O0Ooo / OoOoOO00
      if 93 - 93: i1IIi
      if 53 - 53: OoooooooOO + Oo0Ooo + oO0o
     IIiI1 [ O0o0 , OO00Oo , : ] = oooooOOO000Oo [ : ]
     if 24 - 24: iII111i - IiII - iII111i * I1ii11iIi11i . OoooooooOO / IiII
 print " Done"
 if 66 - 66: Oo0Ooo
 if 97 - 97: i1IIi - OoooooooOO / I1Ii111 * I1IiiI
 OO0 = 0
 i1OO0oOOoo = 0
 print "Calculating storage" ,
 if 55 - 55: o0oOOo0O0Ooo . iII111i
 oOo00o00oO = ii + iIIiIiI1I1 * Iii ;
 o0000 = OO0O0Ooo + ooO * Iii ;
 Ooo00OoOOO = numpy . zeros ( 5 , dtype = arrayType )
 OOOOoO000 = numpy . zeros ( 7 , dtype = arrayType )
 if 42 - 42: I1Ii111 + I1Ii111 * II111iiii
 for O0o0 in range ( xsz ) :
  if ( i1OO0oOOoo % int ( xsz / 10 ) ) == 0 : print "." ,
  i1OO0oOOoo += 1
  for OO00Oo in range ( ysz ) :
   O0OOO0OOoO0O = xll + O0o0 * cellSize
   O00Oo000ooO0 = O0OOO0OOoO0O + cellSize
   if 78 - 78: OoooooooOO
   OoO0O00IIiII = yll + OO00Oo * cellSize
   o0 = OoO0O00IIiII + cellSize
   if 77 - 77: I1ii11iIi11i / i1IIi / Oo0Ooo % OOooOOo
   if O0OOO0OOoO0O < ii or O00Oo000ooO0 > oOo00o00oO or OoO0O00IIiII < OO0O0Ooo or o0 > o0000 :
    ii1 [ O0o0 , OO00Oo , : ] = - 9999.
    continue
    if 48 - 48: I11i - IiII + iIii1I11I1II1 + OoooooooOO
   if storageFunc is not None :
    if 4 - 4: II111iiii . I11i + Ii1I * I1Ii111 . ooOoO0o
    if 87 - 87: OoOoOO00 / OoO0O00 / i11iIiiIii
    if catchmentAreaGridFile is not None :
     storageFunc ( arrayType ( O0OOO0OOoO0O ) , arrayType ( OoO0O00IIiII ) , arrayType ( O00Oo000ooO0 ) , arrayType ( o0 ) , I11iIiII , iIIiIiI1I1 , ooO , arrayType ( ii ) , arrayType ( OO0O0Ooo ) , arrayType ( Iii ) , OOOOoO000 , True , i1ii11 , o00oO0oOo00 , oO0oOo0 , IIi , oo0OO , ii1i , chanMult , chanExp , chanAR , chanMaxD , cellSize )
     if 74 - 74: oO0o / I1ii11iIi11i % o0oOOo0O0Ooo
     if 88 - 88: OoOoOO00 - i11iIiiIii % o0oOOo0O0Ooo * I11i + I1ii11iIi11i
     if 52 - 52: II111iiii . I1IiiI + OoOoOO00 % OoO0O00
     if 62 - 62: o0oOOo0O0Ooo
     if 15 - 15: I11i + Ii1I . OOooOOo * OoO0O00 . OoOoOO00
     ii1 [ O0o0 , OO00Oo , : ] = OOOOoO000
    else :
     storageFunc ( arrayType ( O0OOO0OOoO0O ) , arrayType ( OoO0O00IIiII ) , arrayType ( O00Oo000ooO0 ) , arrayType ( o0 ) , I11iIiII , iIIiIiI1I1 , ooO , arrayType ( ii ) , arrayType ( OO0O0Ooo ) , arrayType ( Iii ) , Ooo00OoOOO , False , I11iIiII , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 )
     if 18 - 18: i1IIi % II111iiii + I1Ii111 % Ii1I
     if 72 - 72: iIii1I11I1II1
     if 45 - 45: Oo0Ooo - o0oOOo0O0Ooo % I1Ii111
     if 38 - 38: I1Ii111 % OOooOOo - OoooooooOO
     if 87 - 87: OoO0O00 % I1IiiI
     ii1 [ O0o0 , OO00Oo , : ] = Ooo00OoOOO
   else :
    if catchmentAreaGridFile is not None :
     ii1 [ O0o0 , OO00Oo , : ] = calcStorageParametersChannel ( [ ( O0OOO0OOoO0O , OoO0O00IIiII ) , ( O00Oo000ooO0 , o0 ) ] , I11iIiII , ii , OO0O0Ooo , Iii , i1ii11 , IIi , oo0OO , ii1i , chanExp , chanMult , chanAR , chanMaxD )
     if 77 - 77: iIii1I11I1II1 - i1IIi . oO0o
     if 26 - 26: o0oOOo0O0Ooo * IiII . i1IIi
     if 59 - 59: O0 + i1IIi - o0oOOo0O0Ooo
    else :
     ii1 [ O0o0 , OO00Oo , : ] = calcStorageParameters ( [ ( O0OOO0OOoO0O , OoO0O00IIiII ) , ( O00Oo000ooO0 , o0 ) ] , I11iIiII , ii , OO0O0Ooo , Iii , plotName = o0o0oOoOO0O , csvOutput = False )
     if 62 - 62: i11iIiiIii % OOooOOo . IiII . OOooOOo
     if 84 - 84: i11iIiiIii * OoO0O00
 print "Done."
 if 18 - 18: OOooOOo - Ii1I - OoOoOO00 / I1Ii111 - O0
 fileIO . saveConveyanceParametersCSV ( iIi1Ii1i1iI , IIiI1 , xll , yll , cellSize , outputPrefix + "conveyanceParams.csv" )
 if 30 - 30: O0 + I1ii11iIi11i + II111iiii
 fileIO . saveStorageParametersCSV ( ii1 , xll , yll , cellSize , outputPrefix + "storageParams.csv" )
 if 14 - 14: o0oOOo0O0Ooo / OOooOOo - iIii1I11I1II1 - oO0o % ooOoO0o
 if 49 - 49: ooOoO0o * oO0o / o0oOOo0O0Ooo / Oo0Ooo * iIii1I11I1II1
 return iIi1Ii1i1iI , IIiI1 , ii1
 if 57 - 57: OoOoOO00 - oO0o / ooOoO0o % i11iIiiIii
def saveMonitoringPoints ( f , pts , t , wlGrid , sp , qX = None , qY = None ) :
 f . write ( "%f" % t )
 for I11 in pts :
  O0o0 = I11 [ 0 ]
  OO00Oo = I11 [ 1 ]
  f . write ( ",%f,%f" % ( wlGrid [ O0o0 , OO00Oo ] , wlGrid [ O0o0 , OO00Oo ] - sp [ O0o0 , OO00Oo , 0 ] ) )
  if 100 - 100: I1ii11iIi11i + i11iIiiIii - i1IIi
  if qX is not None :
   f . write ( ",%f,%f,%f,%f" % ( qX [ O0o0 , OO00Oo ] , qX [ O0o0 + 1 , OO00Oo ] , qY [ O0o0 , OO00Oo ] , qY [ O0o0 , OO00Oo + 1 ] ) )
   if 29 - 29: o0oOOo0O0Ooo / i11iIiiIii / I1IiiI % oO0o % i11iIiiIii
   if 18 - 18: OOooOOo + I1Ii111
 f . write ( "\n" )
 return
 if 80 - 80: oO0o + o0oOOo0O0Ooo * Ii1I + OoO0O00
 if 75 - 75: I11i / o0oOOo0O0Ooo / OOooOOo / IiII % ooOoO0o + II111iiii
def dryCheck ( v , Qx , Qy , sp , dt , verbose = False , sources = None ) :
 if 4 - 4: iII111i - Oo0Ooo - IiII - I11i % i11iIiiIii / OoO0O00
 i1iii11 = 10
 if 92 - 92: OoOoOO00 . OoooooooOO - I1Ii111
 for O0o0 in range ( i1iii11 ) :
  Oo0000o0O0O = IIiIiIIiIi ( v , Qx , Qy , sp , dt , verbose , sources )
  if verbose :
   print "In dryCheck, number of flows changed=" , Oo0000o0O0O
  if Oo0000o0O0O == 0 : break
  if 90 - 90: II111iiii / I1IiiI
 return O0o0 + 1
 if 45 - 45: iIii1I11I1II1
 if 28 - 28: oO0o
def IIiIiIIiIi ( v , Qx , Qy , sp , dt , verbose = False , sources = None ) :
 o0oO0 , oo00 = v . shape
 if 52 - 52: I1IiiI + iIii1I11I1II1
 ooOO = 0
 if 45 - 45: OoOoOO00
 if sources is not None :
  for oOI11Ii111I in sources :
   v [ oOI11Ii111I [ 0 ] , oOI11Ii111I [ 1 ] ] += dt * oOI11Ii111I [ 2 ]
   if 98 - 98: iIii1I11I1II1 + I1Ii111 % OoOoOO00 + I11i % OoOoOO00
 for O0o0 in range ( o0oO0 ) :
  for OO00Oo in range ( oo00 ) :
   iI1I1I11IiII = - Qx [ O0o0 , OO00Oo ] * dt
   iIii11iI1II = Qx [ O0o0 + 1 , OO00Oo ] * dt
   I1II1I1I = - Qy [ O0o0 , OO00Oo ] * dt
   OOo0 = Qy [ O0o0 , OO00Oo + 1 ] * dt
   if 58 - 58: OoOoOO00 - iII111i - OoooooooOO
   if ( iI1I1I11IiII + iIii11iI1II + I1II1I1I + OOo0 ) > v [ O0o0 , OO00Oo ] :
    o00 = v [ O0o0 , OO00Oo ] / ( iI1I1I11IiII + iIii11iI1II + I1II1I1I + OOo0 )
    if 38 - 38: O0 - IiII % I1Ii111
    if 64 - 64: iIii1I11I1II1
    if 15 - 15: I1ii11iIi11i + OOooOOo / I1ii11iIi11i / I1Ii111
    if 31 - 31: ooOoO0o + O0 + ooOoO0o . iIii1I11I1II1 + Oo0Ooo / o0oOOo0O0Ooo
    if 6 - 6: Oo0Ooo % IiII * I11i / I1IiiI + Oo0Ooo
    if 39 - 39: OoOoOO00 - Oo0Ooo / iII111i * OoooooooOO
    if 100 - 100: O0 . I11i . OoO0O00 + O0 * oO0o
    if 42 - 42: oO0o % OoooooooOO + o0oOOo0O0Ooo
    if 56 - 56: OoooooooOO + I1ii11iIi11i - iII111i
    if 24 - 24: o0oOOo0O0Ooo + ooOoO0o + I11i - iIii1I11I1II1
    if iI1I1I11IiII > 0 :
     Qx [ O0o0 , OO00Oo ] *= o00
     ooOO += 1
    if iIii11iI1II > 0 :
     Qx [ O0o0 + 1 , OO00Oo ] *= o00
     ooOO += 1
    if I1II1I1I > 0 :
     Qy [ O0o0 , OO00Oo ] *= o00
     ooOO += 1
    if OOo0 > 0 :
     Qy [ O0o0 , OO00Oo + 1 ] *= o00
     ooOO += 1
     if 49 - 49: I11i . ooOoO0o * OoOoOO00 % IiII . O0
    if verbose :
     if 48 - 48: O0 * Ii1I - O0 / Ii1I + OoOoOO00
     oO00oo000O = Qx [ O0o0 , OO00Oo ]
     ii11 = - Qx [ O0o0 + 1 , OO00Oo ]
     oOoo0 = Qy [ O0o0 , OO00Oo ]
     iIiI = - Qy [ O0o0 , OO00Oo + 1 ]
     if 81 - 81: OoOoOO00 % Ii1I
     oo0 = v [ O0o0 , OO00Oo ] + dt * ( oO00oo000O + ii11 + oOoo0 + iIiI )
     if 16 - 16: Ii1I * OoO0O00 / oO0o
     if oo0 < - 1 :
      print "vNew=" , oo0
      print O0o0 , OO00Oo , o00 , v [ O0o0 , OO00Oo ]
      print iI1I1I11IiII , iIii11iI1II , I1II1I1I , OOo0
      print oO00oo000O * dt , ii11 * dt , oOoo0 * dt , iIiI * dt
      if 22 - 22: oO0o + iIii1I11I1II1 % Oo0Ooo / I11i / Ii1I
      if 54 - 54: OoOoOO00 % IiII . i11iIiiIii
 if sources is not None :
  for oOI11Ii111I in sources :
   v [ oOI11Ii111I [ 0 ] , oOI11Ii111I [ 1 ] ] -= dt * oOI11Ii111I [ 2 ]
   if 93 - 93: ooOoO0o % i11iIiiIii % I1Ii111
 return ooOO
 if 64 - 64: I1Ii111 + I1IiiI * O0 / Oo0Ooo - I11i % I11i
def timeStep ( v , Qx , Qy , dt , sources = None ) :
 o0oO0 , oo00 = v . shape
 if 59 - 59: OOooOOo + OoooooooOO
 if sources is not None :
  for oOI11Ii111I in sources :
   v [ oOI11Ii111I [ 0 ] , oOI11Ii111I [ 1 ] ] += dt * oOI11Ii111I [ 2 ]
   if 55 - 55: i11iIiiIii % iIii1I11I1II1 . i1IIi + OoooooooOO / i11iIiiIii
 for O0o0 in range ( o0oO0 ) :
  for OO00Oo in range ( oo00 ) :
   oO00oo000O = Qx [ O0o0 , OO00Oo ]
   ii11 = - Qx [ O0o0 + 1 , OO00Oo ]
   oOoo0 = Qy [ O0o0 , OO00Oo ]
   iIiI = - Qy [ O0o0 , OO00Oo + 1 ]
   if 10 - 10: iII111i - oO0o * iIii1I11I1II1 % iIii1I11I1II1 * IiII - I1ii11iIi11i
   v [ O0o0 , OO00Oo ] += dt * ( oO00oo000O + ii11 + oOoo0 + iIiI )
   if 97 - 97: II111iiii % I1Ii111 + I1Ii111 - OoO0O00 / Ii1I * I1IiiI
def trappedVol ( sp , cpx , cpy ) :
 o0oO0 , oo00 , iIii1iII1Ii = sp . shape
 if 50 - 50: Ii1I
 I1iiIiI1iI1I = 0.
 if 27 - 27: I1ii11iIi11i * I1Ii111 - OoO0O00 + Ii1I * Ii1I
 for O0o0 in range ( o0oO0 ) :
  for OO00Oo in range ( oo00 ) :
   oOO00O = sp [ O0o0 , OO00Oo , 0 ]
   o0OO0O0OO0oO0 = cpx [ O0o0 , OO00Oo , 0 ]
   iIiiIi11IIi = cpx [ O0o0 + 1 , OO00Oo , 0 ]
   Oo0 = cpy [ O0o0 , OO00Oo , 0 ]
   oOII1ii1ii11I1 = cpy [ O0o0 , OO00Oo + 1 , 0 ]
   if 88 - 88: I1ii11iIi11i
   if oOO00O < min ( o0OO0O0OO0oO0 , iIiiIi11IIi , Oo0 , oOII1ii1ii11I1 ) :
    I1iiIiI1iI1I += volFromWl ( min ( o0OO0O0OO0oO0 , iIiiIi11IIi , Oo0 , oOII1ii1ii11I1 ) , O0o0 , OO00Oo , sp )
    if 93 - 93: iIii1I11I1II1
 return I1iiIiI1iI1I
 if 66 - 66: i11iIiiIii * iIii1I11I1II1 % OoooooooOO
 if 5 - 5: OoOoOO00 % OoooooooOO
def volFromWl ( wl , i , j , sp , dx ) :
 if 60 - 60: OoOoOO00 . i1IIi % OoO0O00 % ooOoO0o % OOooOOo
 if sp . shape [ 2 ] == 5 :
  oOO00O = sp [ i , j , 0 ]
  OOOoo0OO = sp [ i , j , 1 ]
  if 33 - 33: iIii1I11I1II1 - Ii1I * I1ii11iIi11i % iIii1I11I1II1 + OoO0O00 . OOooOOo
  ooo0O0oOoooO0 = sp [ i , j , 2 : ]
  if 42 - 42: OOooOOo % oO0o / OoO0O00 - oO0o * i11iIiiIii
  if wl >= OOOoo0OO :
   iI1IiiiIiI1Ii = dx * dx * ( ooo0O0oOoooO0 [ - 1 ] + ( wl - OOOoo0OO ) )
  else :
   iI1IiiiIiI1Ii = dx * dx * numpy . interp ( wl , [ oOO00O , oOO00O + 1 , oOO00O + 5 , OOOoo0OO ] , numpy . concatenate ( ( [ 0 ] , ooo0O0oOoooO0 ) ) )
   if 78 - 78: OoooooooOO / OOooOOo % OoOoOO00 * OoooooooOO
 else :
  oOO00O = sp [ i , j , 0 ]
  oO0 = sp [ i , j , 1 ]
  OOOoo0OO = sp [ i , j , 2 ]
  if 68 - 68: oO0o
  ooo0O0oOoooO0 = sp [ i , j , 3 : ]
  if 29 - 29: iII111i + i11iIiiIii % I11i
  if wl >= OOOoo0OO :
   iI1IiiiIiI1Ii = ooo0O0oOoooO0 [ - 1 ] + dx * dx * ( wl - OOOoo0OO )
  else :
   iI1IiiiIiI1Ii = dx * dx * numpy . interp ( wl , [ oOO00O , oO0 , oO0 + 1 , oO0 + 5 , OOOoo0OO ] , numpy . concatenate ( ( [ 0 ] , ooo0O0oOoooO0 ) ) )
   if 93 - 93: OoOoOO00 % iIii1I11I1II1
 return iI1IiiiIiI1Ii
 if 90 - 90: I1IiiI - OOooOOo / Ii1I / O0 / I11i
def wlFromVol ( v , i , j , sp , dx ) :
 if 87 - 87: OoOoOO00 / IiII + iIii1I11I1II1
 if len ( sp ) == 5 :
  if 93 - 93: iIii1I11I1II1 + oO0o % ooOoO0o
  oOO00O = sp [ i , j , 0 ]
  OOOoo0OO = sp [ i , j , 1 ]
  iii1IiI1I1 = sp [ i , j , 2 ]
  O00o = sp [ i , j , 3 ]
  oO0o00ooO0OoO = sp [ i , j , 4 ]
  if 1 - 1: OoOoOO00 . i11iIiiIii % OoOoOO00 - iII111i % i1IIi + I1ii11iIi11i
  if v / ( dx * dx ) >= oO0o00ooO0OoO :
   OoOOoOooooOOo = ( v / ( dx * dx ) - oO0o00ooO0OoO ) + OOOoo0OO
  else :
   OoOOoOooooOOo = numpy . interp ( v / ( dx * dx ) , [ 0 , iii1IiI1I1 , O00o , oO0o00ooO0OoO ] , [ oOO00O , oOO00O + 1 , oOO00O + 5 , OOOoo0OO ] )
 else :
  oOO00O = sp [ i , j , 0 ]
  oO0 = sp [ i , j , 1 ]
  OOOoo0OO = sp [ i , j , 2 ]
  if 2 - 2: iIii1I11I1II1 * oO0o / OoOoOO00 . I11i / IiII
  ooo0O0oOoooO0 = sp [ i , j , 3 : ]
  if 75 - 75: OoOoOO00
  if v / ( dx * dx ) >= ooo0O0oOoooO0 [ - 1 ] :
   OoOOoOooooOOo = ( v / ( dx * dx ) - ooo0O0oOoooO0 [ - 1 ] ) + OOOoo0OO
  else :
   OoOOoOooooOOo = numpy . interp ( v / ( dx * dx ) , numpy . concatenate ( ( [ 0 ] , ooo0O0oOoooO0 ) ) , [ oOO00O , oO0 , oO0 + 1 , oO0 + 5 , OOOoo0OO ] )
   if 78 - 78: Ii1I + OoOoOO00 + IiII - IiII . i11iIiiIii / OoO0O00
 return OoOOoOooooOOo
 if 27 - 27: Ii1I - O0 % I11i * I1Ii111 . IiII % iIii1I11I1II1
 if 37 - 37: OoooooooOO + O0 - i1IIi % ooOoO0o
 if 24 - 24: OoOoOO00
 if 94 - 94: i1IIi * i1IIi % II111iiii + OOooOOo
 if 28 - 28: I1IiiI
 if 49 - 49: I11i . o0oOOo0O0Ooo % oO0o / Ii1I
def addVol ( wlg , i , j , sp , dv ) :
 iI1IiiiIiI1Ii = volFromWl ( wlg [ i , j ] , i , j , sp )
 iI1IiiiIiI1Ii += dv
 wlg [ i , j ] = wlFromVol ( iI1IiiiIiI1Ii , i , j , sp )
 if 95 - 95: O0 * OoOoOO00 * IiII . ooOoO0o / iIii1I11I1II1
def resample ( wl , xll , yll , dx , dtmFileName ) :
 o0oO0 , oo00 = wl . shape
 if 28 - 28: IiII + oO0o - ooOoO0o / iIii1I11I1II1 - I1IiiI
 I11iIiII , Iii , iIIiIiI1I1 , ooO , ii , OO0O0Ooo = fileIO . readScalarGrid ( dtmFileName )
 if 45 - 45: O0 / i1IIi * oO0o * OoO0O00
 I11iiI1i1 = numpy . zeros ( ( iIIiIiI1I1 , ooO ) , dtype = arrayType )
 II11I = numpy . zeros ( ( iIIiIiI1I1 , ooO ) , dtype = arrayType ) - 9999.
 if 31 - 31: Ii1I
 for O0o0 in range ( iIIiIiI1I1 ) :
  for OO00Oo in range ( ooO ) :
   i11iIIi = ii + O0o0 * Iii + 0.5 * Iii
   O000O = OO0O0Ooo + OO00Oo * Iii + 0.5 * Iii
   if 16 - 16: I1IiiI . iIii1I11I1II1
   iiiIIIii = int ( ( i11iIIi - xll ) / dx )
   ooOoo00 = int ( ( O000O - yll ) / dx )
   if 34 - 34: II111iiii + ooOoO0o * iIii1I11I1II1 - I1Ii111 - OoO0O00
   if i11iIIi < xll or O000O < yll :
    continue
    if 69 - 69: iIii1I11I1II1 * I1IiiI - iII111i + O0 + O0
   if iiiIIIii >= 0 and iiiIIIii < o0oO0 and ooOoo00 >= 0 and ooOoo00 < oo00 :
    I11iiI1i1 [ O0o0 , OO00Oo ] = max ( 0 , wl [ iiiIIIii , ooOoo00 ] - I11iIiII [ O0o0 , OO00Oo ] )
    if wl [ iiiIIIii , ooOoo00 ] > I11iIiII [ O0o0 , OO00Oo ] :
     II11I [ O0o0 , OO00Oo ] = wl [ iiiIIIii , ooOoo00 ]
   else :
    I11iiI1i1 [ O0o0 , OO00Oo ] = 0.
    if 65 - 65: I1Ii111 / i11iIiiIii / OoO0O00 - OOooOOo
 I11iiI1i1 [ numpy . where ( I11iiI1i1 <= 0 ) ] = - 9999.
 if 9 - 9: I1IiiI / I1Ii111 - Oo0Ooo * iIii1I11I1II1
 return I11iiI1i1 , II11I , ii , OO0O0Ooo , Iii
 if 86 - 86: II111iiii + ooOoO0o + IiII
def nint ( f ) :
 return int ( f + 0.5 )
 if 9 - 9: ooOoO0o + II111iiii % ooOoO0o % IiII + iIii1I11I1II1
def resample2 ( wl , v , xll , yll , dx , dtmFileName , ndv = None , ndr = None ) :
 o0oO0 , oo00 = wl . shape
 if 59 - 59: i1IIi
 I11iIiII , Iii , iIIiIiI1I1 , ooO , ii , OO0O0Ooo = fileIO . readScalarGrid ( dtmFileName )
 if 48 - 48: O0 * Ii1I * OoO0O00 . OoO0O00 * I11i - Ii1I
 if ndv is not None and ndr is not None :
  I11iIiII [ numpy . where ( I11iIiII == ndv ) ] = ndr
  if 14 - 14: I1ii11iIi11i + i11iIiiIii
 I11iIiII [ numpy . where ( numpy . isnan ( I11iIiII ) ) ] = ndr
 if 83 - 83: I1ii11iIi11i / i11iIiiIii + II111iiii . iII111i * OOooOOo + IiII
 I11iiI1i1 = numpy . zeros ( ( iIIiIiI1I1 , ooO ) , dtype = arrayType )
 II11I = numpy . zeros ( ( iIIiIiI1I1 , ooO ) , dtype = arrayType ) - 9999.
 if 42 - 42: i1IIi % II111iiii . ooOoO0o
 IIIII = numpy . zeros ( ( nint ( dx / Iii ) , nint ( dx / Iii ) ) , dtype = arrayType )
 II1II1iI = numpy . zeros ( ( nint ( dx / Iii ) , nint ( dx / Iii ) ) , dtype = arrayType )
 Ooo = numpy . zeros ( ( nint ( dx / Iii ) , nint ( dx / Iii ) ) , dtype = arrayType )
 if 88 - 88: OoooooooOO
 for O0o0 in range ( o0oO0 ) :
  for OO00Oo in range ( oo00 ) :
   if v [ O0o0 , OO00Oo ] < 1e-3 :
    continue
    if 28 - 28: Oo0Ooo * o0oOOo0O0Ooo / I1Ii111
   Oo0OOo00oO0oOO = xll + O0o0 * dx
   iiiii1I1III1 = yll + OO00Oo * dx
   i1oO00O = Oo0OOo00oO0oOO + dx
   oo0o0ooooo = iiiii1I1III1 + dx
   if 70 - 70: Ii1I . i11iIiiIii % Ii1I . O0 - iIii1I11I1II1
   III1Iiii1I11 = int ( ( Oo0OOo00oO0oOO - ii ) / Iii )
   IIII = III1Iiii1I11 + nint ( dx / Iii )
   if 26 - 26: OOooOOo
   o00oooO0Oo = int ( ( iiiii1I1III1 - OO0O0Ooo ) / Iii )
   o0O0OOO0Ooo = o00oooO0Oo + nint ( dx / Iii )
   if 76 - 76: i1IIi * OoooooooOO * O0 + I1Ii111 * I1Ii111
   try :
    IIIII [ : , : ] = I11iIiII [ III1Iiii1I11 : IIII , o00oooO0Oo : o0O0OOO0Ooo ]
   except :
    print Oo0OOo00oO0oOO , i1oO00O
    print iiiii1I1III1 , oo0o0ooooo
    print III1Iiii1I11 , o00oooO0Oo , IIII , o0O0OOO0Ooo
    print
    print ( Oo0OOo00oO0oOO - ii ) / Iii
    print ( i1oO00O - ii ) / Iii
    if 35 - 35: o0oOOo0O0Ooo
    print ( iiiii1I1III1 - OO0O0Ooo ) / Iii
    print ( oo0o0ooooo - OO0O0Ooo ) / Iii
    if 73 - 73: O0 - I1ii11iIi11i
   II1II1iI = wl [ O0o0 , OO00Oo ] - IIIII
   II1II1iI [ numpy . where ( II1II1iI < 0 ) ] = - 9999.
   Ooo [ : , : ] = wl [ O0o0 , OO00Oo ]
   Ooo [ numpy . where ( II1II1iI == - 9999. ) ] = - 9999.
   if 2 - 2: II111iiii / I1Ii111
   I11iiI1i1 [ III1Iiii1I11 : IIII , o00oooO0Oo : o0O0OOO0Ooo ] = II1II1iI
   II11I [ III1Iiii1I11 : IIII , o00oooO0Oo : o0O0OOO0Ooo ] = Ooo
   if 54 - 54: i1IIi . I11i - I1ii11iIi11i + ooOoO0o + Oo0Ooo / Oo0Ooo
   if 22 - 22: ooOoO0o . iIii1I11I1II1
 return I11iiI1i1 , II11I , ii , OO0O0Ooo , Iii
 if 12 - 12: Ii1I
def formatTime ( t ) :
 Oo = int ( t / 86400. )
 ii1IIi1ii = int ( ( t - Oo * 86400. ) / 3600. )
 oo0OoOOooO = int ( ( t - Oo * 86400 - ii1IIi1ii * 3600 ) / 60. )
 o0o0OO0o00o0O = int ( t - Oo * 86400 - ii1IIi1ii * 3600 - oo0OoOOooO * 60. )
 if 28 - 28: OoO0O00 - oO0o + OoOoOO00 + Ii1I / iIii1I11I1II1
 return "%02id:%02ih:%02im:%02is" % ( Oo , ii1IIi1ii , oo0OoOOooO , o0o0OO0o00o0O )
 if 26 - 26: iIii1I11I1II1 - O0 . O0
 if 68 - 68: OOooOOo + oO0o . O0 . Ii1I % i1IIi % OOooOOo
def loadCppLib ( libPath ) :
 if 50 - 50: IiII + o0oOOo0O0Ooo
 if arrayType == numpy . float64 :
  o0OoOOo = ndpointer ( ctypes . c_double )
  O0Oo0 = ctypes . c_double
  iI1II1III1 = libPath
 else :
  o0OoOOo = ndpointer ( ctypes . c_float )
  O0Oo0 = ctypes . c_float
  iI1II1III1 = libPath
  if 62 - 62: I1IiiI * i11iIiiIii . iII111i
 I1iIIIiI = ctypes . cdll . LoadLibrary ( iI1II1III1 )
 if 60 - 60: I1IiiI . i11iIiiIii + OoOoOO00 / I1ii11iIi11i * II111iiii * OOooOOo
 OOO0o0 = I1iIIIiI . sum
 OOO0o0 . restype = O0Oo0
 OOO0o0 . argtypes = [ o0OoOOo , ctypes . c_int , ctypes . c_int ]
 if 34 - 34: I1IiiI % Oo0Ooo - OoOoOO00 + iII111i
 if 79 - 79: II111iiii - ooOoO0o . i1IIi + O0 % O0 * I1IiiI
 Ii1Ii1I = I1iIIIiI . calcFlow
 Ii1Ii1I . restype = O0Oo0
 Ii1Ii1I . argtypes = [ O0Oo0 , O0Oo0 , O0Oo0 , O0Oo0 , O0Oo0 , O0Oo0 , O0Oo0 ]
 if 54 - 54: oO0o + OoOoOO00
 if 77 - 77: I11i . IiII
 o0O0OO0OOOOOo = I1iIIIiI . calcFlowGrid
 o0O0OO0OOOOOo . argtypes = [ o0OoOOo , o0OoOOo , o0OoOOo , O0Oo0 , O0Oo0 , ctypes . c_int , ctypes . c_int ,
 # OOooOOo
 # I1IiiI + iII111i
 # i1IIi
 o0OoOOo , o0OoOOo , o0OoOOo ]
 if 42 - 42: II111iiii - OoO0O00 - OoooooooOO . iII111i / OoOoOO00
 if 56 - 56: i11iIiiIii - iIii1I11I1II1 . II111iiii
 if 81 - 81: IiII / OoOoOO00 * IiII . O0
 OOOOo00oo00O = I1iIIIiI . dryCheck
 OOOOo00oo00O . restype = ctypes . c_int
 OOOOo00oo00O . argtypes = [ o0OoOOo , o0OoOOo , o0OoOOo , O0Oo0 , ctypes . c_int , ctypes . c_int , ndpointer ( ctypes . c_int ) , ndpointer ( ctypes . c_int ) , o0OoOOo , ctypes . c_int ]
 if 83 - 83: II111iiii * i1IIi * iII111i . I1ii11iIi11i / I11i + i1IIi
 if 43 - 43: OoooooooOO
 if 97 - 97: I1ii11iIi11i / Oo0Ooo + I1Ii111
 if 32 - 32: ooOoO0o % I1Ii111 * Oo0Ooo
 if 72 - 72: ooOoO0o . iII111i - I1Ii111 - Ii1I % i1IIi
 oO0o00O0O0oo0 = I1iIIIiI . dryCheckDiagnostic
 oO0o00O0O0oo0 . restype = ctypes . c_int
 oO0o00O0O0oo0 . argtypes = [ o0OoOOo , o0OoOOo , o0OoOOo , O0Oo0 , ctypes . c_int , ctypes . c_int , ndpointer ( ctypes . c_int ) , ndpointer ( ctypes . c_int ) , o0OoOOo , ctypes . c_int , o0OoOOo ]
 if 24 - 24: I1Ii111 * oO0o
 if 88 - 88: i11iIiiIii + iII111i * OoOoOO00 * iII111i + I11i
 if 88 - 88: OOooOOo % Oo0Ooo - iII111i - OoOoOO00 % i11iIiiIii
 if 6 - 6: Ii1I - OoO0O00 . I1IiiI - O0
 if 16 - 16: iII111i * iII111i % Ii1I % I1IiiI
 if 48 - 48: OOooOOo / Ii1I % OoO0O00 / IiII / I1Ii111
 o0OO00o0oOOoo = I1iIIIiI . timeStep
 o0OO00o0oOOoo . argtypes = [ o0OoOOo , o0OoOOo , o0OoOOo ,
 # iIii1I11I1II1 * OoO0O00 % I1Ii111
 O0Oo0 , ctypes . c_int , ctypes . c_int ,
 ndpointer ( ctypes . c_int ) , ndpointer ( ctypes . c_int ) , o0OoOOo , ctypes . c_int ]
 if 46 - 46: I11i . IiII / II111iiii % iIii1I11I1II1 + IiII
 if 61 - 61: OOooOOo / OoO0O00 + II111iiii . oO0o / Oo0Ooo * OOooOOo
 ii1O0ooooo000 = I1iIIIiI . wlFromVolGrid
 ii1O0ooooo000 . argtypes = [ o0OoOOo , o0OoOOo , o0OoOOo , ctypes . c_int , ctypes . c_int , ctypes . c_bool , O0Oo0 ]
 if 97 - 97: i11iIiiIii % oO0o / Oo0Ooo / Oo0Ooo
 if 97 - 97: II111iiii - I1Ii111 - iIii1I11I1II1 * I1IiiI
 if 54 - 54: iIii1I11I1II1
 i111i1I1ii1i = I1iIIIiI . conveyanceParameters
 i111i1I1ii1i . argtypes = [ ctypes . c_int , ctypes . c_int , ctypes . c_int , ctypes . c_int , o0OoOOo , O0Oo0 , ctypes . c_int , ctypes . c_int , O0Oo0 , o0OoOOo , ctypes . c_bool , O0Oo0 , o0OoOOo , ctypes . c_int , ctypes . c_int , ctypes . c_int , ctypes . c_int , ctypes . c_int , ctypes . c_int , O0Oo0 , O0Oo0 , O0Oo0 , O0Oo0 , o0OoOOo ]
 if 100 - 100: IiII . Ii1I - iIii1I11I1II1 . i11iIiiIii / II111iiii
 if 71 - 71: I1Ii111 * Oo0Ooo . I11i
 if 49 - 49: IiII * O0 . IiII
 if 19 - 19: II111iiii - IiII
 if 59 - 59: o0oOOo0O0Ooo * OoO0O00 - Ii1I . OOooOOo
 if 89 - 89: OOooOOo
 if 69 - 69: ooOoO0o - OoooooooOO * O0
 if 84 - 84: ooOoO0o + i11iIiiIii - OOooOOo * ooOoO0o
 I1IiiIiii1 = I1iIIIiI . maxVolGrid
 I1IiiIiii1 . argtypes = [ o0OoOOo , o0OoOOo , o0OoOOo , o0OoOOo , o0OoOOo , o0OoOOo , ctypes . c_int , ctypes . c_int ]
 if 39 - 39: ooOoO0o / O0 * IiII
 if 17 - 17: Ii1I / iIii1I11I1II1 - OoO0O00 + I1IiiI % OOooOOo
 if 14 - 14: o0oOOo0O0Ooo % IiII + I1ii11iIi11i + OoO0O00
 if 76 - 76: OoO0O00 - i11iIiiIii + OoOoOO00 + OOooOOo / OoooooooOO
 if 50 - 50: II111iiii - I1Ii111 + iIii1I11I1II1 + iIii1I11I1II1
 if 91 - 91: II111iiii - O0 . iIii1I11I1II1 . O0 + I1ii11iIi11i - II111iiii
 if 26 - 26: o0oOOo0O0Ooo
 if 12 - 12: OoooooooOO / O0 + II111iiii * I1ii11iIi11i
 Ii11ii1I1 = I1iIIIiI . resample2
 Ii11ii1I1 . argtypes = [ o0OoOOo , o0OoOOo , O0Oo0 , O0Oo0 , O0Oo0 , ctypes . c_int , ctypes . c_int , o0OoOOo , O0Oo0 , ctypes . c_int , ctypes . c_int , O0Oo0 , O0Oo0 , o0OoOOo , o0OoOOo ]
 if 11 - 11: iIii1I11I1II1 . OoOoOO00 / IiII % ooOoO0o
 if 61 - 61: ooOoO0o - OOooOOo + OOooOOo
 if 40 - 40: i11iIiiIii . iIii1I11I1II1
 if 2 - 2: i1IIi * oO0o - oO0o + OoooooooOO % OoOoOO00 / OoOoOO00
 if 3 - 3: OoooooooOO
 if 71 - 71: IiII + i1IIi - iII111i - i11iIiiIii . I11i - ooOoO0o
 if 85 - 85: I1ii11iIi11i - OoOoOO00 / I1ii11iIi11i + OOooOOo - iII111i
 if 49 - 49: OoO0O00 - O0 / OoO0O00 * OoOoOO00 + I1Ii111
 Ii = I1iIIIiI . calcStorageParameters
 Ii . argtypes = [ O0Oo0 , O0Oo0 , O0Oo0 , O0Oo0 , o0OoOOo , ctypes . c_int , ctypes . c_int , O0Oo0 , O0Oo0 , O0Oo0 , o0OoOOo , ctypes . c_bool , o0OoOOo , ctypes . c_int , ctypes . c_int , O0Oo0 , O0Oo0 , O0Oo0 , O0Oo0 , O0Oo0 , O0Oo0 , O0Oo0 , O0Oo0 ]
 if 20 - 20: i1IIi / I1IiiI * oO0o
 if 85 - 85: II111iiii . ooOoO0o % OOooOOo % I11i
 if 80 - 80: oO0o * I11i / iIii1I11I1II1 % oO0o / iIii1I11I1II1
 if 42 - 42: i1IIi / i11iIiiIii . Oo0Ooo * iII111i . i11iIiiIii * O0
 if 44 - 44: i1IIi . I1IiiI / i11iIiiIii + IiII
 if 27 - 27: OOooOOo
 if 52 - 52: I1Ii111 % OoOoOO00 + iIii1I11I1II1 * oO0o . Ii1I
 if 95 - 95: iIii1I11I1II1 . IiII - OoooooooOO * OoO0O00 / o0oOOo0O0Ooo
 oOo0OO0o0 = I1iIIIiI . flowPaths
 oOo0OO0o0 . argtypes = [ O0Oo0 , O0Oo0 , O0Oo0 , ctypes . c_int , ctypes . c_int , o0OoOOo , O0Oo0 , ctypes . c_int , ctypes . c_int , O0Oo0 , O0Oo0 , ctypes . c_char_p , o0OoOOo , o0OoOOo , ctypes . c_int ]
 if 35 - 35: Oo0Ooo . Oo0Ooo % OoooooooOO - Ii1I
 if 43 - 43: OoO0O00 % OoO0O00
 if 46 - 46: Oo0Ooo % iIii1I11I1II1 . iII111i . O0 * ooOoO0o / OoooooooOO
 if 7 - 7: oO0o - O0 * I11i - o0oOOo0O0Ooo - II111iiii
 if 41 - 41: I1IiiI - I1Ii111 % II111iiii . I1Ii111 - I11i
 if 45 - 45: Ii1I - OOooOOo
 if 70 - 70: OoO0O00 % I1IiiI / I1IiiI . I11i % ooOoO0o . II111iiii
 if 10 - 10: Ii1I - i11iIiiIii . I1ii11iIi11i % i1IIi
 OooOOOoOoo0O0 = I1iIIIiI . resample3
 OooOOOoOoo0O0 . argtypes = [ o0OoOOo , o0OoOOo , o0OoOOo , o0OoOOo , o0OoOOo , o0OoOOo , O0Oo0 , O0Oo0 , O0Oo0 , O0Oo0 , ctypes . c_int , ctypes . c_int , o0OoOOo , O0Oo0 , ctypes . c_int , ctypes . c_int , O0Oo0 , O0Oo0 , o0OoOOo , o0OoOOo ]
 if 81 - 81: IiII - o0oOOo0O0Ooo - Oo0Ooo - Ii1I / OOooOOo % I11i
 if 52 - 52: I1ii11iIi11i / iII111i
 if 37 - 37: I11i
 if 83 - 83: O0
 if 89 - 89: Oo0Ooo + I1ii11iIi11i - o0oOOo0O0Ooo
 if 40 - 40: OoO0O00 + OoO0O00
 if 94 - 94: iII111i * iIii1I11I1II1 . I11i
 if 13 - 13: iIii1I11I1II1 * OoOoOO00 / I1Ii111 % ooOoO0o + oO0o
 if 41 - 41: I1ii11iIi11i
 if 5 - 5: Oo0Ooo
 if 100 - 100: Ii1I + iIii1I11I1II1
 if 59 - 59: IiII
 oOoO0OOO00O = I1iIIIiI . lazyFlowPaths
 oOoO0OOO00O . argtypes = [ O0Oo0 , O0Oo0 , O0Oo0 ,
 # Ii1I - OoOoOO00 - OoO0O00 + I1IiiI * Ii1I
 ctypes . c_int , ctypes . c_int ,
 o0OoOOo , O0Oo0 , ctypes . c_int , ctypes . c_int ,
 O0Oo0 , O0Oo0 ,
 o0OoOOo , o0OoOOo , o0OoOOo , O0Oo0 , O0Oo0 , O0Oo0 ]
 if 67 - 67: Oo0Ooo / ooOoO0o - IiII
 O0O00OoOoOOo = I1iIIIiI . fillWlGrid
 O0O00OoOoOOo . argtypes = [ o0OoOOo , o0OoOOo , ctypes . c_int , ctypes . c_int , O0Oo0 , ctypes . c_int , ctypes . c_int , O0Oo0 ]
 if 58 - 58: IiII + iIii1I11I1II1
 if 94 - 94: Ii1I . i1IIi
 if 71 - 71: iII111i + OoO0O00 - IiII . OoO0O00 . IiII + I1IiiI
 iii = I1iIIIiI . burnFlowPaths
 iii . argtypes = [ o0OoOOo , o0OoOOo , ctypes . c_int , ctypes . c_int ]
 if 26 - 26: I1IiiI
 iiiiIiIiI = I1iIIIiI . makeWlGrid
 iiiiIiIiI . argtypes = [ o0OoOOo , o0OoOOo , o0OoOOo , ctypes . c_int , ctypes . c_int , O0Oo0 ]
 if 83 - 83: i1IIi
 if 76 - 76: Ii1I + iIii1I11I1II1 + OoOoOO00 . OoO0O00
 i1i1 = I1iIIIiI . clipZero
 i1i1 . argtypes = [ o0OoOOo , ctypes . c_int , ctypes . c_int ]
 if 68 - 68: Ii1I - I1IiiI
 if 41 - 41: oO0o
 I11II1 = I1iIIIiI . scsAdditionalRunoff
 I11II1 . argtypes = [ o0OoOOo , o0OoOOo , o0OoOOo , o0OoOOo , ctypes . c_int , ctypes . c_int ]
 if 46 - 46: OoOoOO00
 if 83 - 83: i11iIiiIii * I1Ii111
 if 49 - 49: Oo0Ooo * oO0o + o0oOOo0O0Ooo - i11iIiiIii
 OOooO = I1iIIIiI . calcFlowEdges
 OOooO . restype = O0Oo0
 OOooO . argtypes = [ o0OoOOo , o0OoOOo , o0OoOOo , o0OoOOo , O0Oo0 , ctypes . c_int , ctypes . c_int , o0OoOOo , o0OoOOo , o0OoOOo ]
 if 28 - 28: ooOoO0o - OOooOOo / I1IiiI
 if 27 - 27: i1IIi + I1IiiI * I1ii11iIi11i + OOooOOo . oO0o
 i1I111I1Iii1 = I1iIIIiI . checkLicence
 return Ii1Ii1I , o0O0OO0OOOOOo , OOOOo00oo00O , o0OO00o0oOOoo , ii1O0ooooo000 , i111i1I1ii1i , I1IiiIiii1 , Ii11ii1I1 , OooOOOoOoo0O0 , oOo0OO0o0 , OOO0o0 , Ii , oOoO0OOO00O , O0O00OoOoOOo , iii , iiiiIiIiI , i1i1 , oO0o00O0O0oo0 , I11II1 , OOooO , i1I111I1Iii1
 if 68 - 68: iII111i + Oo0Ooo % Ii1I / i11iIiiIii % OoOoOO00
 if 94 - 94: i11iIiiIii / I1Ii111 / Oo0Ooo
 if 9 - 9: I11i / OoOoOO00 / II111iiii + I1Ii111
 if 71 - 71: iII111i / Oo0Ooo
 if 87 - 87: I1ii11iIi11i + I1ii11iIi11i - I1ii11iIi11i % O0
 if 13 - 13: II111iiii
 if 57 - 57: Ii1I - OoooooooOO
 if 68 - 68: o0oOOo0O0Ooo % I1ii11iIi11i / I1Ii111 + I1Ii111 - I1Ii111 . OoO0O00
def Ii1I1I11I ( topoProfile , dxt , width , depth , nChan , nFP ) :
 oOO00ooOOo = int ( width / dxt )
 i11ii1iI = min ( topoProfile )
 if 20 - 20: I1ii11iIi11i
 if 3 - 3: OoO0O00 * i1IIi . I1IiiI . O0 - OoOoOO00
 if 81 - 81: I1IiiI - iIii1I11I1II1 / I1IiiI / O0
 topoProfile . sort ( )
 if 34 - 34: Ii1I * Ii1I - I1ii11iIi11i - O0 . i11iIiiIii
 i11 = [ nFP ] * len ( topoProfile )
 if 32 - 32: iIii1I11I1II1 . OoO0O00 * oO0o / OOooOOo . II111iiii - Oo0Ooo
 if 10 - 10: I1ii11iIi11i / i11iIiiIii - Ii1I + oO0o * I1IiiI
 for O0o0 in range ( oOO00ooOOo ) :
  topoProfile [ O0o0 ] = i11ii1iI - depth
  i11 [ O0o0 ] = nChan
  if 94 - 94: I1IiiI + iIii1I11I1II1 / O0 - OoooooooOO % I1ii11iIi11i
  if 64 - 64: I11i + OoO0O00
  if 25 - 25: I1IiiI . ooOoO0o + I1IiiI % Ii1I * iIii1I11I1II1
  if 31 - 31: i11iIiiIii + OOooOOo - O0
 return i11
 if 51 - 51: OoO0O00 * i1IIi / Ii1I * OOooOOo + ooOoO0o % I1ii11iIi11i
 if 34 - 34: oO0o * OoooooooOO + Ii1I + i11iIiiIii
 if 22 - 22: i1IIi
 if 24 - 24: I11i / I1IiiI * i1IIi % OoooooooOO
def saveResults ( volGrid , wlGrid , flowX , flowY , storagePar , xsz , ysz , cellSize , xll , yll ,
 defaultDepth , flowThreshold , channel , flowPathOutput ,
 dtmFileName , noDataValue , noDataReplacement ,
 outputDirectory , outputPrefix ,
 resampleFunction , lfpFunction ,
 extendWlGrid , cppBurnFlowPaths , cppMakeWlGrid , cppWlFill , cppClipZero ,
 saveCsv = True , zeroPolyList = None ) :
 if 99 - 99: i11iIiiIii . II111iiii . OoooooooOO
 if flowPathOutput is not None :
  defaultDepth = arrayType ( defaultDepth )
  flowThreshold = arrayType ( flowThreshold )
  Ooi1IIii11i1I1 = arrayType ( 0.25 )
  if 12 - 12: i1IIi / OOooOOo % ooOoO0o * IiII * O0 * iIii1I11I1II1
  if 93 - 93: Oo0Ooo / I1ii11iIi11i + i1IIi * oO0o . OoooooooOO
  if 54 - 54: O0 / IiII % ooOoO0o * i1IIi * O0
 if zeroPolyList is not None :
  IIOOOoO00O = 1e20
  iIiii1Ii1I1II = - 1e20
  iIIIIII = 1e20
  IIiiI = - 1e20
  if 36 - 36: iII111i
  for O0ooooooo00 in zeroPolyList :
   I1111ii11IIII = O0ooooooo00 . bounds
   IIOOOoO00O = min ( I1111ii11IIII [ 0 ] , IIOOOoO00O )
   iIIIIII = min ( I1111ii11IIII [ 1 ] , iIIIIII )
   if 48 - 48: iIii1I11I1II1 % i1IIi + OoOoOO00 % o0oOOo0O0Ooo
   iIiii1Ii1I1II = max ( I1111ii11IIII [ 2 ] , iIiii1Ii1I1II )
   IIiiI = max ( I1111ii11IIII [ 3 ] , IIiiI )
   if 79 - 79: OoOoOO00 % I1IiiI % Ii1I / i1IIi % OoO0O00
   oo0o00OO = tempfile . _get_candidate_names ( ) . next ( ) + '.csv'
   if 69 - 69: o0oOOo0O0Ooo % i11iIiiIii / Ii1I
   ooOOO00oOOooO = open ( oo0o00OO , "w" )
   ooOOO00oOOooO . write ( "id;wkt\n" )
   OO0 = 0
   if 46 - 46: iIii1I11I1II1 . i11iIiiIii - OoOoOO00 % O0 / II111iiii * i1IIi
   for O0ooooooo00 in zeroPolyList :
    ooOOO00oOOooO . write ( "%i;%s\n" % ( OO0 , O0ooooooo00 . wkt ) )
    OO0 += 1
    if 66 - 66: O0
   ooOOO00oOOooO . close ( )
   if 52 - 52: OoO0O00 * OoooooooOO
 Ii11iiI = os . path . join ( outputDirectory , outputPrefix )
 if 71 - 71: I1Ii111 - o0oOOo0O0Ooo - OOooOOo
 if 28 - 28: iIii1I11I1II1
 if 7 - 7: o0oOOo0O0Ooo % IiII * OoOoOO00
 if 58 - 58: IiII / I11i + II111iiii % iII111i - OoooooooOO
 if 25 - 25: OoOoOO00 % OoooooooOO * Oo0Ooo - i1IIi * II111iiii * oO0o
 I1iI1I1ii1 = 5 ;
 iIIi1 = 5 ;
 o0Ooo0o0Oo = 0.01 ;
 if 55 - 55: iIii1I11I1II1 * iII111i
 oOO00O = numpy . zeros ( ( xsz , ysz ) , dtype = numpy . float32 , order = 'C' )
 OOOoo0OO = numpy . zeros ( ( xsz , ysz ) , dtype = numpy . float32 , order = 'C' )
 if 85 - 85: iIii1I11I1II1 . II111iiii
 if channel :
  oOO00O [ : , : ] = storagePar [ : , : , 1 ]
  OOOoo0OO [ : , : ] = storagePar [ : , : , 2 ]
 else :
  oOO00O [ : , : ] = storagePar [ : , : , 0 ]
  OOOoo0OO [ : , : ] = storagePar [ : , : , 1 ]
  if 54 - 54: Ii1I . OoooooooOO % Oo0Ooo
 ii111I11Ii , i11IiiI1Ii1 , iIIiIiI1I1 , ooO , I1iiIiiIiiI , oOoO = fileIO . readScalarGridObj ( dtmFileName )
 if 32 - 32: O0 + oO0o % Oo0Ooo
 iI1iI = 5000
 O0O0 = int ( iIIiIiI1I1 / iI1iI ) + 1
 O0oO0o0OOOOOO = int ( ooO / iI1iI ) + 1
 if 24 - 24: i1IIi * IiII - I11i / Ii1I
 for ooooo0 in range ( O0O0 ) :
  for OooO0O0Ooo in range ( O0oO0o0OOOOOO ) :
   if 85 - 85: o0oOOo0O0Ooo / I1Ii111
   print "Processing tile %i/%i,%i/%i" % ( ooooo0 + 1 , O0O0 , OooO0O0Ooo + 1 , O0oO0o0OOOOOO )
   if 67 - 67: I11i % oO0o
   O0OOO0OOoO0O = I1iiIiiIiiI + ooooo0 * iI1iI * i11IiiI1Ii1
   O00Oo000ooO0 = min ( O0OOO0OOoO0O + iI1iI * i11IiiI1Ii1 , O0OOO0OOoO0O + i11IiiI1Ii1 * iIIiIiI1I1 )
   if 39 - 39: i11iIiiIii + IiII
   o0 = oOoO + ooO * i11IiiI1Ii1 - OooO0O0Ooo * iI1iI * i11IiiI1Ii1
   OoO0O00IIiII = max ( o0 - iI1iI * i11IiiI1Ii1 , oOoO )
   if 7 - 7: iIii1I11I1II1 - i1IIi
   I1ii1i1iiii = ooooo0 * iI1iI
   I1i1I = OooO0O0Ooo * iI1iI
   if 17 - 17: I11i - iII111i % Ii1I
   i11Ii1iIIIIi = min ( iI1iI , iIIiIiI1I1 - ooooo0 * iI1iI )
   iiiO000OOO = min ( iI1iI , ooO - OooO0O0Ooo * iI1iI )
   if 59 - 59: I1Ii111 . I1ii11iIi11i + OoooooooOO
   i1II11I11ii1 = ii111I11Ii . ReadAsArray ( xoff = I1ii1i1iiii , yoff = I1i1I , xsize = i11Ii1iIIIIi , ysize = iiiO000OOO ) . transpose ( ) . copy ( )
   i1II11I11ii1 = numpy . array ( i1II11I11ii1 [ : , : : - 1 ] , arrayType )
   if 64 - 64: oO0o % OoOoOO00 / II111iiii % ooOoO0o - iII111i
   if noDataValue is not None and noDataReplacement is not None :
    replaceArrayVals ( i1II11I11ii1 , arrayType ( noDataValue ) , noDataReplacement )
    if 2 - 2: I1Ii111 - I1ii11iIi11i + o0oOOo0O0Ooo * OoO0O00 / iII111i
   iIIiI11iI1Ii1 = numpy . zeros ( ( i11Ii1iIIIIi , iiiO000OOO ) , dtype = arrayType ) - 9999.
   o00oo = numpy . zeros ( ( i11Ii1iIIIIi , iiiO000OOO ) , dtype = arrayType ) - 9999.
   if 70 - 70: I11i - Oo0Ooo / OoooooooOO % OoooooooOO
   resampleFunction ( wlGrid , volGrid , flowX , flowY , oOO00O , OOOoo0OO , flowThreshold , xll , yll , cellSize , xsz , ysz , i1II11I11ii1 , i11IiiI1Ii1 , i11Ii1iIIIIi , iiiO000OOO , O0OOO0OOoO0O , OoO0O00IIiII , iIIiI11iI1Ii1 , o00oo )
   if 95 - 95: OoooooooOO % OoooooooOO . Ii1I
   if 26 - 26: oO0o + IiII - II111iiii . II111iiii + I1ii11iIi11i + OoOoOO00
   if 68 - 68: O0
   if 76 - 76: I1ii11iIi11i
   ooO000OO = "tile%02i%02i" % ( ooooo0 , OooO0O0Ooo )
   if 43 - 43: ooOoO0o * I1Ii111 % OOooOOo
   fileIO . saveScalarGrid ( o00oo , O0OOO0OOoO0O , OoO0O00IIiII , i11IiiI1Ii1 , Ii11iiI + "_depth_" + ooO000OO + ".tif" )
   if 38 - 38: Oo0Ooo
   if 34 - 34: OoOoOO00
   if 70 - 70: iIii1I11I1II1 * IiII - OOooOOo / Oo0Ooo % oO0o
   if 66 - 66: OoooooooOO + ooOoO0o * iII111i
   if zeroPolyList is not None and IIOOOoO00O < O00Oo000ooO0 and iIiii1Ii1I1II > O0OOO0OOoO0O and iIIIIII < o0 and IIiiI > OoO0O00IIiII :
    if 2 - 2: iII111i . OoO0O00 / oO0o
    print "Zeroing polys..."
    reservoirs . maskZeroPoly ( oo0o00OO , os . path . join ( outputDirectory , outputPrefix + "_depth_" + ooO000OO + ".tif" ) )
    if 41 - 41: OoO0O00 . I1Ii111 * IiII * I1Ii111
    if 74 - 74: iIii1I11I1II1 / o0oOOo0O0Ooo
    if 58 - 58: iIii1I11I1II1 - I1IiiI % o0oOOo0O0Ooo % OoooooooOO * iIii1I11I1II1 + OOooOOo
   fileIO . saveScalarGrid ( iIIiI11iI1Ii1 , O0OOO0OOoO0O , OoO0O00IIiII , i11IiiI1Ii1 , Ii11iiI + "_wl_" + ooO000OO + ".tif" )
   if 25 - 25: OOooOOo % O0
   if 44 - 44: I1Ii111 . Ii1I * II111iiii / IiII + iIii1I11I1II1
   Ii1111III1 = numpy . zeros ( ( i11Ii1iIIIIi , iiiO000OOO ) , dtype = arrayType )
   if 74 - 74: I1ii11iIi11i - iII111i * i1IIi
   lfpFunction ( xll , yll , cellSize , xsz , ysz , i1II11I11ii1 , i11IiiI1Ii1 , i11Ii1iIIIIi , iiiO000OOO , O0OOO0OOoO0O , OoO0O00IIiII , Ii1111III1 , flowX , flowY , flowThreshold , defaultDepth , Ooi1IIii11i1I1 )
   if 12 - 12: O0
   if 75 - 75: iIii1I11I1II1 % IiII + I1ii11iIi11i * O0 . iII111i - ooOoO0o
   fileIO . saveScalarGrid ( Ii1111III1 , O0OOO0OOoO0O , OoO0O00IIiII , i11IiiI1Ii1 , Ii11iiI + "_lfp_" + ooO000OO + ".tif" )
   if 32 - 32: Ii1I % oO0o - i1IIi
   if 40 - 40: iIii1I11I1II1 + iII111i * OoOoOO00 + oO0o
   if 15 - 15: I11i % I1IiiI - iIii1I11I1II1 * ooOoO0o
   oO0O0o0o000 = numpy . where ( ( o00oo < defaultDepth ) & ( Ii1111III1 > 0 ) )
   o00oo [ oO0O0o0o000 ] = defaultDepth
   if 6 - 6: OoOoOO00 / ooOoO0o + iII111i - o0oOOo0O0Ooo * OOooOOo + ooOoO0o
   fileIO . saveScalarGrid ( o00oo , O0OOO0OOoO0O , OoO0O00IIiII , i11IiiI1Ii1 , Ii11iiI + "_merge_" + ooO000OO + ".tif" )
   if 76 - 76: II111iiii - OoooooooOO % IiII
   if 40 - 40: Ii1I
   if 59 - 59: I11i * OoooooooOO + OOooOOo . iIii1I11I1II1 / i1IIi
   if extendWlGrid :
    if 75 - 75: I11i . OOooOOo - iIii1I11I1II1 * OoO0O00 * iII111i
    if 93 - 93: ooOoO0o
    if 18 - 18: ooOoO0o
    if 66 - 66: oO0o * i11iIiiIii + OoOoOO00 / OOooOOo
    if 96 - 96: OOooOOo + OOooOOo % IiII % OOooOOo
    if 28 - 28: iIii1I11I1II1 + OoOoOO00 . o0oOOo0O0Ooo % i11iIiiIii
    if Ii1111III1 is not None :
     cppBurnFlowPaths ( o00oo , Ii1111III1 , i11Ii1iIIIIi , iiiO000OOO )
     if 58 - 58: I11i / OoooooooOO % oO0o + OoO0O00
    iIIiI11iI1Ii1 [ : , : ] = 0
    if 58 - 58: O0
    cppMakeWlGrid ( i1II11I11ii1 , o00oo , iIIiI11iI1Ii1 , iIIiIiI1I1 , ooO , 0.1 )
    if 91 - 91: iII111i / I1ii11iIi11i . iII111i - o0oOOo0O0Ooo + I1ii11iIi11i
    o00oo [ : , : ] = 0
    if 72 - 72: Ii1I . IiII * I1ii11iIi11i / I1ii11iIi11i / iII111i
    if 13 - 13: i1IIi
    if I1iI1I1ii1 is None :
     I1iI1I1ii1 = int ( 0.5 * cellSize / i11IiiI1Ii1 )
     if 17 - 17: i11iIiiIii * o0oOOo0O0Ooo * o0oOOo0O0Ooo + OoO0O00
    if iIIi1 is None :
     iIIi1 = 5
     if 95 - 95: I1IiiI
    if o0Ooo0o0Oo is None :
     o0Ooo0o0Oo = 2. * defaultDepth / cellSize
     if 95 - 95: OOooOOo % I1ii11iIi11i + o0oOOo0O0Ooo % ooOoO0o
    Ii1i = numpy . zeros ( ( i11Ii1iIIIIi , iiiO000OOO ) , dtype = arrayType ) - 9999.
    i1IiIiIiiI1 = numpy . zeros ( ( i11Ii1iIIIIi , iiiO000OOO ) , dtype = arrayType ) - 9999.
    if 41 - 41: II111iiii
    cppWlFill ( iIIiI11iI1Ii1 , Ii1i , iIIiIiI1I1 , ooO , i11IiiI1Ii1 , I1iI1I1ii1 , iIIi1 , o0Ooo0o0Oo )
    if 43 - 43: O0 - ooOoO0o % OoooooooOO % OOooOOo + iII111i
    fileIO . saveScalarGrid ( Ii1i , I1iiIiiIiiI , oOoO , i11IiiI1Ii1 , Ii11iiI + "_wl_fill_" + ooO000OO + ".tif" )
    if 61 - 61: ooOoO0o . i11iIiiIii + oO0o
    if 8 - 8: iIii1I11I1II1
    if 55 - 55: oO0o
    i1IiIiIiiI1 [ : , : ] = Ii1i - i1II11I11ii1
    if 37 - 37: IiII / i11iIiiIii / Oo0Ooo
    cppClipZero ( i1IiIiIiiI1 , iIIiIiI1I1 , ooO )
    if 97 - 97: I1Ii111 . I11i / I1IiiI
    fileIO . saveScalarGrid ( i1IiIiIiiI1 , I1iiIiiIiiI , oOoO , i11IiiI1Ii1 , Ii11iiI + "_depth_fill_" + ooO000OO + ".tif" )
    if 83 - 83: I11i - I1ii11iIi11i * oO0o
    if 90 - 90: Oo0Ooo * I1IiiI
    if 75 - 75: I1ii11iIi11i - OoOoOO00 * i11iIiiIii . OoooooooOO - Oo0Ooo . I11i
    if 6 - 6: I11i * oO0o / OoooooooOO % Ii1I * o0oOOo0O0Ooo
 print "Generating VRTs..." ,
 if 28 - 28: IiII * I1IiiI % IiII
 if 95 - 95: O0 / I11i . I1Ii111
 iII11II1II = [ 'gdalbuildvrt' ]
 iII11II1II . append ( Ii11iiI + '_depth.vrt' )
 if 100 - 100: OoO0O00 % I1Ii111 - I11i % I11i % I11i / ooOoO0o
 for ooooo0 in range ( O0O0 ) :
  for OooO0O0Ooo in range ( O0oO0o0OOOOOO ) :
   ooO000OO = "tile%02i%02i" % ( ooooo0 , OooO0O0Ooo )
   iII11II1II . append ( Ii11iiI + "_depth_" + ooO000OO + ".tif" )
   if 83 - 83: oO0o - ooOoO0o - IiII % i1IIi - iII111i . o0oOOo0O0Ooo
 call ( iII11II1II )
 if 96 - 96: Oo0Ooo + I1Ii111 . i1IIi
 if 54 - 54: II111iiii . i1IIi / I1ii11iIi11i % I1IiiI / I1Ii111
 iII11II1II = [ 'gdalbuildvrt' ]
 iII11II1II . append ( Ii11iiI + '_wl.vrt' )
 if 65 - 65: OoOoOO00 . OoOoOO00 - oO0o + Oo0Ooo / i11iIiiIii
 for ooooo0 in range ( O0O0 ) :
  for OooO0O0Ooo in range ( O0oO0o0OOOOOO ) :
   ooO000OO = "tile%02i%02i" % ( ooooo0 , OooO0O0Ooo )
   iII11II1II . append ( Ii11iiI + "_wl_" + ooO000OO + ".tif" )
   if 90 - 90: iIii1I11I1II1 + OoOoOO00
 call ( iII11II1II )
 if 9 - 9: iIii1I11I1II1 . OoooooooOO + i1IIi - Oo0Ooo
 if 30 - 30: iII111i / OoO0O00 . iII111i
 iII11II1II = [ 'gdalbuildvrt' ]
 iII11II1II . append ( Ii11iiI + '_lfp.vrt' )
 if 17 - 17: Oo0Ooo + OoooooooOO * OoooooooOO
 for ooooo0 in range ( O0O0 ) :
  for OooO0O0Ooo in range ( O0oO0o0OOOOOO ) :
   ooO000OO = "tile%02i%02i" % ( ooooo0 , OooO0O0Ooo )
   iII11II1II . append ( Ii11iiI + "_lfp_" + ooO000OO + ".tif" )
   if 5 - 5: I1Ii111 % OoooooooOO . OoOoOO00
 call ( iII11II1II )
 if 67 - 67: I1ii11iIi11i + Ii1I
 if 72 - 72: IiII % o0oOOo0O0Ooo
 iII11II1II = [ 'gdalbuildvrt' ]
 iII11II1II . append ( Ii11iiI + '_merge.vrt' )
 if 93 - 93: iIii1I11I1II1 + i11iIiiIii . o0oOOo0O0Ooo . i1IIi % I1IiiI % ooOoO0o
 for ooooo0 in range ( O0O0 ) :
  for OooO0O0Ooo in range ( O0oO0o0OOOOOO ) :
   ooO000OO = "tile%02i%02i" % ( ooooo0 , OooO0O0Ooo )
   iII11II1II . append ( Ii11iiI + "_merge_" + ooO000OO + ".tif" )
   if 74 - 74: OoOoOO00 / i1IIi % OoooooooOO
 call ( iII11II1II )
 if 52 - 52: IiII % ooOoO0o
 if 25 - 25: I11i / I11i % OoooooooOO - I1ii11iIi11i * oO0o
 if extendWlGrid :
  iII11II1II = [ 'gdalbuildvrt' ]
  iII11II1II . append ( Ii11iiI + '_depth_fill.vrt' )
  if 23 - 23: i11iIiiIii
  for ooooo0 in range ( O0O0 ) :
   for OooO0O0Ooo in range ( O0oO0o0OOOOOO ) :
    ooO000OO = "tile%02i%02i" % ( ooooo0 , OooO0O0Ooo )
    iII11II1II . append ( Ii11iiI + "_depth_fill_" + ooO000OO + ".tif" )
    if 100 - 100: oO0o + O0 . I1IiiI + i1IIi - OoOoOO00 + o0oOOo0O0Ooo
  call ( iII11II1II )
  if 65 - 65: II111iiii / Oo0Ooo
  iII11II1II = [ 'gdalbuildvrt' ]
  iII11II1II . append ( Ii11iiI + '_wl_fill.vrt' )
  if 42 - 42: i11iIiiIii . O0
  for ooooo0 in range ( O0O0 ) :
   for OooO0O0Ooo in range ( O0oO0o0OOOOOO ) :
    ooO000OO = "tile%02i%02i" % ( ooooo0 , OooO0O0Ooo )
    iII11II1II . append ( Ii11iiI + "_wl_fill_" + ooO000OO + ".tif" )
    if 75 - 75: I1Ii111 + iIii1I11I1II1
  call ( iII11II1II )
  if 19 - 19: I1IiiI + i11iIiiIii . IiII - I11i / Ii1I + o0oOOo0O0Ooo
 print "done."
 if 38 - 38: Oo0Ooo / iIii1I11I1II1 * iIii1I11I1II1 % I1ii11iIi11i
 if saveCsv :
  fileIO . saveVectorCSV ( flowX , flowY , xll , yll , cellSize , Ii11iiI + "_flow.csv" , thresholdVal = flowThreshold )
  if 92 - 92: I11i / O0 * I1IiiI - I11i
  if 99 - 99: i11iIiiIii % OoooooooOO
  fileIO . saveScalarCSV ( wlGrid , xll , yll , cellSize , Ii11iiI + "_wl.csv" , headerList = [ 'WL' ] )
  if 56 - 56: IiII * I1Ii111
  if 98 - 98: I11i + O0 * I1Ii111 + i11iIiiIii - OOooOOo - iIii1I11I1II1
  if 5 - 5: OOooOOo % Oo0Ooo % IiII % ooOoO0o
 if zeroPolyList is not None :
  os . remove ( oo0o00OO )
  if 17 - 17: Ii1I + II111iiii + OoooooooOO / OOooOOo / IiII
  if 80 - 80: o0oOOo0O0Ooo % i1IIi / I11i
  if 56 - 56: i1IIi . i11iIiiIii
def deleteFileComplete ( progName , str ) :
 Ii1Ii1IiIIIi1 = os . path . splitext ( os . path . basename ( progName ) ) [ 0 ]
 if 55 - 55: oO0o + O0 / iII111i % ooOoO0o / OoooooooOO
 Ii1Ii1IiIIIi1 += '_' + str + '.txt'
 if 98 - 98: Ii1I * iIii1I11I1II1 % Oo0Ooo % OOooOOo
 try :
  os . remove ( Ii1Ii1IiIIIi1 )
 except :
  pass
  if 88 - 88: iII111i - II111iiii / iII111i - Ii1I
 return
 if 16 - 16: Oo0Ooo % I1Ii111
def notifyFileComplete ( progName , str ) :
 Ii1Ii1IiIIIi1 = os . path . splitext ( os . path . basename ( progName ) ) [ 0 ]
 if 10 - 10: IiII / OoooooooOO
 Ii1Ii1IiIIIi1 += '_' + str + '.txt'
 if 50 - 50: i11iIiiIii - OoooooooOO . oO0o + O0 . i1IIi
 ooo = open ( Ii1Ii1IiIIIi1 , 'w' )
 ooo . close ( )
 if 91 - 91: o0oOOo0O0Ooo . iII111i % Oo0Ooo - iII111i . oO0o % i11iIiiIii
 return
 if 25 - 25: iIii1I11I1II1
 if 63 - 63: ooOoO0o
 if 96 - 96: I11i
def uncompressGeoTiff ( fName , tiled = False , extents = None ) :
 IIIIii1 = tempfile . _get_candidate_names ( ) . next ( ) + '.tiff'
 if 63 - 63: iII111i
 i1i1iIiI = [ 'gdal_translate' ]
 i1i1iIiI += [ '-co' , 'BIGTIFF=YES' ]
 if 23 - 23: IiII + iIii1I11I1II1 % iIii1I11I1II1 / ooOoO0o . oO0o + iIii1I11I1II1
 if 93 - 93: oO0o * o0oOOo0O0Ooo / OOooOOo - OOooOOo . iII111i / I1IiiI
 if 11 - 11: I1Ii111 - I11i % i11iIiiIii . iIii1I11I1II1 * I1IiiI - Oo0Ooo
 if tiled :
  i1i1iIiI += [ '-co' , 'TILED=YES' ]
  if 73 - 73: O0 + ooOoO0o - O0 / OoooooooOO * Oo0Ooo
 if extents is not None :
  i1i1iIiI += [ '-projwin' ]
  i1i1iIiI += [ "%f" % extents [ 0 ] ]
  i1i1iIiI += [ "%f" % extents [ 3 ] ]
  i1i1iIiI += [ "%f" % extents [ 2 ] ]
  i1i1iIiI += [ "%f" % extents [ 1 ] ]
  if 32 - 32: OoO0O00 % I1IiiI % iII111i
 i1i1iIiI += [ fName ]
 i1i1iIiI += [ IIIIii1 ]
 if 66 - 66: OoOoOO00 + o0oOOo0O0Ooo
 call ( i1i1iIiI )
 if 54 - 54: I1ii11iIi11i + I1ii11iIi11i + I11i % i1IIi % i11iIiiIii
 return IIIIii1
 if 100 - 100: I1ii11iIi11i
 if 96 - 96: I1IiiI . IiII * II111iiii % IiII . I1Ii111 * i1IIi
 if 83 - 83: iIii1I11I1II1
def clipRasterPoly ( rasterFileName , polyFileName ) :
 i1i1iIiI = [ 'gdal_rasterize' ]
 i1i1iIiI += [ '-burn' , '-9999' ]
 i1i1iIiI += [ '-i' ]
 i1i1iIiI += [ polyFileName ]
 i1i1iIiI += [ rasterFileName ]
 if 97 - 97: i11iIiiIii + Oo0Ooo * OOooOOo % iII111i . IiII
 call ( i1i1iIiI )
 if 4 - 4: O0 . iII111i - iIii1I11I1II1
 return
 if 19 - 19: OOooOOo % OoO0O00 / Ii1I + II111iiii % OoooooooOO
 if 89 - 89: Ii1I
def calcTimeStep ( wlGrid , zGrid , dx ) :
 I11iiI1i1 = wlGrid - zGrid
 I11iiI1i1 [ numpy . where ( zGrid == - 9999 ) ] = 0.
 if 51 - 51: iII111i
 return 0.7 * dx / numpy . sqrt ( 9.81 * I11iiI1i1 . max ( ) )
# dd678faae9ac167bc83abf78e5cb2f3f0688d3a3
