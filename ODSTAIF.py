#! /usr/bin/python2.7
# coding=utf-8

########################################################

import matplotlib.pyplot as plt
import numpy as np
from numpy import *
import pylab
import os
import poudre		# Contient la fonction der_auto pour caler automatiquement 
from poudre import * 	# la droite de pente nulle de la dérivée à l'IARF

########################################################
print ''
print '################################################################'
print '################################################################'
print '###                                                          ###'
print '###   Outil d analyse des essais en bout de tige pour les    ###'
print '###   fluides compressibles et incompressibles. Methodes     ###'
print '###   tirees de Ezekwe - 2010 - Petroleum reservoir engin-   ###'
print '###   eering practice - Prentice Hall                        ###'
print '###                                                          ###'
print '################################################################'
print '################################################################'
print ''

print 'Récuperation des données de pression et temps'
print''

#
# Tests logiques pour entrer les unites de pression et temps
#
# Pression 
#

q_p=float(raw_input('Unité de pression? psi = 0, kPa = 1, Mpa = 2\n'))
if q_p==2:
	unit_p='MPa'
	print 'Conversion MPa en psi'
	print ''
	p_MPa=loadtxt("pression.txt")
	p=p_MPa/0.00689475729
elif q_p==1:
	unit_p='kPa'
	print 'Conversion kPa en psi'
	print ''
	p_kPa=loadtxt("pression.txt")
	p=p_kPa/6.89475729
else:
	unit_p='psi'
	p = loadtxt("pression.txt")

print 'Pressions',p,'psi'
print ''

#
# Temps de test
#

q_t=float(raw_input('Unité de temps? min = 0, h = 1\n'))
if q_t==0:
	unit_t='min'
	t = loadtxt("temps.txt")
else:
	unit_t='heure'
	print 'Conversion h en min'
	print ''
	t=60*loadtxt("temps.txt")

if t[0]==0:
	t=t[1:] # on veut exclure t=0 car
	p=p[1:]	# un calcul de log de temps est effectué juste après
 
	
print 'Temps',t,'min'
print ''

#
# Temps de production 
#

tp=float(raw_input('Temps de production de fluides? (en min)\n'))
print ''



#############################################################
###                                                       ###
###         Calcul de la dérivée de la pression           ###
###                                                       ###
#############################################################


te,lt=agarwal(t,tp)

pf, dp, dt =deriv(lt,te,p)

#
# Plot loglog
#

popp='Pression'
plot_bourd(te,pf,dt,dp,popp)

hap='n'
while hap=='n':
	pf, dp, dt = deriv(lt,te,p)	
	
	x_IARF, reg_dp, ind_IARF_deb, ind_IARF_fin, y_dp_0 = der_auto(dt,dp)

	### Plotty time!	
	plot_bourd_s(te,pf,dt,dp,x_IARF,reg_dp,popp)
	hap=raw_input('Alors heureux? (n = retrace Bourdet plot)')
	print''

# y_dp_0=float(raw_input('Valeur de dp relevee avec pente de 0?'))



#############################################
###                                       ###
###             Horner plot               ###
###                                       ###
#############################################

#
# Calcul du temps de Horner
#

tH=(tp+t)/t

#
# Calcul de la pente de l'IARF sur un Horner plot
#

slope = penteH(p,tH,ind_IARF_deb,ind_IARF_fin)


#
# Calcul de p0
#

p0=slope*log10(tH[ind_IARF_fin]/1)+p[ind_IARF_fin]

print ''
print 'Pente de Horner', slope, 'psi/cycle log'
print ''
print 'Pseudo-pression initiale du réservoir', p0, 'psi'
print ''

#
# Tracé de la droite correspondant à la pente de l'IARF
#

tH_pente=1
tH_pente=np.append(tH_pente,tH)
p_pente=-slope*log10(tH_pente)+p0

#
# Horner plot de la pseudo-pression
#

hap='n'
while hap=='n':

	ylab='pression (psi)'
	plot_horner(tH,p,tH_pente,p_pente,ylab,popp)
	hap=raw_input('Alors heureux? (n = retrace Horner plot)\n')




###################################################################
###                                                             ###
###     Récupération des données de débit et détermination      ###
###     de la méthode d'interprétation                          ###
###                                                             ###
###################################################################

unit_h=float(raw_input('Unité de hauteur de fluides et longueur des tiges? Mètres= 0 et pieds = 1\n'))
print''
h_flu=float(raw_input('Hauteur des fluides?\n'))
print''
l_MT=float(raw_input('Longueur des masses tiges?\n'))
print''
d_MT_FU=float(raw_input('Diamètre interne des masses tiges (en pouces) ?\n'))
print''
d_T_FU=float(raw_input('Diamètre interne des tiges (en pouces) ?\n'))
print''

A_MT=(d_MT_FU*0.0254/2)**2*np.pi
A_T=(d_T_FU*0.0254/2)**2*np.pi

print"Aire interne des masses tiges",A_MT,'m'
print''
print"Aire interne des tiges",A_T,'m'
print''

#
# Conversion des longueurs en pieds vers les mètres
#

if unit_h==1:
	h_flu=h_flu*0.3048
	l_MT=l_MT*0.3048

#
# Calcul du volume produit dans les tiges en unités SI
#
	
if h_flu>l_MT:
	v_flu=l_MT*A_MT+(h_flu-l_MT)*A_T
else:
	v_flu=h_flu*A_MT

v_flu_FU=v_flu*6.2898

#
#Calcul du débit en fonction du volume produit et du temps de production
#

q=v_flu_FU/(tp/1440)

###################################################################
###                                                             ###
###             Calcul des propriétés des fluides               ###
###                                                             ###
###################################################################

temp=raw_input('Temperature en Farenheit (0 si inconnue)?\n')
T_FU=float(temp)
print ''
if T_FU==0:
	print 'Guestimation de la temperature des fluides au niveau du testeur'
	print ''
	prof=raw_input('Profondeur du testeur en pieds?\n')
	print ''
	z_ft=float(prof)
	temp_moy_an=5
	grad_G=0.02 # moyenne de Lefebvre - 1982 et Tran Ngoc et al - 2011
	T_USI=temp_moy_an+grad_G*(z_ft*0.3048)
	T_FU=T_USI*9/5+32
	print 'T =',T_USI,'deg C ou',T_FU,'deg F'
print ''
T_USI=(T_FU-32)*5/9

S=float(raw_input('Salinité (en pourcentage de masse) ?\n'))

#
# Densité des eaux de formation en lb/ft**3
#

rau_w=62.368+0.438603*S+1.60074e-3*S**2 

# 
# Viscosité des eaux de formation d'après McCain - 1991
#

mu_w1=(109.574-8.40564*S+0.313314*S**2+8.72213e-3*S**3)*T_FU**(-(1.12166-2.63951e-2*S+6.79461e-4*S**2+5.47119e-5*S**3-1.55586e-6*S**4))
mu_w=mu_w1*(0.9994+4.0295e-5*p0+3.1062e-9*p0**2)

#
# Formation Volume Factor des eaux de formations d'après McCain - 1991
#

B_w=(1+(-1.0001e-2+1.3391e-4*T_FU+5.50654e-7*T_FU**2))*(1+(-1.95301e-9*p0*T_FU-1.72834e-13*p0**2*T_FU-3.58922e-7*p0-2.25341e-10*p0**2))


###################################################################
###                                                             ###
###           Calcul de la perméabilité du réservoir            ###
###                                                             ###
###################################################################

z_inf=float(raw_input('Profondeur du bas de la zone investiguée ? (en pieds)\n'))
print ''
z_sup=float(raw_input('Profondeur du haut de la zone investiguée ? (en pieds)\n'))
print ''
h=z_inf-z_sup
print('Épaisseur de la formation testée:',h,'pieds')
print''


#
# Calcul de la perméabilité d'après la pente de la droite de Horner Ezekwe - 11.30
#
 
m=slope
k=162.6*q*B_w*mu_w/m/h

#
# Calcul de la perméabilité d'après la hauteur de la dérivée à l'IARF
#

k_Bourdet=162.6*q*B_w*mu_w*0.5/h/y_dp_0

print ''	
print "Perméabilité de la formation d'après la méthode de Horner k =",k,'md'
print ''
print "Perméabilité de la formation, d'après la dérivée, k =",k_Bourdet,'md'
print ''

###################################################################
###                                                             ###
###                    Sortie des résultats                     ###
###                                                             ###
###################################################################

#
# Plots de Horner et Bourdet avec les k et p0
#

plot_horner_f(tH,p,tH_pente,p_pente,ylab,popp,k,p0)
plot_bourd_f(te,pf,dt,dp,x_IARF,reg_dp,popp,k_Bourdet)

#
# Définition d'une matrice de résultats inéressants modifiable
#

results=array([k,k_Bourdet,p0,T_USI])

#
# Enregistrement des résultats dans un fichier *.csv
#
# Critique programmation: on ne sait pas à quoi sert "format='%.6f'"
#


results.tofile(file='resultats.csv', format='%.6f', sep=';')

#
# Ouverture d'un éditeur de texte pour visualiser/copier les résultats
#

os.system('leafpad resultats.csv')

#
# Message de fin
#

print''
print " - - - - - Tout est bien qui fini bien. - - - - - "
print"                          -"
print"                          -"
print"                          -"
print ''


