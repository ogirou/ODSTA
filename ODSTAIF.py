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

q_p=float(raw_input('Unité de pression? psi = 0, kPa = 1, MPa = 2, Pa = 3\n'))
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
elif q_p==3:
	unit_p='Pa'
	print 'Conversion Pa en psi'
	print ''
	p_Pa=loadtxt("pression.txt")
	p=p_Pa/6894.75729
else:
	unit_p='psi'
	p = loadtxt("pression.txt")

print 'Pressions',p,'psi'
print ''

#
# Temps de test
#

q_t=float(raw_input('Unité de temps? s = 0, min = 1, h = 2\n'))

if q_t==0:
	unit_t='sec'
	print 'Conversion sec en min'
	print ''
	t=loadtxt("temps.txt")/60

elif q_t==1:
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

mano=raw_input('Recalcul manuel de la pente ? o = oui, retrace Horner plot.\n')
if mano=='o':
	hap='n'
	while hap=='n':
		p0, slope, tH_pente, p_pente = penteM(tH,popp)
		ylab='pression (psi)'
		plot_horner(tH,p,tH_pente,p_pente,ylab,popp)
		hap=raw_input('Alors heureux? (n = retrace Horner plot)\n')


###################################################################
###                                                             ###
###     Récupération des données de débit et détermination      ###
###     de la méthode d'interprétation                          ###
###                                                             ###
###################################################################
Vprod=raw_input('Volume produit connu directement? o = oui\n')

if Vprod=='o':
	v_flu=float(raw_input('Volume produit? (m3)\n')) 
else:
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

	print"Aire interne des masses tiges",A_MT,'m2'
	print''
	print"Aire interne des tiges",A_T,'m2'
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

T_USI, T_FU =Temp()

S=float(raw_input('Salinité (en pourcentage de masse) ?\n'))
print''

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
rau_w=rau_w*B_w

###################################################################
###                                                             ###
###           Calcul de la perméabilité du réservoir            ###
###                                                             ###
###################################################################

unit_z=raw_input("Unité de profondeur des obturateurs et d'élévation du KB? (m)ètres ou (f)eet\n")
print''

z_inf=float(raw_input('Profondeur du bas de la zone investiguée ?\n'))
print ''
z_sup=float(raw_input('Profondeur du haut de la zone investiguée ?\n'))
print ''

if unit_z=='m':
	z_inf_USI=z_inf
	z_inf=z_inf/0.3048
	z_sup_USI=z_sup
	z_sup=z_sup/0.3048

if unit_z=='f':
	z_inf_USI=z_inf*0.3048
	z_sup_USI=z_sup*0.3048

h=z_inf-z_sup

print'Épaisseur de la formation testée:',h,'pieds'
print''

h_USI=h*0.3048

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

#######################################
###                                 ###
###           Rendu sexy            ###
###                                 ###
#######################################

#
# Plots de Horner et Bourdet avec les k et p0
#

pouits=raw_input("Cote du puits?")
print ''
num=raw_input("Numéro d'essai?")
print ''
typ=raw_input("ISI ou FSI?")
print ''
Essai=pouits+'DST'+num+typ

hkb=float(raw_input('Élévation du Kelly Bushing (même unité que z obturateurs)'))
print''

if unit_z=='f':
	hkb=hkb*0.3048


plot_horner_f(tH,p,tH_pente,p_pente,ylab,popp,k,p0,Essai)
plot_bourd_f(te,pf,dt,dp,x_IARF,reg_dp,popp,k_Bourdet,Essai)

#
# Définition d'une matrice de résultats inéressants modifiable
#

if Vprod=='o':
	results=array([Essai,'k Horner (md)',k,'k Bourdet (md)',k_Bourdet,'pression intiale (psi)',p0,'Température °C',T_USI,'Température °F',T_FU,'Débit (bbl/j)',q,'Bas de la zone investiguée (m)',z_inf_USI,'Haut de la zone investiguée (m)',z_sup_USI,'Épaisseur (m)',h_USI,'Élévation du KB (m)',hkb,'FVF (-)',B_w,'Viscosité des fluides (cp)',mu_w,'Temps de production (min)',tp,'Ordonnée de la droite de pente nulle (psi)',y_dp_0,'Pente de la droite de Horner de pp (psi/cycle)',slope])
else:
	results=array([Essai,'k Horner (md)',k,'k Bourdet (md)',k_Bourdet,'pression intiale (psi)',p0,'Température °C',T_USI,'Température °F',T_FU,'Débit (bbl/j)',q,'Bas de la zone investiguée (m)',z_inf_USI,'Haut de la zone investiguée (m)',z_sup_USI,'Épaisseur (m)',h_USI,'Élévation du KB (m)',hkb,'FVF (-)',B_w,'Viscosité des fluides (cp)',mu_w,'Temps de production (min)',tp,'Ordonnée de la droite de pente nulle (psi)',y_dp_0,'Pente de la droite de Horner de pp (psi/cycle)',slope,'Hauteur de fluide récupérée (m)',h_flu,'Diamètre interne des masses tiges (in)',d_MT_FU,'Diamètre interne des tiges (in)',d_T_FU,'Longueur des masses tiges (m)',l_MT])

#
# Enregistrement des résultats dans un fichier *.csv
#

results.tofile(file='resultatsfi.csv', sep='\n')

#
# Ouverture d'un éditeur de texte pour visualiser/copier les résultats
#

#os.system('leafpad resultatsfi.csv')

#
# Fun avec LaTeX
#

l1='\\documentclass[10pt]{article} \n\\usepackage[utf8x]{inputenc} \n\\usepackage[frenchb,english]{babel} \n\\usepackage[T1]{fontenc} \n\\usepackage{lmodern} \n \\usepackage{graphicx} \n\n\\begin{document} \n\n\\begin{table}[!h]'
l2='%s %s %s' % ("\\caption{Résultats de l'essai",Essai,'}')
l3='\\begin{center}\n\\begin{tabular}{|c|c|c|}\\hline'
l4='\\bfseries Perméabilité Horner & %.4f & md\\\\ ' % k
l5='\\bfseries Perméabilité Bourdet & %.4f & md\\\\ ' % k_Bourdet
l6='\\bfseries Pression intiale & %.2f & psi\\\\ ' % p0
l7='\\bfseries Température & %.1f & ° C\\\\ ' % T_USI
l8='\\bfseries Température & %.1f & ° F\\\\ ' % T_FU
l9='\\bfseries Bas de la zone investiguée & %s & m\\\\ ' % z_inf_USI
l10='\\bfseries Haut de la zone investiguée & %s & m\\\\ ' % z_sup_USI
l11='\\bfseries Épaisseur & %s & m\\\\ ' % h_USI
l12='\\bfseries Élévation du KB & %s & m\\\\ ' % hkb
l13='\\bfseries Temps de production & %.2f & min\\\\ ' % tp
l14='\\bfseries Débits & %.2f & bbl/j\\\\ ' % q
l15='\\bfseries FVF & %.2f & -\\\\ ' % B_w
l16='\\bfseries Viscosité & %.3f & cp\\\\ ' % mu_w
l17='\\bfseries Ordonnée de la droite de pente nulle & %.1f & psi$^2$/cp\\\\' % y_dp_0
l18='\\bfseries Pente de la droite de Horner & %.1f & psi/cycle log\\\\' % slope
llastab='\\hline\n\\end{tabular}\n\\end{center}\n\\end{table}\n'
lfig1="\\begin{figure}[!h]\\centering\\includegraphics[scale=0.55]{%sHp.eps}\\caption{Horner plot de l'essai %s}\\end{figure}" % (Essai, Essai)
lfig2="\\begin{figure}[!h]\\centering\\includegraphics[scale=0.55]{%sBp.eps}\\caption{Derivative plot de l'essai %s}\\end{figure}" % (Essai, Essai)
llast='\n\\end{document}'

if Vprod !='o':
	l19='\\bfseries Hauteur de fluides récupérés & %s & m\\\\ ' % h_flu
	l20='\\bfseries Longueur des masses tiges & %s & m\\\\ ' % l_MT
	l21='\\bfseries Diamètre interne des MT & %s & in\\\\ ' % d_MT_FU
	l22='\\bfseries Diamètre interne  des tiges & %s & in\\\\ ' % d_T_FU

	crlatex=array([l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16,l17,l18,l19,l20,l21,l22,llastab,lfig1,lfig2,llast])

else:

	crlatex=array([l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16,l17,l18,llastab,lfig1,lfig2,llast])

	


crlatex.tofile(file='woupfi.tex', sep='\n')

#
# Partie de compilation du fichier tex, affichage du dvi et effacage des déchets
# sert surtout à se faire mousser
#

os.system('latex woupfi.tex')
os.system('xdvi woupfi.dvi')
os.system('rm woupfi.log')
os.system('rm woupfi.aux')

#
# Message de fin
#

print''
print"                          -"
print"                          -"
print"                          -"
print " - - - - - Tout est bien qui fini bien. - - - - - "
print"                          -"
print"                          -"
print"                          -"
print ''

