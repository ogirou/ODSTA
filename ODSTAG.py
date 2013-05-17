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

q_p=int(raw_input('Unité de pression? psi = 0, kPa = 1, Mpa = 2\n'))
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

q_t=int(raw_input('Unité de temps? min = 0, h = 1'))
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

tp=int(raw_input('Temps de production de fluides? (en min)\n'))
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

# y_dp_0=int(raw_input('Valeur de dp relevee avec pente de 0?'))

np.savetxt('dt_dp.csv', (dt,dp), delimiter=',')

###################################################################
###                                                             ###
###     Récupération des données de débit et détermination      ###
###     de la méthode d'interprétation                          ###
###                                                             ###
###################################################################



var_q=int(raw_input('Le débit est il variable? Non = 0, Oui un peu (test court) = 1 et Beaucoup (test long) = reste'))
print ''

#
# 3 cas: q cst = méthode de Horner classique, q variable = temps apparent de 
# Horner, q très variable = déconvolution à intégrer (pas fait)
#


if var_q==0: # Méthode: Horner plot classique (cas 1)
	metho_q=0
	q=int(raw_input('Valeur du débit en Mcf/d')) # Il suffit d'entrer le débit
	print ''
	tq=tp
	tH=(tp+t)/t 		# Fonction de temps de Horner
	tpH='null'

elif var_q==1: # Méthode: Horner plot avec un temps apparent de Horner qui dépend de q
	# (cas 2)
	metho_q=1	
	q_q=int(raw_input('Unité de débit? Mcf/d = 0, m3/s = 1, m3/j = 2')) # conversion du temps
	print ''	
	if q_q==0:
		unit_p='Mcf/d'
		q=loadtxt("debits.txt")
	elif q_p==1:
		unit_p='m^3/s'
		print 'Conversion m^3/s en Mcf/d'
		print ''
		q_m3s=loadtxt("debits.txt")
		q=3051187.20*q_m3s
	elif q_p==2:
		print 'Conversion m^3/j en Mcf/d\n'
		q_m3j=loadtxt("debits.txt")
		q=35.3146667*q_m3j
	else:
		print 'Nouvelle conversion à intégrer au code'
		print ''

	unit_tq=raw_input('Unité de temps du débit de production? Jour = 0, heure = 1, minute = 2')

	if unit_tq==1:
		print 'Conversion h en min'
		print ''
		tq=60*loadtxt("tflow.txt")
	elif unit_t==0:
		print 'Conversion j en min'
		print ''
		tq=1440*loadtxt("tflow.txt")		
	else:
		tq=loadtxt("tflow.txt")
	
	#
	# Calcul du temps apparent de Horner
	#
	
	Gp=sum(q*tq)
	qglast=q[-1]
	tpH=24*Gp/qglast
	tH=(tpH+t)/t
	print 'tpH', tpH, 'min'
	print ''


else:
	print 'Déconvolution à intégrer au code, respire, va boire un café, va manger un chien chaud'
	

print 'Temps',tq,'min'
print ''
print 'Débits',q,'Mcf/d', tq, 'min'
print ''

np.savetxt('t_p_te_pf_tH.csv', (t,p,te,pf,tH), delimiter=',')
np.savetxt('tq_q.csv', (tq,q), delimiter=',')

##################################################
###                                            ###
###                Température                 ###
###                                            ###
##################################################



print 'Calcul des variables des gaz'
print ''
print('En l absence de donnees sur la nature du gaz on choisi les proprietes PVT du methane pur')
print ''

T_FU=float(raw_input('Temperature en Farenheit (0 si inconnue)?'))
print ''

#
# Test logique pour estimer la température si inconnue
#

if T_FU==0:
	print 'Guestimation de la temperature des fluides au niveau du testeur'
	print ''
	z_ft=int(raw_input('Profondeur du testeur en pieds?'))
	print ''
	temp_moy_an=5
	grad_G=0.02 # moyenne de Lefebvre - 1982 et Tran Ngoc et al - 2011
	T_USI=temp_moy_an+grad_G*(z_ft*0.3048)
	T_FU=T_USI*9/5+32
	print 'T =',T_USI,'deg C ou',T_FU,'deg F'
	print ''

T_USI=(T_FU-32)*5/9

###################################################################
###                                                             ###
###             Calcul des variables des gaz                    ###
###                                                             ###
###################################################################

print"Calcul des variables des gaz, attention aux gaz choisis"
print''

R=10.732		# Constante des gaz parfaits en unités de terrain

#
# Calcul de T et p pseudo-réduites
#

# Methane
#M_g=16.042		# Masse molaire du méthane
#T_pc=-116.66		# T pseudo-critique du méthane
#p_pc=667.0		# p pseudo-critique du méthane

# Ethane
#M_g=30.069
#T_pc=89.92
#p_pc=706.6


mix=float(raw_input('Connaissance exacte de la mixture de gaz? Oui = 1, Non = 0'))
print''

if mix==1:
	#
	# Calcul de p_pr et T_pr par la méthode de Kay - 1936
	#
	
	print'Méthode Kay - 1936, Valable quand la portion de C7+ est faible.\n'

	table=array([[16.042,667.0,-116.66],
	[30.069,706.6,89.92],
	[44.096,615.5,205.92],
	[58.122,527.9,274.41],
	[58.122,550.9,305.55],
	[72.149,490.4,369],
	[72.149,488.8,385.8],
	[86.175,436.9,453.8],
	[100.202,396.8,512.9],
	[114.229,360.7,564.2],
	[128.255,330.7,610.8],
	[142.282,304.6,652.2],
	[28.01,506.7,-220.63],
	[44.01,1070.0,87.76],
	[34.082,1306.5,212.81],
	[28.959,551.9,-220.97],
	[2.0159,190.7,-399.9],
	[31.9988,731.4,-181.43],
	[28.0135,492.5,-232.53],
	[18.0153,3200.1,705.1]])

	fichier = open("compo.csv", "r")

	compo=array([])
	for n in arange(1,21,1):
		l=fichier.readline().rstrip('\n\r').split(",")
		compo=np.append(compo,float(l[1]))

	p_pc=sum(table[0,1]*compo)
	T_pc=sum(table[0,2]*compo)

else:
	#
	# Calcul de p_pr et T_pr selon la méthode de Sutton - 2005
	#

	gamma_g=float(raw_input('Gravité du gaz? (-)'))
	print''
	M_g=gamma_g*28.9586

	#
	# 4.45 et 4.46 - Ezekwe (Pour les gaz condensés)
	#
	# Existe aussi pour les hydrocarbures associées mais ne fonctionne pas 
	# avec les gravités faibles observées
	#

	p_pc=744-125.4*gamma_g+5.9*gamma_g**2
	T_pc=164.3+357.7*gamma_g-67.7*gamma_g**2	# en deg R
	
#
# 4.17 et 4.18 - Ezekwe
#

p_pr=p/p_pc

T_pr=(T_FU+460)/(T_pc)


print 'Temperature pseudo-reduite',T_pr,'R'
print ''
print 'Pression pseudo-reduite', p_pr,'psig (',len(p_pr),')'
print ''

#
# Calcul de z d'après la méthode de Newton-Raphson tirée de Sutton?
#
# Voir Ezekwe 4.????
#

A1=.3265
A2=-1.07
A3=-.5339
A4=.01569
A5=-.05165
A6=.5475
A7=-.7361
A8=.1844
A9=.1056
A10=.6134
A11=.721
c1=A1+A2/T_pr+A3/T_pr**3+A4/T_pr**4+A5/T_pr**5
c2=A6+A7/T_pr+A8/T_pr**2
c3=A9*(A7/T_pr+A8/T_pr**2)

z=1
for n in arange (1,50,1):		
	rho_r=.27*p_pr/(z*T_pr)
	c4=A10*(1+A11*rho_r**2)*(rho_r**2/T_pr**3)*exp(-A11*rho_r**2)
	z=z-(z-(1+c1*rho_r+c2*rho_r**2-c3*rho_r**5+c4))/(1+c1*rho_r/z+2*c2*rho_r**2/z-5*c3*rho_r**5/z+2*A10*rho_r**2/(z*T_pr**3)*(1+A11*rho_r**2-(A11*rho_r**2)**2)*exp(-A11*rho_r**2))

print 'Valeurs de z',z,'(-)','(',len(z),')'
print ''

#
# Calcul des variables manquantes
#
# 4.80-84 - Ezekwe from Lee et al - 1966
#
# Densité
#

rau_g=M_g*p/z/R/(T_FU+460) # résultat en lbm/ft/ft/ft

rau_gUSI=rau_g*0.01601846

print 'Densite du gaz', rau_gUSI, 'kg/m3','(',len(rau_gUSI),')'
print ''

#
# Viscosité
#

K=((9.379+0.01607*M_g)*(T_FU+460)**1.5)/(209.2+19.26*M_g+(T_FU+460));
X=3.448+986.4/(T_FU+460)+0.01009*M_g;
Y=2.447-0.2224*X;
mu_g=0.0001*K*exp(X*rau_gUSI**Y)

print 'Viscosite du gaz', mu_g, 'cp','(',len(mu_g),')'
print ''

#
# Formation Volume Factor
#
# 4.75 - Ezekwe
#
 
B_g=0.02819*z*(T_FU+459.67)/p
print 'FVF', B_g,'(-)','(',len(B_g),')'
print ''

#
# Il est encore possible de calculer c_t mais on ne s'en sert pas
#


###################################################################
###                                                             ###
###             Calcul de la pseudo-pression                    ###
###                                                             ###
###################################################################

#
# Détermination de la pseudo pression d'apres Al-Hussainy et al - 1966
#
# Calcul de l'intégrale de la fonction p²/(mu_g*z) par la méthode des rectangles 
# 
# (l'intégrale correspond à la surface des rectangles sous la courbe)
#
# Voir dans Ezekwe - 2010, Exemple 10.8 p 324-326 et pseudop.py
#

f=2*p/(mu_g*z)				# Définition de la fonction p²/(mu_g*z)

H=(p[1:]-p[:-1])			# Début du calcul de l'intégrale 
A=(p[1:]-p[:-1])/2*(f[1:]+f[:-1])	# (base des rectanles)
pp=H[0]/2*f[0]				# calcul de la première pseudo-pression
p_p=pp				# on imprime les valeurs de pseudo-pression dans p_p
for A in A:
	pp=A+pp				# nouvelle pp = A + ancienne pp
	p_p=np.append(p_p,pp)	# on imprime les valeurs de pseudo-pression dans p_p
print 'pseudo-pression',p_p,'psi²/cp','(',len(p_p),')'
print''


np.savetxt('ppr_z_rhog_mug_FVF_pp.csv', (p_pr,z,rau_g,mu_g,B_g,p_p), delimiter=',')

###################################################################
###                                                             ###
###           Bourdet plot de la pseudo-pression                ###
###                                                             ###
###################################################################

#
# Calcul de la dérivé de la pseudo pression 
#

pf, dp, dt=deriv(lt,te,p_p)

#
# Plot loglog
#

popp='Pseudo-pression'
plot_bourd(te,pf,dt,dp,popp)

#
# Boucle pour caler la droite de pente nulle de la dérivée
#
# sur le plot loglog
#

hap='n'
while hap=='n':

	pf, dp, dt = deriv(lt,te,p_p)	
	
	### on reprend de nouvelles limites et on recalcule la dérivée à l'IARF

	x_IARF, reg_dp, ind_IARF_deb, ind_IARF_fin, y_dp_0 = der_auto(dt,dp)

	### Plotty time!	

	plot_bourd_s(te,pf,dt,dp,x_IARF,reg_dp,popp)
	hap=raw_input('Alors heureux? (n = retrace Bourdet plot)')

np.savetxt('dtpp_dpp.csv', (dt,dp), delimiter=',')

###################################################################
###                                                             ###
###           Horner plot de la pseudo-pression                 ###
###                                                             ###
###################################################################


#
# Calcul de la pente de l'IARF
#

slope = penteH(p_p,tH,ind_IARF_deb,ind_IARF_fin)

#
# Calcul de pp0
#

pp0=slope*log10(tH[ind_IARF_fin]/1)+p_p[ind_IARF_fin]
#pp0=int(raw_input('Pseudo-pression initiale du réservoir pp0? (psi2/cp)'))

print ''
print 'Pente de Horner', slope, '(psi²/cp)/cycle log'
print ''
print 'Pseudo-pression initiale du réservoir', pp0, '(psi²/cp)/cycle'
print ''

#
# Tracé de la droite correspondant à la pente de l'IARF
#

tH_pente=1
tH_pente=np.append(tH_pente,tH)
p_pente=-slope*log10(tH_pente)+pp0

#
# Horner plot de la pseudo-pression
#

hap='n'
while hap=='n':

	ylab='pseudo-pression (psi2/cp)'
	plot_horner(tH,p_p,tH_pente,p_pente,ylab,popp)
	hap=raw_input('Alors heureux? (n = retrace Horner plot)')

#
# Horner plot avec pente déterminée graphiquement
#

mano=int(raw_input('Recalcul manuel de la pente ? 0 = Non, 1 = Oui'))
if mano==1:
	hap='n'
	while hap=='n':
		pp0, slope, tH_pente, p_pente = penteM(tH)
		ylab='pseudo-pression (psi2/cp)'
		plot_horner(tH,p_p,tH_pente,p_pente,ylab,popp)
		hap=raw_input('Alors heureux? (n = retrace Horner plot)')	

np.savetxt('tHpp.csv', (tH), delimiter=',')

###################################################################
###                                                             ###
###           Extrapolation de la pseudo-pression               ###
###           initiale pour déterminer la pression              ###
###           initiale du réservoir par régression              ###
###           polynomiale                                       ###
###                                                             ###
###################################################################

#
# Régression polynomiale de pp / p de degré 2
#
# polyfit retourne les valeurs a, b et c du polynome pp = a*p^2 + b*p + c
#

cor_pp_p=np.polyfit(p_p,p,2)		

#
# p_calc est la fonction pp = f(p) définie par polyfit
#

p_calc=np.poly1d(cor_pp_p)		

#
# Critique science: calcul de R pour qualité corrélation
#
# graph pour visualiser la corrélation
#

#pylab.plot(p_p,p,'s',p_p,p_calc(p_p),'r-')
#pylab.xlabel('pseudo-presion (psi2/cp)')
#pylab.ylabel('pression (psi)')
#leg = plt.legend(('p mesure','p calcule'),'upper left', shadow=False)
#plt.grid(True)
#pylab.show()

p0=p_calc(pp0)

print 'Pression initiale dans le réservoir p0=',p0, 'psi'
print ''

###################################################################
###                                                             ###
###           Calcul de la perméabilité du réservoir            ###
###                                                             ###
###################################################################

z_inf=int(raw_input('Profondeur du bas de la zone investiguée ? (en pieds)'))
print ''
z_sup=int(raw_input('Profondeur du haut de la zone investiguée ? (en pieds)'))
print ''
h=z_inf-z_sup
print 'Épaisseur de la formation testée:',h,'pieds'
print''

#
# Calcul de k
#

m=slope

if var_q==1:
	k=(1637*qglast*(T_FU+460))/(m*h) # calcul de k utilisant le temps apparent pour q variable
else:
	k=(1637*q*(T_FU+460))/(m*h) # calcul de k utilisant le temps apparent pour q variable
 	
print 'Perméabilité de la formation k=',k,'md'

#
# Horner plot avec k et p0
#

plot_horner_f(tH,p_p,tH_pente,p_pente,ylab,popp,k,p0)

#
# Définition d'une matrice de résultats intéressants modifiable 
#

pouits=raw_input("Cote du puits?")
print ''
num=raw_input("Numéro d'essai?")
print ''
typ=raw_input("ISI ou FSI?")
print ''
Essai=pouits+num+typ
hkb=raw_input('Élévation du Kelly Bushing (même unité que z packer)')
print''

results=array([Essai,'k (md)',k,'pression intiale (psi)',p0,'Température °C',T_USI,'Température °F',T_FU,'Bas de la zone investiguée (ft)',z_inf,'Haut de la zone investiguée (ft)',z_sup,'Épaisseur (ft)',h,'Élévation du KB (ft)',hkb,'Temps de production (min)',tp,"Temps de production de Horner (min)",tpH,'Gravité du gaz (-)',gamma_g,'Température pseudo réduite (-)',T_pr,'Ordonnée de la droite de pente nulle (psi2/cp)',y_dp_0,'Pente de la droite de Horner de pp (psi2/cp)',slope,'Pseudo-pression initiale (psi2/cp)',pp0])

#
# Enregistrement des résultats dans un fichier *.csv
#

results.tofile(file='resultats.csv', sep='\n')

#
# Ouverture d'un éditeur de texte pour visualiser/copier les résultats
#

os.system('leafpad resultats.csv')

#
# Fun avec LaTeX
#

l1='\\documentclass[10pt]{article} \n\\usepackage[utf8x]{inputenc} \n\\usepackage[frenchb,english]{babel} \n\\usepackage[T1]{fontenc} \n\\usepackage{lmodern} \n\n\\begin{document} \n\n\\begin{table}'
l2='%s %s %s' % ('\\caption{Resultats de l essai',Essai,'}')
l3='\\begin{center}\n\\begin{tabular}{|c|c|c|}\n\\hline'
l4='\n\\bfseries Perméabilité & %s & md\\\\ \n' % k
l5='\n\\bfseries Pression intiale & %s & psi\\\\ \n' % p0
l6='\n\\bfseries Température & %s & ° C\\\\ \n' % T_USI
l7='\n\\bfseries Température & %s & ° F\\\\ \n' % T_FU
l8='\n\\bfseries Bas de la zone investiguée & %s & ft\\\\ \n' % z_inf
l9='\n\\bfseries Haut de la zone investiguée & %s & ft\\\\ \n' % z_sup
l10='\n\\bfseries Épaisseur & %s & ft\\\\ \n' % h
l11='\n\\bfseries Élévation du KB & %s & ft\\\\ \n' % hkb
l12='\n\\bfseries Temps de production & %s & min\\\\ \n' % tp
l13='\n\\bfseries Temps de production de Horner & %s & min\\\\ \n' % tpH
l14='\n\\bfseries Gravité du gaz & %s & -\\\\ \n' % gamma_g
l15='\n\\bfseries Température pseudo-réduite & %s & -\\\\ \n' % T_pr
l16='\n\\bfseries Ordonnée de la droite de pente nulle & %s & psi$^2$/cp\\\\ \n' % y_dp_0
l17='\n\\bfseries Pente de la droite de Horner & %s & psi$^2$/cp/cycle log\\\\ \n' % slope
l18='\n\\bfseries Pseudo-pression intiale & %s & psi$^2$/cp\\\\ \n' % pp0

llast='\\hline\n\\end{tabular}\n\\end{center}\n\\end{table}\n\n\\end{document}'

crlatex=array([l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16,l17,l18,llast])
crlatex.tofile(file='woup.tex', sep='\n')

os.system('latex woup.tex')

os.system('xdvi woup.dvi')

#
# Message de fin
#


print ''
print" - - - - - Tout est bien qui fini bien. - - - - - "
print"                         -"
print"                         -"
print"                         -"
print ''


#>>> a='%s %d %s' % ('\\bfseries k &',k,'& md \\ \\')
#>>> results=array([a,b,c])
#>>> results.tofile(file='resultats.txt', sep='\n')




