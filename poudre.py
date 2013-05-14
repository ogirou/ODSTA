# coding=utf-8
#

import numpy as np
from numpy import *
import matplotlib.pyplot as plt

#
# Calcul du log du temps equivalent Agarwal - 1980
#
# Ezekwe 12.16
#

def agarwal(t,tp):
	te=(tp*t)/(tp+t)
	lt=log(te)
	return te, lt

#
# Algorithme de dérivation d'après Bourdet et al - 1989 
#
# Différences centrées à pas variable
#
# Voir Ezekwe - 2010 - 12C.1 pour le calcul
#

def deriv(lt,te,p):
	pf=p-p[0] # Calcul de l evolution de la pression (rabattement ou remontee)
	vois=int(raw_input('Quel voisin aller chercher pour les différences centrées?(>0)')) 
	print ''	
	dX1=lt[1*vois:-1*vois]-lt[:-2*vois] 
	dX2=lt[2*vois:]-lt[1*vois:-1*vois]
	dp1=pf[1*vois:-1*vois]-pf[:-2*vois]
	dp2=pf[2*vois:]-pf[1*vois:-1*vois]
	dp=(dp1/dX1*dX2+dp2/dX2*dX1)/(dX1+dX2)
	dt=te[1*vois:-1*vois]
	return pf, dp, dt

#
# Dans cette section on veut pouvoir tracer automatiquement une droite de régression
# de pente nulle des valeurs de dp durant l IARF.
#
# Elle devra donner automatiquement la valeur de dp à l'IARF
#
# On chercher à définir les valeurs de temps Agarwal qui correspondent au IARF
#

def der_auto(dt,dp):

	x_dp_0_beg=float(raw_input('Valeur de temps du premier point pour lequel la pente est 0?'))
	print ''
	print "Si l'IARF se poursuit jusqu'à la fin, entrer une valeur après le dernier temps"
	print''
	x_dp_0_end=float(raw_input('Valeur de temps du dernier point pour lequel la pente est 0?'))
	print ''
	
	#
	# On défini deux tuples des index de 'dt' > ou = au début et < ou = à la fin de IARF
	#
	
	x_IARF_deb=np.where(dt>=x_dp_0_beg)
	x_IARF_fin=np.where(dt<=x_dp_0_end)

	#
	# On convertit ces tuples en matrices
	#
	
	xa_IARF_deb=array(x_IARF_deb)
	xa_IARF_fin=array(x_IARF_fin)
	
	#
	# On affiche les 1ères valeurs de chacunes des matrices: ce seront les bornes de IARF
	#
	
	print "Valeur de l'index du début de l'IARF",xa_IARF_deb[0,0]
	print''
	print "Valeur de l'index de fin de l'IARF",xa_IARF_fin[0,-1]
	print''

	#
	# On défini et on affiche les valeurs de te à l'IARF
	#
	# Les '+1' à la fin des deux expressions suivantes permettent de sélectionner toute la longueur
	# de l'index. Ex: si a = array([1,3,5,7,9]), a[4]=9 et a[1:4]=array([3,5,7])
	#
	# On peut supprimer les +1 en cas d'instabilité de calcul des différences centrées
	#
	
	x_IARF=dt[xa_IARF_deb[0,0]:xa_IARF_fin[0,-1]+1]
	y_IARF=dp[xa_IARF_deb[0,0]:xa_IARF_fin[0,-1]+1]
	print "Valeurs de temps de l'IARF (x_IARF)",x_IARF
	print''
	print'Longueurs de x_IARF et y_IARF',len(x_IARF),len(y_IARF)
	#print "Valeurs de dérivée de pression de l'IARF (y_IARF)",y_IARF,"(Contrôle visuel de l'homogénéité)"
	#print''
	
	#
	# polyfit est utilisé trouver l'ordonné à l'origine du polynome de degré 0
	# qui correspond à l'équation de dp à l'IARF 
	#
	
	reg_IARF=np.polyfit(x_IARF,y_IARF,0)
	y_dp_0=reg_IARF[0]
	print "Valeur en ordonnée de la dérivée",y_dp_0
	print''
	
	#
	# On défini la droite qui correspond à l'IARF pour valider l'oào
	#
	
	reg_dp=np.poly1d(reg_IARF)
	#print "Valeurs de la dérivée",reg_dp(x_IARF)
	#print''
	
	return x_IARF, reg_dp(x_IARF), xa_IARF_deb[0,0], xa_IARF_fin[0,-1], y_dp_0

#
# Calcul de la pente de l'IARF sur un Horner plot
#

def penteH(p,tH,ind_IARF_deb,ind_IARF_fin):
	p_1=p[ind_IARF_deb]
	p_2=p[ind_IARF_fin]
	t_1=tH[ind_IARF_deb]
	t_2=tH[ind_IARF_fin]
	slope=(p_2-p_1)/log10(t_1/t_2)
	return slope

#
# Calcul manuel de la pente de l'IARF sur le Horner plot
#

def penteM(tH):
	pp0=int(raw_input('Pseudo-pression initiale du réservoir pp0? (psi2/cp)'))
	slope=int(raw_input('Pente de de la droite de Horner? ((psi2/cp)/cycle_log)'))
	tH_pente=1
	#tH_pente=np.append(tH_pente,tH[12])
	tH_pente=np.append(tH_pente,tH[0])
	p_pente=-slope*log10(tH_pente)+pp0
	return pp0, slope, tH_pente, p_pente

#
# Représentation du graphe log-log de la pression et la dérivée
#

def plot_bourd(te,pf,dt,dp,popp):

	print"Représentation de la",popp,"et de sa dérivé"
	print''	
	plt.loglog(te,pf,'bs-',dt,dp,'gv-')
	if popp=='Pression':
		plt.ylabel('p et dp (psi)')
	else:
		plt.ylabel('pp et dpp (psi^2/cp)')
	plt.xlabel('temps equivalent Agarwal (min)')
	plt.grid(b=True, which='major', linestyle='solid')
	plt.grid(b=True, which='minor', linestyle='dashed')
	#plt.grid(True) 	# maillage simplifié (convient mal au log mais esthétique)
	leg = plt.legend((popp, 'Derivee'),'upper left', shadow=False)
	#plt.savefig('foo.png')
	plt.show()

#
# Représentation du graphe log-log de la pression et la dérivée 
#
# avec la droite de pente nulle à l'IARF
#

def plot_bourd_s(te,pf,dt,dp,x_IARF,reg_dp,popp):

	print"Représentation de la",popp,"et de sa dérivé"
	print''	
	plt.loglog(te,pf,'bs-',dt,dp,'gv-',x_IARF,reg_dp,'r-')
	if popp=='Pression':
		plt.ylabel('p et dp (psi)')
	else:
		plt.ylabel('pp et dpp (psi^2/cp)')
	plt.xlabel('temps equivalent Agarwal (min)')
	plt.grid(b=True, which='major', linestyle='solid')
	plt.grid(b=True, which='minor', linestyle='dashed')
	leg = plt.legend((popp, 'Derivee','Pente nulle'),'upper left', shadow=False)
	plt.show()
	
#
# Représentation du Horner plot avec la droite de même pente que l'IARF
#

def plot_horner(tH,p,tH_pente,p_pente,ylab,popp):

	plt.semilogx(tH,p,'s-',tH_pente,p_pente,'r-')
	plt.ylabel(ylab)
	# pylab.xlim([1,tH[0]])
	plt.xlabel('Fonction du temps de Horner (-)')
	# plt.gca().invert_xaxis()
	plt.grid(b=True, which='major', linestyle='solid')
	plt.grid(b=True, which='minor', linestyle='dashed')
	leg = plt.legend((popp,'Droite de Horner'),'best', shadow=False)
	plt.show()

#
# Horner plot avec k et p0
#

def plot_horner_f(tH,p,tH_pente,p_pente,ylab,popp,k,p0):

	plt.semilogx(tH,p,'s-',tH_pente,p_pente,'r-')
	plt.ylabel(ylab)
	# pylab.xlim([1,tH[0]])
	plt.xlabel('Fonction du temps de Horner (-)')
	# plt.gca().invert_xaxis()
	plt.grid(b=True, which='major', linestyle='solid')
	plt.grid(b=True, which='minor', linestyle='dashed')
	leg = plt.legend((popp,'Droite de Horner'),'upper right', shadow=False)
	plt.text(tH[-1], p[0]/2, "k = %.4f md\np0 = %d psi" % ( k , p0)   ,
		bbox={'facecolor':'white', 'alpha':1, 'pad':15})
	plt.show()

#
# Bourdet plot avec k_bourd et p0
#

def plot_bourd_f(te,pf,dt,dp,x_IARF,reg_dp,popp,k_Bourdet):

	print"Représentation de la",popp,"et de sa dérivé"
	print''	
	plt.loglog(te,pf,'bs-',dt,dp,'gv-',x_IARF,reg_dp,'r-')
	plt.ylabel('p et dp (psi)')
	plt.xlabel('temps equivalent Agarwal (min)')
	plt.grid(b=True, which='major', linestyle='solid')
	plt.grid(b=True, which='minor', linestyle='dashed')
	leg = plt.legend((popp, 'Derivee','Pente nulle'),'best', shadow=False)
	plt.text(te[5], pf[0]/2, "k = %d"%k_Bourdet   ,
		bbox={'facecolor':'white', 'alpha':1, 'pad':15})
	plt.show()

#
# end pour le moment
#

