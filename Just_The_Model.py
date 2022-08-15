###################################################################################################################
###################################################################################################################
##### THIS IS THE SEBM FROM "THE PHYSICS OF HEAT WAVES" PAPER #####################################################
##### First function calculates saturation vapor pressure assuming constant surface pressure ######################
##### Second function accepts radiation and precipitation forcing and returns temperature and soil moisture #######
###################################################################################################################
###################################################################################################################

import numpy as np
from netCDF4 import Dataset

def q_s(TK): 					# ACCEPTS TEMPERAUTRE IN KELVIN
	T = TK - 273.15 			# Now in Celsius
	e_s = 6.11*10**(7.5*T/(237.5+T)) 	# Formula @ https://www.weather.gov/media/epz/wxcalc/vaporPressure.pdf
	q_s = 0.61*e_s/1013			# Converting to kg/kg
	return(q_s) # Given in units of kg H2O/ kg air

def the_model(F_forcing,P_forcing):

	########### F forcing is W/m^2
	########### P forcing is in mm per unit of time

	############# Climatological constants #####################################

	q_clim = 0.005 			# kg/kg 
	TD_clim = 280			# K
	DT = 86400			# in seconds (appropriate for daily forcing)

	########### Tunable parameters #################################################

	r_s = 100			# s/m
	v_H = 10			# W/m^2/K

############ Geometry

	h_s = 0.1			# surface layer depth [m]
	theta_max = 0.4656		# porosity
	N = len(F_forcing)
	
############### Physical Constants

	rho_l = 1000			# denisty of water [kg/m^3]
	rho_s = 1000			# density of dry soil [kg/m^3]
	rho_a = 1.25			# density of air [kg/m^3]
	c_ps  = 1000			# heat capacity of dry soil [J/kg/K]
	L = 2257000			# Latent enthalpy of vaporization [J/kg]

################### Combined parameters

	mu_s = rho_l*h_s*theta_max	# storage capacity of surface layer [kg/m^2]
	C_s = rho_s*h_s*c_ps		# effective heat capacity of land surface [J/m^2/K]

################# State Variables

	T = np.zeros(N)
	m_s = np.zeros(N)
	m_d = np.zeros(N)
	LHF = np.zeros(N)

	T[0] = 280			# Initial condition for surface temperature
	m_s[0] = 1			# ""			surface moisture

	i = 0
	while i < N-1:
		
		######### CALCULATING FLUXES #####################################
			
		H = v_H*(T[i] - TD_clim)		# Sensible Heat Flux [W/m^2]
		VPD = q_s(T[i]) - q_clim		# Vapor Pressure Deficit [kPa]

		if VPD > 0 and T[i] > 273.15:
			ET = rho_a*m_s[i]*VPD/r_s	# Surface Evaporation [kg/m^2/s]
		else:
			ET = 0

		LHF[i] = L*ET

		########## TENDENCY EQUATIONS #####################################

		dT_dt = (F_forcing[i] - H - LHF[i])/C_s
		dms_dt = -(ET/mu_s)

		############# INTEGRATING FORWARD IN TIME #########################

		T[i+1] 		= T[i] + dT_dt*DT
		m_s[i+1]	= m_s[i] + dms_dt*DT + P_forcing[i]/(h_s*theta_max)

		############# If bucket overfills, we take the moisture out #######

		if m_s[i+1] > 1:
			m_s[i+1] = 1

		i+=1

	return(T,m_s,LHF)


