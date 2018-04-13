#!/usr/bin/env python

#/home/trwood/MCEq_dev_mulossesin/MCEq/

import os
import matplotlib.pyplot as plt
import numpy as np
import math
#os.chdir('/home/trwood/MCEq_newDev_april25_all3in/MCEq')
#os.chdir('/home/trwood/MCEq/')
os.chdir('/home/trwood/MCEq_rc1/')
import sys

sys.path.append('/home/trwood/MCEq_rc1/')
sys.path.append('/home/trwood/MCEq_rc1/MCEq')


#sys.path.append('/home/trwood/MCEq_newDev_april25_all3in')
#sys.path.append('/home/trwood/MCEq_newDev_april25_all3in/MCEq')
#sys.path.append('/home/trwood/MCEq_newDev_april25_all3in/MCEq/MCEq')


import argparse

#import solver related modules
from MCEq.core import MCEqRun
from mceq_config import config
#import primary model choices
import CRFluxModels as pm


from MCEq.misc import  cornertext
#from MCEq.misc import set_ticks_y, cornertext
import matplotlib
import cPickle as pickle
from os.path import join
import calendar

#data_dir = '/Users/afedynitch/Documents/KIT/artifacts/matrix_method/data_files'
#data_dir = '/Users/trwood/MCEq/MCEq/data'
#data_dir = '/home/trwood/MCEq_newDev_april25_all3in/MCEq/data'

data_dir = '/home/trwood/MCEq_rc1/data'


parser = argparse.ArgumentParser(description='script for calculating neutrino fluxes')
parser.add_argument('theta', help='the zenith angle to calculate the density profile with')
#parser.add_argument('hadmodel', help='hadronic model to use')
#parser.add_argument('CRFluxModel', help='CRFlux model to use')
#parser.add_argument('ofile', help='name of output file for the table')
#parser.add_argument('-s', '--sim', help='turn on extra processing for sim files', action='store_true')
args = parser.parse_args()


coszenith=str(args.theta)

print coszenith
#translate the input cos(theta) (which was in radians
# 1) find theta ( theta = arccos(imput_value)   0.1 = cos(theta), theta = arccos(0.1)
    # we read in angles interms of cos(theta) so that we get equal spacing in the bins (stradians)
# 2) theta from radians into degrees

zen_to_deg_interm = np.arccos(float(coszenith))
print zen_to_deg_interm
zenith_theta = math.degrees(zen_to_deg_interm)
print zen_to_deg_interm
print zenith_theta
theta_use = float(zenith_theta)
#*************** INITIALization of config.. only do this once per script run  ******************
#WARNING.  I am not setting out config here anymore.  changes made to mceq_config.py will
#propogate into this scirpt.

#note! theta_deg does not have a default in mceq_config anymore.. accourdingly set_theta_deg seems to be gone as well
#now I am setting it in the initialization of the mceq_run

def mceq_config_without(key_list):
    r = dict(config)  # make a copy
    for key in key_list:
        del r[key]
    return r

 # This is the modification function, it can explicitely depend on energy
# Particle production matrices will be multiplied with it
def mod_linear(xmat, e_grid, a,*args):
    return (1. + a)*np.ones_like(xmat)

def apply_mod_p(fac_pi_pl, fac_pi_mi, MCEq_obj):

    # Delay re-initialization until all modifications are set
    r = MCEq_obj.set_mod_pprod(2212,  211, mod_linear, (fac_pi_pl,), delay_init=True)
    r += MCEq_obj.set_mod_pprod(2212, -211, mod_linear, (fac_pi_mi,), delay_init=True)

    r = MCEq_obj.set_mod_pprod(2112,  211, mod_linear, (fac_pi_pl,), delay_init=True)
    r += MCEq_obj.set_mod_pprod(2112, -211, mod_linear, (fac_pi_mi,), delay_init=True)

    if r > 0:
        MCEq_obj._init_default_matrices(skip_D_matrix=True)


def apply_mod(fac_pi_pl, fac_pi_mi, fac_k_pl, fac_k_mi, fac_k0, fac_k0s, fac_k0l, MCEq_obj):
    #month = str(month)
    # Delay re-initialization until all modifications are set
    r = MCEq_obj.set_mod_pprod(2212,  211, mod_linear, (fac_pi_pl,), delay_init=True)
    r += MCEq_obj.set_mod_pprod(2212, -211, mod_linear, (fac_pi_mi,), delay_init=True)
    r += MCEq_obj.set_mod_pprod(2212,  321, mod_linear, (fac_k_pl,), delay_init=True)
    r += MCEq_obj.set_mod_pprod(2212, -321, mod_linear, (fac_k_mi,), delay_init=True)
    r += MCEq_obj.set_mod_pprod(2212, 311, mod_linear, (fac_k0,), delay_init=True)
    r += MCEq_obj.set_mod_pprod(2212, -311, mod_linear, (fac_k0,), delay_init=True)
    r += MCEq_obj.set_mod_pprod(2212, 130, mod_linear, (fac_k0l,), delay_init=True)
    r += MCEq_obj.set_mod_pprod(2212, 310, mod_linear, (fac_k0s,), delay_init=True)




    r = MCEq_obj.set_mod_pprod(2112,  211, mod_linear, (fac_pi_pl,), delay_init=True)
    r += MCEq_obj.set_mod_pprod(2112, -211, mod_linear, (fac_pi_mi,), delay_init=True)
    r += MCEq_obj.set_mod_pprod(2112,  321, mod_linear, (fac_k_pl,), delay_init=True)
    r += MCEq_obj.set_mod_pprod(2112, -321, mod_linear, (fac_k_mi,), delay_init=True)
    r += MCEq_obj.set_mod_pprod(2112, 311, mod_linear, (fac_k0,), delay_init=True)
    r += MCEq_obj.set_mod_pprod(2112, -311, mod_linear, (fac_k0,), delay_init=True)
    r += MCEq_obj.set_mod_pprod(2112, 130, mod_linear, (fac_k0l,), delay_init=True)
    r += MCEq_obj.set_mod_pprod(2112, 310, mod_linear, (fac_k0s,), delay_init=True)


    if r > 0:
        MCEq_obj._init_default_matrices(skip_D_matrix=True)




#pmodel = (pm.HillasGaisser2012, "H3a")  #  (pm.Thunman, None)   # (pm.GaisserStanevTilav, "3-gen") # (PolyGonato, False)
pmodel = (pm.HillasGaisser2012, "H3a")  #  (pm.Thunman, None)   # (pm.GaisserStanevTilav, "3-gen") # (PolyGonato, False)

mceq_run_January = MCEqRun(interaction_model='DPMJET-III', density_model=('MSIS00_IC',('SouthPole','January')),
                       primary_model=pmodel, theta_deg=theta_use, **mceq_config_without(['density_model']) )

mceq_run_January.unset_mod_pprod(dont_fill=False)
#Power of energy to scale the flux
mag = 3

#obtain energy grid (nver changes) of the solution for the x-axis of the plots
e_grid_January = mceq_run_January.e_grid

#e_grid = mceq_run_January.y.e_grid
#mceq_run.unset_mod_pprod()

#mod_0 = (0.0,0.0,0.0,0.0,0.0,0.0,0.0,mceq_run_January )
mod_0 = (-1.0,-1.0, mceq_run_January )
apply_mod_p(*mod_0)
mceq_run_January.y.print_mod_pprod()
#end New stuff
mceq_run_January.solve()

mceq_run_January.unset_mod_pprod(dont_fill=False)


#output during the MCEqRun initialization. 
flux_January = {}
#_conv means conventional (mostly pions and kaons)

#flux_January['numu_total'] = (mceq_run_January.get_solution('total_numu', mag)
#                      + mceq_run_January.get_solution('total_antinumu', mag))


#flux_January['nue_total'] = (mceq_run_January.get_solution('total_nue', mag)
#                     + mceq_run_January.get_solution('total_antinue', mag))



#to use .. 

#part one nu/anti nu ratio
flux_January['numu_only'] = (mceq_run_January.get_solution('total_numu', mag))

flux_January['antinumu_only'] = (mceq_run_January.get_solution('total_antinumu', mag))

flux_January['nue_only'] = (mceq_run_January.get_solution('total_nue', mag))

flux_January['antinue_only'] = (mceq_run_January.get_solution('total_antinue', mag))

#part two  k/pi ratio
flux_January['numu_from_kaon'] = (mceq_run_January.get_solution('k_numu', mag))

flux_January['antinumu_from_kaon'] = (mceq_run_January.get_solution('k_antinumu', mag))

flux_January['numu_from_pion'] = (mceq_run_January.get_solution('pi_numu', mag))

flux_January['antinumu_from_pion'] = (mceq_run_January.get_solution('pi_antinumu', mag))

mceq_run_January.unset_mod_pprod(dont_fill=False)
#mod_0 = (-1.0,-1.0,0.0,0.0,0.0,0.0,0.0)
#apply_mod(*mod_0)
#**********END JANUARY
#**** Febuary
mceq_run_February = MCEqRun(interaction_model='DPMJET-III', density_model=('MSIS00_IC',('SouthPole','February')),
                        primary_model=pmodel, theta_deg=theta_use, **mceq_config_without(['density_model']))


mceq_run_February.unset_mod_pprod(dont_fill=False)
#mceq_run_January.unset_mod_pprod(dont_fill=False)

#mod_0 = (0.0,0.0,0.0,0.0,0.0,0.0,0.0, mceq_run_February )
mod_0 = (-1.0,-1.0,mceq_run_February )
apply_mod_p(*mod_0)
#end New stuff
mceq_run_February.solve()


flux_February = {}


flux_February['numu_only'] = (mceq_run_February.get_solution('total_numu', mag))

flux_February['antinumu_only'] = (mceq_run_February.get_solution('total_antinumu', mag))

flux_February['nue_only'] = (mceq_run_February.get_solution('total_nue', mag))

flux_February['antinue_only'] = (mceq_run_February.get_solution('total_antinue', mag))




flux_February['numu_from_kaon'] = (mceq_run_February.get_solution('k_numu', mag))

flux_February['antinumu_from_kaon'] = (mceq_run_February.get_solution('k_antinumu', mag))

flux_February['numu_from_pion'] = (mceq_run_February.get_solution('pi_numu', mag))

flux_February['antinumu_from_pion'] = (mceq_run_February.get_solution('pi_antinumu', mag))
mceq_run_February.unset_mod_pprod(dont_fill=False)
#end feb
#march
mceq_run_March = MCEqRun(interaction_model='DPMJET-III', density_model=('MSIS00_IC',('SouthPole','March')),
                        primary_model=pmodel, theta_deg=theta_use, **mceq_config_without(['density_model']))

mceq_run_March.unset_mod_pprod(dont_fill=False)
#mceq_run_January.unset_mod_pprod(dont_fill=False)

#mod_0 = (0.0,0.0,0.0,0.0,0.0,0.0,0.0, mceq_run_February )
mod_0 = (-1.0,-1.0,mceq_run_March)
apply_mod_p(*mod_0)

mceq_run_March.solve()

#Power of energy to scale the flux


#Retrieve the flux at the surface by using the aliases which were listed in the
#output during the MCEqRun initialization. 
flux_March = {}
#_conv means conventional (mostly pions and kaons)
# _pr means prompt (the mother of the muon had a critical energy
# higher than a D meson. Includes all charm and direct resonance
# contribution)

#to use .. 


flux_March['numu_only'] = (mceq_run_March.get_solution('total_numu', mag))

flux_March['antinumu_only'] = (mceq_run_March.get_solution('total_antinumu', mag))

flux_March['nue_only'] = (mceq_run_March.get_solution('total_nue', mag))

flux_March['antinue_only'] = (mceq_run_March.get_solution('total_antinue', mag))




flux_March['numu_from_kaon'] = (mceq_run_March.get_solution('k_numu', mag))

flux_March['antinumu_from_kaon'] = (mceq_run_March.get_solution('k_antinumu', mag))

flux_March['numu_from_pion'] = (mceq_run_March.get_solution('pi_numu', mag))

flux_March['antinumu_from_pion'] = (mceq_run_March.get_solution('pi_antinumu', mag))

mceq_run_March.unset_mod_pprod(dont_fill=False)
#end March
#April
mceq_run_April = MCEqRun(interaction_model='DPMJET-III', density_model=('MSIS00_IC',('SouthPole','April')),
                        primary_model=pmodel, theta_deg=theta_use, **mceq_config_without(['density_model']))

mceq_run_April.unset_mod_pprod(dont_fill=False)

mod_0 = (-1.0,-1.0,mceq_run_April)
apply_mod_p(*mod_0)

mceq_run_April.solve()


#Retrieve the flux at the surface by using the aliases which were listed in the
#output during the MCEqRun initialization. 
flux_April = {}


# since there are no conventional tau neutrinos, prompt=total

#to use
flux_April['numu_only'] = (mceq_run_April.get_solution('total_numu', mag))

flux_April['antinumu_only'] = (mceq_run_April.get_solution('total_antinumu', mag))

flux_April['nue_only'] = (mceq_run_April.get_solution('total_nue', mag))

flux_April['antinue_only'] = (mceq_run_April.get_solution('total_antinue', mag))


flux_April['numu_from_kaon'] = (mceq_run_April.get_solution('k_numu', mag))

flux_April['antinumu_from_kaon'] = (mceq_run_April.get_solution('k_antinumu', mag))

flux_April['numu_from_pion'] = (mceq_run_April.get_solution('pi_numu', mag))

flux_April['antinumu_from_pion'] = (mceq_run_April.get_solution('pi_antinumu', mag))
mceq_run_April.unset_mod_pprod(dont_fill=False)
## END APRIL **************
## Start May

mceq_run_May= MCEqRun(interaction_model='DPMJET-III', density_model=('MSIS00_IC',('SouthPole','May')),
                        primary_model=pmodel, theta_deg=theta_use, **mceq_config_without(['density_model']))
mceq_run_May.unset_mod_pprod(dont_fill=False)

mod_0 = (-1.0,-1.0,mceq_run_May)
apply_mod_p(*mod_0)


mceq_run_May.solve()


#Retrieve the flux at the surface by using the aliases which were listed in the
#output during the MCEqRun initialization. 
flux_May = {}
#_conv means conventional (mostly pions and kaons)

#to use .. 


flux_May['numu_only'] = (mceq_run_May.get_solution('total_numu', mag))

flux_May['antinumu_only'] = (mceq_run_May.get_solution('total_antinumu', mag))

flux_May['nue_only'] = (mceq_run_May.get_solution('total_nue', mag))

flux_May['antinue_only'] = (mceq_run_May.get_solution('total_antinue', mag))




flux_May['numu_from_kaon'] = (mceq_run_May.get_solution('k_numu', mag))

flux_May['antinumu_from_kaon'] = (mceq_run_May.get_solution('k_antinumu', mag))

flux_May['numu_from_pion'] = (mceq_run_May.get_solution('pi_numu', mag))

flux_May['antinumu_from_pion'] = (mceq_run_May.get_solution('pi_antinumu', mag))

mceq_run_May.unset_mod_pprod(dont_fill=False)
## END MAY ***************
#  START JUNE
mceq_run_June= MCEqRun(interaction_model='DPMJET-III', density_model=('MSIS00_IC',('SouthPole','June')),
                        primary_model=pmodel, theta_deg=theta_use, **mceq_config_without(['density_model']))


mceq_run_June.unset_mod_pprod(dont_fill=False)

mod_0 = (-1.0,-1.0,mceq_run_June)
apply_mod_p(*mod_0)

mceq_run_June.solve()

#output during the MCEqRun initialization. 
flux_June = {}
#_conv means conventional (mostly pions and kaons)


flux_June['numu_only'] = (mceq_run_June.get_solution('total_numu', mag))

flux_June['antinumu_only'] = (mceq_run_June.get_solution('total_antinumu', mag))

flux_June['nue_only'] = (mceq_run_June.get_solution('total_nue', mag))

flux_June['antinue_only'] = (mceq_run_June.get_solution('total_antinue', mag))




flux_June['numu_from_kaon'] = (mceq_run_June.get_solution('k_numu', mag))

flux_June['antinumu_from_kaon'] = (mceq_run_June.get_solution('k_antinumu', mag))

flux_June['numu_from_pion'] = (mceq_run_June.get_solution('pi_numu', mag))

flux_June['antinumu_from_pion'] = (mceq_run_June.get_solution('pi_antinumu', mag))

mceq_run_June.unset_mod_pprod(dont_fill=False)
#end June **********************
#Start July  ****
mceq_run_July= MCEqRun(interaction_model='DPMJET-III', density_model=('MSIS00_IC',('SouthPole','July')),
                        primary_model=pmodel, theta_deg=theta_use, **mceq_config_without(['density_model']))

mceq_run_July.unset_mod_pprod(dont_fill=False)

mod_0 = (-1.0,-1.0,mceq_run_July)
apply_mod_p(*mod_0)

mceq_run_July.solve()

#Retrieve the flux at the surface by using the aliases which were listed in the
#output during the MCEqRun initialization. 
flux_July = {}

#to use .. 
flux_July['numu_only'] = (mceq_run_July.get_solution('total_numu', mag))

flux_July['antinumu_only'] = (mceq_run_July.get_solution('total_antinumu', mag))

flux_July['nue_only'] = (mceq_run_July.get_solution('total_nue', mag))

flux_July['antinue_only'] = (mceq_run_July.get_solution('total_antinue', mag))




flux_July['numu_from_kaon'] = (mceq_run_July.get_solution('k_numu', mag))

flux_July['antinumu_from_kaon'] = (mceq_run_July.get_solution('k_antinumu', mag))

flux_July['numu_from_pion'] = (mceq_run_July.get_solution('pi_numu', mag))

flux_July['antinumu_from_pion'] = (mceq_run_July.get_solution('pi_antinumu', mag))
mceq_run_July.unset_mod_pprod(dont_fill=False)
### *********END JULY
#*** Start August
mceq_run_August= MCEqRun(interaction_model='DPMJET-III', density_model=('MSIS00_IC',('SouthPole','August')),
                        primary_model=pmodel, theta_deg=theta_use, **mceq_config_without(['density_model']))


mceq_run_August.unset_mod_pprod(dont_fill=False)

mod_0 = (-1.0,-1.0,mceq_run_August)
apply_mod_p(*mod_0)

mceq_run_August.solve()


#Retrieve the flux at the surface by using the aliases which were listed in the
#output during the MCEqRun initialization. 
flux_August = {}

# to use ..

flux_August['numu_only'] = (mceq_run_August.get_solution('total_numu', mag))

flux_August['antinumu_only'] = (mceq_run_August.get_solution('total_antinumu', mag))

flux_August['nue_only'] = (mceq_run_August.get_solution('total_nue', mag))

flux_August['antinue_only'] = (mceq_run_August.get_solution('total_antinue', mag))




flux_August['numu_from_kaon'] = (mceq_run_August.get_solution('k_numu', mag))

flux_August['antinumu_from_kaon'] = (mceq_run_August.get_solution('k_antinumu', mag))

flux_August['numu_from_pion'] = (mceq_run_August.get_solution('pi_numu', mag))

flux_August['antinumu_from_pion'] = (mceq_run_August.get_solution('pi_antinumu', mag))
mceq_run_August.unset_mod_pprod(dont_fill=False)
###*********END OF AUGUST
#*SEPTEMBER
mceq_run_September= MCEqRun(interaction_model='DPMJET-III', density_model=('MSIS00_IC',('SouthPole','September')),
                        primary_model=pmodel, theta_deg=theta_use, **mceq_config_without(['density_model']))

mceq_run_September.unset_mod_pprod(dont_fill=False)

mod_0 = (-1.0,-1.0,mceq_run_September)
apply_mod_p(*mod_0)


mceq_run_September.solve()


#Retrieve the flux at the surface by using the aliases which were listed in the
#output during the MCEqRun initialization. 
flux_September = {}

#to use .. 


flux_September['numu_only'] = (mceq_run_September.get_solution('total_numu', mag))

flux_September['antinumu_only'] = (mceq_run_September.get_solution('total_antinumu', mag))

flux_September['nue_only'] = (mceq_run_September.get_solution('total_nue', mag))

flux_September['antinue_only'] = (mceq_run_September.get_solution('total_antinue', mag))




flux_September['numu_from_kaon'] = (mceq_run_September.get_solution('k_numu', mag))

flux_September['antinumu_from_kaon'] = (mceq_run_September.get_solution('k_antinumu', mag))

flux_September['numu_from_pion'] = (mceq_run_September.get_solution('pi_numu', mag))

flux_September['antinumu_from_pion'] = (mceq_run_September.get_solution('pi_antinumu', mag))
mceq_run_September.unset_mod_pprod(dont_fill=False)

## END September ****
#  start October

mceq_run_October= MCEqRun(interaction_model='DPMJET-III', density_model=('MSIS00_IC',('SouthPole','October')),
                        primary_model=pmodel, theta_deg=theta_use, **mceq_config_without(['density_model']))



mceq_run_October.unset_mod_pprod(dont_fill=False)

mod_0 = (-1.0,-1.0,mceq_run_October)
apply_mod_p(*mod_0)


#output during the MCEqRun initialization. 
flux_October = {}

#to use .. 


flux_October['numu_only'] = (mceq_run_October.get_solution('total_numu', mag))

flux_October['antinumu_only'] = (mceq_run_October.get_solution('total_antinumu', mag))

flux_October['nue_only'] = (mceq_run_October.get_solution('total_nue', mag))

flux_October['antinue_only'] = (mceq_run_October.get_solution('total_antinue', mag))




flux_October['numu_from_kaon'] = (mceq_run_October.get_solution('k_numu', mag))

flux_October['antinumu_from_kaon'] = (mceq_run_October.get_solution('k_antinumu', mag))

flux_October['numu_from_pion'] = (mceq_run_October.get_solution('pi_numu', mag))

flux_October['antinumu_from_pion'] = (mceq_run_October.get_solution('pi_antinumu', mag))
mceq_run_October.unset_mod_pprod(dont_fill=False)

# ENd OCTOBER ******************
# ** start November

mceq_run_November= MCEqRun(interaction_model='DPMJET-III', density_model=('MSIS00_IC',('SouthPole','November')),
                        primary_model=pmodel, theta_deg=theta_use, **mceq_config_without(['density_model']))

mceq_run_November.unset_mod_pprod(dont_fill=False)

mod_0 = (-1.0,-1.0,mceq_run_November)
apply_mod_p(*mod_0)
mceq_run_November.solve()
#Retrieve the flux at the surface by using the aliases which were listed in the
#output during the MCEqRun initialization. 
flux_November = {}
#_conv means conventional (mostly pions and kaons)
#to use .. 


flux_November['numu_only'] = (mceq_run_November.get_solution('total_numu', mag))

flux_November['antinumu_only'] = (mceq_run_November.get_solution('total_antinumu', mag))

flux_November['nue_only'] = (mceq_run_November.get_solution('total_nue', mag))

flux_November['antinue_only'] = (mceq_run_November.get_solution('total_antinue', mag))




flux_November['numu_from_kaon'] = (mceq_run_November.get_solution('k_numu', mag))

flux_November['antinumu_from_kaon'] = (mceq_run_November.get_solution('k_antinumu', mag))

flux_November['numu_from_pion'] = (mceq_run_November.get_solution('pi_numu', mag))

flux_November['antinumu_from_pion'] = (mceq_run_November.get_solution('pi_antinumu', mag))
mceq_run_November.unset_mod_pprod(dont_fill=False)

#end november
#start december
mceq_run_December= MCEqRun(interaction_model='DPMJET-III', density_model=('MSIS00_IC',('SouthPole','December')),
                        primary_model=pmodel, theta_deg=theta_use, **mceq_config_without(['density_model']))

mceq_run_December.unset_mod_pprod(dont_fill=False)

mod_0 = (-1.0,-1.0,mceq_run_December)
apply_mod_p(*mod_0)

mceq_run_December.solve()
#Retrieve the flux at the surface by using the aliases which were listed in the
#output during the MCEqRun initialization. 
flux_December = {}

#to use ..

flux_December['numu_only'] = (mceq_run_December.get_solution('total_numu', mag))

flux_December['antinumu_only'] = (mceq_run_December.get_solution('total_antinumu', mag))

flux_December['nue_only'] = (mceq_run_December.get_solution('total_nue', mag))

flux_December['antinue_only'] = (mceq_run_December.get_solution('total_antinue', mag))




flux_December['numu_from_kaon'] = (mceq_run_December.get_solution('k_numu', mag))

flux_December['antinumu_from_kaon'] = (mceq_run_December.get_solution('k_antinumu', mag))

flux_December['numu_from_pion'] = (mceq_run_December.get_solution('pi_numu', mag))

flux_December['antinumu_from_pion'] = (mceq_run_December.get_solution('pi_antinumu', mag))
mceq_run_December.unset_mod_pprod(dont_fill=False)

#end of December
#****************end of flux maddess, ie
####END BULLSHIT 12 time part.  this can be cleaned up. the rest of the script we are likely stuck with

#get path of the home directory + Desktop
save_pdf = True
desktop = os.path.join(os.path.expanduser("~"),'Desktop')
#shall i make a folder for each type of output i want?  probably. sighhhhhhhhh
#desktop = os.path.join(os.path.expanduser("~"),'Desktop/chris_tables')
# number was 7213
'''mceq_run_January.pdg2pref[7213].inverse_decay_length(mceq_run_January.e_grid_January)

mceq_run_February.pdg2pref[7213].inverse_decay_length(mceq_run_February.e_grid_February)

mceq_run_March.pdg2pref[7213].inverse_decay_length(mceq_run_March.e_grid_March)

mceq_run_April.pdg2pref[7213].inverse_decay_length(mceq_run_April.e_grid_April)

mceq_run_May.pdg2pref[7213].inverse_decay_length(mceq_run_May.e_grid_May)

mceq_run_June.pdg2pref[7213].inverse_decay_length(mceq_run_June.e_grid_June)

mceq_run_July.pdg2pref[7213].inverse_decay_length(mceq_run_July.e_grid_July)

mceq_run_August.pdg2pref[7213].inverse_decay_length(mceq_run_August.e_grid_August)

mceq_run_September.pdg2pref[7213].inverse_decay_length(mceq_run_September.e_grid_September)

mceq_run_October.pdg2pref[7213].inverse_decay_length(mceq_run_October.e_grid_October)

mceq_run_November.pdg2pref[7213].inverse_decay_length(mceq_run_November.e_grid_November)

mceq_run_December.pdg2pref[7213].inverse_decay_length(mceq_run_December.e_grid_December)
'''

#make these folders!!! how to just get them made if they don't exist yet? not sure ... bah. 
#folders for first part (nu/anti nu ratios, both for nu_e and nu_mu types)
desktop_numu_all_source = os.path.join(os.path.expanduser("~"),'Desktop/h3a_DPMJETIII_MCEq_rc1_pions_off/nu_antinu_ratio_allSources/numu_all')
desktop_antinumu_all_source = os.path.join(os.path.expanduser("~"),'Desktop/h3a_DPMJETIII_MCEq_rc1_pions_off/nu_antinu_ratio_allSources/antinumu_all')

desktop_nue_all_source = os.path.join(os.path.expanduser("~"),'Desktop/h3a_DPMJETIII_MCEq_rc1_pions_off/nu_antinu_ratio_allSources/nue_all')
desktop_antinue_all_source = os.path.join(os.path.expanduser("~"),'Desktop/h3a_DPMJETIII_MCEq_rc1_pions_off/nu_antinu_ratio_allSources/antinue_all')


#folders for second part (for_k_pi_ratio)
desktop_numu_from_kaon = os.path.join(os.path.expanduser("~"),'Desktop/h3a_DPMJETIII_MCEq_rc1_pions_off/for_k_pi_ratio/numu_from_kaon')
desktop_antinumu_from_kaon = os.path.join(os.path.expanduser("~"),'Desktop/h3a_DPMJETIII_MCEq_rc1_pions_off/for_k_pi_ratio/antinumu_from_kaon')

desktop_numu_from_pion = os.path.join(os.path.expanduser("~"),'Desktop/h3a_DPMJETIII_MCEq_rc1_pions_off/for_k_pi_ratio/numu_from_pion')
desktop_antinumu_from_pion = os.path.join(os.path.expanduser("~"),'Desktop/h3a_DPMJETIII_MCEq_rc1_pions_off/for_k_pi_ratio/antinumu_from_pion')


#**********************************************************************************************
#A TEST PLOT .. delete this after tests complete

#for pref, lab in [('mu_',r'\mu'), ('numu_',r'\nu_\mu'), ('nue_',r'\nu_e')]:
#    plt.figure(figsize=(10., 6))
#    plt.loglog(e_grid, flux[pref + 'total'], color='k', ls='-', lw=1.5)
#    plt.loglog(e_grid, flux[pref + 'conv'], color='b', ls='-.', lw=1.5,
#               label=r'conventional ${0}$'.format(lab))
#    plt.loglog(e_grid, flux[pref + 'pr'], color='r',ls='--', lw=1.5,
#               label='prompt ${0}$'.format(lab))
#    plt.xlim(1,1e9)
#    plt.ylim(1e-6,1)
#    plt.xlabel(r"$E_{{{0}}}$ [GeV]".format(lab))
#    plt.ylabel(r"$\Phi_{" + lab + "}$ (E/GeV)$^{" + str(mag) +" }$" +
#               "(cm$^{2}$ s sr GeV)$^{-1}$")
#    plt.legend(loc='upper right')
#    plt.tight_layout()
#    if save_pdf: plt.savefig(os.path.join(desktop, pref + str(zenith_theta) + 'fluxtanny.eps'))
#**********************************************************************************************
#test March ONly] 
#desktop_numu_march =  os.path.join(os.path.expanduser("~"),'Desktop/GH_DPMJETIII_MCEq_rc1pions_off/nu_antinu_ratio_allSources/numu_march')

#np.savetxt(open(os.path.join(desktop_numu_march, 'coszenith_' + str(coszenith) + '_numu_avg_march.txt'),'w'),
#zip(flux_March['numu_only']),fmt='%6.5E',header=('lepton flux_name scaled with E**{0}').format(mag))

#desktop_antinumu_march =  os.path.join(os.path.expanduser("~"),'Desktop/GH_DPMJETIII_MCEq_rc1pions_off/nu_antinu_ratio_allSources/antinumu_march')

#np.savetxt(open(os.path.join(desktop_antinumu_march, 'coszenith_' + str(coszenith) + '_antinumu_avg_march.txt'),'w'),
#zip(flux_March['antinumu_only']),fmt='%6.5E',header=('lepton flux_name scaled with E**{0}').format(mag))


#FINAL AVERAGES FOR THIS PARTICULAR COS(ZENITH) *****************************
#likely could make only one dictionary object called 'averages' and give them whataever name i like
flux_avg_year = {}

flux_avg_year['from_kaons_antinumu_only'] = ( np.array(flux_January['antinumu_from_kaon']) +  np.array(flux_February['antinumu_from_kaon']) +  np.array(flux_March['antinumu_from_kaon']) +
         np.array(flux_April['antinumu_from_kaon']) +  np.array(flux_May['antinumu_from_kaon']) +  np.array(flux_June['antinumu_from_kaon']) +
         np.array(flux_July['antinumu_from_kaon']) +  np.array(flux_August['antinumu_from_kaon']) +  np.array(flux_September['antinumu_from_kaon']) +
         np.array(flux_October['antinumu_from_kaon']) +  np.array(flux_November['antinumu_from_kaon']) +  np.array(flux_December['antinumu_from_kaon']) ) / 12.0

flux_avg_year['from_kaons_numu_only'] = ( np.array(flux_January['numu_from_kaon']) +  np.array(flux_February['numu_from_kaon']) +  np.array(flux_March['numu_from_kaon']) +
        np.array(flux_April['numu_from_kaon']) + np.array(flux_May['numu_from_kaon']) +  np.array(flux_June['numu_from_kaon']) +  np.array(flux_July['numu_from_kaon']) +
        np.array(flux_August['numu_from_kaon']) +  np.array(flux_September['numu_from_kaon']) +  np.array(flux_October['numu_from_kaon']) +  np.array(flux_November['numu_from_kaon']) +
        np.array(flux_December['numu_from_kaon']) )/ 12.0

flux_avg_year['from_pions_antinumu_only'] =  ( np.array(flux_January['antinumu_from_pion']) +  np.array(flux_February['antinumu_from_pion']) +  np.array(flux_March['antinumu_from_pion']) +
        np.array(flux_April['antinumu_from_kaon']) +  np.array(flux_May['antinumu_from_pion']) +  np.array(flux_June['antinumu_from_pion']) +  np.array(flux_July['antinumu_from_pion']) +
        np.array(flux_August['antinumu_from_pion']) +  np.array(flux_September['antinumu_from_pion']) +  np.array(flux_October['antinumu_from_pion']) +  np.array(flux_November['antinumu_from_pion']) +
        np.array(flux_December['antinumu_from_pion']) ) / 12.0

flux_avg_year['from_pions_numu_only'] =  (np.array(flux_January['numu_from_pion']) +  np.array(flux_February['numu_from_pion']) +  np.array(flux_March['numu_from_pion']) +
        np.array(flux_April['numu_from_pion']) +  np.array(flux_May['numu_from_pion']) +  np.array(flux_June['numu_from_pion']) +  np.array(flux_July['numu_from_pion']) +
        np.array(flux_August['numu_from_pion']) +  np.array(flux_September['numu_from_pion']) +  np.array(flux_October['numu_from_pion']) +  np.array(flux_November['numu_from_pion']) +
        np.array(flux_December['numu_from_pion']) ) / 12.0

#for nu/antinu ratios
#flux_avg_year['total_numu_from_all_sources'] = flux['numu_only'] + all of the months avg
flux_avg_year['year_avg_numu_from_all_sources'] = (np.array(flux_January['numu_only']) +  np.array(flux_February['numu_only']) +  np.array(flux_March['numu_only']) +
        np.array(flux_April['numu_only']) +  np.array(flux_May['numu_only']) +  np.array(flux_June['numu_only']) +  np.array(flux_July['numu_only']) +
        np.array(flux_August['numu_only']) +  np.array(flux_September['numu_only']) +  np.array(flux_October['numu_only']) +  np.array(flux_November['numu_only']) +
        np.array(flux_December['numu_only']) ) / 12.0
#'antinumu_only'
flux_avg_year['year_avg_antinumu_from_all_sources'] =  (np.array(flux_January['antinumu_only']) +  np.array(flux_February['antinumu_only']) +  np.array(flux_March['antinumu_only']) +
        np.array(flux_April['antinumu_only']) +  np.array(flux_May['antinumu_only']) +  np.array(flux_June['antinumu_only']) +  np.array(flux_July['antinumu_only']) +
        np.array(flux_August['antinumu_only']) +  np.array(flux_September['antinumu_only']) +  np.array(flux_October['antinumu_only']) +  np.array(flux_November['antinumu_only']) +
        np.array(flux_December['antinumu_only']) ) / 12.0

flux_avg_year['year_avg_nue_all_sources'] = (np.array(flux_January['nue_only']) +  np.array(flux_February['nue_only']) +  np.array(flux_March['nue_only']) +
        np.array(flux_April['nue_only']) +  np.array(flux_May['nue_only']) +  np.array(flux_June['nue_only']) +  np.array(flux_July['nue_only']) +
        np.array(flux_August['nue_only']) +  np.array(flux_September['nue_only']) +  np.array(flux_October['nue_only']) +  np.array(flux_November['nue_only']) +
        np.array(flux_December['nue_only']) ) / 12.0


flux_avg_year['year_avg_anti_nue_all_sources'] = (np.array(flux_January['antinue_only']) +  np.array(flux_February['antinue_only']) +  np.array(flux_March['antinue_only']) +
        np.array(flux_April['antinue_only']) +  np.array(flux_May['antinue_only']) +  np.array(flux_June['antinue_only']) +  np.array(flux_July['antinue_only']) +
        np.array(flux_August['antinue_only']) +  np.array(flux_September['antinue_only']) +  np.array(flux_October['antinue_only']) +  np.array(flux_November['antinue_only']) +
        np.array(flux_December['antinue_only']) ) / 12.0



# Now print out the averages as txt files
# ***************************** ***************************** *****************************
#THESE WILL GIVE EAVH MONTH ALONE TOTALS  Egrid, flux
#print each of these
#0)  once, the energy grid (bin size of energy grid)
# **** first part for nu/anti_nu ratios
#1) total anti_numu year year avg
#2) total numu year avg
#3) anti nuE year avg
#4) nuE year avg

# **** second part for mu type neutrios only, k/pi ratio
#1) flux_avg_year_Anumu_from_kaons, 
#2) flux_avg_year_numu_from_k
#3) flux_avg_year_antinumu_from_pion
#4) flux_avg_year_numu_from_pion

 ############  EGRID   ########ie. energy grid (bin size of energy grid)
if zenith_theta == 0.1:
    np.savetxt(open(os.path.join(desktop, 'coszenith_' + str(coszenith) + 'Egrid_totalsn.txt'),'w'),
    zip(e_grid_January),
    fmt='%6.5E',
    header=('lepton flux_name scaled with E**{0}. note only has conv plus prompt Order (E' ).format(mag)
    )

    # 1****** antiNUMU TOTALS FULL YEAR AVG  ******* 
np.savetxt(open(os.path.join(desktop_antinumu_all_source, 'coszenith_' + str(coszenith) + '_antinumu_avg_year_totalsn.txt'),'w'),
zip(flux_avg_year['year_avg_antinumu_from_all_sources']),
fmt='%6.5E',
header=('lepton flux_name scaled with E**{0}. note only has conv plus prompt Order (E, antinumu_total' ).format(mag)
)

    #2 ****** NUMU TOTALS FULL YEAR AVG (For a given cos(zenith) *******
np.savetxt(open(os.path.join(desktop_numu_all_source, 'coszenith_' + str(coszenith) + '_numu_avg_year_totalsn.txt'),'w'),
zip(flux_avg_year['year_avg_numu_from_all_sources']),
fmt='%6.5E',
header=('lepton flux_name scaled with E**{0}. note only has conv plus prompt Order ( numu_total)' ).format(mag)
)

    #3****** antiNUe TOTALS FULL YEAR AVG (For a given cos(zenith) *******
np.savetxt(open(os.path.join(desktop_antinue_all_source, 'coszenith_' + str(coszenith) + '_antinue_avg_year_totalsn.txt'),'w'),
zip(flux_avg_year['year_avg_anti_nue_all_sources']),
fmt='%6.5E',
header=('lepton flux_name scaled with E**{0}. note only has conv plus prompt Order (E, antinue_total' ).format(mag)
)

    #4******NUE TOTALS FULL YEAR AVG (For a given cos(zenith) *******
np.savetxt(open(os.path.join(desktop_nue_all_source, 'coszenith_' + str(coszenith) + 'nue_avg_year_totalsn.txt'),'w'),
zip(flux_avg_year['year_avg_nue_all_sources']),
fmt='%6.5E',
header=('lepton flux_name scaled with E**{0}. note only has conv plus prompt Order ( nue_total)' ).format(mag)
)



#  ********************SECOND PART *****************************************************************

#  ********************SECOND PART *****************************************************************
    #  ***** for the numu's only, numu/antinumu flux from pions, kaons, seprately 
# 5) flux_avg_year['from_kaons_antinumu_only']
# 6) flux_avg_year['from_kaons_numu_only']
# 7) flux_avg_year['from_pions_antinumu_only']
# 8) flux_avg_year['from_pions_numu_only']

# 5) flux_avg_year['from_kaons_antinumu_only']
np.savetxt(open(os.path.join(desktop_antinumu_from_kaon, 'coszenith_' + str(coszenith) + '_year_avg_antinumu_flux_from_kaonn.txt'),'w'),
zip(flux_avg_year['from_kaons_antinumu_only']),
fmt='%6.5E',
header=('lepton flux_name scaled with E**{0}. note only has conv plus prompt Order (year_avg_antinumu_from_kaon' ).format(mag)
)
#this egrid_test failsed.  ia mnot sure why. egrid_januray should exist.  error says 'e_grid' not exists
#6) flux_avg_year['from_kaons_numu_only'] plus an egrid test

#np.savetxt(open(os.path.join(desktop_numu_from_kaon_test, 'coszenith_' + str(coszenith) + 'with_egrid_test_year_avg_numu_flux_from_kaon.txt'),'w'),
#zip(e_grid_January, flux_avg_year['from_kaons_numu_only']),
#fmt='%6.5E',
#header=('lepton flux_name scaled with E**{0}. note only has conv plus prompt Order (E, year_avg_numu_from_kaon' ).format(mag)
#)

np.savetxt(open(os.path.join(desktop_numu_from_kaon, 'coszenith_' + str(coszenith) + '_year_avg_numu_flux_from_kaonn.txt'),'w'),
zip(flux_avg_year['from_kaons_numu_only']),
fmt='%6.5E',
header=('lepton flux_name scaled with E**{0}. note only has conv plus prompt Order (year_avg_numu_from_kaon' ).format(mag)
)

# 7) flux_avg_year['from_pions_antinumu_only']
np.savetxt(open(os.path.join(desktop_antinumu_from_pion, 'coszenith_' + str(coszenith) + '_year_avg_antinumu_flux_from_pionsn.txt'),'w'),
zip(flux_avg_year['from_pions_antinumu_only']),
fmt='%6.5E',
header=('lepton flux_name scaled with E**{0}. note only has conv plus prompt Order (year_avg_antinumu_from_pions' ).format(mag)
)

# 8) flux_avg_year['from_pions_numu_only']
np.savetxt(open(os.path.join(desktop_numu_from_pion, 'coszenith_' + str(coszenith) + '_year_avg_numu_flux_from_pionsn.txt'),'w'),
zip(flux_avg_year['from_pions_numu_only']),
fmt='%6.5E',
header=('lepton flux_name scaled with E**{0}. note only has conv plus prompt Order (year_avg_numu_from_pions' ).format(mag)
)



#******* fix this bit
#NO IDEA what this does or if i need it. bit worried. i guess include it.  sighnnnnnn
'''mceq_run_January.pdg2pref[7213].inverse_decay_length(mceq_run_January.e_grid_January)

mceq_run_February.pdg2pref[7213].inverse_decay_length(mceq_run_February.e_grid_February)

mceq_run_March.pdg2pref[7213].inverse_decay_length(mceq_run_March.e_grid_March)

mceq_run_April.pdg2pref[7213].inverse_decay_length(mceq_run_April.e_grid_April)

mceq_run_May.pdg2pref[7213].inverse_decay_length(mceq_run_May.e_grid_May)

mceq_run_June.pdg2pref[7213].inverse_decay_length(mceq_run_June.e_grid_June)

mceq_run_July.pdg2pref[7213].inverse_decay_length(mceq_run_July.e_grid_July)

mceq_run_August.pdg2pref[7213].inverse_decay_length(mceq_run_August.e_grid_August)

mceq_run_September.pdg2pref[7213].inverse_decay_length(mceq_run_September.e_grid_September)

mceq_run_October.pdg2pref[7213].inverse_decay_length(mceq_run_October.e_grid_October)

mceq_run_November.pdg2pref[7213].inverse_decay_length(mceq_run_November.e_grid_November)

mceq_run_December.pdg2pref[7213].inverse_decay_length(mceq_run_December.e_grid_December)
'''






