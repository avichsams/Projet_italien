## creer le fichier avec les valeur de nkr et et deriv_conv_order
import os
import shutil
from os import listdir
from os.path import isfile, join

def cree_repertoir(nkr,deriv_conv_order,deriv_visc_order):

    if not os.path.exists('output'):
        os.makedirs('output')

    chemin_1 =  'output/rung_kutta_'+str(nkr)+'/'
    if not os.path.exists(chemin_1):
        os.makedirs(chemin_1)

    chemin_2=chemin_1+'deriv_conv_order_'+str(deriv_conv_order)+'_'+'deriv_visc_order_'+str(deriv_visc_order)+'/'
    if not os.path.exists(chemin_2):
        os.makedirs(chemin_2)
    
    return chemin_2