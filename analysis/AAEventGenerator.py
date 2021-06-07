#!/usr/bin/env python

import os
import xml.etree.ElementTree as ET
import argparse

pTHat_list = [5., 10., 20., 40., 60., 100.]
current_path = os.getcwd()

for i, new_pT_hat_min in enumerate(pTHat_list[:-1]): 
    with open(current_path+'/../config/jetscape_user_AA200.xml', 'rb') as xml_file: 
        tree = ET.parse(xml_file)
        root = tree.getroot()
        
        name = root.find('outputFilename') 
        file_name = current_path+'/../../JETSCAPE-output/AuAu200/%.6f' %(new_pT_hat_min)
        name.text = file_name
        name.set('updated', 'yes')

        hard = root.find('Hard')
        pythia = hard.find('PythiaGun')
        pT_hat_min = pythia.find('pTHatMin')
        pT_hat_max = pythia.find('pTHatMax')

        pT_hat_min.text = str(new_pT_hat_min)
        pT_hat_min.set('updated', 'yes')
        new_pT_hat_max = pTHat_list[i+1]
        pT_hat_max.text = str(new_pT_hat_max)
        pT_hat_max.set('updated', 'yes')

        tree.write(current_path+'/../config/jetscape_user_AA200.xml', xml_declaration=True, encoding='utf-8')
    print(current_path+'/../build/runJetscape '+current_path+'/../config/jetscape_user_AA200.xml')
    # os.system(current_path+'/../build/runJetscape '+current_path+'/../config/jetscape_user_AA200.xml')

