# mesh_utils.py
from antares import *
import os
import re
import numpy as np

class Constants:
    def __init__(self, tip_gap, span):
        self.z_mid_span = (tip_gap + span/2)
        self.z_2inch_tip = (tip_gap - 0.0508)
        self.z_1inch_tip = (tip_gap - 0.0254)
        self.z_025inch_tip = (tip_gap - 0.00635)
        self.z_5mm_tip = (tip_gap - 0.005)
        self.z_25mm_tip = (tip_gap - 0.025)
        self.z_tip_gap = tip_gap


# Mappign the cut selection location
def map_cut(data_type:str,cut_style:str,tip_gap:float,span:float,AoA:int):
    """Function to map the cut selection to the location of the cut plane given the tip gap and span constnats
    Args:
        data_type (str): The style of the cut, either 'plane' or 'cylinder'
        cut_style (str): The location of the cut plane, defined either explicitly or with a %_TE location designating the distance from the trailing edge
        tip_gap (float): the tip gap size
        span (float): the span size
        AoA (int): the angle of attack
    Returns:
        origin (list): The origin of the cut plane
        normal (list): The normal of the cut plane
    """
    z_loc = constants(tip_gap,span)
    if data_type == 'EXTRACT':
        if cut_style.find("midspan") != -1:
            origin = [1.225,0.,z_loc.z_mid_span]
            normal = [0.,0.,1.]
        elif cut_style.find("2inch_tip") != -1:
            origin = [1.225,0.,z_loc.z_2inch_tip]
            normal = [0.,0.,1.]
        elif cut_style.find("1inch_tip") != -1:  
            origin= [1.225,0.,z_loc.z_1inch_tip]
            normal = [0.,0.,1.]
        elif cut_style.find("025inch_tip") != -1:  
            origin = [1.225,0.,z_loc.z_025inch_tip]
            normal = [0.,0.,1.]
        elif cut_style.find("25mm_tip") != -1:  
            origin= [1.225,0.,z_loc.z_25mm_tip]
            normal = [0.,0.,1.]
        elif cut_style.find("5mm_tip") != -1:  
            origin= [1.225,0.,z_loc.z_5mm_tip]
            normal = [0.,0.,1.]
        elif cut_style.find('PIV1') != -1:
            x,y,z = 1.42222035, 0, z_loc.z_mid_span
            origin = [x,y,z]
            normal  = [1,0,0]
        elif cut_style.find('PIV2') != -1:
            x,y,z = 1.48172998, 0, z_loc.z_mid_span
            origin = [x,y,z]
            normal  = [1,0,0]
        elif cut_style.find('PIV3') != -1:
            x,y,z = 1.5641908, 0, z_loc.z_mid_span
            origin= [x,y,z]
            normal = [1,0,0]
        elif cut_style.find("TE") != -1:
            Loc = float(re.findall(r"\d+", cut_style)[0])/100
            PIV = 1.25 + np.array(Loc)*0.3048*np.cos(AoA*np.pi/180)
            origin =  [PIV,0.,z_loc.z_mid_span]
            normal = [1.,0.,0.]
    elif data_type == 'CLIP':
        origin = [1.42222035,0.,z_loc.z_tip_gap]
        normal = [1.,0.,0.]
    print('The selected cut is of style: {0}'.format(data_type))
    print('The selected cut origin is: {0}'.format(origin))
    print('The selected cut normal is: {0}'.format(normal))
    return origin,normal

# Extracting the airfoil surface mesh for the base cut_locationof DMD
def extract_surface(mesh_fileDir:str,mesh_fileName:str,output:str, data_type:str, cut_location:str,\
    span:float=-0.2286,tip_gap:float=-0.1034,AoA:int=10,reload:bool=False):
    """ Extracting the surface mesh for the DMD computation
    Args:
        mesh_fileDir (str): The directory where the avbp mesh file is located (post zone merge)
        mesh_fileName (str): The name of the AVBP mesh file of the whole domain
        output(str): The destination directory where the extracted mesh is to be stored
        data_type (str): The type of data to be extracted, either 'FWH', 'OUTPUT', 'EXTRACT', or 'CLIP'
        cut_location (str): The location of the cut plane, defined either explicitly via PIV planes
            or with a %_TE location designating the distance from the trailing edge
        span (float): the span size
        tip_gap (float): the tip gap size
        AoA (int): the angle of attack
        reload (bool): The option to reload the files if they are already extracted
    Returns:
        dest_mesh (str): The destination directory where the extracted mesh is stored
        nodes (int): The number of nodes on the surface
    """
    text = 'Extracting Mesh'
    print(f'\n{text:.^80}\n')
    # The name of the mesh file is the cut location name
    mesh_name = cut_location + '_Mesh'  
    mesh = os.path.join(output,mesh_name+'.h5')
    if os.path.exists(mesh) == True and reload == False:
        print('----> LES Mesh already extracted at {0:s}'.format(mesh))
        # Loading the mesh
        r = Reader('hdf_antares')
        r['filename'] = mesh
        mesh = r.read()
        mesh.show()
        nodes = mesh[0][0]['x'].shape[0]
    elif data_type == 'FWH':
        ## Loading the mesh
        text = '----> Extracting the Surface Mesh from FWH Surface'
        print(f'{text}')  
        ## Loading the Main LES Mesh File
        r = Reader('hdf_avbp')
        r['filename'] = os.path.join(mesh_fileDir,mesh_fileName)
        base  = r.read() # b is the Base object of the Antares API
        airfoil_base = Family()
        airfoil_base['Airfoil_Surface'] = base.families['Patches']['Airfoil_Surface']
        airfoil_base['Airfoil_Trailing_Edge'] = base.families['Patches']['Airfoil_Trailing_Edge']
        airfoil_base['Airfoil_Side_LE'] = base.families['Patches']['Airfoil_Side_LE']
        airfoil_base['Airfoil_Side_Mid'] = base.families['Patches']['Airfoil_Side_Mid']
        airfoil_base['Airfoil_Side_TE']  = base.families['Patches']['Airfoil_Side_TE']
        base.families['SKIN'] = airfoil_base
        skin_base = base[base.families['SKIN']]
        text = '----> The Extracted Base Objects'
        print(f'\n{text}')   
        print(skin_base)
        ## Merging the extracted base objects to the same zone
        text = 'Merging the Base Objects'
        print(f'\n{text}')   
        myt = Treatment('merge')
        myt['base'] = skin_base
        myt['duplicates_detection'] = False
        myt['tolerance_decimals'] = 13
        # Writing the extraced mesh
        text = '----> Writing the Mesh File'
        print(f'\n{text}')   
        merged = myt.execute()
        writer = Writer('hdf_antares')
        writer['base'] = merged
        writer['filename'] = os.path.join(output,mesh_name)
        writer.dump()
        # The data for the original mesh
        text = '----> The Original Mesh Surface'
        print(f'\n{text}')   
        skin_base.show()
        # The data for the extracted and merged mesh
        text = '----> The Post Extraced Mesh Surface'
        print(f'\n{text}')   
        merged.show()
        nodes = merged['0000'][0].shape[0]
        print('\nThe number of nodes on the surface is: {0:d}'.format(nodes))
        mesh = os.path.join(output,mesh_name+'.h5')
        print('\nThe Extracted surface mesh is saved in {0:s}'.format(mesh))
    elif data_type == 'OUTPUT':
        # This is not suggested as requires extraction at run time, better to save the full soltuion and run Cutplanes post processing
        text = '----> Extracting the output mesh from AVBP OUTPUT POSTPROC'
        print(f'{text}')   
        r = Reader('hdf_avbp')
        mesh = os.path.join(mesh_fileDir,mesh_fileName)
        mesh = os.path.join(output,mesh_name+'.h5')
        r['filename'] = mesh
        base  = r.read() # b is the Base object of the Antares API
        shutil.copyfile(mesh, mesh)
        nodes = base['0000'][0].shape[0]
        text = '----> The AVBP OUTPUT POSTPROC Database Mesh Surface'
        base.show()
        print('\nThe number of nodes on the surface is: {0:d}'.format(nodes))
        print('\nThe copied mesh is saved in {0:s}'.format(mesh))
    elif data_type == 'EXTRACT' or data_type == 'CLIP':
        # Loading the mesh
        text = '----> Extracting the Cut Mesh'
        print(f'{text}')  
        # Loading the Main LES Mesh File
        r = Reader('hdf_avbp')
        r['filename'] = os.path.join(mesh_fileDir,mesh_fileName)
        base  = r.read() # b is the Base object of the Antares API
        # Compute the origin and normal of the cut plane location
        if data_type == 'CLIP':
            t= Treatment('clip')
            t['base'] = base
            cut_style = 'cylinder'
            t['type'] = cut_style
            t['axis'] = 'x'
            t['radius'] = 0.1
        elif data_type == 'EXTRACT':
            t= Treatment('cut')
            t['base'] = base
            t['type'] = 'plane'
            cut_style = cut_location
        origin, normal = map_cut(data_type,cut_style,tip_gap,span,AoA)
        t['origin'] = origin
        t['normal'] = normal
        inter = t.execute()     
        # Writing the extraced mesh
        text = '----> Writing the Mesh File'
        print(f'\n{text}')   
        writer = Writer('hdf_antares')
        writer['base'] = inter
        writer['filename'] = os.path.join(output,mesh_name)
        writer.dump()
        # The data for the extracted and merged mesh
        text = '----> The Post Extraced Mesh Surface'
        print(f'\n{text}')   
        inter.show()
        nodes = inter['0000'][0].shape[0]
        print('\nThe number of nodes on the surface is: {0:d}'.format(nodes))
        mesh = os.path.join(output,mesh_name+'.h5')
        print('\nThe Extracted surface mesh is saved in: {0:s}'.format(mesh))
    text = 'Mesh Extraction Complete!'
    print(f'\n{text:.^80}\n')  
    return mesh, nodes
