"""
Module to create dicts for multiple (or single) mooring extractions.
"""
from lo_user_tools import llxyfun as lxf

def get_sta_dict(job_name):

    # specific job definitions

    if job_name == 'willapa_bc': # Willapa Bay Center PCSGA Mooring
        sta_dict = {
        'wbc': (-123.9516, 46.6290)
        }
        
    elif job_name == 'mickett_1':
        sta_dict = {
        'ORCA_Hansville': (-122.6270, 47.9073),
        'ORCA_Hoodsport': (-123.1126, 47.4218),
        'ORCA_Point_Wells': (-122.3972, 47.7612),
        'Central_Main_Stem_Hood_Canal': (-122.989507, 47.574352),
        'North_Central_Main_Basin': (-122.440755, 47.825099)
        }
        
    elif job_name == 'mickett_2':
        sta_dict = {
        'Carr_Inlet_ORCA': (-122 - 43.8/60, 47 + 16.8/60),
        'East_of_Fox_Island': (-122 - 35.158/60, 47 + 13.185/60)
        }
        
    elif job_name == 'stoll_corals':
        sta_dict = {
        'Carson_D01_Lopez': (-122.8728, 48.36816),
        'Carson_D02_Admiralty': (-122.7883, 48.19252),
        'Carson_D04_Admiralty': (-122.8166, 48.19764),
        'Carson_D05_Keystone': (-122.6576, 48.12828),
        'Carson_D07_NorthAdmiralty': (-122.8898, 48.22245),
        'Carson_D08_Canada': (-123.149, 48.36136),
        'USNM_19270_Canada': (-123.233, 48.35),
        'USNM_92626_Admiralty': (-122.80, 48.1917),
        'USNM_19228_Dungeness': (-123.189, 48.225),
        'USNM_19272_Admiralty': (-122.817, 48.20),
        'USNM_92620_Lopez': (-122.85, 48.3667)
        }
        
    elif job_name == 'stoll_obs':
        sta_dict = {
        'DOE_SJF002': (-123.025, 48.25),
        'DOE_ADM002': ( -122.8417151, 48.1875056),
        'DOE_ADM001': ( -122.616715, 48.0300056),
        'WOAC_STN21': (-122.8504, 48.1883),
        'WOAC_STN20': (-122.6848, 48.142),
        'WOAC_STN19': (-122.6318, 48.0915)
        }
            
    elif job_name == 'Kelly':
        # note I pushed two of the locations a bit West to get off the landmask
        sta_dict = {
        'Seal_Rock': (-122.87004, 47.70557),
        'Little_Dewatto': (-123.08612-.005, 47.44489),
        'Red_Bluff': (-123.10438-.007, 47.41625)
        }
        
    elif job_name == 'jazzy':
        sta_dict = {
        'Middle_Bank': (-123.09651, 48.40935),
        'East_Bank': (-122.97376, 48.30042),
        'Upright_Channel': (-122.923005, 48.55410),
        'Blakely_Orcas': (-122.82880, 48.58790),
        'Rosario_Strait': (-122.74001, 48.64631),
        'North_Station': (-123.04166, 48.58330),
        'South_Station': (-122.94330, 48.42000),
        'Hein_Bank': (-123.03940, 48.35825)
        }
        
    elif job_name == 'ooi':
        sta_dict = {
            'CE01':(-124.095, 44.6598), # Oregon Inshore (25 m)
            'CE02':(-124.304, 44.6393), # Oregon Shelf (80 m)
            'CE04':(-124.956, 44.3811), # Oregon Offshore (588 m)
            'PN01A':(-125.3983, 44.5096) # Slope Base (2905 m)
        }
        
    elif job_name == 'erika_esci491w2022':
        sta_dict = {
        'Olympia': (-122.9165, 47.0823),
        'Tacoma': (-122.4758, 47.3015),
        'Seattle_West_Point': (-122.4435, 47.6813),
        'Bellingham': (-122.5519, 48.7348),
        'Central_Hood_Canal': (-122.9895, 47.5744),
        'Skokomish': (-123.1272, 47.3639),
        'Hein_Bank': (-123.0394, 48.35825),
        'Admiralty': (-122.6949, 48.1370),
        'Everett': (-122.2806, 47.9855)
        }
        
    elif job_name == 'scoot':
        sta_dict= {
        'F_090_WEH': (-126.92523957836097, 50.86003686529567),
        'F_087_SIM': (-126.88867577473951, 50.8779837900623),
        'F_083_CEC': (-126.71014703214891, 50.86003686529567),
        'M_083_087': (-126.79852585292794, 50.89607323113631),
        'F_084_CYP': (-126.65795536471262, 50.84223133403236),
        'M_083_084': (-126.69268072600828, 50.86003686529567),
        'F_088_SIR': (-126.60638135698926, 50.84223133403236),
        'M_084_089': (-126.6235045170437, 50.86003686529567),
        'F_082_BRE': (-125.34911668789157, 50.28660059365715),
        'F_089_VEN': (-125.33704071737607, 50.300007509695305),
        'F_086_RAZ': (-125.01765650844104, 50.327117635547935),
        'E_079_TOF': (-126.04793822349058, 49.21393673803513),
        'E_129_STG': (-125.12767093555838, 50.0962109697609),
        'E_004_NOO': (-126.60638135698926, 49.57825722433716),
        'E_002_ESP': (-126.9806322902372, 49.829746612851736),
        'E_016_RIV': (-126.13824776832129, 49.66598688018571),
        'F_004_CON': (-126.47181502357597, 49.65693905574285),
        'F_022_GOR': (-126.43883566700802, 49.64795343794384),
        'F_016_MUC': (-126.3414536354723, 49.63902959909838),
        'M_023_014': (-126.3414536354723, 50.607192072687155),
        'F_008_BAR': (-125.2773744875855, 50.31351066854337),
        'F_005_AHL': (-124.15671501359827, 49.77957267482829),
        'F_013_VAN': (-123.85335983884296, 49.675097341923546),
        'F_011_SAL': (-123.83316138272168, 49.62136556218934),
        'F_006_NEW': (-123.65810809633727, 49.64795343794384),
        'E_006_RIV': (-123.5436501783167, 49.693507914816124)
        }
        
    elif job_name == 'orca_eb':
        sta_dict = {
        'CI': (-122.7300, 47.2800),
        'PW': (-122.3972, 47.7612),
        'NB': (-122.6270, 47.9073),
        'DB': (-122.8029, 47.8034),
        'HP': (-123.1126, 47.4218),
        'TW': (-123.0083, 47.3750)
        }

    elif job_name == 'orca_adcp_eb':
        sta_dict = {
        'NB': (-122.6270, 47.9073)
        }

    elif job_name == 'sill0_spinup_test':
        sta_dict = {
        'basin1': (1.77859, 45),
        'basin2': (1.27042, 45),
        'basin3': (0.76225, 45)
        }

    elif job_name == 'basin_sill':
        sta_dict = {
        'basin1': (1.7815475813532617, 45),
        'sill1L': (1.5524914637506995, 45),
        'sill1M': (1.5270407840170814, 45),
        'sill1S': (1.5015901042834634, 45),
        'basin2': (1.2725339866809011, 45),
        'sill2L': (1.043477869078339, 45),
        'sill2M': (1.018027189344721, 45),
        'sill2S': (0.9925765096111029, 45),
        'basin3': (0.7635203920085407, 45),
        'sill3L': (0.5344642744059785, 45),
        'sill3M': (0.5090135946723605, 45),
        'sill3S': (0.4835629149387424, 45),
        'basin4': (0.25450679733618026, 45)
        }

    elif job_name == 'sill1_test':
        sta_dict = {
        'innerbasin': (0.8144217514757768, 45), #64km
        'innersill': (0.5599149541395965, 45), #44km
        'midsill': (0.5344642744059785, 45), #42km
        'outersill': (0.5090135946723605, 45), #40km
        'outerbasin': (0.25450679733618026, 45) #20km
        }

    elif job_name == 'sill2_center':
        sta_dict = {
        'outer1': (0.12725339866809013, 45), #10km
        'outer2': (0.25450679733618026, 45), #20km
        'outer3': (0.38176019600427036, 45), #30km
        'sill1': (0.5090135946723605, 45), #40km
        'sill2': (0.5344642744059785, 45), #42km
        'sill3': (0.5599149541395965, 45), #44km
        'sill4': (0.5853656338732145, 45), #46km
        'sill5': (0.6108163136068325, 45), #48km
        'inner1': (0.7380697122749227, 45), #58km
        'inner2': (0.8653231109430128, 45), #68km
        'inner3': (0.9925765096111029, 45) #78km
        }

    elif job_name == 'sill2_north':
        sta_dict = {
        'outer1': (0.12725339866809013, 45.017996348225445), #10km
        'outer2': (0.25450679733618026, 45.017996348225445), #20km
        'outer3': (0.38176019600427036, 45.017996348225445), #30km
        'sill1': (0.5090135946723605, 45.017996348225445), #40km
        'sill2': (0.5344642744059785, 45.017996348225445), #42km
        'sill3': (0.5599149541395965, 45.017996348225445), #44km
        'sill4': (0.5853656338732145, 45.017996348225445), #46km
        'sill5': (0.6108163136068325, 45.017996348225445), #48km
        'inner1': (0.7380697122749227, 45.017996348225445), #58km
        'inner2': (0.8653231109430128, 45.017996348225445), #68km
        'inner3': (0.9925765096111029, 45.017996348225445) #78km
        }

    elif job_name == 'sill2_south':
        sta_dict = {
        'outer1': (0.12725339866809013, 44.982003651774555), #10km
        'outer2': (0.25450679733618026, 44.982003651774555), #20km
        'outer3': (0.38176019600427036, 44.982003651774555), #30km
        'sill1': (0.5090135946723605, 44.982003651774555), #40km
        'sill2': (0.5344642744059785, 44.982003651774555), #42km
        'sill3': (0.5599149541395965, 44.982003651774555), #44km
        'sill4': (0.5853656338732145, 44.982003651774555), #46km
        'sill5': (0.6108163136068325, 44.982003651774555), #48km
        'inner1': (0.7380697122749227, 44.982003651774555), #58km
        'inner2': (0.8653231109430128, 44.982003651774555), #68km
        'inner3': (0.9925765096111029, 44.982003651774555) #78km
        }

    elif job_name == 'sill3_center':
        sta_dict = {
        'outer1': (lxf.x2lon(10e3,0,45), 45), #10km
        'outer2': (lxf.x2lon(20e3,0,45), 45), #20km
        'outer3': (lxf.x2lon(30e3,0,45), 45), #30km
        'sill1': (lxf.x2lon(40e3,0,45), 45), #40km
        'sill2': (lxf.x2lon(42e3,0,45), 45), #42km
        'sill3': (lxf.x2lon(44e3,0,45), 45), #44km
        'sill4': (lxf.x2lon(46e3,0,45), 45), #46km
        'sill5': (lxf.x2lon(48e3,0,45), 45), #48km
        'inner1': (lxf.x2lon(58e3,0,45), 45), #58km
        'inner2': (lxf.x2lon(68e3,0,45), 45), #68km
        'inner3': (lxf.x2lon(78e3,0,45), 45) #78km
        }

    elif job_name == 'sill3_north':
        sta_dict = {
        'outer1': (lxf.x2lon(10e3,0,45), lxf.y2lat(2e3,45)), #10km
        'outer2': (lxf.x2lon(20e3,0,45), lxf.y2lat(2e3,45)), #20km
        'outer3': (lxf.x2lon(30e3,0,45), lxf.y2lat(2e3,45)), #30km
        'sill1': (lxf.x2lon(40e3,0,45), lxf.y2lat(2e3,45)), #40km
        'sill2': (lxf.x2lon(42e3,0,45), lxf.y2lat(2e3,45)), #42km
        'sill3': (lxf.x2lon(44e3,0,45), lxf.y2lat(2e3,45)), #44km
        'sill4': (lxf.x2lon(46e3,0,45), lxf.y2lat(2e3,45)), #46km
        'sill5': (lxf.x2lon(48e3,0,45), lxf.y2lat(2e3,45)), #48km
        'inner1': (lxf.x2lon(58e3,0,45), lxf.y2lat(2e3,45)), #58km
        'inner2': (lxf.x2lon(68e3,0,45), lxf.y2lat(2e3,45)), #68km
        'inner3': (lxf.x2lon(78e3,0,45), lxf.y2lat(2e3,45)) #78km
        }

    elif job_name == 'sill3_south':
        sta_dict = {
        'outer1': (lxf.x2lon(10e3,0,45), lxf.y2lat(-2e3,45)), #10km
        'outer2': (lxf.x2lon(20e3,0,45), lxf.y2lat(-2e3,45)), #20km
        'outer3': (lxf.x2lon(30e3,0,45), lxf.y2lat(-2e3,45)), #30km
        'sill1': (lxf.x2lon(40e3,0,45), lxf.y2lat(-2e3,45)), #40km
        'sill2': (lxf.x2lon(42e3,0,45), lxf.y2lat(-2e3,45)), #42km
        'sill3': (lxf.x2lon(44e3,0,45), lxf.y2lat(-2e3,45)), #44km
        'sill4': (lxf.x2lon(46e3,0,45), lxf.y2lat(-2e3,45)), #46km
        'sill5': (lxf.x2lon(48e3,0,45), lxf.y2lat(-2e3,45)), #48km
        'inner1': (lxf.x2lon(58e3,0,45), lxf.y2lat(-2e3,45)), #58km
        'inner2': (lxf.x2lon(68e3,0,45), lxf.y2lat(-2e3,45)), #68km
        'inner3': (lxf.x2lon(78e3,0,45), lxf.y2lat(-2e3,45)) #78km
        }

    elif job_name == 'sill5km_center':
        sta_dict = {
        # 'outer1': (lxf.x2lon(10e3,0,45), 45), #10km
        # 'outer2': (lxf.x2lon(20e3,0,45), 45), #20km
        # 'outer3': (lxf.x2lon(30e3,0,45), 45), #30km
        # 'sill1': (lxf.x2lon(40e3,0,45), 45), #40km
        # 'sill2': (lxf.x2lon(41.25e3,0,45), 45), #41.25km
        # 'sill3': (lxf.x2lon(42.5e3,0,45), 45), #42.5km
        # 'sill4': (lxf.x2lon(43.75e3,0,45), 45), #43.75km
        # 'sill5': (lxf.x2lon(45e3,0,45), 45), #45km
        # 'inner1': (lxf.x2lon(55e3,0,45), 45), #55km
        # 'inner2': (lxf.x2lon(65e3,0,45), 45), #65km
        # 'inner3': (lxf.x2lon(75e3,0,45), 45) #75km
        'mouth': (lxf.x2lon(0,0,45), 45), #0km (mouth)
        'outer': (lxf.x2lon(20e3,0,45), 45), #20km (middle of outer basin)
        'sill': (lxf.x2lon(42.5e3,0,45), 45), #42.5km (middle of sill)
        'inner': (lxf.x2lon(65e3,0,45), 45) #65km (middle of inner basin)
        }

    elif job_name == 'sill20kmdeep_center':
        sta_dict = {
        # 'outer1': (lxf.x2lon(10e3,0,45), 45), #10km
        # 'outer2': (lxf.x2lon(20e3,0,45), 45), #20km
        # 'outer3': (lxf.x2lon(30e3,0,45), 45), #30km
        # 'sill1': (lxf.x2lon(40e3,0,45), 45), #40km
        # 'sill2': (lxf.x2lon(45e3,0,45), 45), #45km
        # 'sill3': (lxf.x2lon(50e3,0,45), 45), #50km
        # 'sill4': (lxf.x2lon(55e3,0,45), 45), #55km
        # 'sill5': (lxf.x2lon(60e3,0,45), 45), #60km
        # 'inner1': (lxf.x2lon(70e3,0,45), 45), #70km
        # 'inner2': (lxf.x2lon(80e3,0,45), 45), #80km
        # 'inner3': (lxf.x2lon(90e3,0,45), 45) #90km
        'mouth': (lxf.x2lon(0,0,45), 45), #0km (mouth)
        'outer': (lxf.x2lon(20e3,0,45), 45), #20km (middle of outer basin)
        'sill': (lxf.x2lon(50e3,0,45), 45), #50km (middle of sill)
        'inner': (lxf.x2lon(80e3,0,45), 45) #80km (middle of inner basin)
        }

    elif job_name == 'sill80km_center':
        sta_dict = {
        # 'outer1': (lxf.x2lon(10e3,0,45), 45), #10km
        # 'outer2': (lxf.x2lon(20e3,0,45), 45), #20km
        # 'outer3': (lxf.x2lon(30e3,0,45), 45), #30km
        # 'sill1': (lxf.x2lon(40e3,0,45), 45), #40km
        # 'sill2': (lxf.x2lon(60e3,0,45), 45), #60km
        # 'sill3': (lxf.x2lon(80e3,0,45), 45), #80km
        # 'sill4': (lxf.x2lon(100e3,0,45), 45), #100km
        # 'sill5': (lxf.x2lon(120e3,0,45), 45), #120km
        # 'inner1': (lxf.x2lon(130e3,0,45), 45), #130km
        # 'inner2': (lxf.x2lon(140e3,0,45), 45), #140km
        # 'inner3': (lxf.x2lon(150e3,0,45), 45) #150km
        'mouth': (lxf.x2lon(0,0,45), 45), #0km (mouth)
        'outer': (lxf.x2lon(20e3,0,45), 45), #20km (middle of outer basin)
        'sill': (lxf.x2lon(80e3,0,45), 45), #80km (middle of sill)
        'inner': (lxf.x2lon(1400e3,0,45), 45) #1400km (middle of inner basin)
        }

    else:
        print('Unsupported job name!')
        a = dict()
        return a
        
    return sta_dict