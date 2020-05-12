from sim_funcs import *

import numpy as np

# Generate fake objects
@jit
def make_fakes(N_obj, a_power, a_min, a_max, H_power, H_min, H_max,
               inc_width, Brown=True, flat_a=False):
    
    
    if flat_a == True:
        a = np.random.uniform(a_min, a_max, N_obj)
    else:
        a = rand_power_law(a_min, a_max, a_power, N_obj)
    
    H = rand_power_law(H_min, H_max, H_power, N_obj)
    #varpi = np.random.uniform(-np.pi, np.pi, N_obj)
    #node = np.random.uniform(-np.pi, np.pi, N_obj)
    varpi = np.random.normal(np.radians(71), np.radians(60), N_obj)
    node = np.random.normal(np.radians(113), np.radians(60), N_obj)
    e = np.random.uniform(0.69, 0.999, N_obj)
    q = a * (1 - e)
    ma = np.radians(np.random.exponential(20/9, N_obj) *
                np.random.choice([-1,1], N_obj))

    if Brown == True:
        inc = np.random.rayleigh(inc_width, N_obj)
    else:
        inc = np.random.uniform(0, np.pi/3, N_obj)
    
    return a, e, q, varpi, node, ma, inc, H

# Draw random values from a power law distribution
@jit
def rand_power_law(x_min, x_max, x_power, N_obj):
    distribution = ((x_max**(x_power+1.0) - x_min**(x_power+1.0))
                    *np.random.uniform(0, 1, N_obj) +
                    x_min**(x_power+1.0))**(1.0/(x_power+1.0))
    return distribution


import sys

run_description = '_P9_6' # Manually set a unique string for each run.
DES_out_name = 'DES{}'.format(str(run_description))
OSSOS_out_name = 'OSSOS{}'.format(str(run_description))
ST_out_name = 'ST{}'.format(str(run_description))

'''
All angle objects are in RADIANS!!!
'''

# Set the number of simulated detections in each survey
N_obj = 10000

# Set the number of fakes to produce per loop
N_per_loop = 10000000

# Specify the semimajor axis distribution (single-slope)
flat = False
a_power = 0.7
a_min = 250
a_max = 1250

# Specify the absolute magnitude distribution (single-slope)
H_power = 0.8
H_min = 3.5
H_max = 10.0

# Specify the center and width of the Brown Distribution
Brown = True
inc_center = 0
inc_width = np.pi/12 # 15 degrees

# Initialize dataframes for storing 'detected' objects
DES = pd.DataFrame()
OSSOS = pd.DataFrame()
ST = pd.DataFrame()

# Set initial number of detections in each survey to zero
N_in_DES, N_in_OSSOS, N_in_ST = 0, 0, 0
kk = 0

while min(N_in_DES, N_in_OSSOS, N_in_ST) < N_obj:
    kk += 1
    
    scipy.random.seed() # Set a different random seed in each loop to avoid redundancy
    a, e, peri, varpi, node, ma, inc, H = make_fakes(N_per_loop, a_power, 
                                                     a_min, a_max, H_power,
                                                     H_min, H_max, inc_width,
                                                     Brown=Brown, flat_a=flat)

    E = calcE(e, ma, 100)
    ra, dec, delta, r, mag, mag2 = calcradec(e, E, a, varpi-node, node, inc, H)
    hpix64 = radec_to_index(ra, dec, NSIDE=64)
    hpix128 = radec_to_index(ra, dec, NSIDE=128)
    hpix1024 = radec_to_index(ra, dec, NSIDE=1024)

    
    # Initiate a dictionary to pass to a pandas dataframe
    data = {'a':a,
            'e':e,
            'inc':inc,
            'node':node,
            'varpi':varpi,
            'H':H,
            'epoch':epoch,
            'mag':mag,
            'mean_anomaly':ma,
            'peri':a * (1 - e),
            'sun_distance':r,
            'earth_distance':delta,
            'hpix64':hpix64,
            'hpix128':hpix128,
            'hpix1024':hpix1024}
    df = pd.DataFrame(data)

    if N_in_DES < N_obj:
        
        DES_wide = check_if_in_survey(df, DES_hpix, 23.5, 64)
        DES_SN_SHALLOW = check_if_in_survey(df, DES_SN_shallow, 23.5, 1024)
        DES_SN_DEEP = check_if_in_survey(df, DES_SN_deep, 24.4, 1024)
        inDES = pd.concat([DES_wide, DES_SN_SHALLOW, DES_SN_DEEP])
        DES = DES.append(inDES, ignore_index=True)
        N_in_DES += len(inDES)
        

    if N_in_OSSOS < N_obj:
        
        in15BS = check_if_in_survey(df, B15BS, 25.12, 64)
        in15BT = check_if_in_survey(df, B15BT, 24.93, 64)
        in13BL = check_if_in_survey(df, B13BL, 24.42, 64)
        in14BH = check_if_in_survey(df, B14BH, 24.67, 64)
        in15BC = check_if_in_survey(df, B15BC, 24.78, 64)
        in15BC_2 = check_if_in_survey(df, B15BC_2, 24.78, 64)
        in15BD = check_if_in_survey(df, B15BD, 25.15, 64)
        in15AP = check_if_in_survey(df, B15AP, 24.80, 64)        
        in13AE = check_if_in_survey(df, B13AE, 24.09, 64)
        in15AM = check_if_in_survey(df, B15AM, 24.87, 64)
        in13AO = check_if_in_survey(df, B13AO, 24.40, 64)
        inOSSOS = pd.concat([in15BS, in15BT, in13BL, in14BH, in15BC, 
                             in15BC_2, in15BD, in15AP, in13AE, in15AM,
                             in13AO]).drop_duplicates()
        OSSOS = OSSOS.append(inOSSOS, ignore_index=True)
        N_in_OSSOS += len(inOSSOS)
        
        
    if N_in_ST < N_obj:
        
        inctio = check_if_in_survey(df, ctio_st, 24.5, 64, ST=True)
        inkpno = check_if_in_survey(df, kpno_st, 24.6, 128, ST=True)
        inmagellan = check_if_in_survey(df, magellan_st, 25.5, 128, ST=True)
        insubaru = check_if_in_survey(df, subaru_st, 25.7, 128, ST=True)
        inlbt = check_if_in_survey(df, lbt_st, 25.2, 128, ST=True)
        inST = pd.concat([inctio, inkpno, inlbt,
                          inmagellan, insubaru]).drop_duplicates()
        
        ST = ST.append(inST, ignore_index=True)
        N_in_ST += len(inST)

    sys.stdout.write('{}% complete after {} loops \r'.format(math.floor(100 * min(N_in_ST, N_in_DES,
                                                  N_in_OSSOS)/N_obj), kk))
        
outdir = '~/Research/TNOs/ClusteringOfETNOs/Data/sim_data/'
DES[:N_obj].to_csv(os.path.join(outdir, '{}.csv'.format(DES_out_name)))
OSSOS[:N_obj].to_csv(os.path.join(outdir, '{}.csv'.format(OSSOS_out_name)))
ST[:N_obj].to_csv(os.path.join(outdir, '{}.csv'.format(ST_out_name)))
