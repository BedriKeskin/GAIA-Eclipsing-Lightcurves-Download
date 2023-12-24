import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from astropy import units as u
from astropy.io.votable import parse_single_table
from astropy.time import Time
from astropy.timeseries import TimeSeries
from astropy.timeseries import aggregate_downsample
from astroquery.gaia import Gaia
import datetime
import os

now1 = datetime.datetime.now()
print("Start time: ", now1)

if not os.path.exists("PNG"):
    os.makedirs("PNG")
if not os.path.exists("VOTable"):
    os.makedirs("VOTable")

job = Gaia.launch_job(
    "SELECT TOP 100 source_id, global_ranking, reference_time, frequency FROM gaiadr3.vari_eclipsing_binary")
results = job.get_results()

for index, result in enumerate(results):
    print(index, result['source_id'])
    url = 'https://gea.esac.esa.int/data-server/data?ID=Gaia+DR3+' + str(
        result['source_id']) + '&RETRIEVAL_TYPE=EPOCH_PHOTOMETRY&VALID_DATA=true'
    byte = requests.get(url).content

    with open('VOTable/' + str(result['source_id']) + '.xml', 'wb') as file:
        file.write(byte)
        file.close()

    table = parse_single_table('VOTable/' + str(result['source_id']) + '.xml').to_table()
    LC = pd.DataFrame.from_dict(table, orient='index').transpose()
    LC = LC.astype({'mag': float})
    LC = LC.astype({'flux': float})
    LC = LC.astype({'flux_error': float})
    LC = LC.astype({'flux_over_error': float})
    LC = LC.astype({'rejected_by_variability': bool})
    LC = LC.astype({'rejected_by_photometry': bool})
    LC = LC[LC['time'] != '--']  # time boş satırları gözardı et
    LC = LC[(LC['rejected_by_variability'] == False) | (
            LC['rejected_by_variability'] == False)]  # rejected olanları gözardı et
    LC = LC[(LC['band'] == 'G')]  # sadece G bandı al

    LC['time'] = LC['time'] + 2455197.5

    LC = LC.drop(['source_id'], axis=1)
    LC = LC.drop(['transit_id'], axis=1)
    LC = LC.drop(['solution_id'], axis=1)
    LC = LC.drop(['other_flags'], axis=1)
    LC = LC.drop(['rejected_by_variability'], axis=1)
    LC = LC.drop(['rejected_by_photometry'], axis=1)
    LC = LC.drop(['band'], axis=1)

    LC.index = pd.to_datetime(LC['time'], origin='julian', unit='D')
    LC.rename(columns={'time': 'bjd'}, inplace=True)

    timeSeries = TimeSeries.from_pandas(LC)

    T0 = result['reference_time'] + 2455197.5
    T0 = pd.to_datetime(T0, origin='julian', unit='D')
    T0 = Time(T0, format='datetime')
    P = 1 / result['frequency']

    ts_folded = timeSeries.fold(period=P * u.day, epoch_time=T0)
    ts_binned = aggregate_downsample(ts_folded, time_bin_size=10 * u.min, aggregate_func=np.nanmedian)

    fig, ax = plt.subplots(1, 1)
    #plt.axis('off')
    ax.plot(ts_binned.time_bin_start.jd, ts_binned['flux'], 'k.')
    plt.savefig('PNG/' + str(result['source_id']) + '.png')
    plt.close()

now2 = datetime.datetime.now()
print("End time: ", now2)
print("Elapsed time: ", now2 - now1)
