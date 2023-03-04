from typing import Any
import datetime
import dateutil.parser

from astropy.time import Time
import numpy as np

def infer_time_format(epoch: list[Any]):
    if isinstance(epoch[0], Time):
        format = Time
    elif isinstance(epoch[0], datetime.datetime):
        format = datetime.datetime
    elif isinstance(epoch[0], str):
        epoch = [dateutil.parser.parse(x, fuzzy_with_tokens=True)[0] for x in epoch]
        format = datetime.datetime
    elif np.isscalar(epoch[0]):
        if np.all(epoch > 100000):
            format = 'jd'
        else:
            format = 'mjd'
    return epoch, format