from utils import infer_time_format
from typing import Any, Optional
import numpy as np

class Time:

    def __init__(self, epoch: list[Any], format: Optional[str] = None, scale: Optional[str] = 'utc'):

        epoch = np.atleast_1d(epoch)

        if format is None:
            self.epoch, self.format = infer_time_format(epoch)
        else:
            self.epoch = epoch
            self.format = format

        self.scale = scale

    def __repr__(self):
        return f'Time(epoch={self.epoch}, format={self.format}, scale={self.scale})'
    
def main():
    epoch = ['28 February 2023']
    time = Time(epoch)
    print(time)

if __name__ == '__main__':
    main()