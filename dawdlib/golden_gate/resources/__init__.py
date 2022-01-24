import glob
import os

ligation_data = {
    csv: os.path.abspath(csv)
    for csv in glob.glob("*.csv", root_dir=os.path.dirname(__file__))
}


__all__ = {ligation_data}
