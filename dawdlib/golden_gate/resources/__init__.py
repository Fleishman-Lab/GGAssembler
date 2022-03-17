import glob
import os

ligation_data = {
    os.path.basename(csv): os.path.abspath(csv)
    for csv in glob.glob(os.path.join(os.path.dirname(__file__), "*.csv"))
}
