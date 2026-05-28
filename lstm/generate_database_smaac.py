from tools_database import generate_data_csv
from simuEF.tools_fea import read_props


props = read_props("simuEF/params_sma_init.txt")
generate_data_csv(props, 10)
