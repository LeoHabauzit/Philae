import sys
from pathlib import Path
import time

sys.path.append(str(Path(__file__).resolve().parent.parent))
from tools_homogeneisation import *
from tools_database import generate_data_csv, read_data
import torch

read_data("lstm/dataset/test_dataset.csv", i=1)
read_data("lstm/dataset/predictions.csv", i=1)
