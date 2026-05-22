#!/usr/bin/env python3
import sys, os
sys.path.insert(0, '/home/dministrator/Lab/clones/molcraft-agent/src')
from config import RECEPTOR_PDBQT, DOCKING_CENTER, DOCKING_SIZE
print(f'receptor: {RECEPTOR_PDBQT}')
print(f'center: {DOCKING_CENTER}')
print(f'size: {DOCKING_SIZE}')
print(f'receptor exists: {os.path.exists(RECEPTOR_PDBQT)}')