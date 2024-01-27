"""
flick wrapper

example:
  import flick
  m = flick.run("radiator planck 6000 9")
  print(m)
"""
import numpy as np
import subprocess

def run(arguments):
    flick_output = _run_os("flick "+arguments)
    if flick_output.stdout:
        return _to_matrix(flick_output)

def _run_os(command):
    c = subprocess.run(command, stdout=subprocess.PIPE, shell=True, check=True)
    return c

def _to_matrix(flick_output):
    try:
        l = flick_output.stdout.decode('utf-8');
        l = l.splitlines()
        l = list(filter(None, l))
        if not l:
            return np.ndarray(shape=(0,0),dtype=float) 
        numbers_first_line = [float(x) for x in l[0].split()]
        rows = len(l);
        cols = len(numbers_first_line)
        matrix = np.ndarray(shape=(rows,cols),dtype=float)
        for i in range(rows):
            floats = [float(x) for x in l[i].split()]
            matrix[i] = floats
        return matrix
    except Exception as e:
        raise Exception(flick_output.stdout.decode('utf-8')) 

def _to_string(numbers):
    if isinstance(numbers, np.ndarray):
        return "\""+" ".join(str(x) for x in numbers)+"\""
    return str(numbers)

def config(file_name, parameter_name, value):
    command = "flick text "+file_name+" set "+parameter_name+" "+ \
        _to_string(value)+" > "+file_name+"_tmp"
    _run_os(command)
    _run_os("mv -f "+file_name+"_tmp "+ file_name)

def to_streams(n_angles):
    n_streams = np.floor(n_angles**(1/1.6))    
    if (n_streams % 2) != 0:
        n_streams += 1
    return str(n_streams).rstrip('0').rstrip('.');





