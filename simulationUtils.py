import sys,os,errno

def make_sure_path_exists(path):
    """ 
    Funciton: Create the output folder for simulations.
      If the folder already exists then do nothing.
      If some other error occurs then raise that error.
    
    Input: (str) path name
    Output: None
    """
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

