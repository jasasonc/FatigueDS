
import numpy as np
import scipy


def convert_Q_damp(self,Q=None,damp=None): 

    if damp is not None:
        self.damp = damp
        self.Q = 1/(2*self.damp)

    elif Q is not None:
        self.Q = Q
        self.damp = 1/(2*self.Q)

def get_f0_range(self,f0_data):
        """
        Natural frequency range [Hz] -> X-axis of MRS/FDS plot.
        """
        if isinstance(f0_data, (tuple)) and len(f0_data)==3:
            f0_start, f0_stop, f0_step = f0_data
            f0_range = np.arange(f0_start, f0_stop + f0_step, f0_step, dtype=float)
        else:
             f0_range = f0_data        
        
        if f0_range[0]==0:
            f0_range[0] = 1e-3    # sets frequency to a small number to avoid dividing by 0
        else:
            pass

        return f0_range

def update_docstring(target_method, doc_method=None, delimiter='---', added_doc=''):
    """
    Update the docstring in target_method with the docstring from doc_method.
    
    :param target_method: The method that waits for the docstring
    :type target_method: method
    :param doc_method: The method that holds the desired docstring
    :type doc_method: method
    :param delimiter: insert the desired docstring between two delimiters, defaults to '---'
    :type delimiter: str, optional
    """
    
    docstring = target_method.__doc__.split(delimiter)
    
    leading_spaces = len(docstring[1].replace('\n', '')) - len(docstring[1].replace('\n', '').lstrip(' '))
    
    if doc_method is not None:
        if doc_method.__doc__:
            docstring[1] = doc_method.__doc__
            
        else:
            docstring[1] = '\n' + ' '*leading_spaces + \
                'The selected method does not have a docstring.\n'
    else:
        docstring[1] = added_doc.replace('\n', '\n' + ' '*leading_spaces)
    
    target_method.__func__.__doc__ = delimiter.join(docstring)








