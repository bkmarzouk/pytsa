import os
import pickle as pk

def log_stats(modelnumber, key, loc):
    
    # Build log path for pk file
    log_path = os.path.join(loc, "{}.stats".format(modelnumber))
    
    # Define possible keys for dict
    keys = ['short', 'kexit', 'feoi', 'back', "end"]
    
    # Assert that key is valid option, i.e. custom keys should not be passed!
    assert key in keys, "Key not vailid for log file: {}".format(key)
    
    # If path exists
    if os.path.exists(log_path):
        
        # Load existing binary file
        f = open(log_path, "rb")
        with f:
            log = pk.load(f)
    
        # Remove existing log file in order to update it
        os.remove(log_path)
    
    # Otherwise build new log file
    else:
        log = {'short': 0, 'kexit': 0, 'feoi': 0, 'back': 0}
        
    # If key corresponds to an iterable update
    if key in ['short', 'kexit', 'feoi', 'back']:
        log[key] += 1
    
    # If key corresponds to a sample being obtained
    elif key == "end":
        # Add key values
        keys.remove("end")
        tot = sum([log[k] for k in keys])
        log["end"] = tot + 1 # Final number of iterations required to obtain successful sample
        
    else: pass
    
    # Write new dictionary
    f = open(log_path, "wb")
    with f:
        pk.dump(log_path, f)