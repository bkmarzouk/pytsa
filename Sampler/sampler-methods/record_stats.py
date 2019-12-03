import os
import pickle as pk
import gc

def log_stats(modelnumber, flag, loc, ttotal=None):
    """ Integration flags:

    int SHORT = -50;          // Inflation too short
    int KEXIT = -49;          // Unable to find Fourier mode
    int FEOI = -48;           // Integration failure in feoi
    int BACK = -47;           // Integration failure in background
    int VIOLATED = -46;       // Field space position violates model
    int ETERNAL = -45;        // Unable to find end of inflation
    int TIMEOUT = -44;        // Integration time exceeded

    """
    
    # Build log path for pk file
    log_path = os.path.join(loc, "{}.stats".format(modelnumber))
    
    # Define integrator flags and dict keys: Note that "end" does not have a flag number, so acts as key if passed
    flags = [-50, -49, -48, -47, -46, -45, -44]
    keys = ['short', 'kexit', 'feoi', 'back', "violated", "eternal", "timeout", "end"]
    if flag == "end":
        key = flag
    else:
        key = keys[flags.index(flag)] # Get dict key from flag
    
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
        log = {'short': 0, 'kexit': 0, 'feoi': 0, 'back': 0, "violated": 0, "eternal": 0, "timeout": 0}
        
    # If key corresponds to an iterable update
    if key in ['short', 'kexit', 'feoi', 'back', "violated", "eternal", "timeout"]:
        log[key] += 1
    
    # If key corresponds to a sample being obtained
    elif key == "end":
        
        assert ttotal is not None, "Final time required"
        
        # Add key values
        keys.remove("end")
        tot = sum([log[k] for k in keys])
        log["end"] = tot + 1 # Final number of iterations required to obtain successful sample
        log["time"] = ttotal

    else: pass
    
    # Write new dictionary
    f = open(log_path, "wb")
    with f:
        pk.dump(log, f)

    # Clean up final bits
    del flags; del log; del log_path; del f; del keys; del modelnumber; del ttotal; del loc; del flag; del key
    if key == "end": del tot
    assert len([item for item in locals() if not item.startswith("__")]), "Unexpected local variable: {}".format(item)
    gc.collect()