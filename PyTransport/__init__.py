from PyTransport import PyTransSetup
import os

PyTransSetup.pathSet()

sym_cache = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sym_cache")  # cache locations for sym calcs
fmet_cache = os.path.join(sym_cache, "fmet")
pot_cache = os.path.join(sym_cache, "pot")

for d in [sym_cache, fmet_cache, pot_cache]:
    if not os.path.exists(d):
        os.mkdir(d)

os.environ["model_cache"] = sym_cache
os.environ["fmet_cache"] = fmet_cache
os.environ["pot_cache"] = pot_cache