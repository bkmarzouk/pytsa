import os
import site
import sys
import platform

sym_cache = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sym_cache")  # cache locations for sym calcs
fmet_cache = os.path.join(sym_cache, "fmet")
pot_cache = os.path.join(sym_cache, "pot")
covd_cache = os.path.join(sym_cache, "covd")

for d in [sym_cache, fmet_cache, pot_cache, covd_cache]:
    if not os.path.exists(d):
        os.mkdir(d)

os.environ["model_cache"] = sym_cache
os.environ["fmet_cache"] = fmet_cache
os.environ["pot_cache"] = pot_cache
os.environ["covd_cache"] = covd_cache

root = os.path.dirname(__file__)
site.addsitedir(root)
install_lib = os.path.join(root, "pyt")
version_str = ".".join(platform.python_version().split(".")[:2])
site_lib = os.path.join(install_lib, "lib", 'python' + version_str, "site-packages")
site.addsitedir(site_lib)
os.environ['PYTHONPATH'] += ":{}".format(install_lib)
os.environ['PYTHONPATH'] += ":{}".format(site_lib)

assert os.path.exists(site_lib)

print(os.environ['PYTHONPATH'])