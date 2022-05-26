import os

try:
    from mpi4py import MPI

    use_mpi = True
    comm = MPI.COMM_WORLD
    barrier = comm.Barrier
    finalize = MPI.Finalize
    rank = comm.Get_rank()
    size = comm.Get_size()

    if rank == 0:
        print("-- MPI import successful")

except ImportError:

    MPI = None

    print("-- Unable to import mpi4py, procs will be serial")
    use_mpi = False
    rank = 0
    size = 1


    def barrier():
        pass


def choose_items(*items):
    """ Split list of objects into sublists of ~uniform length across avail. procs
        then return list item indexed by proc. rank """

    items = list(items)

    r = rank
    s = size
    l = len(items)
    min_items = l // s  # Populates list with min. number of items chose by int. div.
    rem_items = l % s  # Find num. of remaining items with mod. div. and pop. lists

    items_packets = []

    for i in range(s):
        packet = items[i * min_items: (i + 1) * min_items]
        items_packets.append(packet)

    if rem_items > 0:
        for j in range(rem_items):
            extra = items[min_items * s + j]
            items_packets[j].append(extra)

    return items_packets[r]


def single_proc_task(func, *args, **kwargs):
    if rank == 0:
        func(*args, **kwargs)
    barrier()


def make_dir(dir_name):
    if not os.path.exists(dir_name):
        if rank == 0:
            os.makedirs(dir_name)
    barrier()
