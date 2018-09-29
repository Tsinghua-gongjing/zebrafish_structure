
#   http://github.com/pixelb/scripts/commits/master/scripts/ps_mem.py

#
#
#   Estimate Memory Usage of Current Process
#
#
#

import getopt
import time
import errno
import os
import sys


class Proc:
    def __init__(self):
        uname = os.uname()
        if uname[0] == "FreeBSD":
            self.proc = '/compat/linux/proc'
        else:
            self.proc = '/proc'
    def path(self, *args):
        return os.path.join(self.proc, *(str(a) for a in args))
    def open(self, *args):
        try:
            if sys.version_info < (3,):
                #print self.path(*args)
                return open(self.path(*args))
            else:
                #print self.path(*args)
                return open(self.path(*args), errors='ignore')
        except (IOError, OSError):
            val = sys.exc_info()[1]
            if (val.errno == errno.ENOENT or # kernel thread or process gone
                val.errno == errno.EPERM or
                val.errno == errno.EACCES):
                raise LookupError
            raise

proc = Proc()

def kernel_ver():
    kv = proc.open('sys/kernel/osrelease').readline().split(".")[:3]
    last = len(kv)
    if last == 2:
        kv.append('0')
    last -= 1
    while last > 0:
        for char in "-_":
            kv[last] = kv[last].split(char)[0]
        try:
            int(kv[last])
        except:
            kv[last] = 0
        last -= 1
    return (int(kv[0]), int(kv[1]), int(kv[2]))

def getMemStats(pid):
    PAGESIZE = os.sysconf("SC_PAGE_SIZE") / 1024 #KiB
    have_pss = 0
    have_swap_pss = 0
    mem_id = pid #unique
    Private_lines = []
    Shared_lines = []
    Pss_lines = []
    Rss = (int(proc.open(pid, 'statm').readline().split()[1])
           * PAGESIZE)
    Swap_lines = []
    Swap_pss_lines = []
    #
    Swap = 0
    Swap_pss = 0
    #
    if os.path.exists(proc.path(pid, 'smaps')):  # stat
        lines = proc.open(pid, 'smaps').readlines()  # open
        # Note we checksum smaps as maps is usually but
        # not always different for separate processes.
        mem_id = hash(''.join(lines))
        for line in lines:
            if line.startswith("Shared"):
                Shared_lines.append(line)
            elif line.startswith("Private"):
                Private_lines.append(line)
            elif line.startswith("Pss"):
                have_pss = 1
                Pss_lines.append(line)
            elif line.startswith("Swap:"):
                Swap_lines.append(line)
            elif line.startswith("SwapPss:"):
                have_swap_pss = 1
                Swap_pss_lines.append(line)
        Shared = sum([int(line.split()[1]) for line in Shared_lines])
        Private = sum([int(line.split()[1]) for line in Private_lines])
        #Note Shared + Private = Rss above
        #The Rss in smaps includes video card mem etc.
        if have_pss:
            pss_adjust = 0.5 # add 0.5KiB as this avg error due to truncation
            Pss = sum([float(line.split()[1])+pss_adjust for line in Pss_lines])
            Shared = Pss - Private
        # Note that Swap = Private swap + Shared swap.
        Swap = sum([int(line.split()[1]) for line in Swap_lines])
        if have_swap_pss:
            # The kernel supports SwapPss, that shows proportional swap share.
            # Note that Swap - SwapPss is not Private Swap.
            Swap_pss = sum([int(line.split()[1]) for line in Swap_pss_lines])
    elif (2,6,1) <= kernel_ver() <= (2,6,9):
        Shared = 0 #lots of overestimation, but what can we do?
        Private = Rss
    else:
        Shared = int(proc.open(pid, 'statm').readline().split()[2])
        Shared *= PAGESIZE
        Private = Rss - Shared
    return (Private, Shared, mem_id, Swap, Swap_pss)


def human(num, power="Ki", units=None):
    if units is None:
        powers = ["Ki", "Mi", "Gi", "Ti"]
        while num >= 1000: #4 digits
            num /= 1024.0
            power = powers[powers.index(power)+1]
        return "%.1f %sB" % (num, power)
    else:
        return "%.f" % ((num * 1024) / units)


def memory():
    current_pid = os.getpid()
    try:
        (Private, Shared, mem_id, Swap, Swap_pss) = getMemStats(current_pid)
    except LookupError:
        return " unpermitted in this OS "
    else:
        return human(Private)


