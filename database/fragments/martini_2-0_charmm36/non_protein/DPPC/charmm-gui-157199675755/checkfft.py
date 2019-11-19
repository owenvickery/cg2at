from __future__ import print_function
import sys

def usage():
    print("""return FFT grid size for PME

usage:
# checkfft.py A B C
""")

def is_factor(n):
    if (n % 2 != 0): return False  # favors even number
    while n:
        flag = False
        for x in (2,3,5):
            if n % x == 0:
               n = n / x
               flag = True
               break

        if flag: continue
        break

    if n == 1: return True
    return False

def checkfft(n, margin = 5):
    n = int(n) + margin
    while 1:
        if is_factor(n): break
        else: n += 1
    return n

if __name__ == '__main__':
    try:
        xtla = float(sys.argv[1])
        xtlb = float(sys.argv[2])
        xtlc = float(sys.argv[3])
    except:
        usage()
        sys.exit()

    fftx = checkfft(xtla)
    fftz = checkfft(xtlc)

    if xtla == xtlb:
        ffty = fftx
    else: ffty = checkfft(xtlb)

    print("""set fftx %d
set ffty %d
set fftz %d
""" % (fftx, ffty, fftz))

