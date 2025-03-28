import os
import sys
import struct
import glob


def check(fname, fname_ref, data_len=4, threshold=1e-4):

    if not os.path.exists(fname):
        msg = 'Error [check]: File <%s> does NOT exist.' % (fname)
        print(msg)
        return

    if not os.path.exists(fname_ref):
        msg = 'Error [check]: File <%s> does NOT exist.' % (fname_ref)
        print(msg)
        return

    file_size = os.path.getsize(fname)
    file_size_ref = os.path.getsize(fname_ref)

    if file_size != file_size_ref:
        msg = 'Error [check]: File size does NOT match (%d byte vs %d byte).' % (file_size, file_size_ref)
        print(msg)
        return

    Ndata = file_size//data_len

    with open(fname, 'rb') as f, open(fname_ref, 'rb') as f_ref:

        for idata in range(Ndata):

            data = struct.unpack('=f', f.read(data_len))[0]
            data_ref = struct.unpack('<f', f_ref.read(data_len))[0]

            diff = data - data_ref
            if abs(diff)>threshold:
                msg = 'Error [check]: iData=%d, Data=%f, Data_Ref=%f, Diff=%f' % (idata+1, data, data_ref, diff)
                print(msg)
                return

    print('Message [check]: Passed file check for <%s>.' % fname)
    return


if __name__ == '__main__':

    fdir_ref='data/output'

    fnames = sorted(glob.glob('tmp-data/shdom-out*.sHdOm-out*'))

    for fname in fnames:
        filename = os.path.basename(fname)
        fname_ref = os.path.join(fdir_ref, filename)
        check(fname, fname_ref)

    pass
