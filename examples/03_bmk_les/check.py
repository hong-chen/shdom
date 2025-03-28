import os
import sys
import struct
import glob


def check(fname, fdir_ref='data/output', data_len=4):

    filename = os.path.basename(fname)
    fname_ref = os.path.join(fdir_ref, filename)

    file_size = os.path.getsize(fname)
    file_size_ref = os.path.getsize(fname_ref)

    if file_size != file_size_ref:
        msg = 'Error [check]: File size does NOT match (%d byte vs %d byte).' % (file_size, file_size_ref)
        print(msg)
        return

    Ndata = file_size//data_len

    with open(fname, 'rb') as f, open(fname_ref, 'rb') as f_ref:

        for idata in range(Ndata):

            data = struct.unpack('=f', f.read(data_len))
            data_ref = struct.unpack('<f', f_ref.read(data_len))

            if data != data_ref:
                msg = 'Error [check]: iData=%d, Data=%f, Data_Ref=%f, Diff=%f' % (idata+1, data, data_ref, data-data_ref)
                print(msg)
                return

    print('Message [check]: Passed file check for <%s>.' % fname)
    return


if __name__ == '__main__':

    fnames = sorted(glob.glob('tmp-data/shdom-out*.sHdOm-out*'))

    for fname in fnames:
        check(fname)

    pass
