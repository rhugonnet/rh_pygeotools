# -*- coding: utf-8 -*-
"""
Created on Mon Feb 05 11:55:37 2018

@author: hugonnet

SHELL LIBRARY:
Library of Python functions for file, directory and path manipulation

LIST:
"""
from __future__ import print_function
import os, sys
import shutil
import tarfile, zipfile
from contextlib import contextmanager
import errno
import traceback
# from types import ModuleType, FunctionType
# from gc import get_referents
# from pympler import asizeof

def create_tmp_dir_for_outfile(file_out):

    tmp_dir = os.path.join(os.path.dirname(file_out), 'tmp_'+os.path.splitext(os.path.basename(file_out))[0]) + os.path.sep
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)

    return tmp_dir

def remove_tmp_dir_for_outfile(file_out):

    tmp_dir = os.path.join(os.path.dirname(file_out), 'tmp_'+os.path.splitext(os.path.basename(file_out))[0]) + os.path.sep
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir,ignore_errors=True)

def make_pdirs(dir_path):

    outdir = os.path.abspath(dir_path)
    try:
        os.makedirs(outdir)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(outdir):
            pass
        else:
            raise

def add_suf_before_ext(basename,suf):

    filename=os.path.splitext(basename)[0]
    ext=os.path.splitext(basename)[1]

    sufname=filename+suf+ext

    return sufname

def extract_file_from_tar_gz(tar_in,filename_in,file_out):

    #ref: https://stackoverflow.com/questions/37752400/how-do-i-extract-only-the-file-of-a-tar-gz-member
    with tarfile.open(tar_in, "r") as tar:
        counter = 0

        for member in tar:
            if member.isfile():
                filename = os.path.basename(member.name)
                if filename != filename_in:  # do your check
                    continue

                with open(file_out, "wb") as output:
                    print('Extracting '+filename_in + ' from archive '+tar_in+' to '+file_out+'...')
                    shutil.copyfileobj(tar.fileobj, output, member.size)


                break  # got our file

            counter += 1
            if counter % 1000 == 0:
                tar.members = []  # free ram... yes we have to do this manually


def extract_file_from_zip(zip_in,filename_in,file_out):

    #ref: https://stackoverflow.com/questions/4917284/extract-files-from-zip-without-keeping-the-structure-using-python-zipfile
    with zipfile.ZipFile(zip_in) as zip_file:

        for member in zip_file.namelist():
            # if member.isfile():
            filename = os.path.basename(member)
            if filename != filename_in:
                # skip directories
                continue

            # copy file (taken from zipfile's extract)
            source = zip_file.open(member)
            target = open(file_out, "wb")
            with source, target:
                shutil.copyfileobj(source, target)


#redirecting stdout and stderr, source: https://stackoverflow.com/questions/4675728/redirect-stdout-to-a-file-in-python/22434262#22434262
def fileno(file_or_fd):
    fd = getattr(file_or_fd, 'fileno', lambda: file_or_fd)()
    if not isinstance(fd, int):
        raise ValueError("Expected a file (`.fileno()`) or a file descriptor")
    return fd

@contextmanager
def stdout_redirected(to=os.devnull, stdout=None):
    if stdout is None:
       stdout = sys.stdout

    stdout_fd = fileno(stdout)
    # copy stdout_fd before it is overwritten
    #NOTE: `copied` is inheritable on Windows when duplicating a standard stream
    with os.fdopen(os.dup(stdout_fd), 'wb') as copied:
        stdout.flush()  # flush library buffers that dup2 knows nothing about
        try:
            os.dup2(fileno(to), stdout_fd)  # $ exec >&to
        except ValueError:  # filename
            with open(to, 'wb') as to_file:
                os.dup2(to_file.fileno(), stdout_fd)  # $ exec > to
        try:
            yield stdout # allow code to be run with the redirected stdout
        finally:
            # restore stdout to its previous value
            #NOTE: dup2 makes stdout_fd inheritable unconditionally
            stdout.flush()
            os.dup2(copied.fileno(), stdout_fd)  # $ exec >&copied

def merged_stderr_stdout():  # $ exec 2>&1
    return stdout_redirected(to=sys.stdout, stdout=sys.stderr)

# def getsize(obj):
#     BLACKLIST = type, ModuleType, FunctionType
#     #source: https://stackoverflow.com/questions/449560/how-do-i-determine-the-size-of-an-object-in-python
#     # Custom objects know their class.
#     # Function objects seem to know way too much, including modules.
#     # Exclude modules as well.
#     """sum size of object & members."""
#     if isinstance(obj, BLACKLIST):
#         raise TypeError('getsize() does not take argument of type: '+ str(type(obj)))
#     seen_ids = set()
#     size = 0
#     objects = [obj]
#     while objects:
#         need_referents = []
#         for obj in objects:
#             if not isinstance(obj, BLACKLIST) and id(obj) not in seen_ids:
#                 seen_ids.add(id(obj))
#                 size += sys.getsizeof(obj)
#                 need_referents.append(obj)
#         objects = get_referents(*need_referents)
#     return size/1000000.

# def getsizeobj(obj):
#     return asizeof.asizeof(obj)

# def some_function_with_cached_sys_stdout(stdout=sys.stderr):
#     print('cached stdout', file=stdout)

# def some_function_with_stdout_and_error():
#     print('stdout in function')
#     x = 1/0
#
#
# with open('/home/atom/ongoing/test.log','w') as f:
#
#     with stdout_redirected(to=f), merged_stderr_stdout():
#         print('stdout goes to devnull')
#         print('stderr also goes to stdout that goes to devnull', file=sys.stderr)
#         try:
#             some_function_with_stdout_and_error()
#         except Exception:
#             print(traceback.format_exc())
#
#     print('stdout is back')
#     some_function_with_stdout_and_error()
#     print('stderr is back', file=sys.stderr)