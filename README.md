fast5tofastq
============

This program (fast5tofastq) opens a HDF5 file and checks the paths

- /Analyses/Basecall_2D_000/BaseCalled_2D/Fastq
- /Analyses/Basecall_2D_000/BaseCalled_template/Fastq
- /Analyses/Basecall_2D_000/BaseCalled_complement/Fastq

for string objects in the given order. When it finds any such object it
prints the string contents on the standard output channel and quits. If it
finds no such element than it quits without printing anything on standard
output.

Compilation requires a HDF5 installation in version >= 1.8.13 
(see http://www.hdfgroup.org/HDF5/release/obtainsrc.html#src),
a C compiler and a make program.

A non-standard installation path for HDF5 can be given as an argument for
make, e.g.

	make HDF5=/path/to/installed/hdf5-1.8.13

fast5tofastq is distributed under version 3 of the GNU General Public
License (see the file GPL-3).
