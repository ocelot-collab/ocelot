#!/bin/bash
 
# Example job script for GENESIS2
# Christoph Lechner, European XFEL GmbH

# MpiQueueLauncher class uses command line parameters to override some
# of these parameters, if needed.
# Currently applies to partition, nodes, and ntasks-per-node.

#SBATCH --partition=maxwell
#SBATCH --time=0-23:59:00		# Time limit (0 days, 23h59min)
#SBATCH --exclude=max-wn025 # exclude node that doesn't work all the time (2021-02-05)

### BEGIN: requested job parameters
# Note that if the job cannot be scheduled (for instance because
# the requested number of nodes is beyond your permissions), the
# job may wait forever.
#
# number of nodes (=machines) requested
#SBATCH --nodes=2			# Number of nodes
# Requesting fixed number of processors per node (if a machine has more CPUs, some will remain idle)
#SBATCH --ntasks-per-node=32		# Number of threads per node
# memory per node
#SBATCH --mem 200G			# for most cases, 200GB of RAM is ok
### END: requested job parameters

#SBATCH --job-name  G2
#SBATCH --output    jobout/G2-%j-%N.out	# File to which STDOUT will be written
#SBATCH --error     jobout/G2-%j-%N.err	# File to which STDERR will be written
#SBATCH --mail-type END			# Type of email notification- BEGIN,END,FAIL,ALL
# per default email goes to job submitter -> no need to specify email address
# SBATCH --mail-user jane.doe@desy.de   # Email to which notifications will be sent

# Check if mpirun is available. If not, load the module
source /etc/profile.d/modules.sh
if ! [ -x "$(command -v mpirun)" ]; then
	module load mpi/openmpi-x86_64
fi

# these files can be used to signal the status of the simulation to other processes
rm -rf flag.start flag.finish

# when your simulations are crashing, a coredump can be useful
# ulimit -c unlimited

# launch the job
touch flag.start
mpirun /gpfs/exfel/data/group/wp72/_code/products/genesis/genesis < tmp.cmd
touch flag.finish
