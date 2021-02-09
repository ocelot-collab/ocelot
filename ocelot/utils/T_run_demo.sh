#!/bin/bash
 
# Example job script (demo)
# Christoph Lechner, European XFEL GmbH

# MpiQueueLauncher class uses command line parameters to override some
# of these parameters, if needed.
# Currently applies to partition, nodes, and ntasks-per-node.

#SBATCH --partition=maxwell
#SBATCH --time=0-23:59:00		# Time limit (0 days, 23h59min)

### BEGIN: requested job parameters
# Note that if the job cannot be scheduled (for instance because
# the requested number of nodes is beyond your permissions), the
# job may wait forever.
#
# number of nodes (=machines) requested
#SBATCH --nodes=1			# Number of nodes
# Requesting fixed number of processors per node (if a machine has more CPUs, some will remain idle)
#SBATCH --ntasks-per-node=32		# Number of threads per node
# memory per node
#SBATCH --mem 200G			# for most cases, 200GB of RAM is ok
### END: requested job parameters

#SBATCH --job-name  "queue demo"

# simple demo script => don't use jobout sub-directory, but write to cwd
#SBATCH --output    demo-%j-%N.out	# File to which STDOUT will be written
#SBATCH --error     demo-%j-%N.err	# File to which STDERR will be written

#SBATCH --mail-type END			# Type of email notification- BEGIN,END,FAIL,ALL
# per default email goes to job submitter -> no need to specify email address
# SBATCH --mail-user jane.doe@desy.de   # Email to which notifications will be sent

# Check if mpirun is available. If not, load the module
source /etc/profile.d/modules.sh
if ! [ -x "$(command -v mpirun)" ]; then
	module load mpi/openmpi-x86_64
fi

rm flag.finish

NODE=`uname -n`
sleep 60  # some wait time to demonstrate transitions of job state...
echo "Hallo from compute node $NODE" > T_run_demo_out.txt

touch flag.finish
