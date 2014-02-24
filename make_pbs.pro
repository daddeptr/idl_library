pro make_pbs, runtag, localdir
;## runtag = "pix_3-2co_bp_detMap"
;## runtag = "pix_2-1co_bp_detMap"

   restore, localdir+runtag+".par.sav"

   spawn, 'mv '+localdir+'*.tmp ~dpietrob/rubbish'

   openw,1,localdir+'run_cv.pbs'
   printf,1,"#PBS -S /bin/bash"
   printf,1,"#PBS -V"
   printf,1,"#PBS -q usplanck"
   printf,1,"#PBS -l walltime=01:33:00"
;## printf,1,"#PBS -l nodes="+strtrim(string(nchannels),2)+":ppn="+strtrim(string(8/nprocperband),2)+",pvmem="+strtrim(string(20000/(8/nprocperband)),2)+"mb"
   printf,1,"#PBS -l nodes="+strtrim(string(nchannels),2)+":ppn=4,pvmem=5000mb"
   printf,1,"#PBS -A usplanck"
   printf,1,"### #PBS -l advres=usplanck.6428"
   printf,1,"#PBS -j eo"
   printf,1,"#PBS -N dx9_ns128_"+runtag+".job"
   printf,1,"#PBS -o dx9_ns128_"+runtag+".log"
   printf,1,"#PBS -m be"
   printf,1,"#PBS -M davide.pietrobon@jpl.nasa.gov"
   printf,1,""
   printf,1,"cd $PBS_O_WORKDIR"
   printf,1,""
   printf,1," . /project/projectdirs/cmb/modules/carver/hpcports.sh"
   printf,1," hpcports gnuk"
   printf,1," module load hpcp"
   printf,1," module load planck"
   printf,1," MACHINE=hpcports"
   printf,1," CMD_DIR=/project/projectdirs/planck/user/dpietrob/software/commander/install_${MACHINE}/bin/"
   printf,1," ### CMD_DIR=/project/projectdirs/planck/user/dpietrob/software/commander_latest/install_${MACHINE}/bin/"
   printf,1,""
   printf,1," module list"
   printf,1,""
   printf,1,"### cmbenv gnu ;module load cmb ;module load planck ;module load python ;module load cmbdev"
   printf,1,"### MACHINE=carver.gnu"
   printf,1,"### ### MACHINE=carver.gnu.ddt"
   printf,1,"### CMD_DIR=/project/projectdirs/planck/user/dpietrob/software/mycommander/install_${MACHINE}/bin/"
   printf,1,""
   printf,1,"### PARDIR=chains/pix_3co_bp_detMap/"
   printf,1,"### PARFILE=param_dx9_pix_3co.txt"
   printf,1,"PARDIR="+localdir
   printf,1,"PARFILE=param_"+runtag+".txt"
   printf,1,""
   printf,1,"mpirun -np "+strtrim(string(nchannels*nprocperband),2)+" ${CMD_DIR}commander ${PARDIR}${PARFILE}"
   close,1
end
