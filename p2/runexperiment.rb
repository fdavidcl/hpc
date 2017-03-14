#!/usr/bin/env ruby

header = <<EOF
num_procs, num_ints, err_excess, err_defect, err_middle, t_calc, t_comms, t_create, t_destroy
EOF

if ARGV[0] == "atcgrid"
  `scp pi-mpi acap4@atcgrid.ugr.es:/home/acap4`
end

outlines = [1, 2, 3, 6].map do |proc_num|
  (1 .. 9).map do |exp|
    if ARGV[0] == "atcgrid"
      cmd = %(ssh acap4@atcgrid.ugr.es 'echo "/usr/lib64/openmpi/bin/mpiexec -n #{proc_num} pi-mpi #{10 ** exp}" | qsub -q acap')
      puts cmd
      `#{cmd}`
      sleep 13
    else
      `mpirun -n #{proc_num} pi-mpi #{10 ** exp}`
    end
  end
end.flatten

File.write("out.csv", header + outlines.join("\n"))
