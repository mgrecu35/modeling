p1='/PANFS/user/home/sciteam/data/johnk/rh8_orr/combined/20230801/'

import glob

fs1=sorted(glob.glob(p1+'src/*.f90'))
fs2=sorted(glob.glob('src/*.f90'))
