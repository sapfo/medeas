import options

options.TESTING = True  # Shows some debugging info and many plots
options.FST = False  # Also calculates F_ST
options.BOOTRUNS = 50  # How many bootstrap runs we need.
options.BOOTSIZE  = 100
options.SIMULATION = True  # Are we running on simulated datas
options.K_OVERRIDE = 3 #number of population shoule be overwritten? k == 0 -> no; K > 0 -> number of population
options.MATHPLOTLIB_BACKEND = 'agg' #which backend should be used by mathplotlib
